//
// Created by srikanthmaturu on 10/17/2017.
//

#include "core/cluster.h"
#include "core/falconn.h"
#include "core/similar_fasta_sequences.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include "edlib.h"

using namespace std;
using namespace SequencesAnalyzer::core;

void performClustering(string fastaFile){
    ofstream resultsFile("outputfile.txt", ofstream::out);
    vector<string> sequences;
    fasta::getFastaSequences(fastaFile, sequences, 30);
    uint64_t kmerSize = 30;

    FALCONNIndexConfiguration falconnConfig;
    falconnConfig.lshType = 2;
    falconnConfig.ngl = 5;
    falconnConfig.numberOfHashBits = 11;
    falconnConfig.numberOfProbes = 2032;
    falconnConfig.threshold = 1.5;
    falconnConfig.numberOfHashTables = 32;
    falconnConfig.dataset_type = 0;
    ClusterConfiguration clusterConfig;
    clusterConfig.editDistanceThreshold = 14;
    clusterConfig.kmerSize = kmerSize;
    ClustersGenerator generator(sequences, falconnConfig, clusterConfig);
    /*for(string sequence:sequences){
        cout << sequence << endl;
    }*/

    generator.generateClusters();
    cout << "Clusters generation complete." << endl;
    for(string sequence:sequences){
        vector<uint64_t>& similarSequences = generator.getSequencesCluster(sequence, kmerSize);
        resultsFile << ">" << sequence << endl;
        for(uint64_t sequenceId: similarSequences){
            //cout << sequenceId << endl;
            //cout << sequences[sequenceId] << endl;
            resultsFile << sequences[sequenceId] << endl;
        }
    }
}

void findSimilarSequences(string databaseFile, string queryFile, uint8_t dataset_type, uint64_t maxED){
    EdlibAlignConfig edlibConfig = edlibNewAlignConfig(maxED, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
    vector<string> databaseFastaSequences, queryFastaSequences;
    fasta::getFastaSequences(databaseFile, databaseFastaSequences);
    fasta::getFastaSequences(queryFile, queryFastaSequences);

    FALCONNIndexConfiguration falconnConfig;
    falconnConfig.lshType = 2;
    falconnConfig.ngl = 5;
    falconnConfig.numberOfHashBits = 11;
    falconnConfig.numberOfProbes = 300;
    falconnConfig.threshold = 1.5;
    falconnConfig.numberOfHashTables = 32;
    falconnConfig.dataset_type = dataset_type;
    SimilarFastaSequencesFinder sequencesFinder(databaseFastaSequences, falconnConfig);
    sequencesFinder.initialize();
    ofstream resultsFile("results.txt");
    for(uint64_t i = 0; i < queryFastaSequences.size(); i++){
        auto res = sequencesFinder.match(queryFastaSequences[i]);
        cout << "Query: " << i << " Candidates: " << res.first << endl;
        vector<int32_t>& result_matches = *res.second;
        resultsFile << ">Sequence " << i << " matches" << endl;
        if(result_matches.size() == 0){
            continue;
        }
        for(uint64_t j = 0; j < result_matches.size(); j++){
            EdlibAlignResult ed_result = edlibAlign(queryFastaSequences[i].c_str(), queryFastaSequences[i].size(),
                                                    databaseFastaSequences[result_matches[j]].c_str(), databaseFastaSequences[result_matches[j]].size(), edlibConfig);
            if(ed_result.editDistance >= 0){
                resultsFile << result_matches[j] ;
                if(j < (result_matches.size() - 1)){
                    resultsFile << " , ";
                }
                else {
                    resultsFile << std::endl;
                }
            }
        }
    }
}

int main(int argc, char** argv){

    if(argc < 2){
        cout << "Usage ./executable [type_of_task] [additional_options]" << endl;
        return 1;
    }
    switch(stoi(argv[1])){
        case 0:
            if(argc < 3) {
                cout << "Usage ./executable 0 [fasta_file]" << endl;
                return 2;
            }
            performClustering(argv[2]);
            break;
        case 1:
            if(argc < 6) {
                cout << "Usage ./executable 1 [database_file] [query_file] [dataset_type] [maxED]" << endl;
                return 3;
            }
            findSimilarSequences(argv[2], argv[3], stoi(argv[4]), stoi(argv[5]));
            break;
        default:
            cout << "Invalid task. Ex task: 0 for clustering or 1 for similar fasta seq finder" << endl;
            return 100;
    }
}
