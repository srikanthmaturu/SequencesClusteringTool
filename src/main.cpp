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

void findSimilarSequences(string databaseFile, string queryFile){
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
        for(uint64_t j = 0; j < result_matches.size() - 1; j++){
            resultsFile << result_matches[j] << ", ";
        }
        resultsFile << result_matches[result_matches.size() - 1] << endl;
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
            if(argc < 4) {
                cout << "Usage ./executable 1 [database_file] [query_file]" << endl;
                return 3;
            }
            findSimilarSequences(argv[2], argv[3]);
            break;
        default:
            cout << "Invalid task. Ex task: 0 for clustering or 1 for similar fasta seq finder" << endl;
            return 100;
    }
}
