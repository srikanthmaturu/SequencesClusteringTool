//
// Created by srikanthmaturu on 10/17/2017.
//

#include "core/cluster.h"
#include "core/falconn.h"
#include "core/similar_fasta_sequences.h"
#include "core/parameters_tuner.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include <unistd.h>
#include <cstdlib>
#include "edlib.h"
#include <chrono>
#include <algorithm>

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

using namespace std;
using namespace SequencesAnalyzer::core;

void load_sequences(string sequences_file, vector<string>& sequences){
    ifstream input_file(sequences_file, ifstream::in);

    for(string sequence; getline(input_file, sequence);){
        uint64_t pos;
        if((pos=sequence.find('\n')) != string::npos){
            sequence.erase(pos);
        }
        if((pos=sequence.find('\r')) != string::npos){
            sequence.erase(pos);
        }
        sequences.push_back(sequence);
    }
}

void performClustering(string fastaFile){
    ofstream resultsFile("outputfile.txt", ofstream::out);
    vector<string> sequences;
    fasta::getFastaSequences(fastaFile, sequences, 30);
    uint64_t kmerSize = 30;

    FALCONNIndexConfiguration falconnConfig;
    falconnConfig.lshType = LSH_HASH_TYPE;
    falconnConfig.ngl = NGRAM_LENGTH;
    falconnConfig.numberOfHashBits = NUMBER_OF_HASH_BITS;
    falconnConfig.numberOfProbes = NUMBER_OF_PROBES;
    falconnConfig.threshold = THRESHOLD/100.0;
    falconnConfig.numberOfHashTables = NUMBER_OF_HASH_TABLES;
    falconnConfig.dataset_type = 0;
    falconnConfig.data_type = 0;
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

void loadFileByType(string file, string fileType, vector<string>& sequences){
    cout << "Loading " << file << " of type " << fileType << endl;
    if(fileType == "kmers"){
        load_sequences(file, sequences);
    }
    else if(fileType == "fasta"){
        fasta::getFastaSequences(file, sequences);
    }
    else{
        cout << "Wrong file type." << endl;
        exit(1);
    }
}

void findSimilarSequences(string databaseFile, string databaseFileType, string queryFile, string queryFileType, uint8_t dataset_type, uint8_t data_type, uint64_t maxED, bool parallel){
    if(!parallel){
        omp_set_num_threads(1);
    }
    EdlibAlignConfig edlibConfig = edlibNewAlignConfig(maxED, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
    vector<string> databaseFastaSequences, queryFastaSequences;

    loadFileByType(databaseFile, databaseFileType, databaseFastaSequences);
    loadFileByType(queryFile, queryFileType, queryFastaSequences);

    FALCONNIndexConfiguration falconnConfig;
    falconnConfig.lshType = LSH_HASH_TYPE;
    falconnConfig.ngl = NGRAM_LENGTH;
    falconnConfig.numberOfHashBits = NUMBER_OF_HASH_BITS;
    falconnConfig.numberOfProbes = NUMBER_OF_PROBES;
    falconnConfig.threshold = THRESHOLD/100.0;
    falconnConfig.numberOfHashTables = NUMBER_OF_HASH_TABLES;
    falconnConfig.dataset_type = dataset_type;
    falconnConfig.data_type = data_type;
    SimilarFastaSequencesFinder sequencesFinder(databaseFastaSequences, falconnConfig);
    sequencesFinder.initialize();
    ofstream resultsFile(queryFile + "_falconn_results.txt");

    uint64_t batchSize = 1000;
    uint64_t numberOfQueries = queryFastaSequences.size();
    uint64_t lastBatchSize = numberOfQueries % batchSize;
    uint64_t numberOfBatches = numberOfQueries / batchSize + ( (lastBatchSize > 0) ? 1 : 0);
    for(uint64_t i = 0; i < numberOfBatches; i++) {
        uint64_t batchBegin = i * batchSize;
        uint64_t batchEnd = (i + 1) * batchSize;
        if (i == (numberOfBatches - 1) && lastBatchSize > 0) {
            batchEnd = i * batchSize + lastBatchSize;
        }

        uint64_t curretBatchSize = (batchEnd - batchBegin);
        std::vector<std::vector<int32_t>> &blockRes = sequencesFinder.processBatch(queryFastaSequences, batchBegin, curretBatchSize);
        for (uint64_t j = 0; j < blockRes.size(); j++) {
            vector<int32_t> &result_matches = blockRes[j];
            cout << "Query: " << batchBegin + j  << " Candidates: " << result_matches.size() << endl;
            resultsFile << ">Sequence: " << batchBegin + j << "-" << queryFastaSequences[batchBegin + j].size() << std::endl;
            if (result_matches.size() == 0) {
                continue;
            }
            bool firstTruePositivesFound = false;
            for (uint64_t k = 0; k < result_matches.size(); k++) {
                EdlibAlignResult ed_result = edlibAlign(queryFastaSequences[batchBegin + j].c_str(), queryFastaSequences[batchBegin + j].size(), databaseFastaSequences[result_matches[k]].c_str(),
                                                        databaseFastaSequences[result_matches[k]].size(), edlibConfig);
                if (ed_result.editDistance >= 0) {
                    if (firstTruePositivesFound) {
                        resultsFile << ",";
                    } else {
                        firstTruePositivesFound = true;
                    }
                    resultsFile << result_matches[k] << "-" << databaseFastaSequences[result_matches[k]].size() << "-" << ed_result.editDistance;
                }
                edlibFreeAlignResult(ed_result);
            }
            if (firstTruePositivesFound) {
                resultsFile << std::endl;
            }
        }
    }
}

std::vector<std::vector<pair<int32_t, int16_t>>>& processBatchByEdlib(std::vector<std::string>& referenceSequences, std::vector<std::string>& querySequences, uint64_t offset, uint64_t batchSize, EdlibAlignConfig edlibConfig, bool parallel){
    std::vector<std::vector<pair<int32_t, int16_t>>> * batchResults = new std::vector<std::vector<pair<int32_t, int16_t>>>(batchSize, std::vector<pair<int32_t, int16_t>>());
    if(parallel){
#pragma omp parallel for
        for(uint64_t i = offset; i < offset + batchSize; i++){
            std::vector<pair<int32_t,int16_t>>* nearestNeighbours = new std::vector<pair<int32_t,int16_t>>();
            for(int32_t j = 0; j <(int32_t)referenceSequences.size(); j++){
                EdlibAlignResult ed_result = edlibAlign(querySequences[i].c_str(), querySequences[i].size(),
                                                        referenceSequences[j].c_str(), referenceSequences[j].size(), edlibConfig);
                if(ed_result.editDistance >= 0){
                    nearestNeighbours->push_back(make_pair((int32_t)j, (int16_t)ed_result.editDistance));
                    //std::cout << j << " - " << ed_result.editDistance << std::endl;
                }
                edlibFreeAlignResult(ed_result);
            }
            (*batchResults)[i - offset] = *nearestNeighbours;
        }
    }
    else {
        for(uint64_t i = offset; i < offset + batchSize; i++){
            std::vector<pair<int32_t,int16_t>>* nearestNeighbours = new std::vector<pair<int32_t,int16_t>>();
            for(int32_t j = 0; j <(int32_t)referenceSequences.size(); j++){
                EdlibAlignResult ed_result = edlibAlign(querySequences[i].c_str(), querySequences[i].size(),
                                                        referenceSequences[j].c_str(), referenceSequences[j].size(), edlibConfig);
                if(ed_result.editDistance >= 0){
                    nearestNeighbours->push_back(make_pair((int32_t)j, (int16_t)ed_result.editDistance));
                    //std::cout << j << " - " << ed_result.editDistance << std::endl;
                }
                edlibFreeAlignResult(ed_result);
            }
            (*batchResults)[i - offset] = *nearestNeighbours;
        }
    }

    return *batchResults;
}

void findSimilarSequencesByEdlib(string databaseFile, string databaseFileType, string queryFile, string queryFileType, uint64_t maxED, bool parallel){
    EdlibAlignConfig edlibConfig = edlibNewAlignConfig(maxED, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
    vector<string> databaseFastaSequences, queryFastaSequences;
    loadFileByType(databaseFile, databaseFileType, databaseFastaSequences);
    loadFileByType(queryFile, queryFileType, queryFastaSequences);

    ofstream resultsFile(queryFile + "_edlib_results.txt");
    uint64_t batchSize = 1000;
    uint64_t numberOfQueries = queryFastaSequences.size();
    uint64_t lastBatchSize = numberOfQueries % batchSize;
    uint64_t numberOfBatches = numberOfQueries / batchSize + ( (lastBatchSize > 0) ? 1 : 0);
    for(uint64_t i = 0; i < numberOfBatches; i++){
        uint64_t batchBegin = i * batchSize;
        uint64_t batchEnd = (i + 1) * batchSize;
        if(i == (numberOfBatches - 1) && lastBatchSize > 0){
            batchEnd =  i * batchSize + lastBatchSize;
        }

        uint64_t curretBatchSize = (batchEnd - batchBegin);
        std::vector<std::vector<pair<int32_t, int16_t>>>& res = processBatchByEdlib(databaseFastaSequences, queryFastaSequences, batchBegin, curretBatchSize, edlibConfig, parallel);
        cout << "Processed batch: " << batchBegin << " - " << batchEnd << endl;
        for(uint64_t j = 0; j < res.size(); j++){
            resultsFile << ">Sequence: " << batchBegin + j << "-" << queryFastaSequences[batchBegin + j].size() << std::endl;
            if(res[j].size() == 0 ){
                continue;
            }
            for(uint64_t k = 0; k < res[j].size() - 1; k++){
                resultsFile << res[j][k].first << "-" << databaseFastaSequences[res[j][k].first].size() << "-" << res[j][k].second << ",";
            }
            if(res[j].size() > 0){
                resultsFile << res[j][res[j].size() - 1].first << "-" << databaseFastaSequences[res[j][res[j].size() - 1].first].size() << "-" << res[j][res[j].size() - 1].second << std::endl;
            }
        }
    }
}

void summarizeResults(string resultsFile) {
    ifstream file(resultsFile);
    vector<vector<pair<int32_t, int16_t>>> results;

    regex e("^>");
    smatch m;

    uint64_t j = 0;
    while(!file.eof()){
        std::string line;
        std::getline(file, line);
        if(!regex_search(line, e)){
            uint64_t pos;
            if((pos=line.find('\n')) != string::npos){
                line.erase(pos);
            }
            if((pos=line.find('\r')) != string::npos){
                line.erase(pos);
            }
            stringstream ss(line);
            string value;
            while(getline(ss, value, ',')){
                int32_t ind;
                int16_t ed;
                uint8_t pos = value.find('-');
                ind = stoi(value.substr(0, pos));
                pos = value.find('-', pos + 1);
                ed = stoi(value.substr(pos+1));
                results[j-1].push_back(make_pair(ind, ed));
            }
        } else {
            j++;
            results.push_back(vector<pair<int32_t,int16_t>>());
        }
    }

    vector<uint64_t> matches;
    vector<bool> matchesFound;
    for(uint64_t i = 0 ; i < results.size(); i++){
        for(uint64_t j = 0; j < results[i].size(); j++){
            if(results[i][j].second >= (int16_t)matches.size()){
                matches.resize(results[i][j].second + 1, 0);
                matchesFound.resize(results[i][j].second + 1, false);
            }
            if(!matchesFound[results[i][j].second]){
                matchesFound[results[i][j].second] = true;
                matches[results[i][j].second]++;
            }
        }
        std::fill(matchesFound.begin(), matchesFound.end(), false);
    }

    vector<int64_t> categoryCounts;

    for(uint64_t j = 0; j < results.size(); j++){
        for(pair<int32_t, int16_t> p: results[j]){
            uint64_t category = (p).second / 5;
            if(categoryCounts.size() < (category + 1)){
                categoryCounts.resize(category + 1, 0);
            }
            categoryCounts[category]++;
        }
    }

    cout << "Printing categorized results: " << endl;
    for(uint64_t i = 0; i < categoryCounts.size(); i++){
        cout << "Ed: [" << i*5 << " - " << (i + 1) * 5 << ") " <<  "  Count: " << categoryCounts[i] << endl;
    }
    cout << "Printing detailed results: " << endl;
    for(uint64_t i = 0; i < matches.size(); i++){
        cout << "Ed-" << i << " Count: " << matches[i] << endl;
    }
}

void compareFalconnAndEdlibResults(string falconnResults, string edlibResults){
    ifstream files[2];
    files[0] = ifstream(falconnResults);
    files[1] = ifstream(edlibResults);
    vector<vector<pair<int32_t, int16_t>>> results[2];

    regex e("^>");
    smatch m;

    for(uint8_t i = 0; i < 2; i++){
        uint64_t j = 0;
        while(!files[i].eof()){
            std::string line;
            std::getline(files[i], line);
            if(!regex_search(line, e)){
                uint64_t pos;
                if((pos=line.find('\n')) != string::npos){
                    line.erase(pos);
                }
                if((pos=line.find('\r')) != string::npos){
                    line.erase(pos);
                }
                stringstream ss(line);
                string value;
                while(getline(ss, value, ',')){
                    int32_t ind;
                    int16_t ed;
                    uint8_t pos = value.find('-');
                    ind = stoi(value.substr(0, pos));
                    pos = value.find('-', pos + 1);
                    ed = stoi(value.substr(pos+1));
                    results[i][j-1].push_back(make_pair(ind, ed));
                }
            } else {
                j++;
                results[i].push_back(vector<pair<int32_t,int16_t>>());
            }
        }
    }

    vector<int64_t> categoryCounts[2];

    for(uint64_t i = 0; i < 2; i++){
        for(uint64_t j = 0; j < results[i].size(); j++){
            for(pair<int32_t, int16_t> p: results[i][j]){
                uint64_t category = (p).second / 5;
                if(categoryCounts[i].size() < (category + 1)){
                    categoryCounts[i].resize(category + 1, 0);
                }
                categoryCounts[i][category]++;
            }
        }
    }

    if(categoryCounts[0].size() > categoryCounts[1].size()){
        categoryCounts[1].resize(categoryCounts[0].size(), 0);
    }
    else if(categoryCounts[1].size() > categoryCounts[0].size()) {
        categoryCounts[0].resize(categoryCounts[1].size(), 0);
    }

    for(uint64_t i = 0; i < categoryCounts[0].size(); i++){
        cout << "Ed: [" << i*5 << " - " << (i + 1) * 5 << ") " <<  "  FR: " << categoryCounts[0][i] << " ER: " << categoryCounts[1][i]  << " diff: " << (categoryCounts[1][i] - categoryCounts[0][i]) << endl;
    }
}

void generateSampleFiles(string falconnResults, string edlibResults, string databaseFile, string databaseFileType, string queryFile, string queryFileType, uint64_t numberOfSequences){
    ifstream files[2];
    files[0] = ifstream(falconnResults);
    files[1] = ifstream(edlibResults);
    vector<vector<pair<int32_t, int16_t>>> results[2];
    vector<string> databaseSequences, querySequences;
    loadFileByType(databaseFile, databaseFileType, databaseSequences);
    loadFileByType(queryFile, queryFileType, querySequences);

    regex e("^>");
    smatch m;

    for(uint8_t i = 0; i < 2; i++){
        uint64_t j = 0;
        while(!files[i].eof()){
            std::string line;
            std::getline(files[i], line);
            if(!regex_search(line, e)){
                uint64_t pos;
                if((pos=line.find('\n')) != string::npos){
                    line.erase(pos);
                }
                if((pos=line.find('\r')) != string::npos){
                    line.erase(pos);
                }
                stringstream ss(line);
                string value;
                while(getline(ss, value, ',')){
                    int32_t ind;
                    int16_t ed;
                    uint8_t pos = value.find('-');
                    ind = stoi(value.substr(0, pos));
                    pos = value.find('-', pos + 1);
                    ed = stoi(value.substr(pos+1));
                    results[i][j-1].push_back(make_pair(ind, ed));
                }
            } else {
                j++;
                results[i].push_back(vector<pair<int32_t,int16_t>>());
            }
        }
    }
    cout << "Read edlib and falconn results file" << endl;
    ofstream sampleQueryFile("query.txt"), sampleDatabaseFile("database.txt");
    uint64_t selectedSequencesSize = 0;
    cout << "Generating " << numberOfSequences << " sequences" << endl;
    srand (time(NULL));
    vector<int32_t> selectedDatabaseSequences;
    while(selectedSequencesSize < numberOfSequences){
        uint64_t queryIndex = rand() % (querySequences.size());
        int32_t databaseIndex = 0;
        if(results[1][queryIndex].size() > 1){
            uint64_t resultSelection = rand() % results[1][queryIndex].size();
            databaseIndex = results[1][queryIndex][resultSelection].first;
        }
        if(find(selectedDatabaseSequences.begin(), selectedDatabaseSequences.end(), databaseIndex) != selectedDatabaseSequences.end()){
            continue;
        }
        sampleQueryFile << querySequences[queryIndex] << std::endl;
        sampleDatabaseFile << databaseSequences[databaseIndex] << std::endl;
        selectedSequencesSize++;
    }
}

void performParameterTuningtests(string databaseFile, string databaseFileType, string queryFile, string queryFileType, uint8_t dataset_type, uint8_t data_type, uint64_t maxED, uint64_t testType, double_t threshold = 0.0) {
    vector<string> databaseSequences, querySequences;

    loadFileByType(databaseFile, databaseFileType, databaseSequences);
    loadFileByType(queryFile, queryFileType, querySequences);

    FALCONNIndexConfiguration falconnConfig;
    falconnConfig.lshType = LSH_HASH_TYPE;
    falconnConfig.ngl = NGRAM_LENGTH;
    falconnConfig.numberOfHashBits = NUMBER_OF_HASH_BITS;
    falconnConfig.numberOfProbes = NUMBER_OF_PROBES;
    falconnConfig.threshold = THRESHOLD / 100.0;
    falconnConfig.numberOfHashTables = NUMBER_OF_HASH_TABLES;
    falconnConfig.dataset_type = dataset_type;
    falconnConfig.data_type = data_type;
    ParametersTuner parametersTuner(databaseSequences, falconnConfig);
    parametersTuner.initialize();

    switch(testType) {
        case 0:
            cout << "Test type: Box test" << endl;
            parametersTuner.process_queries_box_test(querySequences,threshold,maxED);
            break;
        case 1:
            cout << "Test type: Thresholds test" << endl;
            parametersTuner.process_queries_thresholds_test(querySequences);
            break;
        case 2:
            cout << "Test type: Linear test" << endl;
            parametersTuner.process_queries_linear_test(querySequences);
            break;
        default:
            cout << "Invalid test type." << endl;
            break;
    }
}

int main(int argc, char** argv){

    if(argc < 2){
        cout << "Usage ./executable [type_of_task] [additional_options]" << endl;
        return 1;
    }
    auto start = timer::now();
    switch(stoi(argv[1])){
        case 0:
            if(argc < 3) {
                cout << "Usage ./executable 0 [fasta_file]" << endl;
                return 2;
            }
            performClustering(argv[2]);
            break;
        case 1:
            if(argc < 10) {
                cout << "Usage ./executable 1 [database_file] [file_type] [query_file] [file_type] [dataset_type] [data_type] [maxED] [parallel]" << endl;
                return 3;
            }

            findSimilarSequences(argv[2], argv[3], argv[4], argv[5], stoi(argv[6]), stoi(argv[7]), stoi(argv[8]), stoi(argv[9]));
            break;
        case 2:
            if(argc < 6) {
                cout << "Usage ./executable 2 [database_file] [file_type] [query_file] [file_type] [maxED] [parallel]" << endl;
                return 4;
            }
            findSimilarSequencesByEdlib(argv[2], argv[3], argv[4], argv[5], stoi(argv[6]), stoi(argv[7]));
            break;
        case 3:
            if(argc < 3) {
                cout << "Usage ./executable 3 results_file" << endl;
            }
            summarizeResults(argv[2]);
            break;
        case 4:
            if(argc < 4) {
                cout << "Usage ./executable 4 [falconn_results_file] [edlib_results_file]" << endl;
                return 5;
            }
            compareFalconnAndEdlibResults(argv[2], argv[3]);
            break;
        case 5:
            if(argc < 9) {
                cout << "Usagee ./executable 5 falconn_results_file edlib_results_file database_file file_type query_file file_type numberOfSequences " << endl;
                return 6;
            }
            generateSampleFiles(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], stoi(argv[8]));
            break;
        case 6:
            if(argc < 10) {
                cout << "Usage ./executable 6 [database_file] [file_type] [query_file] [file_type] [dataset_type] [data_type] [maxED] test_type threshold " << endl;
                return 3;
            }
            performParameterTuningtests(argv[2], argv[3], argv[4], argv[5], stoi(argv[6]), stoi(argv[7]), stoi(argv[8]), stoi(argv[9]), (argc == 11)?stoi(argv[10])/100.0:0);
            break;
        default:
            cout << "Invalid task. Ex task: 0 for clustering or 1 for similar fasta seq finder" << endl;
            return 100;
    }
    auto stop = timer::now();
    cout << "Total duration in seconds: " << (duration_cast<chrono::seconds>(stop-start).count()) << endl;
}
