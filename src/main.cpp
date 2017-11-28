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
#include <tuple>

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

void performClustering(string sequencesFile, string sequencesFileType, uint64_t dataset_type, uint64_t data_type, uint64_t percentIdentityThreshold, uint64_t stepSize){
    ofstream resultsFile(sequencesFile + ".clusters.txt", ofstream::out);
    vector<string> sequences;
    loadFileByType(sequencesFile, sequencesFileType, sequences);
    cout << "Size of sequences: " << sequences.size() << endl;
    FALCONNIndexConfiguration falconnConfig;
    falconnConfig.lshType = LSH_HASH_TYPE;
    falconnConfig.ngl = NGRAM_LENGTH;
    falconnConfig.numberOfHashBits = NUMBER_OF_HASH_BITS;
    falconnConfig.numberOfProbes = 32;
    //falconnConfig.numberOfProbes = NUMBER_OF_PROBES;
    falconnConfig.threshold = THRESHOLD/100.0;
    falconnConfig.numberOfHashTables = NUMBER_OF_HASH_TABLES;
    falconnConfig.dataset_type = dataset_type;
    falconnConfig.data_type = data_type;
    ClusterConfiguration clusterConfig;
    clusterConfig.percentIdentityThreshold = percentIdentityThreshold;
    ClustersGenerator generator(sequences, falconnConfig, clusterConfig);
    /*for(string sequence:sequences){
        cout << sequence << endl;
    }*/
    generator.initialize();
    generator.generateClusters(stepSize);
    cout << "Clusters generation complete." << endl;
    generator.printClusters(resultsFile);
}

std::vector<pair<int32_t, int16_t>>& processFALCONNCandidatesByEdlib(std::vector<std::string>& referenceSequences, std::string querySequence, std::vector<int32_t> &candidates, uint64_t minPercentIdentity, bool parallel){
    std::vector<pair<int32_t, int16_t>> * candidatesResults = new std::vector<pair<int32_t, int16_t>>(candidates.size(), pair<int32_t, int16_t>());
    std::vector<pair<int32_t, int16_t>> * finalCandidatesResults = new std::vector<pair<int32_t, int16_t>>();
    bool candidatesSelection[candidates.size()] = {false};
    if(parallel){
#pragma omp parallel for
        for(uint64_t i = 0; i < candidates.size(); i++){
            auto percentIdentity = fastPercentIdentity(querySequence, referenceSequences[candidates[i]], minPercentIdentity);
            if(percentIdentity >= minPercentIdentity){
                (*candidatesResults)[i].first = (int32_t)(candidates[i]);
                (*candidatesResults)[i].second = (int16_t)percentIdentity;
                candidatesSelection[i] = true;
                //std::cout << j << " - " << ed_result.editDistance << std::endl;
            }
        }
    }
    else {
        for(uint64_t i = 0; i < candidates.size(); i++){
            auto percentIdentity = fastPercentIdentity(querySequence, referenceSequences[candidates[i]], minPercentIdentity);
            //auto pair = getSequencesComparison(querySequence, referenceSequences[candidates[i]]);

            if(percentIdentity >= minPercentIdentity){
                (*candidatesResults)[i].first = (int32_t)(candidates[i]);
                (*candidatesResults)[i].second = (int16_t)percentIdentity;
                candidatesSelection[i] = true;
                //cout << querySequence << endl;
                //cout << referenceSequences[candidates[i]] << endl;
                //cout << pair.first << "\t" << pair.second << "\t" << (*candidatesResults)[i].second << endl;
                //std::cout << j << " - " << ed_result.editDistance << std::endl;
            }
        }
    }
    for(uint64_t i = 0; i < candidatesResults->size(); i++) {
        if(candidatesSelection[i]){
            finalCandidatesResults->push_back((*candidatesResults)[i]);
        }
    }
    return *finalCandidatesResults;
}

void findSimilarSequences(string databaseFile, string databaseFileType, string queryFile, string queryFileType, uint8_t dataset_type, uint8_t data_type, uint64_t minPercentIdentity, bool parallel){
    if(!parallel){
        omp_set_num_threads(1);
    }
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
    auto start = timer::now();
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
            std::vector<pair<int32_t, int16_t>>& finalCandidates = processFALCONNCandidatesByEdlib(databaseFastaSequences, queryFastaSequences[batchBegin + j], result_matches, minPercentIdentity, parallel);
            if(finalCandidates.size() == 0) {
                continue;
            }
            for(uint64_t k = 0; k < finalCandidates.size() - 1; k++){
                //auto sequencesSimilarity = getSequencesComparison(databaseFastaSequences[finalCandidates[k].first], queryFastaSequences[batchBegin + j]);
                resultsFile << finalCandidates[k].first << "-" << databaseFastaSequences[finalCandidates[k].first].size() << "-" << finalCandidates[k].second << ",";
            }
            if(finalCandidates.size() > 0){
                //auto sequencesSimilarity = getSequencesComparison(databaseFastaSequences[finalCandidates[finalCandidates.size() - 1].first], queryFastaSequences[batchBegin + j]);
                resultsFile << finalCandidates[finalCandidates.size() - 1].first << "-" << databaseFastaSequences[finalCandidates[finalCandidates.size() - 1].first].size() << "-" << finalCandidates[finalCandidates.size() - 1].second << std::endl;
            }
        }
    }
    auto stop = timer::now();
    cout << "Total query time duration in seconds: " << (duration_cast<chrono::seconds>(stop-start).count()) << endl;
}

std::vector<std::vector<pair<int32_t, int16_t>>>& processBatchByEdlib(std::vector<std::string>& referenceSequences, std::vector<std::string>& querySequences, uint64_t offset, uint64_t batchSize, uint64_t minPercentIdentity, bool parallel){
    std::vector<std::vector<pair<int32_t, int16_t>>> * batchResults = new std::vector<std::vector<pair<int32_t, int16_t>>>(batchSize, std::vector<pair<int32_t, int16_t>>());
    if(parallel){
#pragma omp parallel for
        for(uint64_t i = offset; i < offset + batchSize; i++){
            std::vector<pair<int32_t,int16_t>>* nearestNeighbours = new std::vector<pair<int32_t,int16_t>>();
            for(int32_t j = 0; j <(int32_t)referenceSequences.size(); j++){
                auto percentIdentity = fastPercentIdentity(querySequences[i], referenceSequences[j], minPercentIdentity);
                if(percentIdentity >= minPercentIdentity){
                    nearestNeighbours->push_back(make_pair((int32_t)j, (int16_t)percentIdentity));
                    //std::cout << j << " - " << ed_result.editDistance << std::endl;
                }
            }
            (*batchResults)[i - offset] = *nearestNeighbours;
        }
    }
    else {
        for(uint64_t i = offset; i < offset + batchSize; i++){
            std::vector<pair<int32_t,int16_t>>* nearestNeighbours = new std::vector<pair<int32_t,int16_t>>();
            for(int32_t j = 0; j <(int32_t)referenceSequences.size(); j++){
                auto percentIdentity = fastPercentIdentity(querySequences[i], referenceSequences[j], minPercentIdentity);
                if(percentIdentity >= minPercentIdentity){
                    nearestNeighbours->push_back(make_pair((int32_t)j, (int16_t)percentIdentity));
                    //std::cout << j << " - " << ed_result.editDistance << std::endl;
                }
            }
            (*batchResults)[i - offset] = *nearestNeighbours;
        }
    }

    return *batchResults;
}

void findSimilarSequencesByEdlib(string databaseFile, string databaseFileType, string queryFile, string queryFileType, uint64_t minPercentIdentity, bool parallel){
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
        std::vector<std::vector<pair<int32_t, int16_t>>>& res = processBatchByEdlib(databaseFastaSequences, queryFastaSequences, batchBegin, curretBatchSize, minPercentIdentity, parallel);
        cout << "Processed batch: " << batchBegin << " - " << batchEnd << endl;
        for(uint64_t j = 0; j < res.size(); j++){
            resultsFile << ">Sequence: " << batchBegin + j << "-" << queryFastaSequences[batchBegin + j].size() << std::endl;
            if(res[j].size() == 0 ){
                continue;
            }
            for(uint64_t k = 0; k < res[j].size() - 1; k++){
                //auto sequencesSimilarity = getSequencesComparison(databaseFastaSequences[res[j][k].first], queryFastaSequences[batchBegin + j]);
                resultsFile << res[j][k].first << "-" << databaseFastaSequences[res[j][k].first].size() << "-" << res[j][k].second << ",";
            }
            if(res[j].size() > 0){
                //auto sequencesSimilarity = getSequencesComparison(databaseFastaSequences[res[j][res[j].size() - 1].first], queryFastaSequences[batchBegin + j]);
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
                //int32_t sl;
                int16_t pi;
                //uint64_t ml, al;
                uint64_t pos = value.find('-'), prevPos;
                ind = stoi(value.substr(0, pos));
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                //sl = stoi(value.substr(prevPos, pos - prevPos));
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                pi = stoi(value.substr(prevPos, pos - prevPos));
                results[j-1].push_back(make_pair(ind, pi));
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
            if(results[i][j].second >= (int32_t)matches.size()){
                matches.resize(results[i][j].second + 1, 0);
                matchesFound.resize(results[i][j].second + 1, false);
            }
            if(!matchesFound[results[i][j].second]){
                matchesFound[results[i][j].second] = true;
                //matches[results[i][j].second]++;
            }
        }
        for(uint64_t k = 0; k < matchesFound.size(); k++) {
            if(matchesFound[k]) {
                matches[k]++;
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

    vector<int64_t> foundSequencesMinPI;
    vector<bool> foundSequences;
    for(uint64_t i = 0 ; i < results.size(); i++){
        for(uint64_t j = 0; j < results[i].size(); j++){
            if(results[i][j].first >= (int32_t)foundSequencesMinPI.size()){
                foundSequencesMinPI.resize(results[i][j].first + 1, -1);
                foundSequences.resize(results[i][j].first + 1, false);
            }
            if(!foundSequences[results[i][j].first]){
                foundSequences[results[i][j].first] = true;
            }
            if(foundSequencesMinPI[results[i][j].first] == -1 || (foundSequencesMinPI[results[i][j].first] > -1 && results[i][j].second > foundSequencesMinPI[results[i][j].first])) {
                foundSequencesMinPI[results[i][j].first] = results[i][j].second;
            }
        }
    }

    vector<int64_t> PIList;

    for(uint64_t i = 0; i < foundSequences.size(); i++){
        if(foundSequences[i]) {
            if(foundSequencesMinPI[i] >= (int16_t)PIList.size()){
                PIList.resize(foundSequencesMinPI[i] + 1, 0);
            }
            PIList[foundSequencesMinPI[i]]++;
        }
    }

    cout << "Printing categorized results: " << endl;
    for(uint64_t i = 0; i < categoryCounts.size(); i++){
        cout << "PI: [" << i*5 << " - " << (i + 1) * 5 << ") " <<  "  Count: " << categoryCounts[i] << endl;
    }

    cout << "Printing detailed results: " << endl;
    cout << "PI,count,cumcount" << endl;
    uint64_t totalCount = 0;
    for(uint64_t i = 0; i < matches.size(); i++){
        totalCount += matches[i];
        cout << i << "," << matches[i] << "," << totalCount << endl;
    }

    cout << "Printing summary of found Sequences: " << endl;
    totalCount = 0;
    for(uint64_t i = 0; i < PIList.size(); i++){
        totalCount += PIList[i];
        cout << "PI-" << i << " Count: " << PIList[i] << "  CumCount: "<< totalCount << endl;
    }

    cout << "Printing summary of found Sequences in csv format: " << endl;
    totalCount = 0;
    cout << "PI,Count,CumCount" << endl;
    for(uint64_t i = 0; i < PIList.size(); i++){
        totalCount += PIList[i];
        cout << i << "," << PIList[i] << ","<< totalCount << endl;
    }
}

void getClustersFromResults(string resultsFile) {
    cout << "Results file " << resultsFile << endl;
    ifstream file(resultsFile);
    ofstream type2ResultsFile(resultsFile + "_type2_results.txt"), type3ResultsFile(resultsFile + "_type3_results.txt");
    vector<vector<tuple<int32_t, uint32_t, int16_t>>> results;

    regex e("^>");
    smatch m;

    uint64_t j = 0;
    vector<uint64_t> queriesLengths;
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
                uint32_t sl;
                int16_t PI;
                uint64_t pos = value.find('-'), prevPos;
                ind = stoi(value.substr(0, pos));
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                sl = stoi(value.substr(prevPos, pos - prevPos));
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                PI = stoi(value.substr(prevPos, pos - prevPos));
                results[j-1].push_back(make_tuple(ind, sl, PI));
                //cout << ind << " " << sl << " " << PI << " " << ml << " " << al << endl;
            }
        } else {
            j++;
            queriesLengths.push_back(stoi(line.substr(line.find("-")+1)));
            results.push_back(vector<tuple<int32_t, uint32_t, int16_t>>());
        }
    }

    cout << "Queries size: " << queriesLengths.size() << endl;
    vector<bool> foundSequences(queriesLengths.size(), false);
    for(uint64_t i = 0 ; i < results.size(); i++){
        if(foundSequences[i]) {
            continue;
        }
        type2ResultsFile << ">" << "Sequence: " << i << "-" << queriesLengths[i] << endl;
        type3ResultsFile << ">" << "Sequence: " << i << "-" << queriesLengths[i] << endl;
        bool firstItemFoundType2 = false, firstItemFoundType3 = false;
        for(uint64_t j = 0; j < results[i].size(); j++){
            if(firstItemFoundType2) {
                type2ResultsFile << ",";
            }
            else {
                firstItemFoundType2 = true;
            }
            type2ResultsFile << get<0>(results[i][j]) << "-" << get<1>(results[i][j]) << "-" << get<2>(results[i][j]);
            if(!foundSequences[get<0>(results[i][j])]) {
                if(firstItemFoundType3) {
                    type3ResultsFile << ",";
                }
                else {
                    firstItemFoundType3 = true;
                }
                type3ResultsFile << get<0>(results[i][j]) << "-" << get<1>(results[i][j]) << "-" << get<2>(results[i][j]);
            }
            foundSequences[get<0>(results[i][j])] = true;
        }
        if(firstItemFoundType2) {
            type2ResultsFile << endl;
        }
        if(firstItemFoundType3) {
            type3ResultsFile << endl;
        }
        foundSequences[i] = true;
    }

    type2ResultsFile.close();
    type3ResultsFile.close();
}

void printClustersInformation(string resultsFile) {
    cout << "Results file " << resultsFile << endl;
    ifstream file(resultsFile);
    vector<vector<tuple<int32_t, uint32_t, int16_t>>> results;
    vector<uint64_t> uniqueSequences;
    regex e("^>");
    smatch m;

    uint64_t j = 0;
    vector<uint64_t> queriesLengths;
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
                uint32_t sl;
                int16_t PI;
                uint64_t pos = value.find('-'), prevPos;
                ind = stoi(value.substr(0, pos));
                if(find(uniqueSequences.begin(), uniqueSequences.end(), ind) == uniqueSequences.end() ) {
                    uniqueSequences.push_back((uint64_t)ind);
                }
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                sl = stoi(value.substr(prevPos, pos - prevPos));
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                PI = stoi(value.substr(prevPos, pos - prevPos));
                results[j-1].push_back(make_tuple(ind, sl, PI));
                //cout << ind << " " << sl << " " << PI << " " << ml << " " << al << endl;
            }
        } else {
            j++;
            queriesLengths.push_back(stoi(line.substr(line.find("-")+1)));
            results.push_back(vector<tuple<int32_t, uint32_t, int16_t>>());
        }
    }
    cout << "Unique indices size: " << uniqueSequences.size() << endl;
    uint64_t clusterId = 0, seqCount = 0;
    for(uint64_t i = 0; i < results.size(); i++){
        if(results[i].size() > 0) {
            cout << "ClusterId: " << clusterId << " Size: " << results[i].size() << endl;
            clusterId++;
            seqCount += results[i].size();
        }
    }
    cout << "Seq Count " << seqCount << endl;
}

void correctResultsFile(string databaseFile, string databaseFileType, string queryFile, string queryFileType, string resultsFile) {
    vector<string> databaseSequences, querySequences;
    loadFileByType(databaseFile, databaseFileType, databaseSequences);
    loadFileByType(queryFile, queryFileType, querySequences);

    ifstream file(resultsFile);
    vector<vector<pair<int32_t, int16_t>>> results;

    regex e("^>");
    smatch m;

    uint64_t j = 0;
    uint64_t query = 0;
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
                //int32_t sl;
                int16_t pi;
                //uint64_t ml, al;
                uint64_t pos = value.find('-'), prevPos;
                ind = stoi(value.substr(0, pos));
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                //sl = stoi(value.substr(prevPos, pos - prevPos));
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                pi = stoi(value.substr(prevPos, pos - prevPos));
                auto p = getSequencesComparison(databaseSequences[ind], querySequences[query]);
                pi = floor((p.first * 100.0)/(p.second * 1.0));
                results[j-1].push_back(make_pair(ind, pi));
            }
        } else {
            j++;
            uint64_t dashPos = line.find('-');
            query = stoi(line.substr(11, dashPos - 11));
            results.push_back(vector<pair<int32_t,int16_t>>());
        }
    }

    vector<int64_t> foundSequencesMinPI;
    vector<bool> foundSequences;
    for(uint64_t i = 0 ; i < results.size(); i++){
        for(uint64_t j = 0; j < results[i].size(); j++){
            if(results[i][j].first >= (int32_t)foundSequencesMinPI.size()){
                foundSequencesMinPI.resize(results[i][j].first + 1, -1);
                foundSequences.resize(results[i][j].first + 1, false);
            }
            if(!foundSequences[results[i][j].first]){
                foundSequences[results[i][j].first] = true;
            }
            if(foundSequencesMinPI[results[i][j].first] == -1 || (foundSequencesMinPI[results[i][j].first] > -1 && results[i][j].second > foundSequencesMinPI[results[i][j].first])) {
                foundSequencesMinPI[results[i][j].first] = results[i][j].second;
            }
        }
    }

    vector<int64_t> PIList;

    for(uint64_t i = 0; i < foundSequences.size(); i++){
        if(foundSequences[i]) {
            if(foundSequencesMinPI[i] >= (int16_t)PIList.size()){
                PIList.resize(foundSequencesMinPI[i] + 1, 0);
            }
            PIList[foundSequencesMinPI[i]]++;
        }
    }

    cout << "Printing summary of found Sequences: " << endl;
    uint64_t totalCount = 0;
    for(uint64_t i = 0; i < PIList.size(); i++){
        totalCount += PIList[i];
        cout << "PI-" << i << " Count: " << PIList[i] << "  CumCount: "<< totalCount << endl;
    }

}


vector<vector<tuple<int32_t, int16_t, int32_t>>> parseResultsFile(string resultsFile){
    fstream file(resultsFile);
    vector<vector<tuple<int32_t, int16_t, int32_t>>> results;

    regex e("^>");
    smatch m;

    uint64_t j = 0;
    int32_t query = 0;
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
                //int32_t sl;
                int16_t pi;
                //uint64_t ml, al;
                uint64_t pos = value.find('-'), prevPos;
                ind = stoi(value.substr(0, pos));
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                //sl = stoi(value.substr(prevPos, pos - prevPos));
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                pi = stoi(value.substr(prevPos, pos - prevPos));
                results[j-1].push_back(make_tuple(ind, pi, query));
            }
        } else {
            j++;
            uint64_t dashPos = line.find('-');
            query = stoi(line.substr(11, dashPos - 11));
            //cout << query  << "  "  << line << endl;
            results.push_back(vector<tuple<int32_t,int16_t, int32_t>>());
        }
    }
    return results;
}

void processConsensusResults(int argc, char** argv ) {
    int8_t numberOfResultsFiles = argc;
    vector<vector<vector<tuple<int32_t, int16_t, int32_t>>>> results;
    for(int8_t i = 0; i < numberOfResultsFiles; i++) {
        results.push_back(parseResultsFile(argv[i]));
    }

    vector<int64_t> foundSequencesMinPI;
    vector<bool> foundSequences;
    vector<string> foundSequencesQuery;
    for(uint64_t i = 0 ; i < results.size(); i++){
        for(uint64_t j = 0; j < results[i].size(); j++){
            for(uint64_t k = 0; k < results[i][j].size(); k++) {
                if(get<0>(results[i][j][k]) >= (int32_t)foundSequencesMinPI.size()){
                    foundSequencesMinPI.resize(get<0>(results[i][j][k]) + 1, -1);
                    foundSequences.resize(get<0>(results[i][j][k]) + 1, false);
                    foundSequencesQuery.resize(get<0>(results[i][j][k]) + 1, "");
                }
                if(!foundSequences[get<0>(results[i][j][k])]){
                    foundSequences[get<0>(results[i][j][k])] = true;
                }
                if(foundSequencesMinPI[get<0>(results[i][j][k])] == -1 || (foundSequencesMinPI[get<0>(results[i][j][k])] > -1 && get<1>(results[i][j][k]) > foundSequencesMinPI[get<0>(results[i][j][k])])) {
                    foundSequencesMinPI[get<0>(results[i][j][k])] = get<1>(results[i][j][k]);
                    foundSequencesQuery[get<0>(results[i][j][k])] = to_string(i) + "-"+ to_string(get<2>(results[i][j][k]));
                }
            }
        }
    }

    vector<int64_t> PIList;

    for(uint64_t i = 0; i < foundSequences.size(); i++){
        if(foundSequences[i]) {
            if(foundSequencesMinPI[i] >= (int16_t)PIList.size()){
                PIList.resize(foundSequencesMinPI[i] + 1, 0);
            }
            PIList[foundSequencesMinPI[i]]++;
        }
    }

    cout << "Printing summary of found Sequences: " << endl;
    uint64_t totalCount = 0;
    for(uint64_t i = 0; i < PIList.size(); i++){
        totalCount += PIList[i];
        cout << "PI-" << i << " Count: " << PIList[i] << "  CumCount: "<< totalCount << endl;
    }

    cout << "Printing summary of found Sequences in csv format: " << endl;
    totalCount = 0;
    cout << "PI,Count,CumCount" << endl;
    for(uint64_t i = 0; i < PIList.size(); i++){
        totalCount += PIList[i];
        cout << i << "," << PIList[i] << ","<< totalCount << endl;
    }

    cout << "Storing consensus results file" << endl;
    ofstream consensusResultsFile("consensusResultsFile.txt");
    for(uint64_t i = 0; i < foundSequences.size(); i++){
        if(foundSequences[i] && foundSequencesMinPI[i] > 70) {
            consensusResultsFile << ">Sequence: " << i  << "-" << foundSequencesMinPI[i] << endl;
            consensusResultsFile << foundSequencesQuery[i] << endl;
        }
    }
}

void generateConsensusFile(int argc, char** argv ){
    vector<vector<string>> querySets((argc - 1)/2, vector<string>());
    for(int8_t i = 0; i < argc - 1; i = i + 2) {
        loadFileByType(argv[i], argv[i+1], querySets[i/2]);
    }
    ofstream consensusQuerySequences("consensusQuerySequences.txt");
    fstream file(argv[argc-1]);

    regex e("^>");
    smatch m;

    uint64_t j = 0;

    while(!file.eof()){
        std::string line;
        std::getline(file, line);
        uint64_t pos;
        if((pos=line.find('\n')) != string::npos){
            line.erase(pos);
        }
        if((pos=line.find('\r')) != string::npos){
            line.erase(pos);
        }
        if(!regex_search(line, e)){
            stringstream ss(line);
            string value;
            while(getline(ss, value, ',')){
                int32_t ind;
                int32_t seqInd;
                //uint64_t ml, al;
                uint64_t pos = value.find('-'), prevPos;
                ind = stoi(value.substr(0, pos));
                prevPos = pos + 1;
                seqInd = stoi(value.substr(prevPos));
                consensusQuerySequences << querySets[ind][seqInd] << endl;
            }
        } else {
            j++;
            consensusQuerySequences << line << std::endl;
        }
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
                    int16_t PI;
                    uint8_t pos = value.find('-');
                    ind = stoi(value.substr(0, pos));
                    pos = value.find('-', pos + 1);
                    PI = stoi(value.substr(pos+1));
                    results[i][j-1].push_back(make_pair(ind, PI));
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
        cout << "PI: [" << i*5 << " - " << (i + 1) * 5 << ") " <<  "  FR: " << categoryCounts[0][i] << " ER: " << categoryCounts[1][i]  << " diff: " << (categoryCounts[1][i] - categoryCounts[0][i]) << endl;
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
                    int16_t PI;
                    uint8_t pos = value.find('-');
                    ind = stoi(value.substr(0, pos));
                    pos = value.find('-', pos + 1);
                    PI = stoi(value.substr(pos+1));
                    results[i][j-1].push_back(make_pair(ind, PI));
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

void performParameterTuningtests(string databaseFile, string databaseFileType, string queryFile, string queryFileType, uint8_t dataset_type, uint8_t data_type, uint64_t minPercentIdentity, uint64_t testType, double_t threshold = 0.0) {
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
            parametersTuner.process_queries_box_test(querySequences,threshold,minPercentIdentity);
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

void printFASTAInfor(string fastaFile){
    vector<string> sequences;
    loadFileByType(fastaFile, "fasta", sequences);
    int64_t min = -1, max = -1;
    for(uint64_t i = 0; i < sequences.size(); i++) {
        int64_t sequencesSize = sequences[i].size();
        if(min < 0) {
            min = sequencesSize;
        }
        if(max < 0) {
            max = sequencesSize;
        }

        if(min > sequencesSize) {
            min = sequencesSize;
        }
        if(max < sequencesSize) {
            max = sequencesSize;
        }
    }
    cout << "Minimum length of the sequences: " << min << endl;
    cout << "Maximum length of the sequences: " << max << endl;
}

void extractClusterRepresentatives(string clusterResultsFile, string sequencesFile, string sequencesFileType){
    vector<string> sequences;
    loadFileByType(sequencesFile, sequencesFileType, sequences);

    ofstream sequencesRepresentativesFile(clusterResultsFile+"_cluster_representatives.txt");
    fstream file(clusterResultsFile);

    regex e("^>");
    smatch m;

    uint64_t j = 0;

    while(!file.eof()){
        std::string line;
        std::getline(file, line);
        uint64_t pos;
        if((pos=line.find('\n')) != string::npos){
            line.erase(pos);
        }
        if((pos=line.find('\r')) != string::npos){
            line.erase(pos);
        }
        if(!regex_search(line, e)){
            int32_t ind;
            uint64_t pos;
            if(line.find('*') != string::npos){
                pos = line.find('*');
                pos--;
                ind = stoi(line.substr(0, pos));
                sequencesRepresentativesFile << sequences[ind] << endl;
            }
        } else {
            sequencesRepresentativesFile << line << endl;
            j++;
        }
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
            if(argc < 8) {
                cout << "Usage ./executable 0 [sequences_file] [file_type] [dataset_type] [data_type] [minPercentIdentity] stepSize" << endl;
                return 2;
            }
            performClustering(argv[2], argv[3], stoi(argv[4]), stoi(argv[5]), stoi(argv[6]), stoi(argv[7]));
            break;
        case 1:
            if(argc < 10) {
                cout << "Usage ./executable 1 [database_file] [file_type] [query_file] [file_type] [dataset_type] [data_type] [minPercentIdentity] [parallel]" << endl;
                return 3;
            }

            findSimilarSequences(argv[2], argv[3], argv[4], argv[5], stoi(argv[6]), stoi(argv[7]), stoi(argv[8]), stoi(argv[9]));
            break;
        case 2:
            if(argc < 6) {
                cout << "Usage ./executable 2 [database_file] [file_type] [query_file] [file_type] [minPercentIdentity] [parallel]" << endl;
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
                cout << "Usage ./executable 6 [database_file] [file_type] [query_file] [file_type] [dataset_type] [data_type] [minPercentIdentity] test_type threshold " << endl;
                return 7;
            }
            performParameterTuningtests(argv[2], argv[3], argv[4], argv[5], stoi(argv[6]), stoi(argv[7]), stoi(argv[8]), stoi(argv[9]), (argc == 11)?stoi(argv[10])/100.0:0);
            break;
        case 7:
            if(argc < 3) {
                cout << "Usage ./executable 7 results_file" << endl;
                return 8;
            }
            getClustersFromResults(argv[2]);
            break;
        case 8:
            if(argc < 3) {
                cout << "Usage ./executable 8 results_file" << endl;
                return 8;
            }
            printClustersInformation(argv[2]);
            break;
        case 9:
            if(argc < 3) {
                cout << "Usage ./executable 9 fasta_file" << endl;
                return 8;
            }
            printFASTAInfor(argv[2]);
            break;
        case 10:
            if(argc < 7) {
                cout << "Usage ./executable 10 database_file file_type query_file file_type results_file" << endl;
            }
            correctResultsFile(argv[2], argv[3], argv[4], argv[5], argv[6]);
            break;
        case 11:
            if(argc < 3) {
                cout << "Usage ./executable 11 results1 results2 results3" << endl;
            }
            processConsensusResults(argc - 2, &(argv[2]));
            break;
        case 12:
            if(argc < 4) {
                cout << "Usage ./executable 12 queryset1 fileType queryset2 fileType" << endl;
            }
            generateConsensusFile(argc -2, &(argv[2]));
            break;
        case 13:
            if(argc < 5) {
                cout << "Usage ./executable 13 clusterResultsFile sequencesFile fileType" << endl;
            }
            extractClusterRepresentatives(argv[2], argv[3], argv[4]);
            break;
        default:
            cout << "Invalid task. Ex task: 0 for clustering or 1 for similar fasta seq finder" << endl;
            return 100;
    }
    auto stop = timer::now();
    cout << "Total duration in seconds: " << (duration_cast<chrono::seconds>(stop-start).count()) << endl;
}
