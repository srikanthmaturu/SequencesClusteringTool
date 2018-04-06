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

void loadFileByType(string file, string fileType, vector<string>& sequences, vector<string>& descriptionLines){
    cout << "Loading " << file << " of type " << fileType << endl;
    if(fileType == "kmers"){
        load_sequences(file, sequences);
    }
    else if(fileType == "fasta"){
        fasta::getFastaSequences(file, sequences, descriptionLines);
    }
    else{
        cout << "Wrong file type." << endl;
        exit(1);
    }
}

void performClustering(string sequencesFile, string sequencesFileType, uint64_t dataset_type, uint64_t data_type, uint64_t percentIdentityThreshold, uint8_t typeOfAlignment, uint64_t stepSize){
    ofstream resultsFile(sequencesFile + ".clusters.txt", ofstream::out);
    vector<string> sequences;
    vector<string> descriptionLines;
    loadFileByType(sequencesFile, sequencesFileType, sequences, descriptionLines);
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
    clusterConfig.typeOfAlignmentForPercentIdentityThreshold = typeOfAlignment;
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
    vector<string> databaseDescriptionLines, queryDescriptionLines;
    loadFileByType(databaseFile, databaseFileType, databaseFastaSequences, databaseDescriptionLines);
    loadFileByType(queryFile, queryFileType, queryFastaSequences, queryDescriptionLines);

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
    uint64_t totalNumberOfCandidates = 0;
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
            totalNumberOfCandidates += result_matches.size();
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
    cout << "Total number of candidates: " << totalNumberOfCandidates << endl;
    cout << "Average number of candidates: " << (totalNumberOfCandidates * 1.0)/queryFastaSequences.size() << endl;
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
    vector<string> databaseDescriptionLines, queryDescriptionLines;
    loadFileByType(databaseFile, databaseFileType, databaseFastaSequences, databaseDescriptionLines);
    loadFileByType(queryFile, queryFileType, queryFastaSequences, queryDescriptionLines);

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
    vector<uint64_t> bestMatches(results.size(), 0);
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
                if(bestMatches[i] < k) {
                    bestMatches[i] = k;
                }
            }
        }
        std::fill(matchesFound.begin(), matchesFound.end(), false);
    }

    vector<uint64_t> bestMatchesDistribution;
    for(uint64_t i = 0; i < bestMatches.size(); i++) {
        if(bestMatches[i] >= bestMatchesDistribution.size()){
            bestMatchesDistribution.resize(bestMatches[i] + 1, 0);
        }
        bestMatchesDistribution[bestMatches[i]]++;
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

    cout << "Printing best matches distribution in csv format: " << endl;
    totalCount = 0;
    cout << "PI,Count,CumCount" << endl;
    bestMatchesDistribution[0] = 0;
    for(uint64_t i = 0; i < bestMatchesDistribution.size(); i++){
        totalCount += bestMatchesDistribution[i];
        cout << i << "," << bestMatchesDistribution[i] << ","<< totalCount << endl;
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
    vector<string> databaseDescriptionLines, queryDescriptionLines;
    loadFileByType(databaseFile, databaseFileType, databaseSequences, databaseDescriptionLines);
    loadFileByType(queryFile, queryFileType, querySequences, queryDescriptionLines);

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
                //cout << value << endl;
                uint64_t pos = value.find('-'), prevPos;
                ind = stoi(value.substr(0, pos));
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                //sl = stoi(value.substr(prevPos, pos - prevPos));
                prevPos = pos + 1;
                pos = value.find('-', prevPos);
                pi = stoi(value.substr(prevPos, pos - prevPos));
                //cout << ind << " - " << pi << endl;
                //cout << databaseSequences[ind] << endl;
                //cout << querySequences[query] << endl;
                auto p = getSequencesComparison(databaseSequences[ind], querySequences[query]);
                //cout << "Done sequence comparision " << endl;
                pi = floor((p.first * 100.0)/(p.second * 1.0));
                //cout << "Done calculating pi " << endl;
                results[j-1].push_back(make_pair(ind, pi));
                //cout << "Done one loop " << endl;
            }
            //cout << "While exited!" << endl;
        } else {
            j++;
            uint64_t dashPos = line.find('-');
            //cout << line << "  " << line.substr(11, dashPos - 11) << endl;
            query = stoi(line.substr(11, dashPos - 11));
            //cout << query << " - " << endl;
            results.push_back(vector<pair<int32_t,int16_t>>());
        }
    }
    cout << "Read results file." << endl;
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
    vector<vector<string>> descriptionLineSets((argc - 1)/2, vector<string>());
    for(int8_t i = 0; i < argc - 1; i = i + 2) {
        loadFileByType(argv[i], argv[i+1], querySets[i/2], descriptionLineSets[i/2]);
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

uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed) {
    uint32_t h = seed;
    if (len > 3) {
        const uint32_t* key_x4 = (const uint32_t*) key;
        size_t i = len >> 2;
        do {
            uint32_t k = *key_x4++;
            k *= 0xcc9e2d51;
            k = (k << 15) | (k >> 17);
            k *= 0x1b873593;
            h ^= k;
            h = (h << 13) | (h >> 19);
            h = (h * 5) + 0xe6546b64;
        } while (--i);
        key = (const uint8_t*) key_x4;
    }
    if (len & 3) {
        size_t i = len & 3;
        uint32_t k = 0;
        key = &key[i - 1];
        do {
            k <<= 8;
            k |= *key--;
        } while (--i);
        k *= 0xcc9e2d51;
        k = (k << 15) | (k >> 17);
        k *= 0x1b873593;
        h ^= k;
    }
    h ^= len;
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

void simulateAllCombinedResults(uint8_t argc, char** argv) {
    uint8_t numberOfResultsFiles = argc/3;
    vector<vector<string>> queryFiles(numberOfResultsFiles, vector<string>());
    vector<vector<string>> descriptionLineSets(numberOfResultsFiles, vector<string>());
    vector<vector<vector<tuple<int32_t, int16_t, int32_t>>>> results;

    for(uint8_t i = 0; i < numberOfResultsFiles * 2; i = i + 2) {
        loadFileByType(argv[i], argv[i+1], queryFiles[i/2], descriptionLineSets[i/2]);
    }
    cout << "Read all query files." << endl;
    for(uint8_t i = numberOfResultsFiles*2; i < argc; i++) {
        results.push_back(parseResultsFile(argv[i]));
    }
    cout << "Read all results file." << endl;

    map<uint64_t, int32_t> sequencesBestMatches;
    for(uint64_t i = 0 ; i < results.size(); i++){
        cout << "Processing results " << i << endl;
        for(uint64_t j = 0; j < results[i].size(); j++){
            uint64_t hash = murmur3_32((uint8_t*)queryFiles[i][j].c_str(), queryFiles[i][j].size(), 0);
            //0xcc9e2d51
            if(sequencesBestMatches.find(hash) == sequencesBestMatches.end()) {
                sequencesBestMatches[hash] = 0;
            }
            for(uint64_t k = 0; k < results[i][j].size(); k++) {
                if(sequencesBestMatches[hash] < get<1>(results[i][j][k])) {
                    sequencesBestMatches[hash] = get<1>(results[i][j][k]);
                }
            }
        }
    }
    cout << "Processed all results." << endl;
    vector<uint64_t> bestMatchesDistribution;
    for(auto it:sequencesBestMatches) {
        if(it.second >= (int32_t)bestMatchesDistribution.size()){
            bestMatchesDistribution.resize(it.second + 1, 0);
        }
        bestMatchesDistribution[it.second]++;
    }
    cout << "Total number of all combined sequences: " << sequencesBestMatches.size() << endl;
    cout << "Printing best matches distribution in csv format: " << endl;
    uint64_t totalCount = 0;
    cout << "PI,Count,CumCount" << endl;
    bestMatchesDistribution[0] = 0;
    for(uint64_t i = 0; i < bestMatchesDistribution.size(); i++){
        totalCount += bestMatchesDistribution[i];
        cout << i << "," << bestMatchesDistribution[i] << ","<< totalCount << endl;
    }
}

void compareFalconnAndEdlibResults(uint8_t argc, char** argv){
    vector<vector<vector<tuple<int32_t, int16_t, int32_t>>>> results;

    for(uint8_t i = 0; i < argc; i++) {
        results.push_back(parseResultsFile(argv[i]));
    }
    cout << "Read all results file." << endl;

    vector<vector<uint64_t>> matches(2,vector<uint64_t>());

    for(uint64_t i = 0; i < 2; i++) {
        for(uint64_t j = 0 ; j < results[i].size(); j++){
            for(uint64_t k = 0; k < results[i][j].size(); k++){
                if(get<1>(results[i][j][k]) >= (int32_t)matches[i].size()){
                    matches[i].resize(get<1>(results[i][j][k]) + 1, 0);
                }
                matches[i][get<1>(results[i][j][k])]++;
            }
        }
    }

    if(matches[0].size() > matches[1].size()){
        matches[1].resize(matches[0].size(), 0);
    }
    if(matches[1].size() > matches[0].size()){
        matches[0].resize(matches[1].size(), 0);
    }
    cout << "PI,matches" << endl;
    for(uint64_t i = 0; i < matches[0].size(); i++){
        cout << i << "," << matches[0][i] << "," << matches[1][i] << endl;
    }
}

void generateSampleFiles(string falconnResults, string edlibResults, string databaseFile, string databaseFileType, string queryFile, string queryFileType, uint64_t numberOfSequences){
    ifstream files[2];
    files[0] = ifstream(falconnResults);
    files[1] = ifstream(edlibResults);
    vector<vector<pair<int32_t, int16_t>>> results[2];
    vector<string> databaseSequences, querySequences;
    vector<string> databaseDescriptionLines, queryDescriptionLines;
    loadFileByType(databaseFile, databaseFileType, databaseSequences, databaseDescriptionLines);
    loadFileByType(queryFile, queryFileType, querySequences, queryDescriptionLines);

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
    vector<string> databaseDescriptionLines, queryDescriptionLines;
    loadFileByType(databaseFile, databaseFileType, databaseSequences, databaseDescriptionLines);
    loadFileByType(queryFile, queryFileType, querySequences, queryDescriptionLines);

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
    vector<string> descriptionLines;
    loadFileByType(fastaFile, "fasta", sequences, descriptionLines);
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
    vector<string> descriptionLines;
    loadFileByType(sequencesFile, sequencesFileType, sequences, descriptionLines);
    vector<bool> representativeSequencesSelection(sequences.size(), false);
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
                representativeSequencesSelection[ind] = true;
            }
        } else {
            j++;
        }
    }

    for(uint64_t i = 0; i < sequences.size(); i++) {
        if(representativeSequencesSelection[i]) {
            sequencesRepresentativesFile << descriptionLines[i] << endl;
            j = 0;
            for(j = 0; j < (sequences[i].size()/80); j++) {
                sequencesRepresentativesFile << sequences[i].substr(j * 80, 80) << endl;
            }
            sequencesRepresentativesFile << sequences[i].substr(j * 80) << endl;
        }
    }

}

vector<vector<int32_t>> extractFALCONNClusteringResults(string falconnClusteringResults) {
    fstream file(falconnClusteringResults, std::ifstream::in);
    vector<vector<int32_t>> results;
    regex e("^>");
    smatch m;
    uint64_t numberOfClusters = 0;

    while(!file.eof()) {
        std::string line;
        std::getline(file, line);
        if(file.eof()) {
            break;
        }
        int32_t sequenceId;
        uint64_t pos;
        if ((pos = line.find('\n')) != string::npos) {
            line.erase(pos);
        }
        if ((pos = line.find('\r')) != string::npos) {
            line.erase(pos);
        }
        if(!regex_search(line, e)) {
            pos = line.find(' ');
            //cout << line << endl << " -";
            sequenceId = stoi(line.substr(0, pos));
            results[numberOfClusters-1].push_back(sequenceId);
        }
        else{
            numberOfClusters++;
            results.push_back(vector<int32_t>());
        }
    }
    return results;
}

int32_t findFastaSequenceIndexByDescriptionLine(vector<string> descriptionLines, string descriptionLineSegment){
    for(uint64_t i = 0; i < descriptionLines.size(); i++) {
        if(descriptionLines[i].find(descriptionLineSegment) != string::npos) {
            return i;
        }
    }
    return -1;
}

map<string, int32_t> getDescriptionLinesKeyIndexMap(vector<string>& descriptionLines){
    map<string, int32_t> keyIndexMap;
    for(uint64_t i = 0; i < descriptionLines.size(); i++) {
        int32_t pos = descriptionLines[i].find(" ");
        string key = descriptionLines[i].substr(0, pos);
        keyIndexMap[key] = i;
    }
    return keyIndexMap;
}

vector<vector<int32_t>> extractCDHITClusteringResults(string cdhitClusteringResults, string sequencesFile) {
    fstream file(cdhitClusteringResults);
    vector<string> sequences;
    vector<string> descriptionLines;
    loadFileByType(sequencesFile, "fasta", sequences, descriptionLines);
    vector<vector<int32_t>> results;
    map<string, int32_t> descriptionLinesKeyIndexMap = getDescriptionLinesKeyIndexMap(descriptionLines);

    regex e("^>"), e1("([0-9]+)\t(.+), (>.+)... (\\*|at)(.*)");
    smatch m;

    uint64_t j = 0;
    int32_t sequenceIndex = 0;
    while(!file.eof()){
        std::string line;
        std::getline(file, line);
        if(file.eof()) {
            break;
        }
        uint64_t pos;
        if((pos=line.find('\n')) != string::npos){
            line.erase(pos);
        }
        if((pos=line.find('\r')) != string::npos){
            line.erase(pos);
        }
        //cout << line << endl;
        if(!regex_search(line, e)){
            regex_match(line, m, e1);
            //cout << descriptionLinesKeyIndexMap[m[3]] << " " << m[3] << endl;
            sequenceIndex = descriptionLinesKeyIndexMap[m[3]];
            //cout << "Matching results size: " << m.size() << endl;
            //cout << "m0: " << m[0] << " m1: " << m[1] << " m2: " << m[2] << " m3: " << m[3] << " m4: " << m[4] << " m5: " << m[5] << endl;
            if(m[4] == "*") {
                results[j-1].insert(results[j-1].begin(), sequenceIndex);
            }
            else if(m[4] == "at") {
                results[j-1].push_back(sequenceIndex);
            }
            else {
                cout << "Invalid string. st: " << line << endl;
                exit(1);
            }
        } else {
            j++;
            results.push_back(vector<int32_t>());
        }
    }
    //cout << endl;
    return results;
}

vector<vector<int32_t>> extractClusteringResults(string resultsFile, string algorithm, string sequencesFile){
    if(algorithm == "falconn") {
        cout << "Reading falconn clustering results file." << endl;
        return extractFALCONNClusteringResults(resultsFile);
    }
    else if(algorithm == "cd-hit") {
        cout << "Reading cd-hit clustering results file." << endl;
        return extractCDHITClusteringResults(resultsFile, sequencesFile);
    }
    else {
        cout << "Uknown algorithm "<< algorithm << ". Exiting.." << endl;
        exit(1);
    }
}

void compareFALCONNAndCDHITClusteringResults(int8_t argc, char** argv) {
    int numberOfAlgorithms = (argc - 1) / 2;
    cout << "Number of algorithms: " << numberOfAlgorithms << endl;
    vector<vector<vector<int32_t>>> algorithmsClusteringResults(numberOfAlgorithms);
    for(int8_t i = 0; i < numberOfAlgorithms; i++) {
        algorithmsClusteringResults[i] = extractClusteringResults(argv[i*2], argv[i*2+1], argv[argc - 1]);
    }
    cout << "Algorithms results are loaded." << endl;
    vector<string> sequences, descriptionLines;
    loadFileByType(argv[argc - 1], "fasta", sequences, descriptionLines);
    vector<vector<uint64_t>> algorithmsPICounts(numberOfAlgorithms, vector<uint64_t>(101, 0));
    vector<uint64_t> singleClustersSize(numberOfAlgorithms, 0);
    for(int8_t i = 0; i < numberOfAlgorithms; i++) {
        cout << "Algorithm: " << i+1 << " Number Of Clusters: " << algorithmsClusteringResults[i].size() << endl;
        for(uint64_t j = 0; j < algorithmsClusteringResults[i].size(); j++) {
            uint64_t clusterSize = algorithmsClusteringResults[i][j].size();
            if(clusterSize == 1) {
                singleClustersSize[i]++;
            }
            for(uint64_t k = 0; k < clusterSize - 1; k++) {
                for(uint64_t l = k + 1; l < clusterSize; l++) {
                    int8_t percentIdentity = getInfixPercentIdentity(sequences[algorithmsClusteringResults[i][j][k]], sequences[algorithmsClusteringResults[i][j][l]]);
                    if(percentIdentity < 10) {
                        cout << "Algorithm: " << (int)i << endl;
                        cout << sequences[algorithmsClusteringResults[i][j][k]] << endl;
                        cout << sequences[algorithmsClusteringResults[i][j][l]] << endl;
                        cout << "PI: " << (int)percentIdentity << endl;
                    }
                    algorithmsPICounts[i][percentIdentity]++;
                }
            }
        }
    }
    cout << "Number of single clusters: " << endl;
    for(int8_t i = 0; i < numberOfAlgorithms; i++) {
        cout << argv[i*2 + 1] << " : " << singleClustersSize[i] << endl;
    }
    cout << "PI";
    for(int8_t i = 0; i < numberOfAlgorithms; i++) {
        cout << "," << argv[i*2 + 1];
    }
    cout << endl;
    for(int8_t j = 0; j < 101; j++) {
        cout << (int)j << ",";
        for(int8_t i = 0; i < numberOfAlgorithms; i++) {
            cout << algorithmsPICounts[i][j];
            if(i < numberOfAlgorithms - 1) {
                cout << "," ;
            }
        }
        cout << endl;
    }
}

void generateFastaFile(string sequencesFile, string fileType, uint64_t numberOfsequences) {
    vector<string> sequences;
    vector<string> descriptionLines;
    loadFileByType(sequencesFile, fileType, sequences, descriptionLines);
    ofstream outputSequencesFile(sequencesFile + "_" + to_string(numberOfsequences) + "_sequences");
    uint64_t sequencesSize = sequences.size();

    vector<bool> selectedSequences(sequencesSize, false);
    srand (time(NULL));
    for(uint64_t i = 0; i < numberOfsequences;) {
        uint64_t sequenceId = rand() % sequencesSize;
        if(!selectedSequences[sequenceId]) {
            selectedSequences[sequenceId] = true;
            outputSequencesFile << ">Sequence" << to_string(i) << endl;
            outputSequencesFile << sequences[i] << endl;
            i++;
        }
    }

    outputSequencesFile.close();
}



int main(int argc, char** argv){

    if(argc < 2){
        cout << "Usage ./executable [type_of_task] [additional_options]" << endl;
        return 1;
    }
    auto start = timer::now();
    switch(stoi(argv[1])){
        case 0:
            if(argc < 9) {
                cout << "Usage ./executable 0 [sequences_file] [file_type] [dataset_type] [data_type] [minPercentIdentity] typeOfAlignment stepSize" << endl;
                return 2;
            }
            performClustering(argv[2], argv[3], stoi(argv[4]), stoi(argv[5]), stoi(argv[6]), stoi(argv[7]),  stoi(argv[8]));
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
            compareFalconnAndEdlibResults(argc -2, &(argv[2]));
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
                return 9;
            }
            printFASTAInfor(argv[2]);
            break;
        case 10:
            if(argc < 7) {
                cout << "Usage ./executable 10 database_file file_type query_file file_type results_file" << endl;
                return 10;
            }
            correctResultsFile(argv[2], argv[3], argv[4], argv[5], argv[6]);
            break;
        case 11:
            if(argc < 3) {
                cout << "Usage ./executable 11 results1 results2 results3" << endl;
                return 11;
            }
            processConsensusResults(argc - 2, &(argv[2]));
            break;
        case 12:
            if(argc < 4) {
                cout << "Usage ./executable 12 queryset1 fileType queryset2 fileType" << endl;
                return 12;
            }
            generateConsensusFile(argc -2, &(argv[2]));
            break;
        case 13:
            if(argc < 5) {
                cout << "Usage ./executable 13 clusterResultsFile sequencesFile fileType" << endl;
                return 13;
            }
            extractClusterRepresentatives(argv[2], argv[3], argv[4]);
            break;
        case 14:
            if(argc < 5) {
                cout << "Usage ./executable 14 sequencesFile1 fileType .. resultsFile1" << endl;
                return 14;
            }
            simulateAllCombinedResults(argc - 2, &argv[2]);
            break;
        case 15:
            if(argc < 4) {
                cout << "Usage ./executable 15 sequencesFile1 numberOfSequences" << endl;
                return 15;
            }
            generateFastaFile(argv[2], argv[3], stoi(argv[4]));
            break;
        case 16:
            if(argc < 5) {
                cout << "Usage ./executable 16 clusteringResultsFile clusteringAlgorithm sequencesFile" << endl;
                return 16;
            }
            compareFALCONNAndCDHITClusteringResults(argc - 2, &argv[2]);
            break;
        default:
            cout << "Invalid task. Ex task: 0 for clustering or 1 for similar fasta seq finder" << endl;
            return 100;
    }
    auto stop = timer::now();
    cout << "Total duration in seconds: " << (duration_cast<chrono::seconds>(stop-start).count()) << endl;
}
