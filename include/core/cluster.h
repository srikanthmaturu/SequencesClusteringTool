//
// Created by srikanthmaturu on 10/17/2017.
//

#pragma once

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <numeric>
#include <omp.h>
#include <iterator>
#include "fasta.hpp"
#include "falconn.h"
#include "tf_idf_falconn_idx.hpp"

namespace SequencesAnalyzer{
    namespace core{
        struct ClusterConfiguration{
            uint64_t percentIdentityThreshold;
            uint8_t typeOfAlignmentForPercentIdentityThreshold;
            std::string typeOfSimilarityAlgorithm;
            ClusterConfiguration(){

            }
            ClusterConfiguration(uint64_t percentIdentityThreshold, uint8_t typeOfAlignmentForPercentIdentityThreshold, std::string typeOfSimilarityAlgorithm):percentIdentityThreshold(percentIdentityThreshold),typeOfAlignmentForPercentIdentityThreshold(typeOfAlignmentForPercentIdentityThreshold),typeOfSimilarityAlgorithm(typeOfSimilarityAlgorithm){

            }
        };

        struct Cluster{
            uint64_t clusterId;
            int32_t clusterRepresentative;
            std::vector<std::pair<int32_t,double>> clusterItems;
            Cluster(uint64_t clusterId, int32_t clusterRepresentative) : clusterId(clusterId),clusterRepresentative(clusterRepresentative) {

            }

            uint64_t clusterSize() {
                return clusterItems.size();
            }

            void printClusterInfo(ofstream& outputFileStream) {
                outputFileStream << clusterRepresentative << " *" << std::endl;
                for(uint64_t i = 0; i < clusterItems.size(); i++) {
                    outputFileStream << clusterItems[i].first << " at " << clusterItems[i].second << "%" << std::endl;
                }
            }
        };

        class ClustersGenerator{
        public:
            ClustersGenerator() = default;

            ClustersGenerator(const ClustersGenerator &) = default;

            ClustersGenerator(ClustersGenerator &&) = default;

            ClustersGenerator &operator=(const ClustersGenerator &) = default;

            ClustersGenerator &operator=(ClustersGenerator &&) = default;

            ClustersGenerator(std::vector<std::string>& sequences, FALCONNIndexConfiguration config, ClusterConfiguration clusterConfig){
                this->sequences = sequences;
                originalSequences.insert(originalSequences.end(), sequences.begin(), sequences.end());
                sequencesLocations.resize(sequences.size());
                for(int32_t i = 0; i < (int32_t)sequences.size(); i++) {
                    sequencesLocations[i] = i;
                }
                std::cout << "Sorting sequences.... " << std::endl;
                sortSequences(this->sequences);
                std::cout << "Sorting complete. " << std::endl;
                FALCONNConfig = config;
                this->clusterConfig = clusterConfig;
                clusteredSequences = std::vector<bool>(sequences.size(), false);
            }

        public:
            FALCONNIndexConfiguration FALCONNConfig;
            ClusterConfiguration clusterConfig;

            std::vector<std::string> sequences;
            std::vector<std::string> originalSequences;
            std::vector<int32_t> sequencesLocations;
            std::vector<bool> clusteredSequences;
            std::vector<Cluster> clusters;


            tf_idf_falconn_index::tf_idf_falconn_idx<5, false, false, false, 2, 32, 11, 2032, 150, 0, 0, POINT_TYPE> idx;

            void sortSequences(std::vector<std::string>& sequences){
                std::sort(sequencesLocations.begin(), sequencesLocations.end(), [&](int32_t ai, int32_t bi){ return sequences[ai].size() > sequences[bi].size(); } );
                std::sort(sequences.begin(), sequences.end(), [](std::string a, std::string b){ return a.size() > b.size(); });
            }

            void initialize(){
                std::cout << "Initializing.... " << std::endl;
                falconn::LSHConstructionParameters lshParams;

                switch(FALCONNConfig.lshType){
                    case 1:
                        lshParams.lsh_family = falconn::LSHFamily::Hyperplane;
                        break;
                    case 2:
                        lshParams.lsh_family = falconn::LSHFamily::CrossPolytope;
                        break;
                }
                lshParams.l = FALCONNConfig.numberOfHashTables;
                lshParams.num_rotations = 1;
                switch(FALCONNConfig.data_type) {
                    case 0:
                        lshParams.feature_hashing_dimension = pow(4, FALCONNConfig.ngl);
                        break;
                    case 1:
                        lshParams.feature_hashing_dimension = 1024;
                        break;
                    default:
                        std::cerr<< "Invalid datatype. Exiting..." << std::endl;
                        exit(1);
                }
                idx.setlshParams(lshParams);
                idx.setThreshold(FALCONNConfig.threshold);
                idx.setNGL(FALCONNConfig.ngl);
                idx.setNumberOfProbes(FALCONNConfig.numberOfProbes);
                idx.setDatasetType(FALCONNConfig.dataset_type);
                idx.setDataType(FALCONNConfig.data_type);
                idx.initialize(sequences);
                std::cout << "Constructing initial table: " << std::endl;
                idx.construct_table();
            }

            void reinitialize(){
                idx.reconstruct_table_by_removing_data(clusteredSequences);
                int32_t i = 0;
                while(i < (int32_t)clusteredSequences.size()) {
                    if(clusteredSequences[i]) {
                        clusteredSequences.erase(clusteredSequences.begin() + i);
                        sequencesLocations.erase(sequencesLocations.begin() + i);
                        sequences.erase(sequences.begin() + i);
                    }
                    else {
                        i++;
                    }
                }
            }

            std::pair<uint64_t, std::vector<int32_t>*> match(std::string sequence){
                return idx.getNearestNeighbours(sequence);
            }

            std::pair<uint64_t, std::vector<int32_t>*> bruteForceMatch(std::string sequence) {
                std::pair<uint64_t, std::vector<int32_t>*> matchesPair;
                matchesPair.first = sequences.size();
                std::vector<int32_t> * matches = new std::vector<int32_t>(sequences.size(),0);
                std::iota(matches->begin(),matches->end(),0);
                matchesPair.second = matches;
                return matchesPair;
            }

            bool cluster(int32_t sequenceIndex, uint8_t clusteringType) {
                if(clusteredSequences[sequenceIndex]) {
                    return false;
                }
                std::pair<uint64_t, std::vector<int32_t>*> matchesPair;

                if(clusterConfig.typeOfSimilarityAlgorithm == "falconn") {
                    matchesPair = match(sequences[sequenceIndex]);
                }
                else if(clusterConfig.typeOfSimilarityAlgorithm == "brute-force") {
                    matchesPair = bruteForceMatch(sequences[sequenceIndex]);
                }
                else {
                    cout << "Uknown algorithm "<< clusterConfig.typeOfSimilarityAlgorithm << ". Exiting.." << endl;
                    exit(1);
                }

                std::cout << "Matches found: " <<  matchesPair.second->size() << "\t";
                bool candidatesSelection[matchesPair.second->size()] = {false};
                std::vector<int32_t>* unclusteredMatches = new std::vector<int32_t>();
#pragma omp parallel for
                for(uint64_t i = 0; i < matchesPair.second->size(); i++) {
                    int32_t matchIndex = (*(matchesPair.second))[i];
                    if(!clusteredSequences[matchIndex]){
                        double pi = getFastPI(sequences[matchIndex], sequences[sequenceIndex]);
                        if( pi >= clusterConfig.percentIdentityThreshold) {
                            candidatesSelection[i] = true;
                        }
                    }
                }
                for(uint64_t i = 0; i < matchesPair.second->size(); i++) {
                    int32_t matchIndex = (*(matchesPair.second))[i];
                    if(candidatesSelection[i]) {
                        unclusteredMatches->push_back(matchIndex);
                    }
                }
                std::cout << "Clusterable matches: " << unclusteredMatches->size() << std::endl;
                std::sort(unclusteredMatches->begin(), unclusteredMatches->end(), std::less<int32_t>());
                Cluster* cluster = new Cluster(clusters.size(), sequencesLocations[sequenceIndex]);
                clusteredSequences[sequenceIndex] = true;
                for(uint64_t i = 0; i < unclusteredMatches->size(); i++) {
                    if(clusteredSequences[(*unclusteredMatches)[i]]) {
                        continue;
                    }
                    double pi = getFastPI(sequences[sequenceIndex], sequences[(*unclusteredMatches)[i]]);
                    if( pi < clusterConfig.percentIdentityThreshold) {
                        continue;
                    }
                    bool clusterable = true;
                    if(clusteringType == 0) {
                        for(uint64_t j = 0; j < cluster->clusterItems.size(); j++) {
                            if(getFastPI(originalSequences[cluster->clusterItems[j].first], sequences[(*unclusteredMatches)[i]]) < clusterConfig.percentIdentityThreshold) {
                                clusterable = false;
                                break;
                            }
                        }
                    }

                    if(clusterable) {
                        cluster->clusterItems.push_back(make_pair(sequencesLocations[(*unclusteredMatches)[i]], pi));
                        clusteredSequences[(*unclusteredMatches)[i]] = true;
                    }
                }
                clusters.push_back(*cluster);
                std::cout << "Cluster size " << cluster->clusterItems.size() << std::endl;
                return true;
            }

            double getFastPI(string firstSequence, string lastSequence){
                switch(clusterConfig.typeOfAlignmentForPercentIdentityThreshold){
                    case 0:
                        return fastPercentIdentity(firstSequence, lastSequence, clusterConfig.percentIdentityThreshold);
                        break;
                    case 1:
                        return getInfixPercentIdentity(firstSequence, lastSequence);
                        break;
                    default:
                        std::cerr << "Invalid alignment type specified. Valid options are 0 and 1" << std::endl;
                        exit(1);
                }
            }

            void generateClusters(uint64_t stepSize){
                int32_t j = 0;
                int32_t sequencesSize = sequences.size();
                for(int32_t i = 0; i < sequencesSize && j < (int32_t)sequences.size(); i++,j++) {
                    if( i > 0 && i % stepSize == 0) {
                        reinitialize();
                        std::cout << "Reinitialized index " << std::endl;
                        j = 0;
                    }
                    std::cout << "Clustering " << j << std::endl;
                    cluster(j, 1);
                    if( i > 0 && i%stepSize == 0) {
                        std::cout << "Processed sequences: " << i - stepSize << " - " << i << std::endl;
                    }
                }
            }

            void printClusters(ofstream& clustersFile) {
                for(uint64_t i = 0; i < clusters.size(); i++) {
                    clustersFile << ">ClusterId: " << i << endl;
                    clusters[i].printClusterInfo(clustersFile);
                }
            }
        };
    }
}
