//
// Created by srikanthmaturu on 10/17/2017.
//

#pragma once

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <omp.h>
#include <iterator>
#include "fasta.hpp"
#include "falconn.h"
#include "tf_idf_falconn_idx.hpp"

namespace SequencesAnalyzer{
    namespace core{
        struct ClusterConfiguration{
            uint64_t percentIdentityThreshold;
            ClusterConfiguration(){

            }
            ClusterConfiguration(uint64_t percentIdentityThreshold):percentIdentityThreshold(percentIdentityThreshold){

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

            void printClusterInfo(ofstream& outputFileStream, std::vector<int32_t> sequencesLocations) {
                outputFileStream << sequencesLocations[clusterRepresentative] << " *" << std::endl;
                for(uint64_t i = 0; i < clusterItems.size(); i++) {
                    outputFileStream << sequencesLocations[clusterItems[i].first] << " at " << clusterItems[i].second << "%" << std::endl;
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
                offset = 0;
            }

        public:
            FALCONNIndexConfiguration FALCONNConfig;
            ClusterConfiguration clusterConfig;

            std::vector<std::string> sequences;
            std::vector<int32_t> sequencesLocations;
            std::vector<bool> clusteredSequences;
            std::vector<Cluster> clusters;

            int32_t offset = 0;

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
                lshParams.feature_hashing_dimension = pow(4, FALCONNConfig.ngl);
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

            void reinitialize(int32_t dataOffset){
                idx.reconstruct_table_by_offset_data(dataOffset);
            }

            std::pair<uint64_t, std::vector<int32_t>*> match(std::string sequence){
                return idx.getNearestNeighbours(sequence);
            }

            bool cluster(int32_t sequenceIndex) {
                if(clusteredSequences[sequenceIndex]) {
                    return false;
                }
                std::cout << "Clustering " << sequenceIndex << std::endl;
                auto matchesPair = match(sequences[sequenceIndex]);
                std::cout << "Matches found: " <<  matchesPair.second->size() << "\t";
                bool candidatesSelection[matchesPair.second->size()] = {false};
                std::vector<int32_t>* unclusteredMatches = new std::vector<int32_t>();
#pragma omp parallel for
                for(uint64_t i = 0; i < matchesPair.second->size(); i++) {
                    int32_t matchIndex = offset + (*(matchesPair.second))[i];
                    if(!clusteredSequences[matchIndex]){
                        double pi = getFastPI(matchIndex, sequenceIndex);
                        if( pi >= clusterConfig.percentIdentityThreshold) {
                            candidatesSelection[i] = true;
                        }
                    }
                }
                for(uint64_t i = 0; i < matchesPair.second->size(); i++) {
                    int32_t matchIndex = offset + (*(matchesPair.second))[i];
                    if(candidatesSelection[i]) {
                        unclusteredMatches->push_back(matchIndex);
                    }
                }
                std::cout << "Clusterable matches: " << unclusteredMatches->size() << std::endl;
                std::sort(unclusteredMatches->begin(), unclusteredMatches->end(), std::less<int32_t>());
                Cluster* cluster = new Cluster(clusters.size(), sequenceIndex);
                clusteredSequences[sequenceIndex] = true;
                for(uint64_t i = 0; i < unclusteredMatches->size(); i++) {
                    if(clusteredSequences[(*unclusteredMatches)[i]]) {
                        continue;
                    }
                    double pi = getFastPI(sequenceIndex, (*unclusteredMatches)[i]);
                    if( pi < clusterConfig.percentIdentityThreshold) {
                        continue;
                    }
                    bool clusterable = true;
                    for(uint64_t j = 0; j < cluster->clusterItems.size(); j++) {
                        if(getFastPI(cluster->clusterItems[j].first, (*unclusteredMatches)[i]) < clusterConfig.percentIdentityThreshold) {
                            clusterable = false;
                            break;
                        }
                    }
                    if(clusterable) {
                        cluster->clusterItems.push_back(make_pair((*unclusteredMatches)[i], pi));
                        clusteredSequences[(*unclusteredMatches)[i]] = true;
                    }
                }
                clusters.push_back(*cluster);
                std::cout << "Cluster size " << cluster->clusterItems.size() << std::endl;
                return true;
            }

            double getFastPI(uint64_t firstSequenceIndex, uint64_t secondsequenceIndex){
                return fastPercentIdentity(sequences[firstSequenceIndex], sequences[secondsequenceIndex], clusterConfig.percentIdentityThreshold);
            }

            void generateClusters(){
                for(int32_t i = 0; i < (int32_t)sequences.size(); i++) {
                    if( i > 0 && i%2000 == 0) {
                        int32_t dataOffset = i - offset;
                        offset = i;
                        reinitialize(dataOffset);
                        std::cout << "Reinitialized index using offset " << offset << std::endl;
                        std::cout << "Processed sequences: " << i - 1000 << " - " << i << std::endl;
                    }
                    cluster(i);
                    if( i > 0 && i%1000 == 0) {
                        std::cout << "Processed sequences: " << i - 1000 << " - " << i << std::endl;
                    }
                }
            }

            void printClusters(ofstream& clustersFile) {
                for(uint64_t i = 0; i < clusters.size(); i++) {
                    clustersFile << ">ClusterId: " << i << endl;
                    clusters[i].printClusterInfo(clustersFile, sequencesLocations);
                }
            }

        };
    }
}
