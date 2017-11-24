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

            void printClusterInfo(ofstream outputFileStream) {
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
                sortSequences(this->sequences);
                FALCONNConfig = config;
                this->clusterConfig = clusterConfig;
                similarityMatrix = std::vector<std::vector<double>>(this->sequences.size(), std::vector<double>(this->sequences.size(),-1));
                for(uint64_t i = 0; i < similarityMatrix.size(); i++) {
                    similarityMatrix[i][i] = 100;
                }
                clusteredSequences = std::vector<bool>(sequences.size(), false);
            }

        public:
            FALCONNIndexConfiguration FALCONNConfig;
            ClusterConfiguration clusterConfig;

            std::vector<std::string> sequences;
            std::vector<std::vector<double>> similarityMatrix;
            std::vector<bool> clusteredSequences;
            std::vector<Cluster> clusters;
            uint64_t minimumPercentIdentity;

            tf_idf_falconn_index::tf_idf_falconn_idx<5, false, false, false, 2, 32, 11, 2032, 150, 0, 0, POINT_TYPE> idx;

            void sortSequences(std::vector<std::string>& sequences){
                std::sort(sequences.begin(), sequences.end(), [](std::string a, std::string b) { return a.size() >= b.size(); });
            }

            void initialize(){
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
                idx.construct_table();
            }

            std::pair<uint64_t, std::vector<int32_t>*> match(std::string sequence){
                return idx.getNearestNeighbours(sequence);
            }

            bool cluster(int32_t sequenceIndex) {
                if(clusteredSequences[sequenceIndex]) {
                    return false;
                }
                auto matchesPair = match(sequences[sequenceIndex]);
                std::vector<int32_t>* unclusteredMatches = new std::vector<int32_t>();
                for(uint64_t i = 0; i < matchesPair.second->size(); i++) {
                    int32_t matchIndex = (*(matchesPair.second))[i];
                    if(!clusteredSequences[matchIndex]){
                        unclusteredMatches->push_back(matchIndex);
                    }
                }
                std::sort(unclusteredMatches->begin(), unclusteredMatches->end(), std::less<int32_t>());
                Cluster* cluster = new Cluster(clusters.size(), sequenceIndex);
                clusteredSequences[sequenceIndex] = true;
                for(uint64_t i = 0; i < unclusteredMatches->size(); i++) {
                    if(clusteredSequences[(*unclusteredMatches)[i]]) {
                        continue;
                    }
                    if(!isSimilar(sequenceIndex, (*unclusteredMatches)[i])) {
                        continue;
                    }
                    bool clusterable = true;
                    for(uint64_t j = 0; j < cluster->clusterItems.size(); j++) {
                        if(!isSimilar(cluster->clusterItems[j].first, (*unclusteredMatches)[i])) {
                            clusterable = false;
                            break;
                        }
                    }
                    if(clusterable) {
                        cluster->clusterItems.push_back(make_pair((*unclusteredMatches)[i], similarityMatrix[(*unclusteredMatches)[i]][sequenceIndex]));
                        clusteredSequences[(*unclusteredMatches)[i]] = true;
                    }
                }
                clusters.push_back(*cluster);
                return true;
            }

            bool isSimilar(uint64_t firstSequenceIndex, uint64_t secondsequenceIndex){
                double percentIdentity = similarityMatrix[firstSequenceIndex][secondsequenceIndex];
                if(percentIdentity < 0) {
                    percentIdentity = fastPercentIdentity(sequences[firstSequenceIndex], sequences[secondsequenceIndex], clusterConfig.percentIdentityThreshold);
                    similarityMatrix[firstSequenceIndex][secondsequenceIndex] = percentIdentity;
                    similarityMatrix[secondsequenceIndex][firstSequenceIndex] = percentIdentity;
                }
                if(percentIdentity < clusterConfig.percentIdentityThreshold) {
                    return false;
                } else {
                    return true;
                }
            }

        };
    }
}
