//
// Created by srikanthmaturu on 10/17/2017.
//

#pragma once

#include <cstdlib>
#include <vector>
#include <string>
#include <omp.h>
#include <iterator>

#include "fasta.hpp"
#include "falconn.h"
#include "tf_idf_falconn_idx.hpp"

namespace SequencesAnalyzer{
    namespace core{
        struct ClusterConfiguration{
            uint64_t kmerSize, editDistanceThreshold;
            ClusterConfiguration(){

            }
            ClusterConfiguration(uint64_t kmerSize):kmerSize(kmerSize){

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
                FALCONNConfig = config;
                this->clusterConfig = clusterConfig;
                fasta::getSubSequencesMap(sequences, kmers, kmersToSequencesMap, clusterConfig.kmerSize);
            }

        public:
            FALCONNIndexConfiguration FALCONNConfig;
            ClusterConfiguration clusterConfig;

            std::vector<std::string> kmers;
            std::map<uint64_t,std::vector<uint64_t>> kmersToSequencesMap;

            void generateClusters(){
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
                tf_idf_falconn_index::tf_idf_falconn_idx<5, false, false, false, 2, 32, 11, 2032, 150, DenseVectorFloat>  idx(lshParams);
                idx.setThreshold(FALCONNConfig.threshold);
                idx.setNGL(FALCONNConfig.ngl);
                idx.setNumberOfProbes(FALCONNConfig.numberOfProbes);
                idx.initialize(kmers);
                idx.construct_table();

                std::map<uint64_t, unique_ptr<falconn::LSHNearestNeighborQuery<DenseVectorFloat>>> queryObjects;

                std::cout << "Index Creation Complete." << std::endl;
#pragma omp parallel
                {
#pragma omp single
                    for (int64_t threadId = 0; threadId < omp_get_num_threads(); threadId++) {
                        queryObjects[threadId] = idx.createQueryObject();
                    }

#pragma omp for
                    for (uint64_t i = 0; i < kmers.size(); i++) {
                        std::vector<int32_t> nearestNeighbours;
                        idx.getNearestNeighboursByEditDistance(queryObjects[omp_get_thread_num()], kmers[i], nearestNeighbours,
                                                               clusterConfig.editDistanceThreshold);
                        auto bucket = kmersToSequencesMap[getHash(kmers[i])];
                        bucket.insert(bucket.end(), nearestNeighbours.begin(), nearestNeighbours.end());
                    }
                }
                std::cout << "Buckets Creation Complete." << std::endl;
#pragma omp parallel for
                for(uint64_t i = 0; i < kmersToSequencesMap.size(); i++){
                    std::map<long unsigned int, std::vector<long unsigned int> >::iterator it(kmersToSequencesMap.begin());
                    std::advance(it, i);
                    std::sort((*it).second.begin(),(*it).second.end());
                    auto p = std::unique((*it).second.begin(),(*it).second.end());
                    (*it).second.resize(std::distance((*it).second.begin(),p));
                }
            }

            std::vector<uint64_t>& getSequencesCluster(std::string sequence, uint64_t windowSize){
                std::vector<uint64_t> * candidates = new std::vector<uint64_t>();
                for(uint16_t i = 0; i < sequence.size() - windowSize + 1; i++) {
                    uint64_t hash = getHash(sequence.substr(i, windowSize));
                    if(kmersToSequencesMap.find(hash) == kmersToSequencesMap.end()){
                        continue;
                    }
                    candidates->insert(candidates->end(), kmersToSequencesMap[hash].begin(), kmersToSequencesMap[hash].end());
                }
                return *candidates;
            }

        };
    }
}
