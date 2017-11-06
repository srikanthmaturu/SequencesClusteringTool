//
// Created by srikanthmaturu on 11/3/2017.
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

namespace SequencesAnalyzer {
    namespace core  {
        class SimilarFastaSequencesFinder   {
        public:
            SimilarFastaSequencesFinder() = default;

            SimilarFastaSequencesFinder(const SimilarFastaSequencesFinder &) = default;

            SimilarFastaSequencesFinder(SimilarFastaSequencesFinder &&) = default;

            SimilarFastaSequencesFinder &operator=(const SimilarFastaSequencesFinder &) = default;

            SimilarFastaSequencesFinder &operator=(SimilarFastaSequencesFinder &&) = default;

            SimilarFastaSequencesFinder(std::vector<std::string>& fastaSequences, FALCONNIndexConfiguration config){
                FALCONNConfig = config;
                this->fastaSequences = fastaSequences;
            }

        public:
            FALCONNIndexConfiguration FALCONNConfig;
            std::vector<std::string> fastaSequences;
            tf_idf_falconn_index::tf_idf_falconn_idx<5, false, false, false, 2, 32, 11, 2032, 150, 0, DenseVectorFloat> idx;

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
                idx.initialize(fastaSequences);
                idx.construct_table();
            }

            std::pair<uint64_t, std::vector<int32_t>*> match(std::string fastaSequence){
                return idx.getNearestNeighbours(fastaSequence);
            }

            std::vector<std::vector<int32_t>>& processBatch(std::vector<std::string>& fastaSequences, uint64_t offset, uint64_t batchSize){
                std::vector<std::vector<int32_t>> * batchResults = new std::vector<std::vector<int32_t>>(batchSize, std::vector<int32_t>());
                std::map<uint64_t, unique_ptr<falconn::LSHNearestNeighborQuery<DenseVectorFloat>>> queryObjects;
#pragma omp parallel
                {
#pragma omp single
                    for (int64_t threadId = 0; threadId < omp_get_num_threads(); threadId++) {
                        queryObjects[threadId] = idx.createQueryObject();
                    }

#pragma omp for
                    for (uint64_t i = 0; i < batchSize; i++) {
                        auto result = idx.getNearestNeighbours(queryObjects[omp_get_thread_num()], fastaSequences[offset + i]);
                        (*batchResults)[offset + i] = *(result.second);
                    }
                }
                return *batchResults;
            }
        };
    }
}
