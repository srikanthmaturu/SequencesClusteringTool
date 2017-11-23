//
// Created by Srikanth Maturu on 11/11/17.
//

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include "similar_fasta_sequences.h"
#include <chrono>
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

namespace SequencesAnalyzer {
    namespace core {
        class ParametersTuner: public SimilarFastaSequencesFinder {
            using SimilarFastaSequencesFinder::SimilarFastaSequencesFinder;
        public:
            double process_queries_parallely(vector<string>& queries, uint64_t number_of_blocks, uint64_t block_size, vector<vector<string>> linear_results){
                uint64_t realMatchesCount = 0, actualMatchesCount = 0, queries_size = queries.size();
//#pragma omp parallel for
                for(uint64_t bi = 0; bi < number_of_blocks; bi++){
                    uint64_t block_end = (bi == (number_of_blocks-1))? queries_size : (bi + 1)*block_size;
                    for(uint64_t i= bi * block_size, j = 0; i< block_end; i++, j++){
                        vector<string>& linear_res = linear_results[i];
                        auto res = idx.match(queries[i]);
                        auto cs_fp_fn_pair = get_comparison(linear_res, res.second);
                        realMatchesCount += (linear_res.size());
                        actualMatchesCount += (get<0>(cs_fp_fn_pair) - get<1>(cs_fp_fn_pair));
                        //cout << "Linear_res: " << linear_res.size() << " TruePositives: " << (get<0>(cs_fp_fn_pair) - get<1>(cs_fp_fn_pair)) <<  " FalsePositives: " << get<1>(cs_fp_fn_pair) << " FalseNegatives: " << (get<2>(cs_fp_fn_pair)) << endl;
                    }
                }

                double recall = (actualMatchesCount * 1.0) /(realMatchesCount * 1.0);
                return recall;
            }
            
            void process_queries_box_test(vector<string>& queries, double threshold, uint8_t maxPercentIdentity){
                ofstream box_test_results_file("box_test_results_NGL" + to_string(FALCONNConfig.ngl) + ".csv");
                uint64_t block_size = 100000;
                uint64_t queries_size = queries.size();
                std::cout << "Queries Size: " << queries_size << std::endl;
                if(queries_size < block_size){
                    block_size = queries_size;
                }
                uint64_t extra_block = queries_size % block_size;
                uint64_t number_of_blocks =  queries_size / block_size;
                auto start = timer::now();
                if(extra_block > 0) {
                    number_of_blocks++;
                }
                idx.setThreshold(threshold);
                vector<vector<uint64_t>> queries_linear_results(queries.size(),vector<uint64_t>(3,0));
                vector<vector<string>> linear_results;

                for(uint64_t bi = 0; bi < number_of_blocks; bi++){
                    uint64_t block_end = (bi == (number_of_blocks-1))? queries_size : (bi + 1)*block_size;
                    for(uint64_t i= bi * block_size, j = 0; i< block_end; i++, j++){
                        linear_results.push_back(idx.get_nearest_neighbours_by_linear_method_by_percent_identity(queries[i], maxPercentIdentity));
                    }
                }

                uint64_t polytope_vertices = FALCONNConfig.ngl * 2 + 1;

                for(uint64_t l = 32; l <= 256; l += 32){
                    for(uint64_t nhb = polytope_vertices; nhb <= 32; nhb += polytope_vertices){
                        uint64_t np_max = 0;
                        if(nhb == 7){
                            np_max = l * 3;
                        }
                        else {
                            np_max = l * 12;
                        }
                        uint64_t step = 100;
                        for(uint64_t np = l; np < np_max; np = np + step){
                            idx.updateParmeters(l, nhb, np);
                            cout << "Current LSH Parameters: " << endl;
                            idx.printLSHConstructionParameters();
                            idx.construct_table();

                            double recall = process_queries_parallely(queries, number_of_blocks, block_size, linear_results);
                            auto stop = timer::now();
                            box_test_results_file << recall << "," << (duration_cast<chrono::microseconds>(stop-start).count()/1000000.0)/(double)queries.size() << "," << l << "," << nhb << "," << np << endl;
                            if(np_max <= (np + step) && (recall < 0.95)){
                                step *= 2;
                                np_max = np + step * 2;
                            }
                            if(recall >= 0.99 ||  np_max > 1000000){
                                break;
                            }
                        }
                    }
                }
            }

            void process_queries_thresholds_test(vector<string>& queries){
                ofstream thresholds_test_results_file("thresholds_test_results");
                vector< vector< pair<string, uint64_t > > > query_results_vector;
                uint64_t block_size = 100000;
                uint64_t queries_size = queries.size();
                std::cout << queries_size << std::endl;
                if(queries_size < block_size){
                    block_size = queries_size;
                }

                cout << "Current LSH Parameters: " << endl;
                idx.printLSHConstructionParameters();
                idx.construct_table();

                vector<vector<uint64_t>> queries_linear_results(queries.size(),vector<uint64_t>(3,0));
                uint64_t extra_block = queries_size % block_size;
                uint64_t number_of_blocks =  queries_size / block_size;

                if(extra_block > 0) {
                    number_of_blocks++;
                }
                thresholds_test_results_file << "Query Index, tp, (ED <= 15), (15 < ED <= 20), (20 < ED <= 25), (ED > 25), ";
                for(double th = 10; th <= 150; th += 10) {
                    thresholds_test_results_file << "Threshold-" << th / 100.0 << " candidates, (ED <= 15), (15 < ED <= 20), (20 < ED <= 25), (ED > 25) , fp, fn" ;
                    if(th < 150){
                        thresholds_test_results_file << ",";
                    }
                }
                thresholds_test_results_file << endl;
                cout << "Number of blocks: " << number_of_blocks << endl;
                for(uint64_t bi = 0; bi < number_of_blocks; bi++){
                    uint64_t block_end = (bi == (number_of_blocks-1))? queries_size : (bi + 1)*block_size;
                    query_results_vector.resize(block_size);
                    cout << "Block Index: " << bi << endl;
                    cout << "Current Block Size: " << block_end - bi * block_size << endl;
                    //#pragma omp parallel for
                    for(uint64_t i= bi * block_size, j = 0; i< block_end; i++, j++){
                        auto linear_res = idx.get_nearest_neighbours_by_linear_method(queries[i], 30);
                        //std::cout << "Getting Linear Category Counts" << std::endl;
                        auto edCategoryCounts = idx.getCategoryCounts(queries[i], linear_res);
                        //std::cout << "Done" << std::endl;
                        thresholds_test_results_file << i << "," << linear_res.size() << ",";
                        for(auto it : edCategoryCounts){
                            thresholds_test_results_file << it.second << ",";
                        }

                        for(double th = 10; th <= 150; th += 10) {
                            idx.setThreshold(th/100.0);
                            auto res = idx.match(queries[i]);
                            //std::cout << "Falconnn Threshold: " << th/100.0 << std::endl;
                            //std::cout << "Getting Falconn Category Counts"<< std::endl;
                            auto categoryCounts = idx.getCategoryCounts(queries[i], res.second);
                            //std::cout << "Done" << std::endl;
                            auto cs_fp_fn_pair = get_comparison(linear_res, res.second);
                            thresholds_test_results_file << std::get<0>(cs_fp_fn_pair);
                            for(auto it : categoryCounts){
                                thresholds_test_results_file << "," << it.second;
                            }
                            thresholds_test_results_file << "," << std::get<1>(cs_fp_fn_pair) << "," << std::get<2>(cs_fp_fn_pair);
                            if(th < 150){
                                thresholds_test_results_file << ",";
                            }
                        }
                        //cout << "Processed " << i << endl;
                        thresholds_test_results_file << endl;
                    }
                    query_results_vector.clear();
                }
            }

            void process_queries_linear_test(vector<string>& queries){
                ofstream results_file("linear_test_results");
                for(string query: queries){
                    idx.linear_test(query, results_file);
                }
            }
        };
    }
}

