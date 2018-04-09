//
// Created by srikanthmaturu on 4/9/2018.
//

#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include <seqan/align.h>
#include <seqan/graph_msa.h>


/**This function will accept clustering results of falconn or cd-hit and returns average percentIdentity between MSA pairs and also maximum percent identity

 */
std::tuple<int, int,double> performMSA(std::vector<int>& clusteringResult, std::vector<std::string>& sequences){
    seqan::Align<seqan::String<seqan::AminoAcid>> align;
    seqan::resize(seqan::rows(align), clusteringResult.size());
    for (int i = 0; i < clusteringResult.size(); i++) {
        seqan::assignSource(seqan::row(align, i), sequences[clusteringResult[i]]);
    }
    seqan::globalMsaAlignment(align, seqan::Blosum62(-1, -11));
    //std::cout << align << "\n";
    std::vector<int> pairwisePercentIdentities;
    uint64_t numberOfPairs = 0;
    for(int i = 0; i < clusteringResult.size(); i++) {
        for(int j = i+1; j < clusteringResult.size(); j++ ) {
            uint64_t matches = 0;
            for(int k = 0; k < seqan::length(seqan::row(align, 0)); k++) {
                if(seqan::ordValue(seqan::getValue(seqan::row(align, i), k)) == seqan::ordValue(seqan::getValue(seqan::row(align, j), k))) {
                    matches++;
                }
            }
            int percentIdentity = ((matches*1.0)/seqan::length(seqan::row(align, 0))) * 100;
            /*if(percentIdentity > 100) {
                std::cout << " pi: " << percentIdentity << " matches: " << matches << " al: " << seqan::length(seqan::row(align, 0)) << std::endl;
            }*/

            //std::cout << " pi: " << percentIdentity << std::endl;
            pairwisePercentIdentities.push_back(percentIdentity);
            numberOfPairs++;
        }
    }
    int sum = std::accumulate(pairwisePercentIdentities.begin(), pairwisePercentIdentities.end(),0);
    double average = (sum*1.0)/numberOfPairs;
    //std::cout << "Average: " << average << " s: " << sum << " c: " << clusteringResult.size() << std::endl;
    return std::make_tuple(clusteringResult.size(),*std::min_element(pairwisePercentIdentities.begin(), pairwisePercentIdentities.end()), average);
}

bool compareByVectorSize(std::vector<int>& v1, std::vector<int>& v2) {
    return v1.size() > v2.size() ? true : false;
}

std::vector<std::tuple<int, int,double>> getClusteringAlgorithmsMSAResults(std::vector<std::vector<int>>& clusterAlgorithmResults, std::vector<std::string>& sequences, int& numberOfSingleClusters){
    std::vector<std::tuple<int, int,double>> msaResults;
    std::sort(clusterAlgorithmResults.begin(), clusterAlgorithmResults.end(), compareByVectorSize);
    for(int i = 0; i < clusterAlgorithmResults.size(); i++) {
        if(clusterAlgorithmResults[i].size() > 1) {
            auto p = performMSA(clusterAlgorithmResults[i], sequences);
            //std::cout << std::get<0>(p) << " " << std::get<1>(p) << " " << std::get<2>(p) << std::endl;
            msaResults.push_back(p);
        }
        else {
            numberOfSingleClusters++;
        }
    }
    return msaResults;
}