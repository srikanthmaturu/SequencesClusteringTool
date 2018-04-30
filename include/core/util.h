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
#include <seqan/graph_types.h>



/**This function will accept clustering results of falconn or cd-hit and returns average percentIdentity between MSA pairs and also maximum percent identity

 */
std::tuple<int, int,double> performMSA(std::vector<int>& clusteringResult, std::vector<std::string>& sequences){
    seqan::Align<seqan::String<seqan::AminoAcid>> align;
    seqan::resize(seqan::rows(align), clusteringResult.size());
    for (uint64_t i = 0; i < clusteringResult.size(); i++) {
        seqan::assignSource(seqan::row(align, i), sequences[clusteringResult[i]]);
    }
    seqan::globalMsaAlignment(align, seqan::Blosum62(-1, -11));
    //std::cout << align << "\n";

    std::vector<int> pairwisePercentIdentities;
    uint64_t numberOfPairs = 0;
    for(uint64_t i = 0; i < clusteringResult.size(); i++) {
        for(uint64_t j = i+1; j < clusteringResult.size(); j++ ) {
            uint64_t matches = 0;
            for(uint64_t k = 0; k < seqan::length(seqan::row(align, 0)); k++) {
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

double computeUPGMATree(std::vector<int>& clusteringResult, std::vector<std::string>& sequences) {
    seqan::String<double> distanceMatrix;
    seqan::resize(distanceMatrix,clusteringResult.size()*clusteringResult.size(),0);
    //seqan::fill(distanceMatrix, clusteringResult.size()*clusteringResult.size(), 0);

    for(uint64_t i = 0; i < clusteringResult.size(); i++) {
        double pi = getPercentIdentity(sequences[clusteringResult[i]], sequences[clusteringResult[i]]);
        assignValue(distanceMatrix, i * clusteringResult.size() + i, (100 - pi)/100);
        for(uint64_t j = 0; j < i; j++) {
            double pi = getPercentIdentity(sequences[clusteringResult[i]], sequences[clusteringResult[j]]);
            assignValue(distanceMatrix, i * clusteringResult.size() + j, (100 - pi)/100);
            assignValue(distanceMatrix, j * clusteringResult.size() + i, (100 - pi)/100);
        }
    }

    typedef seqan::Graph<seqan::Tree<double> > TGraph;
    typedef seqan::VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef seqan::EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

    TGraph upgmaTreeOut;
    seqan::upgmaTree(distanceMatrix,upgmaTreeOut);

    double root_leaf_distance = 0.0;
    /*
    typedef seqan::Iterator<TGraph, DfsPreorder>::Type TDfsPreorder;
    TDfsPreorder dfsIt(upgmaTreeOut, seqan::getRoot(upgmaTreeOut));
     */
    TVertexDescriptor currentNode = seqan::getRoot(upgmaTreeOut);
    typedef seqan::Iterator<TGraph, seqan::AdjacencyIterator>::Type TAdjacencyIterator;

    while(seqan::outDegree(upgmaTreeOut, currentNode) > 0) {
        TAdjacencyIterator adjacencyIterator(upgmaTreeOut,currentNode);
        TEdgeDescriptor ed = seqan::findEdge(upgmaTreeOut, currentNode, value(adjacencyIterator));
        root_leaf_distance += getCargo(ed);
        currentNode = targetVertex(upgmaTreeOut, ed);
    }

    //std::cout << "Root to leaf distance:- " << root_leaf_distance << std::endl;

    typedef seqan::Iterator<TGraph, seqan::EdgeIterator>::Type TEdgeIter;
    TEdgeIter edIt(upgmaTreeOut);
    //std::cout << "Edges: " << std::endl;
    /*
    for(;!atEnd(edIt);goNext(edIt)) {
        std::cout << sourceVertex(edIt) << " -- " << targetVertex(edIt) << ": Weight = " << getCargo(*edIt) << std::endl;
    }
    */
    return root_leaf_distance;

    //std::cout << upgmaTreeOut << std::endl;
}

double computeJaccardIndex(std::vector<int> firstCluster,std::vector<int> secondCluster) {
    uint64_t firstClusterIndex = 0, secondClusterIndex = 0;
    uint64_t in = 0, un = 0;
    while(firstClusterIndex < firstCluster.size() && secondClusterIndex < secondCluster.size()) {
        if(firstCluster[firstClusterIndex] == secondCluster[secondClusterIndex]) {
            in++;
            un++;
            firstClusterIndex++;
            secondClusterIndex++;
        }
        else if(firstCluster[firstClusterIndex] < secondCluster[secondClusterIndex]) {
            while(firstCluster[firstClusterIndex] < secondCluster[secondClusterIndex]) {
                un++;
                firstClusterIndex++;
            }
        }
        else {
            while(firstCluster[firstClusterIndex] > secondCluster[secondClusterIndex]) {
                un++;
                secondClusterIndex++;
            }
        }
    }

    return (in * 1.0) / (un * 1.0);
}

double computeMaximumJaccardIndex(std::vector<int>& firstClusteringAlgorithmCluster, std::vector<std::vector<int>>& secondClusteringAlgorithmClusters) {
    double maximumJaccardIndex = 0.0;

    for(uint64_t i = 0; i < secondClusteringAlgorithmClusters.size(); i++) {
        double maxJaccardIndexPossible = 1.0*(std::min(firstClusteringAlgorithmCluster.size(), secondClusteringAlgorithmClusters[i].size()))/1.0*(std::max(firstClusteringAlgorithmCluster.size(), secondClusteringAlgorithmClusters[i].size()));
        if(maxJaccardIndexPossible < maximumJaccardIndex) continue;

        double ji = computeJaccardIndex(firstClusteringAlgorithmCluster, secondClusteringAlgorithmClusters[i]);
        if(ji > maximumJaccardIndex) {
            maximumJaccardIndex = ji;
        }
    }

    return maximumJaccardIndex;
}

double computeAverageMaximumJaccardIndex(std::vector<std::vector<int>>& firstClusteringAlgorithmResults, std::vector<std::vector<int>>& secondClusteringAlgorithmResults, uint64_t maxNumberOfClusters) {
    double average = 0.0;

    for(uint64_t i = 0; i < maxNumberOfClusters && i < firstClusteringAlgorithmResults.size(); i++) {
        double mji = computeMaximumJaccardIndex(firstClusteringAlgorithmResults[i], secondClusteringAlgorithmResults);
        average += mji;
    }

    return average/maxNumberOfClusters;
}

bool compareByVectorSize(std::vector<int>& v1, std::vector<int>& v2) {
    return v1.size() > v2.size() ? true : false;
}

std::vector<std::tuple<int, int,double,double>> getClusteringAlgorithmsMSAResults(std::vector<std::vector<int>>& clusterAlgorithmResults, std::vector<std::string>& sequences, int& numberOfSingleClusters){
    std::vector<std::tuple<int, int,double,double>> msaResults;
    std::sort(clusterAlgorithmResults.begin(), clusterAlgorithmResults.end(), compareByVectorSize);
    for(uint64_t i = 0; i < clusterAlgorithmResults.size(); i++) {
        if(clusterAlgorithmResults[i].size() > 1) {
            auto p = performMSA(clusterAlgorithmResults[i], sequences);
            //std::cout << std::get<0>(p) << " " << std::get<1>(p) << " " << std::get<2>(p) << std::endl;
            double rootToLeafDistance = computeUPGMATree(clusterAlgorithmResults[i], sequences);
            std::tuple<int, int,double,double> n_p = std::tuple_cat(p,std::make_tuple(rootToLeafDistance));
            msaResults.push_back(n_p);
        }
        else {
            numberOfSingleClusters++;
        }
    }
    return msaResults;
}

int findClusterBySequenceId(int sequenceId, std::vector<std::vector<int>>& clusterAlgorithmResults) {
    for(uint64_t i = 0; i < clusterAlgorithmResults.size(); i++) {
        for(uint64_t j = 0; j < clusterAlgorithmResults[i].size(); j++) {
            if(clusterAlgorithmResults[i][j] == sequenceId) {
                return i;
            }
        }
    }
    //std::cout << "Couldn't find a cluster with the input sequence id. Critical error. " << std::endl;
    //std::cout << "Sequence Id searched " << sequenceId << std::endl;
    return -1;
}

std::vector<std::vector<double>> compareClusteringAlgorithmsByBruteForceMethod(std::vector<std::vector<std::vector<int>>>& algorithmsClusteringResults, std::vector<std::string>& sequences){
    std::vector<std::vector<double>> comparisonResults;
    for(uint64_t i = 0; i < algorithmsClusteringResults.size(); i++) {
        std::sort(algorithmsClusteringResults[i].begin(), algorithmsClusteringResults[i].end(), compareByVectorSize);
    }

    for(uint64_t i = 0; i < algorithmsClusteringResults[0].size(); i++) {
        bool failed = false;
        std::vector<double> clusterComparisonResult;
        std::tuple<int, int,double,double> groundTruthClusterMSAResult;
        if(algorithmsClusteringResults[0][i].size() > 1) {
            double rootToLeafDistance = computeUPGMATree(algorithmsClusteringResults[0][i], sequences);
            groundTruthClusterMSAResult = std::tuple_cat(performMSA(algorithmsClusteringResults[0][i], sequences),std::make_tuple(rootToLeafDistance));
        }
        else {
            groundTruthClusterMSAResult = std::make_tuple(1,0,0,1);
        }

        clusterComparisonResult.push_back(get<0>(groundTruthClusterMSAResult));
        clusterComparisonResult.push_back(get<1>(groundTruthClusterMSAResult));
        clusterComparisonResult.push_back(get<2>(groundTruthClusterMSAResult));
        clusterComparisonResult.push_back(get<3>(groundTruthClusterMSAResult));
        for(uint64_t j = 1; j < algorithmsClusteringResults.size(); j++) {
            double clusterId = findClusterBySequenceId(algorithmsClusteringResults[0][i][0],algorithmsClusteringResults[j]);
            /*double clusterSize = algorithmsClusteringResults[j][clusterId].size();
            clusterComparisonResult.push_back(clusterSize);*/
            if(clusterId == -1) {
                failed = true;
            }
            std::tuple<int, int,double,double> p;
            if(!failed && algorithmsClusteringResults[j][clusterId].size() > 1) {
                double rootToLeafDistance = computeUPGMATree(algorithmsClusteringResults[j][clusterId], sequences);
                p = std::tuple_cat(performMSA(algorithmsClusteringResults[j][clusterId], sequences),std::make_tuple(rootToLeafDistance));
            }
            else {
                p = std::make_tuple((failed?0:1),0,0,1);
            }
            clusterComparisonResult.push_back(get<0>(p));
            clusterComparisonResult.push_back(get<1>(p));
            clusterComparisonResult.push_back(get<2>(p));
            clusterComparisonResult.push_back(get<3>(p));
        }
        comparisonResults.push_back(clusterComparisonResult);
    }
    return comparisonResults;
}
