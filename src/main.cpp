//
// Created by srikanthmaturu on 10/17/2017.
//

#include "core/cluster.h"
#include <iostream>
#include <string>

using namespace std;
using namespace core;

int main(int argc, char** argv){

    if(argc < 2 ){
        cout << "Usage: " << argv[0] << " [fasta_file] " << endl;
        return 1;
    }

    uint64_t kmerSize = 100;
    FALCONNIndexConfiguration falconnConfig;
    falconnConfig.lshType = 2;
    falconnConfig.ngl = 5;
    falconnConfig.numberOfHashBits = 11;
    falconnConfig.numberOfProbes = 2032;
    falconnConfig.threshold = 1.5;
    falconnConfig.numberOfHashTables = 32;
    ClusterConfiguration clusterConfig;
    clusterConfig.editDistanceThreshold = 30;
    clusterConfig.kmerSize = kmerSize;

    vector<string> sequences;
    getFastaSequences(argv[1], sequences, kmerSize);
    ClustersGenerator generator(sequences, falconnConfig, clusterConfig);
}
