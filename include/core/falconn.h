//
// Created by srikanthmaturu on 11/3/2017.
//

#pragma once

#include <cstdlib>

namespace SequencesAnalyzer{
    namespace core {
        struct FALCONNIndexConfiguration{
            uint8_t numberOfHashTables, numberOfHashBits, ngl, lshType, dataset_type, data_type;
            double threshold;
            uint64_t numberOfProbes;
        };
    }
}