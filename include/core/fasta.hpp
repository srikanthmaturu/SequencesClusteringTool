#pragma  once

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <regex>
#include "xxhash.h"

uint64_t getHash(std::string sequence){
    return XXH64(sequence.c_str(), 100, 0xcc9e2d51);
}

namespace fasta {

    void getFastaSequences(std::string fastaFileName, std::vector<std::string>& sequences, uint64_t subSequenceWidth = 0){
        std::ifstream fastaFile(fastaFileName, std::ifstream::in);

        std::regex e("^>");
        std::smatch m;

        while(!fastaFile.eof()){
            std::string line;
            std::getline(fastaFile, line);

            if(!regex_search(line, e) && line.size() >= subSequenceWidth){
                sequences.push_back(line);
            }
        }

        fastaFile.close();
        return;
    }

    void getSubSequencesMap(std::vector<std::string>& sequences, std::vector<std::string>& subSequences, std::map<uint64_t, std::vector<uint64_t>> subSequencesMap, uint8_t windowSize){
        for(uint64_t i = 0; i < sequences.size(); i++){
            for(uint16_t j = 0; j < sequences[i].size() - windowSize + 1; j++) {
                std::string subSeqeunce = sequences[i].substr(j, windowSize);
                subSequences.push_back(subSeqeunce);
                uint64_t hash = getHash(subSeqeunce);
                auto element = subSequencesMap.find(hash);
                if(element == subSequencesMap.end()){
                    std::vector<uint64_t> * vec = new std::vector<uint64_t>();
                    vec->push_back(i);
                    subSequencesMap[hash] = *vec;
                }
                else {
                    (*element).second.push_back(i);
                }
            }
        }
    }
}