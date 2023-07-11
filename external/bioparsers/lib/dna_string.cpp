//
// Created by Diaz, Diego on 15.10.2021.
//
#include "dna_string.h"

uint8_t dna_string::comp[170] ={0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,84,0,71,0,0,0,67,0,0,0,0,0,0,0,0,
                                0,0,0,0,65,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,168,0,155,0,0,0,151,0,0,0,0,
                                0,0,0,0,0,0,0,0,149,0};

std::ostream& operator<<(std::ostream& os, const dna_string& seq) {
    size_t i=0;
    for(auto const& sym : seq){
        if(sym<=dna_string::dna_offset){
            os << char(sym);
        }else{
            os << "["<<char(sym-dna_string::dna_offset) <<","<<seq.homopolymer_length(i++)<<"]";
        }
    }
    return os;
}
