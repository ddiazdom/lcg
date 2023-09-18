//
// Created by Diaz, Diego on 26.4.2023.
//

#ifndef SIMPLE_EXAMPLE_PARSING_STRATEGIES_H
#define SIMPLE_EXAMPLE_PARSING_STRATEGIES_H

#include <sys/types.h>
#include "cds/int_array.h"
#include <vector>
#include "hashing.h"

struct parsing_opts{
    off_t chunk_size{};
    size_t active_chunks{};
    size_t n_threads{};
    off_t page_cache_limit = 262144000;
    size_t max_mt_sym_in_buff{};
    std::vector<hashing> p_functions;

    //vbyte compression paramters
    float vbyte_threshold{};//if the vbyte compression achieves this threshold, then we use vbyte encoding for the parse
    size_t vbyte_alphabet_threshold{};// if the alphabet size is below this value, then we apply vbyte compression

    //this is the state of the parse
    size_t sep_sym{};
    size_t n_sym;
    size_t max_sym;
    size_t tot_phrases; // number of phrases in T^{i} (i.e., the number of symbols in T^{i+1}
    size_t vbyte_size;
    size_t p_round=0;
    uint8_t p_alph_bytes;// number of bytes required by the parse's alphabet
    std::vector<uint32_t> vbyte_sym_perm;
    std::vector<uint32_t> new_vbyte_sym_perm;
    std::vector<off_t> *str_ptr= nullptr;
    float vbyte_comp_ratio;
    bool parse_compressible=false;
    hashing p_func;
    //

    void next_p_function(){
        if(p_round<p_functions.size()){
            p_func = p_functions[p_round];
        }else{
            p_func = hashing();
        }
    }
};

#include "semi-external-strategy/local_minima_parser.hpp"
#include "semi-external-strategy/create_dict.h"
#include "semi-external-strategy/create_dict_full_scan.h"
#include "semi-external-strategy/parse_text.h"

#endif //SIMPLE_EXAMPLE_PARSING_STRATEGIES_H
