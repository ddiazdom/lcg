//
// Created by Diaz, Diego on 26.4.2023.
//

#ifndef SIMPLE_EXAMPLE_PARSING_STRATEGIES_H
#define SIMPLE_EXAMPLE_PARSING_STRATEGIES_H

#include <sys/types.h>
#include <cds/int_array.h>
#include <vector>

struct parsing_opts{
    off_t chunk_size{};
    size_t active_chunks{};
    size_t n_threads{};
    off_t page_cache_limit = 262144000;
    size_t max_mt_sym_in_buff{};

    //vbyte compression paramters
    float vbyte_threshold{};//if the vbyte compression achieves this threshold, then we use vbyte encoding for the parse
    size_t vbyte_alphabet_threshold{};// if the alphabet size is below this value, then we apply vbyte compression

    //this is the state of the parse
    size_t n_sym;
    size_t max_sym;
    size_t tot_phrases; // number of phrases in T^{i} (i.e., the number of symbols in T^{i+1}
    size_t vbyte_size;
    size_t p_round=0;
    uint8_t p_alph_bytes;// number of bytes required by the parse's alphabet
    std::vector<uint32_t> sym_perm;
    std::vector<uint32_t> new_sym_perm;
    std::vector<off_t> *str_ptr= nullptr;
    float vbyte_comp_ratio;
    bool parse_compressible=false;
    //
};

#include "parsers/local_minima_parser.hpp"
#include "create_dict.h"
#include "create_dict_full_scan.h"
#include "parse_text.h"

#endif //SIMPLE_EXAMPLE_PARSING_STRATEGIES_H
