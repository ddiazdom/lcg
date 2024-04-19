//
// Created by Diaz, Diego on 19.4.2024.
//

#ifndef SE_STRAT_PARSING_OPTS_H
#define SE_STRAT_PARSING_OPTS_H

struct parsing_opts{
    off_t chunk_size{};
    size_t active_chunks{};
    size_t n_threads{};
    off_t page_cache_limit = 262144000;
    size_t max_mt_sym_in_buff{};
    std::vector<uint64_t> p_seeds;

    //vbyte compression parameters
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
    //hashing p_func;
    uint64_t p_seed{};
    //

    void next_seed(){
        p_seed = p_seeds[p_round];
    }
};
#endif //SE_STRAT_PARSING_OPTS_H
