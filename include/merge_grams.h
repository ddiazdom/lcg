//
// Created by Diaz, Diego on 22.8.2024.
//

#ifndef LCG_MERGE_GRAMS_H
#define LCG_MERGE_GRAMS_H

#include "partial_gram.h"
#include "ext_gram.h"
#include "cds/memory_handler.hpp"

#ifdef DEBUG_MODE
#include "malloc_count.h"
#endif

#include <random>

struct merge_data_t{

    std::vector<uint64_t> fps;
    std::vector<uint32_t> map_a;
    std::vector<uint32_t> map_b;
    size_t lvl_sigma=0;
    size_t longest_rule=0;
    bitstream<size_t> buffer;

    void initialize(size_t lvl_sigma_, size_t sym_bytes, size_t longest_rule_){
        lvl_sigma = lvl_sigma_;
        map_a.resize(lvl_sigma);
        map_b.resize(lvl_sigma);
        fps.resize(lvl_sigma);
        for(size_t i=0;i<fps.size();i++){
            fps[i] = XXH3_64bits(&i, sym_bytes);
            map_a[i] = i;
            map_b[i] = i;
        }
        longest_rule = longest_rule_;
    }

    [[nodiscard]] size_t space_usage() const {
        size_t bytes = fps.capacity()*sizeof(uint64_t);
        bytes += map_a.capacity()*sizeof(uint32_t);
        bytes += map_b.capacity()*sizeof(uint32_t);
        return bytes;
    }

    /*void space_breakdown() const {
        std::cout<<report_space(fps.size()*sizeof(uint64_t))<<" "<<report_space(fps.capacity()*sizeof(uint64_t))<<std::endl;
        std::cout<<report_space(map_a.size()*sizeof(uint32_t))<<" "<<report_space(map_a.capacity()*sizeof(uint32_t))<<std::endl;
        std::cout<<report_space(map_b.size()*sizeof(uint32_t))<<" "<<report_space(map_b.capacity()*sizeof(uint32_t))<<std::endl;
    }*/

    ~merge_data_t(){
        destroy(fps);
        destroy(map_a);
        destroy(map_b);
        buffer.destroy();

#ifdef __linux__
        malloc_trim(0);
#endif
    }
};

/*
void merge_many_grams_in_serial(std::string& concat_grammars, size_t n_threads) {
    size_t rem = file_size(concat_grammars);
    int fd_r = open(concat_grammars.c_str(), O_RDONLY);
    malloc_count_reset_peak();
    extensible_gram sink_gram(n_threads);
    size_t read_bytes = sink_gram.load_from_file(fd_r);
    rem-=read_bytes;
    std::cout<<"Esto leimos en la primera: "<<report_space((off_t)read_bytes)<<" con un peak de "<<report_space((off_t)malloc_count_peak())<<std::endl;
    malloc_count_reset_peak();

    size_t n_mrg=0;
    while(rem>0){
        i_gram_stream new_gram(fd_r);
        read_bytes = sink_gram.merge_grams(new_gram, n_threads);
        std::cout<<"Esto leimos en la siguiente "<<report_space((off_t)read_bytes)<<" con un peak de "<<report_space((off_t)malloc_count_peak())<<" and size "<<report_space((off_t)sink_gram.ext_gram_size_in_bytes())<<" merge_round: "<<n_mrg<<std::endl;
        malloc_count_reset_peak();
        rem-=read_bytes;
        n_mrg++;
    }
}*/

template<class p_gram_type>
void merge_two_partial_grammars_in_memory(p_gram_type& p_gram_a, p_gram_type& p_gram_b, std::vector<uint64_t>& fp_seeds) {

    size_t longest_rule = std::max(p_gram_a.longest_rule, p_gram_b.longest_rule);
    assert(p_gram_a.par_seed == p_gram_b.par_seed);

    merge_data_t mg_data;

    //in the first level, tot_symbols represents the alphabet of terminals
    mg_data.initialize(p_gram_a.metadata[0].tot_symbols, sizeof(typename p_gram_type::sym_type), fp_seeds[0], longest_rule);

    //we subtract one because the last level contains the compressed strings,
    // and we do not merge but concatenate
    size_t max_lvl = std::max(p_gram_a.lvl, p_gram_b.lvl)-1;
    size_t min_lvl = std::min(p_gram_a.lvl, p_gram_b.lvl)-1;

    //pointer to the grammar with the least number of levels
    p_gram_type *st_gram = nullptr;
    std::vector<uint32_t> *st_gram_map = nullptr;

    if(p_gram_a.lvl<p_gram_b.lvl){
        st_gram = &p_gram_a;
        st_gram_map = &mg_data.map_a;
    } else if(p_gram_b.lvl<p_gram_a.lvl) {
        st_gram = &p_gram_b;
        st_gram_map = &mg_data.map_b;
    }

    //we will move the compressed string to the back
    if(st_gram!=nullptr){
        st_gram->rules.resize(max_lvl+1);
        st_gram->metadata.resize(max_lvl+2);
        st_gram->rules[min_lvl].swap(st_gram->rules[max_lvl]);
        st_gram->lvl = st_gram->rules.size();
        //+1 because the first element is the terminals' metadata
        std::swap(st_gram->metadata[min_lvl+1], st_gram->metadata[max_lvl+1]);
    }

    lvl_metadata_type buffer_metadata;
    size_t i=0;
    while(i<max_lvl) {
        buffer_metadata = merge_level(p_gram_a.rules[i], p_gram_a.metadata[i+1],
                                      p_gram_b.rules[i], p_gram_b.metadata[i+1],
                                      fp_seeds[i+1], mg_data);
        i++;
        if(i>=min_lvl && i<max_lvl){
            assert(st_gram!=nullptr && st_gram!= nullptr);
            create_fake_level(*st_gram, i, mg_data.fps, fp_seeds[i+1], *st_gram_map);
        }

        p_gram_a.metadata[i] = buffer_metadata;
        mg_data.buffer.copy(p_gram_a.metadata[i].n_bits(), p_gram_a.rules[i-1]);
    }
    p_gram_a.metadata[max_lvl+1] = concatenate_strings(p_gram_a.rules[max_lvl], p_gram_a.metadata[max_lvl+1],
                                                       p_gram_b.rules[max_lvl], p_gram_b.metadata[max_lvl+1],
                                                       mg_data);
    //p_gram_a.rules[max_lvl].swap(mg_data.buffer);
    mg_data.buffer.copy(p_gram_a.metadata[max_lvl+1].n_bits(), p_gram_a.rules[max_lvl]);

    p_gram_a.lvl = p_gram_a.metadata.size()-1;
    p_gram_a.text_size += p_gram_b.text_size;
    p_gram_a.longest_rule = longest_rule;
}

template<class p_gram_type>
void merge_two_partial_grammars_se(std::string& p_gram_a_file, std::string& p_gram_b_file,
                                   std::string& p_gram_c_file) {

    p_gram_type p_gram_a, p_gram_b;
    std::ifstream ifs_a(p_gram_a_file, std::ios::binary);
    p_gram_a.load_metadata(ifs_a);

    std::ifstream ifs_b(p_gram_b_file, std::ios::binary);
    p_gram_b.load_metadata(ifs_b);
    assert(p_gram_a.par_seed == p_gram_b.par_seed);

    size_t longest_rule = std::max(p_gram_a.longest_rule, p_gram_b.longest_rule);
    merge_data_t mg_data;
    //in the first level, tot_symbols represents the alphabet of terminals
    mg_data.initialize(p_gram_a.metadata[0].tot_symbols, sizeof(typename p_gram_type::sym_type),  longest_rule);

    //we subtract one because the last level contains the compressed strings,
    // and we do not merge but concatenate
    size_t max_lvl = std::max(p_gram_a.lvl, p_gram_b.lvl)-1;
    size_t min_lvl = std::min(p_gram_a.lvl, p_gram_b.lvl)-1;

    bitstream<size_t> rules_buffer_a;
    bitstream<size_t> rules_buffer_b;
    lvl_metadata_type buffer_metadata;
    p_gram_a.rules.resize(max_lvl+1);

    size_t i=0;
    while(i<min_lvl) {
        p_gram_a.load_next_rule_set(ifs_a, i, rules_buffer_a);
        p_gram_b.load_next_rule_set(ifs_b, i, rules_buffer_b);

        buffer_metadata = merge_level(rules_buffer_a, p_gram_a.metadata[i+1],
                                      rules_buffer_b, p_gram_b.metadata[i+1], mg_data);
        i++;
        p_gram_a.metadata[i] = buffer_metadata;
        mg_data.buffer.copy(p_gram_a.metadata[i].n_bits(), p_gram_a.rules[i-1]);

        //todo testing space
        size_t tmp_bytes = rules_buffer_a.capacity_in_bytes();
        tmp_bytes+=rules_buffer_b.capacity_in_bytes();
        tmp_bytes+=mg_data.map_a.size()*sizeof(uint32_t);
        tmp_bytes+=mg_data.map_b.size()*sizeof(uint32_t);
        tmp_bytes+=mg_data.buffer.capacity_in_bytes();
        for(size_t u=0;u<i;u++){
            tmp_bytes+=p_gram_a.rules[u].capacity_in_bytes();
        }
        tmp_bytes+=mg_data.fps.size()*sizeof(uint64_t);
        //
        //std::cout<<report_space((off_t)tmp_bytes)<<" "<<malloc_count_peak()<<std::endl;
    }

    //NOTE: now p_gram_a.rules[min_lvl] or p_gram_b.rules[min_lvl] contains the concatenated strings, and we do not merge them

    if(min_lvl==max_lvl){

        //grammars have the same height, not need to deal with different heights.
        //just read load the concatenated strings from disk, update their symbols and concatenate them
        p_gram_a.load_next_rule_set(ifs_a, i, rules_buffer_a);
        p_gram_b.load_next_rule_set(ifs_b, i, rules_buffer_b);

        p_gram_a.metadata[max_lvl+1] = concatenate_strings(rules_buffer_a, p_gram_a.metadata[max_lvl+1],
                                                           rules_buffer_b, p_gram_b.metadata[max_lvl+1],
                                                           mg_data);
        mg_data.buffer.copy(p_gram_a.metadata[max_lvl+1].n_bits(), p_gram_a.rules[max_lvl]);
    } else {
        //grammars have different heights, we need to deal with this

        //pointers for the shortest grammar
        p_gram_type *st_gram;
        std::vector<uint32_t> *st_gram_map;
        bitstream<size_t> *st_rules_buffer;
        std::ifstream * st_ifs;

        //pointers for the tallest grammar
        p_gram_type *ht_gram;
        bitstream<size_t> *ht_rules_buffer;
        std::ifstream * ht_ifs;

        if(p_gram_a.lvl<p_gram_b.lvl){
            st_gram = &p_gram_a;
            st_gram_map = &mg_data.map_a;
            st_rules_buffer = &rules_buffer_a;
            st_ifs = &ifs_a;

            ht_gram = &p_gram_b;
            ht_rules_buffer = &rules_buffer_b;
            ht_ifs = &ifs_b;
        } else {
            st_gram = &p_gram_b;
            st_gram_map = &mg_data.map_b;
            st_rules_buffer = &rules_buffer_b;
            st_ifs = &ifs_b;

            ht_gram = &p_gram_a;
            ht_rules_buffer = &rules_buffer_a;
            ht_ifs = &ifs_a;
        }

        //we will move the concatenated strings to the last level, because we need to update them in each round
        st_gram->rules.resize(max_lvl+1);//+1 because max_lvl does not consider the concatenated strings
        st_gram->metadata.resize(max_lvl+2);//+2 because max_lvl does not consider the concat. strings and the terminals metadata
        st_gram->load_next_rule_set(*st_ifs, i, st_gram->rules[max_lvl]);
        //+1 because the first element is the terminals' metadata
        std::swap(st_gram->metadata[min_lvl+1], st_gram->metadata[max_lvl+1]);

        while(i<max_lvl){
            create_fake_level_se(*st_gram, *st_rules_buffer, i, mg_data.fps, *st_gram_map);
            ht_gram->load_next_rule_set(*ht_ifs, i, *ht_rules_buffer);

            buffer_metadata = merge_level(rules_buffer_a, p_gram_a.metadata[i+1],
                                          rules_buffer_b, p_gram_b.metadata[i+1], mg_data);
            i++;
            p_gram_a.metadata[i] = buffer_metadata;
            mg_data.buffer.copy(p_gram_a.metadata[i].n_bits(), p_gram_a.rules[i-1]);
        }

        ht_gram->load_next_rule_set(*ht_ifs, i, *ht_rules_buffer);

        if(st_gram==&p_gram_b){
            p_gram_a.metadata[max_lvl+1] = concatenate_strings(rules_buffer_a, p_gram_a.metadata[max_lvl+1],
                                                               p_gram_b.rules[max_lvl], p_gram_b.metadata[max_lvl+1],
                                                               mg_data);
        }else{
            p_gram_a.metadata[max_lvl+1] = concatenate_strings(p_gram_a.rules[max_lvl], p_gram_a.metadata[max_lvl+1],
                                                               rules_buffer_b, p_gram_b.metadata[max_lvl+1],
                                                               mg_data);
        }

        mg_data.buffer.copy(p_gram_a.metadata[max_lvl+1].n_bits(), p_gram_a.rules[max_lvl]);
    }

    p_gram_a.lvl = p_gram_a.metadata.size()-1;
    p_gram_a.text_size += p_gram_b.text_size;
    p_gram_a.longest_rule = longest_rule;
    store_to_file(p_gram_c_file, p_gram_a);

    //destroy buffers
    rules_buffer_a.destroy();
    rules_buffer_b.destroy();
#ifdef __linux__
    malloc_trim(0);
#endif
}

void merge_gramms(std::vector<std::string>& grams_to_merge, std::string& merged_grammar){
    std::vector<uint64_t> par_seeds;
    uint64_t orig_seed = get_par_seed_par_gram(grams_to_merge[0]);

    std::mt19937 gen(orig_seed); //Standard mersenne_twister_engine seeded with a fixed value
    std::uniform_int_distribution<uint64_t> distrib(1, std::numeric_limits<uint64_t>::max());
    par_seeds.resize(32);
    for(size_t i=0;i<32;i++){
        par_seeds[i] = distrib(gen);
    }

    std::cout << "Merged grammars: " <<0<<"/"<<grams_to_merge.size()<<std::flush;
    merge_two_partial_grammars_se<partial_gram<uint8_t>>(grams_to_merge[0], grams_to_merge[1], merged_grammar);
    std::cout << "\r" << "Merged grammars: " <<2<<"/"<<grams_to_merge.size()<<std::flush;
    for(size_t i=2;i<grams_to_merge.size();i++){
        merge_two_partial_grammars_se<partial_gram<uint8_t>>(merged_grammar, grams_to_merge[i], merged_grammar);
        std::cout << "\r" << "Merged grammars: " <<(i+1)<<"/"<<grams_to_merge.size()<<std::flush;
    }
    std::cout<<""<<std::endl;
    get_breakdown(merged_grammar);
    std::cout<<"The resulting grammar was stored in "<<merged_grammar<<std::endl;
}

template<class stream_type>
lvl_metadata_type concatenate_strings(stream_type &stream_a, lvl_metadata_type &lvl_met_a,
                                      stream_type &stream_b, lvl_metadata_type &lvl_met_b,
                                      merge_data_t& mg_data){

    lvl_metadata_type c_string_lvl{};
    c_string_lvl.sym_width = sym_width(mg_data.lvl_sigma)+1;//+1 is to mark the end of each phrase in the stream of rules
    c_string_lvl.tot_symbols = lvl_met_a.tot_symbols + lvl_met_b.tot_symbols;
    c_string_lvl.n_rules = 1;
    c_string_lvl.terminals = false;

    mg_data.buffer.reserve_in_bits(c_string_lvl.n_bits());

    size_t d_pos=0;
    size_t d_width = c_string_lvl.sym_width;

    size_t s_pos =0;
    size_t s_width = lvl_met_a.sym_width;
    size_t n_bits = lvl_met_a.n_bits();
    size_t mt_sym=0;
    while(s_pos<n_bits){
        mt_sym = stream_a.read(s_pos, s_pos+s_width-1);
        mt_sym>>=1UL;
        mt_sym = mg_data.map_a[mt_sym];
        mg_data.buffer.write(d_pos, d_pos+d_width-1, mt_sym<<1UL);
        s_pos+=s_width;
        d_pos+=d_width;
    }

    s_pos =0;
    s_width = lvl_met_b.sym_width;
    n_bits = lvl_met_b.n_bits();
    while(s_pos<n_bits){
        mt_sym = stream_b.read(s_pos, s_pos+s_width-1);
        mt_sym>>=1UL;
        mt_sym = mg_data.map_b[mt_sym];
        mg_data.buffer.write(d_pos, d_pos+d_width-1, mt_sym<<1UL);
        s_pos+=s_width;
        d_pos+=d_width;
    }

    //return one position back to mark the last symbol with a bit to indicate it is the end of the stream
    d_pos-=d_width;
    mg_data.buffer.write(d_pos, d_pos+d_width-1, ((mt_sym<<1UL)| 1UL) );
    d_pos+=d_width;
    assert(d_pos==c_string_lvl.n_bits());
    return c_string_lvl;
}

void print_merge_stats(std::vector<lvl_metadata_type>& met_a, std::vector<lvl_metadata_type>& met_b,
                       std::vector<lvl_metadata_type>& met_c){

    for(size_t i=1;i<met_c.size()-1;i++){
        std::cout<<"Level "<<i<<std::endl;
        size_t n_rules_a = i<met_a.size() ? met_a[i].n_rules : 0;
        size_t n_rules_b = i<met_a.size() ? met_b[i].n_rules : 0;
        size_t n_rules_c = met_c[i].n_rules;
        size_t n_sym_a = i<met_a.size() ? met_a[i].tot_symbols : 0;
        size_t n_sym_b = i<met_a.size() ? met_b[i].tot_symbols : 0;
        size_t n_sym_c = met_c[i].tot_symbols;
        std::cout<<"  Number of rules:   A:"<<n_rules_a<<", B:"<<n_rules_b<<" -> C:"<<n_rules_c<<std::endl;
        std::cout<<"  Number of symbols: A:"<<n_sym_a<<", B:"<<n_sym_b<<" -> C:"<<n_sym_c<<std::endl;
    }
    std::cout<<"Compressed string"<<std::endl;
    std::cout<<"  Number of strings: A:"<<met_a.back().tot_symbols<<", B:"<<met_b.back().tot_symbols<<" -> C:"<<met_c.back().tot_symbols<<std::endl;
}

template<class stream_type>
lvl_metadata_type merge_level(stream_type &stream_a, lvl_metadata_type &lvl_met_a,
                              stream_type &stream_b, lvl_metadata_type &lvl_met_b,
                              merge_data_t& mg_data) {
    uint64_t fp_a, fp_b;
    size_t curr_rule_a=0, curr_rule_b=0;
    size_t curr_pos_a=0, curr_pos_b=0;
    size_t len_a, len_b;

    lvl_metadata_type lvl_met_c{};
    lvl_met_c.sym_width = sym_width(mg_data.lvl_sigma)+1;//we use the extra bit to mark the end of each phrase in the stream of rules
    lvl_met_c.terminals = lvl_met_a.terminals;

    uint8_t m_width = lvl_met_c.sym_width;
    size_t curr_pos_m=0;

    mg_data.buffer.reserve_in_bits(std::max(lvl_met_a.n_bits(), lvl_met_b.n_bits()));

    std::vector<uint8_t> merge_marks;
    merge_marks.reserve(lvl_met_a.n_rules+lvl_met_b.n_rules+1);

    std::vector<uint64_t> new_fps;
    new_fps.reserve(lvl_met_a.n_rules+lvl_met_b.n_rules+1);
    new_fps.push_back(0);//fake

    uint8_t a_width = lvl_met_a.sym_width;
    uint8_t b_width = lvl_met_b.sym_width;

    std::vector<uint64_t> fp_seq(mg_data.longest_rule, 0);
    std::vector<uint32_t> phrase_a(mg_data.longest_rule, 0);
    std::vector<uint32_t> phrase_b(mg_data.longest_rule, 0);

    get_rule_info(stream_a, curr_pos_a, a_width, mg_data.fps, mg_data.map_a,  fp_seq,
                  phrase_a, fp_a, len_a);
    assert(len_a<=mg_data.longest_rule);
    get_rule_info(stream_b, curr_pos_b, b_width, mg_data.fps, mg_data.map_b, fp_seq,
                  phrase_b, fp_b, len_b);
    assert(len_b<=mg_data.longest_rule);

    uint64_t prev_fp_a = fp_a;
    uint64_t prev_fp_b = fp_b;

    while(curr_rule_a<lvl_met_a.n_rules && curr_rule_b<lvl_met_b.n_rules){

        if(fp_a<fp_b){
            //write rule from A
            append_rule(phrase_a, len_a, curr_pos_m, m_width, mg_data.buffer);
            curr_rule_a++;

            new_fps.push_back(fp_a);
            lvl_met_c.tot_symbols+=len_a;
            merge_marks.push_back(1);
        } else if(fp_b<fp_a) {
            //write rule from B
            append_rule(phrase_b, len_b, curr_pos_m, m_width, mg_data.buffer);
            curr_rule_b++;

            new_fps.push_back(fp_b);
            lvl_met_c.tot_symbols+=len_b;
            merge_marks.push_back(2);
        } else {

            bool eq = len_a==len_b && (memcmp(phrase_a.data(), phrase_b.data(), len_a*sizeof(uint32_t))==0);

            if(eq){
                //write rule from A
                append_rule(phrase_a, len_a, curr_pos_m, m_width, mg_data.buffer);
                curr_rule_a++;
                curr_rule_b++;

                new_fps.push_back(fp_a);
                lvl_met_c.tot_symbols+=len_a;
                merge_marks.push_back(3);
            }else {

                //a collision occurred (extremely unlikely, but not impossible)
                std::cout << "Collision warning:  " << fp_a << " " << fp_b << " " << curr_pos_a << " " << curr_pos_b<< std::endl;
                // break ties by lex. rank
                bool a_is_lex_smaller = rules_lex_comp(phrase_a, len_a, phrase_b, len_b);

                if(a_is_lex_smaller) {
                    //write rule from A
                    append_rule(phrase_a, len_a, curr_pos_m, m_width, mg_data.buffer);
                    curr_rule_a++;

                    new_fps.push_back(fp_a);
                    lvl_met_c.tot_symbols += len_a;
                    merge_marks.push_back(1);
                } else {
                    //write rule from B
                    append_rule(phrase_b, len_b, curr_pos_m, m_width, mg_data.buffer);
                    curr_rule_b++;

                    new_fps.push_back(fp_b);
                    lvl_met_c.tot_symbols += len_b;
                    merge_marks.push_back(2);
                }
            }
            //
        }

        if((merge_marks.back() & 1) && curr_rule_a<lvl_met_a.n_rules){
            get_rule_info(stream_a, curr_pos_a, a_width, mg_data.fps, mg_data.map_a, fp_seq,
                          phrase_a, fp_a, len_a);
            assert(fp_a>=prev_fp_a);
            prev_fp_a = fp_a;
        }

        if((merge_marks.back() & 2) && curr_rule_b<lvl_met_b.n_rules){
            get_rule_info(stream_b, curr_pos_b, b_width, mg_data.fps, mg_data.map_b, fp_seq,
                          phrase_b, fp_b, len_b);
            assert(fp_b>=prev_fp_b);
            prev_fp_b = fp_b;
        }
    }

    while(curr_rule_a<lvl_met_a.n_rules){
        //write rule from A
        append_rule(phrase_a, len_a, curr_pos_m, m_width, mg_data.buffer);
        curr_rule_a++;

        new_fps.push_back(fp_a);
        lvl_met_c.tot_symbols+=len_a;
        merge_marks.push_back(1);

        if(curr_rule_a<lvl_met_a.n_rules){
            get_rule_info(stream_a, curr_pos_a, a_width, mg_data.fps, mg_data.map_a, fp_seq, phrase_a, fp_a, len_a);
            assert(fp_a>=prev_fp_a);
            prev_fp_a = fp_a;
        }
    }

    while(curr_rule_b<lvl_met_b.n_rules){
        //write rule from B
        append_rule(phrase_b, len_b, curr_pos_m, m_width, mg_data.buffer);
        curr_rule_b++;

        new_fps.push_back(fp_b);
        lvl_met_c.tot_symbols+=len_b;
        merge_marks.push_back(2);

        if(curr_rule_b<lvl_met_b.n_rules){
            get_rule_info(stream_b, curr_pos_b, b_width, mg_data.fps, mg_data.map_b, fp_seq,
                          phrase_b, fp_b, len_b);
            assert(fp_b>=prev_fp_b);
            prev_fp_b = fp_b;
        }
    }

    mg_data.fps.swap(new_fps);
    destroy(new_fps);

    //update mapping values
    //the new metasymbols are one-based to differentiate them from the separator symbol (0) is the parse
    size_t mt_sym_a=1, mt_sym_b=1, mg_mt_sym=1;
    mg_data.map_a.resize(lvl_met_a.n_rules+1);
    mg_data.map_b.resize(lvl_met_b.n_rules+1);
    mg_data.map_a[0] = 0;
    mg_data.map_b[0] = 0;
    for(unsigned char merge_mark : merge_marks){
        if(merge_mark==1){
            mg_data.map_a[mt_sym_a++] = mg_mt_sym;
        } else if(merge_mark==2){
            mg_data.map_b[mt_sym_b++] = mg_mt_sym;
        } else {
            mg_data.map_a[mt_sym_a++] = mg_mt_sym;
            mg_data.map_b[mt_sym_b++] = mg_mt_sym;
        }
        mg_mt_sym++;
    }

    lvl_met_c.n_rules = merge_marks.size();
    mg_data.lvl_sigma = lvl_met_c.n_rules;//store this information for the next round

    mg_data.map_a.shrink_to_fit();
    mg_data.map_b.shrink_to_fit();
    destroy(merge_marks);

#ifdef __linux__
    malloc_trim(0);
#endif

    return lvl_met_c;
}

template<class gram_type>
void create_fake_level(gram_type& p_gram, size_t new_lvl, std::vector<uint64_t>& prev_fps, std::vector<uint32_t>& mt_map){

    // new_level is to the previous level in the metadata because
    // the first element of the metadata vector has the terminal alphabet
    p_gram.metadata[new_lvl+1].sym_width = sym_width(p_gram.metadata[new_lvl].n_rules)+1;
    p_gram.metadata[new_lvl+1].n_rules = p_gram.metadata[new_lvl].n_rules;
    p_gram.metadata[new_lvl+1].tot_symbols = p_gram.metadata[new_lvl].n_rules;
    p_gram.metadata[new_lvl+1].terminals = false;

    assert((p_gram.metadata[new_lvl+1].n_rules+1) == mt_map.size());

    p_gram.rules[new_lvl].reserve_in_bits(p_gram.metadata[new_lvl+1].n_bits());

    std::vector<std::tuple<uint64_t, uint64_t, uint32_t>> perm(mt_map.size());
    perm[0] = {0, 0, 0};
    for(size_t i=1;i<mt_map.size();i++){
        std::get<0>(perm[i]) = i;
        uint64_t fp = prev_fps[mt_map[i]];
        std::get<1>(perm[i]) = XXH3_64bits(&fp, sizeof(uint64_t));
        std::get<2>(perm[i]) = mt_map[i];
    }

    std::sort(perm.begin(), perm.end(), [&](auto const& a, auto const &b) -> bool{
        if(std::get<1>(a)!=std::get<1>(b)){
            return std::get<1>(a) < std::get<1>(b);//break ties using the level fingerprint
        }
        assert(std::get<0>(a)==std::get<0>(b) ||
               prev_fps[std::get<2>(a)]!=prev_fps[std::get<2>(b)]);
        return prev_fps[std::get<2>(a)]<prev_fps[std::get<2>(b)];
    });

    size_t pos=0, width=p_gram.metadata[new_lvl+1].sym_width;
    std::vector<uint32_t> inv_perm(perm.size());
    for(size_t mt_sym=1;mt_sym<perm.size();mt_sym++){
        p_gram.rules[new_lvl].write(pos, pos+width-1, (std::get<0>(perm[mt_sym])<<1UL | 1));
        mt_map[std::get<0>(perm[mt_sym])] = std::get<2>(perm[mt_sym]);
        pos+=width;
        inv_perm[std::get<0>(perm[mt_sym])] = mt_sym;
    }
    assert(pos==p_gram.metadata[new_lvl+1].n_bits());

    //update the compressed string
    size_t last_lvl = p_gram.rules.size()-1;
    pos = 0;
    width = p_gram.metadata[last_lvl+1].sym_width;
    size_t n_bits = p_gram.metadata[last_lvl+1].n_bits();
    size_t mt_sym;
    bool last_sym;
    while(pos<n_bits){
        mt_sym = p_gram.rules[last_lvl].read(pos, pos+width-1);
        last_sym = mt_sym & 1UL;
        mt_sym>>=1UL;
        mt_sym = inv_perm[mt_sym];
        p_gram.rules[last_lvl].write(pos, pos+width-1, (mt_sym<<1UL | last_sym));
        pos+=width;
    }
    assert(pos==p_gram.metadata[last_lvl+1].n_bits());
}

template<class gram_type>
void create_fake_level_se(gram_type& p_gram, bitstream<size_t>& buffer, size_t new_lvl, std::vector<uint64_t>& prev_fps,
                          std::vector<uint32_t>& mt_map){

    // new_level is new_lvl+1 in the metadata because
    // the first element of the metadata vector has the terminal alphabet
    p_gram.metadata[new_lvl+1].sym_width = sym_width(p_gram.metadata[new_lvl].n_rules)+1;
    p_gram.metadata[new_lvl+1].n_rules = p_gram.metadata[new_lvl].n_rules;
    p_gram.metadata[new_lvl+1].tot_symbols = p_gram.metadata[new_lvl].n_rules;
    p_gram.metadata[new_lvl+1].terminals = false;

    assert((p_gram.metadata[new_lvl+1].n_rules+1) == mt_map.size());
    buffer.reserve_in_bits(p_gram.metadata[new_lvl+1].n_bits());

    std::vector<std::tuple<uint64_t, uint64_t, uint32_t>> perm(mt_map.size());
    perm[0] = {0, 0, 0};
    for(size_t i=1;i<mt_map.size();i++){
        std::get<0>(perm[i]) = i;
        uint64_t fp = prev_fps[mt_map[i]];
        std::get<1>(perm[i]) = XXH3_64bits(&fp, sizeof(uint64_t));
        std::get<2>(perm[i]) = mt_map[i];
    }

    std::sort(perm.begin(), perm.end(), [&](auto const& a, auto const &b) -> bool{
        if(std::get<1>(a)!=std::get<1>(b)){
            return std::get<1>(a) < std::get<1>(b);//break ties using the level fingerprint
        }
        assert(std::get<0>(a)==std::get<0>(b) ||
               prev_fps[std::get<2>(a)]!=prev_fps[std::get<2>(b)]);
        return prev_fps[std::get<2>(a)]<prev_fps[std::get<2>(b)];
    });

    size_t pos=0, width=p_gram.metadata[new_lvl+1].sym_width;
    std::vector<uint32_t> inv_perm(perm.size());
    for(size_t mt_sym=1;mt_sym<perm.size();mt_sym++){
        buffer.write(pos, pos+width-1, (std::get<0>(perm[mt_sym])<<1UL | 1));
        mt_map[std::get<0>(perm[mt_sym])] = std::get<2>(perm[mt_sym]);
        pos+=width;
        inv_perm[std::get<0>(perm[mt_sym])] = mt_sym;
    }
    assert(pos==p_gram.metadata[new_lvl+1].n_bits());

    //update the compressed strings
    size_t last_lvl = p_gram.rules.size()-1;
    pos = 0;
    width = p_gram.metadata[last_lvl+1].sym_width;
    size_t n_bits = p_gram.metadata[last_lvl+1].n_bits();
    size_t mt_sym;
    bool last_sym;
    while(pos<n_bits){
        mt_sym = p_gram.rules[last_lvl].read(pos, pos+width-1);
        last_sym = mt_sym & 1UL;
        mt_sym>>=1UL;
        mt_sym = inv_perm[mt_sym];
        p_gram.rules[last_lvl].write(pos, pos+width-1, (mt_sym<<1UL | last_sym));
        pos+=width;
    }
    assert(pos==p_gram.metadata[last_lvl+1].n_bits());
}
#endif //LCG_MERGE_GRAMS_H
