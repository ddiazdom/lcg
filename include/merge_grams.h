//
// Created by Diaz, Diego on 22.8.2024.
//

#ifndef LCG_MERGE_GRAMS_H
#define LCG_MERGE_GRAMS_H

#include "partial_gram.h"

template<class p_gram_type>
void merge_two_partial_grammars_in_memory(p_gram_type& p_gram_a, p_gram_type& p_gram_b, std::vector<uint64_t>& fp_seeds) {

    size_t longest_rule = std::max(p_gram_a.longest_rule, p_gram_b.longest_rule);

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
void merge_two_partial_grammars_semi_external(std::string& p_gram_a_file, std::string& p_gram_b_file, std::string& p_gram_c_file,
                                              std::vector<uint64_t>& fp_seeds) {

    p_gram_type p_gram_a, p_gram_b;
    std::ifstream ifs_a(p_gram_a_file, std::ios::binary);
    p_gram_a.load_metadata(ifs_a);

    std::ifstream ifs_b(p_gram_b_file, std::ios::binary);
    p_gram_b.load_metadata(ifs_b);

    size_t longest_rule = std::max(p_gram_a.longest_rule, p_gram_b.longest_rule);

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
    bitstream<size_t> rules_buffer_a;
    bitstream<size_t> rules_buffer_b;
    p_gram_a.rules.resize(max_lvl+1);

    size_t i=0;
    while(i<max_lvl) {

        p_gram_a.load_next_rule_set(ifs_a, i, rules_buffer_a);
        p_gram_b.load_next_rule_set(ifs_b, i, rules_buffer_b);

        buffer_metadata = merge_level(rules_buffer_a, p_gram_a.metadata[i+1],
                                      rules_buffer_b, p_gram_b.metadata[i+1],
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
    store_to_file(p_gram_c_file, p_gram_a);
}

template<class stream_type>
lvl_metadata_type concatenate_strings(stream_type &stream_a, lvl_metadata_type &lvl_met_a,
                                      stream_type &stream_b, lvl_metadata_type &lvl_met_b,
                                      merge_data_t& mg_data){

    //std::cout<<"Buffer: "<<report_space(mg_data.buffer.capacity_in_bytes())<<" A:"<<report_space(stream_a.capacity_in_bytes())<<" B:"<<report_space(stream_b.capacity_in_bytes())<<std::endl;

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
                              uint64_t &fp_seed, merge_data_t& mg_data) {
    uint64_t fp_a, fp_b;
    size_t curr_rule_a=0, curr_rule_b=0;
    size_t curr_pos_a=0, curr_pos_b=0;
    size_t len_a, len_b;

    lvl_metadata_type lvl_met_c{};
    lvl_met_c.sym_width = sym_width(mg_data.lvl_sigma)+1;//we use the extra bit to mark the end of each phrase in the stream of rules
    lvl_met_c.terminals = lvl_met_a.terminals;

    uint8_t m_width = lvl_met_c.sym_width;
    size_t curr_pos_m=0;

    //std::cout<<"Buffer: "<<report_space(mg_data.buffer.capacity_in_bytes())<<" A:"<<report_space(stream_a.capacity_in_bytes())<<" B:"<<report_space(stream_b.capacity_in_bytes())<<std::endl;

    mg_data.buffer.reserve_in_bits(std::max(lvl_met_a.n_bits(), lvl_met_b.n_bits()));

    std::vector<uint8_t> merge_marks;
    merge_marks.reserve(lvl_met_a.n_rules+lvl_met_b.n_rules+1);

    std::vector<uint64_t> new_fps;
    new_fps.reserve(lvl_met_a.n_rules+lvl_met_b.n_rules+1);
    new_fps.push_back(0);//fake

    uint8_t a_width = lvl_met_a.sym_width;
    uint8_t b_width = lvl_met_b.sym_width;

    std::vector<uint64_t> phrase_fp_buff(mg_data.longest_rule, 0);
    std::vector<uint64_t> phrase_a(mg_data.longest_rule, 0);
    std::vector<uint64_t> phrase_b(mg_data.longest_rule, 0);

    get_rule_info(stream_a, curr_pos_a, a_width, mg_data.fps, mg_data.map_a, fp_seed, phrase_fp_buff,
                  phrase_a, fp_a, len_a);
    assert(len_a<=mg_data.longest_rule);
    get_rule_info(stream_b, curr_pos_b, b_width, mg_data.fps, mg_data.map_b, fp_seed, phrase_fp_buff,
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

            bool eq = len_a==len_b && (memcmp(phrase_a.data(), phrase_b.data(), len_a*sizeof(uint64_t))==0);

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
            get_rule_info(stream_a, curr_pos_a, a_width, mg_data.fps, mg_data.map_a, fp_seed, phrase_fp_buff,
                          phrase_a, fp_a, len_a);
            assert(fp_a>=prev_fp_a);
            prev_fp_a = fp_a;
        }

        if((merge_marks.back() & 2) && curr_rule_b<lvl_met_b.n_rules){
            get_rule_info(stream_b, curr_pos_b, b_width, mg_data.fps, mg_data.map_b, fp_seed, phrase_fp_buff,
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
            get_rule_info(stream_a, curr_pos_a, a_width, mg_data.fps, mg_data.map_a, fp_seed,
                          phrase_fp_buff, phrase_a, fp_a, len_a);
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
            get_rule_info(stream_b, curr_pos_b, b_width, mg_data.fps, mg_data.map_b, fp_seed, phrase_fp_buff,
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
void create_fake_level(gram_type& p_gram, size_t new_lvl, std::vector<uint64_t>& prev_fps,
                       uint64_t fp_seed, std::vector<uint32_t>& mt_map){

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
        std::get<1>(perm[i]) = XXH64(&fp, sizeof(uint64_t), fp_seed);
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
#endif //LCG_MERGE_GRAMS_H
