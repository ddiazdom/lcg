//
// Created by Diaz, Diego on 23.10.2023.
//

#ifndef LCG_PARTIAL_GRAM_H
#define LCG_PARTIAL_GRAM_H
#include "cds/bitstream.h"
#include "lz_map.h"

struct lvl_metadata_type{
    size_t n_rules;
    size_t tot_symbols;
    uint8_t sym_width;
    bool terminals;
};

template<class ter_type>
struct partial_gram {

    bitstream<size_t> rules;
    std::vector<lvl_metadata_type> rules_metadata;
    size_t curr_lvl = 0;
    size_t acc_bits = 0;

    typedef string_map_iterator<partial_gram<ter_type>> iter_type;

    partial_gram(){
        lvl_metadata_type ter{};
        ter.sym_width = std::numeric_limits<ter_type>::digits;
        ter.tot_symbols = std::numeric_limits<ter_type>::max()+1;
        ter.n_rules = ter.tot_symbols;
        ter.terminals = true;
        rules_metadata.push_back(ter);
    }


    off_t serialize(std::ofstream &ofs){
        off_t written_bytes = rules.serialize(ofs);
        written_bytes+= serialize_plain_vector(ofs, rules_metadata);
        return written_bytes;
    }

    size_t serialize_to_fd(int fd){
        size_t g_words = INT_CEIL(acc_bits, (sizeof(size_t)*8));
        write(fd, &g_words, sizeof(size_t));
        write(fd, rules.stream, g_words*sizeof(size_t));
        size_t written_bytes = g_words*sizeof(size_t) + sizeof(size_t);

        size_t len = rules_metadata.size();
        write(fd, &len, sizeof(size_t));
        written_bytes += sizeof(size_t);
        write(fd, rules_metadata.data(), rules_metadata.size()*sizeof(lvl_metadata_type));
        written_bytes += rules_metadata.size()*sizeof(lvl_metadata_type);

        return written_bytes;
    }

    /*void load_metadata(std::ifstream &ifs){
        load_plain_vector(ifs, rules_metadata);
    }
    void load_next_gram_level(std::ifstream  &ifs, lvl_metadata_type& lvl_met, bitstream<size_t>& lvl_rules){
        assert(curr_lvl<rules_metadata.size());
        off_t space = rules_metadata[curr_lvl].sym_width * rules_metadata[curr_lvl].tot_symbols; //number of bits
        space = INT_CEIL(space, 8); //number of bytes
        rules.stream_size = INT_CEIL(space, sizeof(size_t)); //number of size_t words
        rules.stream = (size_t *)realloc(rules.stream, rules.stream_size*sizeof(size_t));
        ifs.read((char *)rules.stream, rules.stream_size*sizeof(size_t));
    }*/

    template<class sym_type>
    void append_new_lvl(sym_type* text, const lz_map::phrase_list_t &phrase_set, size_t tot_symbols,
                        std::vector<std::pair<uint32_t, uint64_t>>& perm){

        lvl_metadata_type lvl_met{};
        lvl_met.n_rules = phrase_set.size();

        //the extra bit is to mark the end of each rule
        lvl_met.sym_width = sym_width(rules_metadata.back().n_rules)+1;
        lvl_met.tot_symbols = tot_symbols;
        lvl_met.terminals = false;

        size_t new_bits = (lvl_met.sym_width * lvl_met.tot_symbols);
        off_t tot_words = INT_CEIL((acc_bits + new_bits), (sizeof(size_t)*8));

        if(rules.stream_size==0){
            rules.stream = (size_t *)malloc(tot_words*sizeof(size_t));
            rules.stream_size = tot_words;
        }else if(rules.stream_size<tot_words){
            rules.stream = (size_t *)realloc(rules.stream, tot_words*sizeof(size_t));
            rules.stream_size = tot_words;
        }

        for(size_t i=0;i<phrase_set.size();i++){

            size_t idx = perm[i].first;
            uint32_t source = phrase_set[idx].source/sizeof(sym_type);
            uint32_t len = (phrase_set[idx].len/sizeof(sym_type));
            uint32_t last = source + len-1;

            for(size_t j=source;j<last;j++){
                rules.write(acc_bits, acc_bits+lvl_met.sym_width-1, text[j]<<1UL);
                acc_bits+=lvl_met.sym_width;
            }
            rules.write(acc_bits, acc_bits+lvl_met.sym_width-1, ((text[last]<<1UL) | 1UL));
            acc_bits+=lvl_met.sym_width;
        }
        rules_metadata.push_back(lvl_met);
    }


    template<class sym_type>
    void add_compressed_string(sym_type* text, size_t size){

        lvl_metadata_type lvl_met{};
        lvl_met.n_rules = 1;
        //the extra bit is to mark the end of each rule
        lvl_met.sym_width = sym_width(rules_metadata.back().n_rules)+1;
        lvl_met.tot_symbols = size/2;
        lvl_met.terminals = false;

        size_t new_bits = (lvl_met.sym_width * lvl_met.tot_symbols);
        off_t tot_words = INT_CEIL((acc_bits + new_bits), (sizeof(size_t)*8));

        if(rules.stream_size<tot_words){
            assert(rules.stream_size>0);
            rules.stream = (size_t *)realloc(rules.stream, tot_words*sizeof(size_t));
            rules.stream_size = tot_words;
        }
        //std::cout<<"  "<<(int)lvl_met.sym_width<<" "<<lvl_met.tot_symbols<<" | "<<acc_bits<<" "<<tot_words<<" "<<new_bits<<std::endl;

        //we skip the separator symbols
        size_t last = size-2;
        for(size_t j=0;j<last;j+=2){
            rules.write(acc_bits, acc_bits+lvl_met.sym_width-1, text[j]<<1UL);
            acc_bits+=lvl_met.sym_width;
        }
        rules.write(acc_bits, acc_bits+lvl_met.sym_width-1, ((text[last]<<1UL) | 1UL));
        acc_bits+=lvl_met.sym_width;

        rules_metadata.push_back(lvl_met);
    }

    void reset_grammar(){
        acc_bits = 0;
        memset(rules.stream, 0,  rules.stream_size*sizeof(size_t));
        rules_metadata.resize(1);
    }
};

std::pair<uint64_t,size_t> get_rule_info(bitstream<size_t>& rule_stream, size_t pos, size_t width,
                                         std::vector<uint64_t>& prev_fps, std::vector<uint32_t>& sym_map, size_t lvl,
                                         size_t seed){
    size_t len=0, sym;
    uint64_t fingerprint;
    std::vector<uint64_t> tmp_seq;

    if(lvl==0){
        do{
            sym = rule_stream.read(pos, pos+width-1);
            size_t tmp = sym>>1UL;
            tmp_seq.push_back(XXH64(&tmp, sizeof(size_t), seed));
            pos+=width;
            len++;
        }while(!(sym & 1UL));
    }else{
        do{
            sym = (rule_stream.read(pos, pos+width-1))>>1UL;
            sym = sym_map[sym-1];// marks the position in the merged grammar
            tmp_seq.push_back(prev_fps[sym]);
            pos+=width;
            len++;
        }while(!(sym & 1UL));
    }

    fingerprint = XXH64(tmp_seq.data(), sizeof(uint64_t)*len, seed);
    return {fingerprint, len};
}

int compare_rules(size_t pos_a, bitstream<size_t>& rule_a, uint8_t width_a,
                  size_t pos_b, bitstream<size_t>& rule_b, uint8_t width_b,
                  std::vector<uint32_t>& mt_map_a, std::vector<uint32_t>& mt_map_b){

    bool last_a=false, last_b, equal;
    size_t sym_a, sym_b;
    while(!last_a){
        sym_a = rule_a.read(pos_a, pos_a+width_a-1);
        last_a = sym_a & 1UL;
        sym_a = mt_map_a[sym_a>>1UL];

        sym_b = rule_b.read(pos_b, pos_b+width_b-1);
        last_b = sym_b & 1UL;
        sym_b = mt_map_b[sym_b>>1UL];
        equal = (sym_a == sym_b) && (last_a == last_b);

        if(!equal){
            if(sym_a<sym_b || (sym_a==sym_b && last_a && !last_b)){
                return -1;
            }else{
                return 1;
            }
        }

        pos_a+=width_a;
        pos_b+=width_b;
    }
    return 0;
}

void append_rule(size_t& s_pos, size_t s_width, size_t s_len, bitstream<size_t>& source, std::vector<uint32_t>& mt_map,
                 size_t& d_pos, size_t d_width, bitstream<size_t>& dest){

    //check the appended rule fits the buffer of merged rules
    size_t min_size = INT_CEIL((d_pos+(s_len*d_width)), (sizeof(size_t)*8));
    if(dest.stream_size<min_size){
        dest.stream_size = (min_size * 120)/100;
        dest.stream = (size_t *)realloc(dest.stream, dest.stream_size);
    }

    bool last;
    size_t sym;
    for(size_t i=0;i<s_len;i++){
        sym = source.read(s_pos, s_pos+s_width-1);
        last = sym & 1UL;
        sym>>=1UL;
        sym = ((mt_map[sym]<<1UL) | last);
        dest.write(d_pos, d_pos+d_width-1, sym);
        d_pos+=d_width;
        s_pos+=s_width;
    }
    assert(last & 1UL);
}

void merge_level(bitstream<size_t> &rl_stream_a, lvl_metadata_type &lvl_met_a, std::vector<uint32_t> &mt_map_a,
                 bitstream<size_t> &rl_stream_b, lvl_metadata_type &lvl_met_b, std::vector<uint32_t> &mt_map_b,
                 bitstream<size_t> &rl_stream_m, lvl_metadata_type &prev_lvl_met,
                 uint64_t &fp_seed, std::vector<uint64_t> &prev_fps, size_t lvl) {

    assert(prev_lvl_met.n_rules>0);

    uint64_t fp_a, fp_b;
    size_t curr_rule_a=0, curr_rule_b=0;
    size_t curr_pos_a=0, curr_pos_b=0;
    size_t len_a, len_b;

    lvl_metadata_type lvl_met_m{};
    lvl_met_m.sym_width = sym_width(prev_lvl_met.n_rules)+1;
    uint8_t m_width = lvl_met_m.sym_width;
    size_t curr_pos_m=0;

    size_t new_size = std::max(lvl_met_a.sym_width*lvl_met_a.tot_symbols, lvl_met_b.sym_width*lvl_met_b.tot_symbols);
    new_size = INT_CEIL(new_size, (sizeof(size_t)*8));

    //initialize the merged grammar
    if(rl_stream_m.stream== nullptr) {
        rl_stream_m.stream = (size_t *)malloc(new_size*sizeof(size_t));
    }else {
        rl_stream_m.stream = (size_t *)realloc(rl_stream_m.stream, new_size*sizeof(size_t));
    }
    rl_stream_m.stream_size = new_size;
    //

    std::vector<uint8_t> merge_marks;
    std::vector<uint64_t> new_fps;

    size_t lvl_rules = std::min(lvl_met_a.n_rules, lvl_met_b.n_rules);

    uint8_t a_width = lvl_met_a.sym_width;
    uint8_t b_width = lvl_met_b.sym_width;

    std::tie(fp_a, len_a) = get_rule_info(rl_stream_a, curr_pos_a, a_width, prev_fps, mt_map_a, lvl, fp_seed);
    std::tie(fp_b, len_b) = get_rule_info(rl_stream_b, curr_pos_b, b_width, prev_fps, mt_map_b, lvl, fp_seed);

    while(curr_rule_a<lvl_rules && curr_rule_b<lvl_rules){

        if(fp_a<fp_b){
            //write rule from A
            append_rule(curr_pos_a, a_width, len_a, rl_stream_a, mt_map_a, curr_pos_m, m_width, rl_stream_m);
            new_fps.push_back(fp_a);
            curr_rule_a++;
            lvl_met_m.tot_symbols+=len_a;

            merge_marks.push_back(1);
        } else if(fp_b<fp_a) {
            //write rule from B
            append_rule(curr_pos_b, b_width, len_b, rl_stream_b, mt_map_b, curr_pos_m, m_width, rl_stream_m);
            new_fps.push_back(fp_b);
            curr_rule_b++;
            lvl_met_m.tot_symbols+=len_b;
            merge_marks.push_back(2);
        } else {

            int eq_seq = compare_rules(curr_pos_a, rl_stream_a, a_width,
                                       curr_pos_b, rl_stream_b, b_width,
                                       mt_map_a, mt_map_b);

            if(eq_seq==0){
                //write rule from A
                append_rule(curr_pos_a, a_width, len_a, rl_stream_a, mt_map_a, curr_pos_m, m_width, rl_stream_m);
                new_fps.push_back(fp_a);
                curr_rule_a++;
                lvl_met_m.tot_symbols+=len_a;

                //ignore rule from B
                curr_pos_b +=len_b*b_width;
                curr_rule_b++;

                merge_marks.push_back(3);
            } else if(eq_seq<0){ // collision: break ties by lex. rank
                //TODO implement this
            } else{
                //TODO implement this
            }
        }

        if((curr_rule_a < lvl_rules) &&
           (merge_marks.back() & 1UL)){
            std::tie(fp_a, len_a) = get_rule_info(rl_stream_a, curr_pos_a, a_width, prev_fps, mt_map_a, lvl, fp_seed);
        }

        if((curr_rule_b < lvl_rules) &&
           (merge_marks.back() & 2UL)){
            std::tie(fp_b, len_b) = get_rule_info(rl_stream_b, curr_pos_b, b_width, prev_fps, mt_map_b, lvl, fp_seed);
        }
    }

    while(curr_rule_a<lvl_met_a.n_rules){
        //write rule from A
        append_rule(curr_pos_a, a_width, len_a, rl_stream_a, mt_map_a, curr_pos_m, m_width, rl_stream_m);
        new_fps.push_back(fp_a);
        curr_rule_a++;
        lvl_met_m.tot_symbols+=len_a;
        merge_marks.push_back(1);

        if(curr_rule_a<lvl_met_a.n_rules){
            std::tie(fp_a, len_a) = get_rule_info(rl_stream_a, curr_pos_a, a_width, prev_fps, mt_map_a, lvl, fp_seed);
        }
    }

    while(curr_rule_b<lvl_met_b.n_rules){
        //write rule from B
        append_rule(curr_pos_b, b_width, len_b, rl_stream_b, mt_map_b, curr_pos_m, m_width, rl_stream_m);
        new_fps.push_back(fp_b);
        curr_rule_b++;
        lvl_met_m.tot_symbols+=len_b;
        merge_marks.push_back(2);

        if(curr_rule_b<lvl_met_b.n_rules){
            std::tie(fp_b, len_b) = get_rule_info(rl_stream_b, curr_pos_b, b_width, prev_fps, mt_map_b, lvl, fp_seed);
        }
    }

    //update mapping values
    size_t pos_a=0, pos_b=0, rank=0;
    mt_map_a.resize(lvl_met_a.n_rules);
    mt_map_b.resize(lvl_met_b.n_rules);
    for(unsigned char merge_mark : merge_marks){
        if(merge_mark==1){
            mt_map_a[pos_a++] = rank;
        } else if(merge_mark==2){
            mt_map_a[pos_b++] = rank;
        } else {
            mt_map_a[pos_a++] = rank;
            mt_map_b[pos_b++] = rank;
        }
        rank++;
    }
    prev_fps.swap(new_fps);
    lvl_met_m.n_rules = rank;

    prev_lvl_met = lvl_met_m;

    //TODO resize the stream with the merged rules
}

template<class gram_type>
void merge_grammars(gram_type& p_gram_a, gram_type& p_gram_b, std::vector<uint64_t>& fp_seeds) {

    //get the information about the terminal symbols
    lvl_metadata_type prev_lvl_met = p_gram_a.rules_metadata[0];

    std::vector<uint32_t> mt_map_a;
    mt_map_a.resize(prev_lvl_met.tot_symbols);

    std::vector<uint32_t> mt_map_b;
    mt_map_b.resize(prev_lvl_met.tot_symbols);

    std::vector<uint64_t> lvl_fps;
    lvl_fps.resize(prev_lvl_met.tot_symbols);

    for(size_t i=0;i<lvl_fps.size();i++){
        lvl_fps[i] = XXH64(&i, sizeof(size_t), fp_seeds[0]);
        mt_map_a[i] = i;
        mt_map_b[i] = i;
    }

    size_t shared_levels = std::min(p_gram_a.rules_metadata.size(), p_gram_b.rules_metadata.size());
    for(size_t i=0;i<shared_levels;i++) {
        bitstream<size_t> merged_rules;
        merge_level(p_gram_a.rules, p_gram_a.rules_metadata[i], mt_map_a,
                    p_gram_b.rules, p_gram_b.rules_metadata[i], mt_map_b,
                    merged_rules, prev_lvl_met,
                    fp_seeds[i+1], lvl_fps, i);
    }

    //TODO concatenate the compressed strings
}

#endif //LCG_PARTIAL_GRAM_H
