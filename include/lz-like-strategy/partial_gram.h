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

    std::vector<bitstream<size_t>> rules;
    std::vector<lvl_metadata_type> rules_metadata;
    size_t lvl = 0;
    //size_t acc_bits = 0;
    //typedef string_map_iterator<partial_gram<ter_type>> iter_type;

    partial_gram(){
        lvl_metadata_type ter{};
        ter.sym_width = std::numeric_limits<ter_type>::digits;
        ter.tot_symbols = std::numeric_limits<ter_type>::max()+1;
        ter.n_rules = ter.tot_symbols;
        ter.terminals = true;
        rules_metadata.push_back(ter);
    }


    size_t serialize(std::ofstream &ofs){
        size_t n_levels = rules.size();
        size_t written_bytes = serialize_elm(ofs, n_levels);
        for(const auto & lvl_rules : rules){
            written_bytes += lvl_rules.serialize(ofs);
        }
        written_bytes+= serialize_plain_vector(ofs, rules_metadata);
        return written_bytes;
    }

    void load(std::ifstream &ifs){
        load_elm(ifs, lvl);
        rules.resize(lvl);
        for(auto & lvl_rules : rules){
            lvl_rules.load(ifs);
        }
        load_plain_vector(ifs, rules_metadata);
    }

    size_t serialize_to_fd(int fd){
        assert(lvl==(rules_metadata.size()-1));

        write(fd, &lvl, sizeof(size_t));
        size_t written_bytes = sizeof(size_t);

        for(size_t i=0;i<lvl;i++){
            write(fd, &rules[i].stream_size, sizeof(size_t));
            written_bytes += sizeof(size_t);
            write(fd, rules[i].stream, rules[i].stream_size*sizeof(size_t));
            written_bytes += rules[i].stream_size*sizeof(size_t);
        }

        write(fd, rules_metadata.data(), rules_metadata.size()*sizeof(lvl_metadata_type));
        written_bytes += rules_metadata.size()*sizeof(lvl_metadata_type);
        return written_bytes;
    }

    size_t load_from_fs(int fd){
        read(fd, &lvl, sizeof(size_t));
        size_t read_bytes = sizeof(size_t);

        rules.resize(lvl);
        for(size_t i=0;i<lvl;i++){
            size_t tot_words;
            read(fd, &tot_words, sizeof(size_t));
            read_bytes+=sizeof(size_t);

            if(rules[i].stream_size==0){
                assert(rules[i].stream== nullptr);
                rules[i].stream = (size_t *)malloc(tot_words*sizeof(size_t));
            }else if(rules[i].stream_size<tot_words){
                rules[i].stream = (size_t *)realloc(rules[i].stream, tot_words*sizeof(size_t));
            }
            rules[i].stream_size = tot_words;

            assert(rules[i].stream!= nullptr);
            read(fd, rules[i].stream, rules[i].stream_size*sizeof(size_t));
            read_bytes+=rules[i].stream_size*sizeof(size_t);
        }

        rules_metadata.resize(lvl+1);
        read(fd, rules_metadata.data(), rules_metadata.size()*sizeof(lvl_metadata_type));
        read_bytes+=rules_metadata.size()*sizeof(lvl_metadata_type);
        return read_bytes;
    }

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
        size_t tot_words = INT_CEIL(new_bits, (sizeof(size_t)*8));

        if(rules.size()<lvl+1){
            rules.resize(lvl+1);
        }

        if(rules[lvl].stream_size==0){
            assert(rules[lvl].stream== nullptr);
            rules[lvl].stream = (size_t *)malloc(tot_words*sizeof(size_t));
        }else if(rules[lvl].stream_size<tot_words){
            rules[lvl].stream = (size_t *)realloc(rules[lvl].stream, tot_words*sizeof(size_t));
        }
        rules[lvl].stream_size = tot_words;

        size_t acc_bits=0;
        for(size_t i=0;i<phrase_set.size();i++){
            size_t idx = perm[i].first;
            uint32_t source = phrase_set[idx].source/sizeof(sym_type);
            uint32_t len = (phrase_set[idx].len/sizeof(sym_type));
            uint32_t last = source + len-1;
            for(size_t j=source;j<last;j++){
                rules[lvl].write(acc_bits, acc_bits+lvl_met.sym_width-1, text[j]<<1UL);
                acc_bits+=lvl_met.sym_width;
            }
            rules[lvl].write(acc_bits, acc_bits+lvl_met.sym_width-1, ((text[last]<<1UL) | 1UL));
            acc_bits+=lvl_met.sym_width;
        }
        assert(INT_CEIL(acc_bits, (sizeof(size_t)*8))==rules[lvl].stream_size);
        rules_metadata.push_back(lvl_met);

        lvl++;
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
        size_t tot_words = INT_CEIL(new_bits, (sizeof(size_t)*8));

        if(rules[lvl].stream_size==0){
            assert(rules[lvl].stream== nullptr);
            rules[lvl].stream = (size_t *)malloc(tot_words*sizeof(size_t));
        }else if(rules[lvl].stream_size<tot_words){
            rules[lvl].stream = (size_t *)realloc(rules[lvl].stream, tot_words*sizeof(size_t));
        }
        rules[lvl].stream_size = tot_words;

        //we skip the separator symbols
        size_t last = size-2;
        size_t acc_bits=0;
        for(size_t j=0;j<last;j+=2){
            rules[lvl].write(acc_bits, acc_bits+lvl_met.sym_width-1, text[j]<<1UL);
            acc_bits+=lvl_met.sym_width;
        }
        rules[lvl].write(acc_bits, acc_bits+lvl_met.sym_width-1, ((text[last]<<1UL) | 1UL));
        acc_bits+=lvl_met.sym_width;

        assert(INT_CEIL(acc_bits, (sizeof(size_t)*8))==rules[lvl].stream_size);
        rules_metadata.push_back(lvl_met);
        lvl++;
    }

    void reset_grammar(){
        for(auto & lvl_rules : rules){
            memset(lvl_rules.stream, 0,  lvl_rules.stream_size*sizeof(size_t));
        }
        lvl=0;
        rules_metadata.resize(1);
    }

    ~partial_gram(){
        for(auto & lvl_rules : rules){
            if(lvl_rules.stream!= nullptr){
                free(lvl_rules.stream);
            }
        }
        std::vector<lvl_metadata_type>().swap(rules_metadata);
    }

    void print_stats(){
        for(size_t i=0;i<rules_metadata.size();i++){
            std::cout<<"Level "<<i+1<<std::endl;
            std::cout<<"  Number of rules? "<<rules_metadata[i].n_rules<<std::endl;
            std::cout<<"  Number of symbols? "<<rules_metadata[i].tot_symbols<<std::endl;
            std::cout<<"  Is terminal? "<<rules_metadata[i].terminals<<std::endl;
        }
    }
};

std::pair<uint64_t,size_t> get_rule_info(bitstream<size_t>& rule_stream, size_t pos, size_t width,
                                         std::vector<uint64_t>& prev_fps, std::vector<uint32_t>& mt_map, size_t lvl,
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
            sym = rule_stream.read(pos, pos+width-1);
            size_t tmp = sym>>1UL;
            assert(tmp>0 && tmp<mt_map.size());
            tmp = mt_map[tmp];
            assert(tmp>0 && tmp<prev_fps.size());
            tmp_seq.push_back(prev_fps[tmp]);
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

        assert(pos_a+width_a-1<rule_a.stream_size*64);
        sym_a = rule_a.read(pos_a, pos_a+width_a-1);
        last_a = sym_a & 1UL;
        assert((sym_a>>1UL)>0);
        assert((sym_a>>1UL)<mt_map_a.size());
        sym_a = mt_map_a[sym_a>>1UL];

        assert((pos_b+width_b-1)<(rule_b.stream_size*64));
        sym_b = rule_b.read(pos_b, pos_b+width_b-1);
        last_b = sym_b & 1UL;

        assert((sym_b>>1UL)>0);
        assert((sym_b>>1UL)<mt_map_b.size());
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
        std::cout<<"No creo que pase por aca , o si "<<dest.stream_size<<" "<<min_size<<" "<<(min_size*120)/100<<std::endl;
        dest.stream_size = (min_size * 120)/100;
        dest.stream = (size_t *)realloc(dest.stream, dest.stream_size);
    }

    bool last;
    size_t sym;
    for(size_t i=0;i<s_len;i++){
        sym = source.read(s_pos, s_pos+s_width-1);
        last = sym & 1UL;
        sym>>=1UL;
        assert(sym>0);
        sym = ((mt_map[sym]<<1UL) | last);
        assert(sym>0);
        dest.write(d_pos, d_pos+d_width-1, sym);
        d_pos+=d_width;
        s_pos+=s_width;
    }
    assert(last & 1UL);
}

void merge_level(bitstream<size_t> &stream_a, lvl_metadata_type &lvl_met_a, std::vector<uint32_t> &mt_map_a,
                 bitstream<size_t> &stream_b, lvl_metadata_type &lvl_met_b, std::vector<uint32_t> &mt_map_b,
                 bitstream<size_t> &stream_m, lvl_metadata_type &prev_lvl_met,
                 uint64_t &fp_seed, std::vector<uint64_t> &prev_fps, size_t lvl) {

    assert(prev_lvl_met.n_rules>0);

    uint64_t fp_a, fp_b;
    size_t curr_rule_a=0, curr_rule_b=0;
    size_t curr_pos_a=0, curr_pos_b=0;
    size_t len_a, len_b;

    lvl_metadata_type lvl_met_m{};
    lvl_met_m.sym_width = sym_width(prev_lvl_met.n_rules)+1;
    lvl_met_m.terminals = lvl_met_a.terminals;

    uint8_t m_width = lvl_met_m.sym_width;
    size_t curr_pos_m=0;

    //size_t new_size = std::max(lvl_met_a.sym_width*lvl_met_a.tot_symbols, lvl_met_b.sym_width*lvl_met_b.tot_symbols);
    size_t new_size = lvl_met_a.sym_width*lvl_met_a.tot_symbols + lvl_met_b.sym_width*lvl_met_b.tot_symbols;
    new_size = INT_CEIL(new_size, (sizeof(size_t)*8));

    //initialize the merged grammar
    if(stream_m.stream== nullptr) {
        stream_m.stream = (size_t *)malloc(new_size*sizeof(size_t));
    }else {
        stream_m.stream = (size_t *)realloc(stream_m.stream, new_size*sizeof(size_t));
    }
    stream_m.stream_size = new_size;
    //

    std::vector<uint8_t> merge_marks;
    std::vector<uint64_t> new_fps(1);

    uint8_t a_width = lvl_met_a.sym_width;
    uint8_t b_width = lvl_met_b.sym_width;

    std::tie(fp_a, len_a) = get_rule_info(stream_a, curr_pos_a, a_width, prev_fps, mt_map_a, lvl, fp_seed);
    std::tie(fp_b, len_b) = get_rule_info(stream_b, curr_pos_b, b_width, prev_fps, mt_map_b, lvl, fp_seed);

    while(curr_rule_a<lvl_met_a.n_rules && curr_rule_b<lvl_met_b.n_rules){

        if(fp_a<fp_b){
            //write rule from A
            append_rule(curr_pos_a, a_width, len_a, stream_a, mt_map_a, curr_pos_m, m_width, stream_m);
            curr_rule_a++;

            new_fps.push_back(fp_a);
            lvl_met_m.tot_symbols+=len_a;
            merge_marks.push_back(1);
        } else if(fp_b<fp_a) {
            //write rule from B
            append_rule(curr_pos_b, b_width, len_b, stream_b, mt_map_b, curr_pos_m, m_width, stream_m);
            curr_rule_b++;

            new_fps.push_back(fp_b);
            lvl_met_m.tot_symbols+=len_b;
            merge_marks.push_back(2);
        } else {

            int eq_seq = compare_rules(curr_pos_a, stream_a, a_width,
                                       curr_pos_b, stream_b, b_width,
                                       mt_map_a, mt_map_b);

            //TODO asdasd
            if(eq_seq!=0){
                std::cout<<"Colision? "<<fp_a<<" "<<fp_b<<std::endl;
                bool last_a=false, last_b, equal;
                size_t sym_a, sym_b;
                size_t pos_a = curr_pos_a;
                size_t pos_b = curr_pos_b;
                while(!last_a){
                    sym_a = stream_a.read(pos_a, pos_a+a_width-1);
                    last_a = sym_a & 1UL;
                    sym_a = mt_map_a[sym_a>>1UL];

                    sym_b = stream_b.read(pos_b, pos_b+a_width-1);
                    last_b = sym_b & 1UL;
                    sym_b = mt_map_b[sym_b>>1UL];

                    equal = (sym_a == sym_b) && (last_a == last_b);
                    std::cout<<sym_a<<" "<<sym_b<<" "<<last_a<<" "<<last_b<<std::endl;
                    if(!equal){
                        break;
                    }
                    pos_a+=a_width;
                    pos_b+=b_width;
                }
            }
            //

            if(eq_seq==0){
                //write rule from A
                append_rule(curr_pos_a, a_width, len_a, stream_a, mt_map_a, curr_pos_m, m_width, stream_m);
                curr_rule_a++;
                curr_rule_b++;
                curr_pos_b +=len_b*b_width;

                new_fps.push_back(fp_a);
                lvl_met_m.tot_symbols+=len_a;
                merge_marks.push_back(3);
            } else if(eq_seq<0) { // collision: break ties by lex. rank
                //write rule from A
                append_rule(curr_pos_a, a_width, len_a, stream_a, mt_map_a, curr_pos_m, m_width, stream_m);
                curr_rule_a++;

                new_fps.push_back(fp_a);
                lvl_met_m.tot_symbols+=len_a;
                merge_marks.push_back(1);
            } else{
                //write rule from B
                append_rule(curr_pos_b, b_width, len_b, stream_b, mt_map_b, curr_pos_m, m_width, stream_m);
                curr_rule_b++;

                new_fps.push_back(fp_b);
                lvl_met_m.tot_symbols+=len_b;
                merge_marks.push_back(2);
            }
        }

        if((merge_marks.back() & 1) && curr_rule_a<lvl_met_a.n_rules){
            std::tie(fp_a, len_a) = get_rule_info(stream_a, curr_pos_a, a_width, prev_fps, mt_map_a, lvl, fp_seed);
        }

        if((merge_marks.back() & 2) && curr_rule_b<lvl_met_b.n_rules){
            std::tie(fp_b, len_b) = get_rule_info(stream_b, curr_pos_b, b_width, prev_fps, mt_map_b, lvl, fp_seed);
        }
    }

    while(curr_rule_a<lvl_met_a.n_rules){
        //write rule from A
        append_rule(curr_pos_a, a_width, len_a, stream_a, mt_map_a, curr_pos_m, m_width, stream_m);
        curr_rule_a++;

        new_fps.push_back(fp_a);
        lvl_met_m.tot_symbols+=len_a;
        merge_marks.push_back(1);

        if(curr_rule_a<lvl_met_a.n_rules){
            std::tie(fp_a, len_a) = get_rule_info(stream_a, curr_pos_a, a_width, prev_fps, mt_map_a, lvl, fp_seed);
        }
    }

    while(curr_rule_b<lvl_met_b.n_rules){
        //write rule from B
        append_rule(curr_pos_b, b_width, len_b, stream_b, mt_map_b, curr_pos_m, m_width, stream_m);
        curr_rule_b++;

        new_fps.push_back(fp_b);
        lvl_met_m.tot_symbols+=len_b;
        merge_marks.push_back(2);

        if(curr_rule_b<lvl_met_b.n_rules){
            std::tie(fp_b, len_b) = get_rule_info(stream_b, curr_pos_b, b_width, prev_fps, mt_map_b, lvl, fp_seed);
        }
    }

    //update mapping values
    //it is one-based is the metasymbols are also one-based because of the separator symbol
    size_t pos_a=1, pos_b=1, rank=1;
    mt_map_a.resize(lvl_met_a.n_rules+1);
    mt_map_b.resize(lvl_met_b.n_rules+1);
    for(unsigned char merge_mark : merge_marks){
        if(merge_mark==1){
            mt_map_a[pos_a++] = rank;
        } else if(merge_mark==2){
            mt_map_b[pos_b++] = rank;
        } else {
            mt_map_a[pos_a++] = rank;
            mt_map_b[pos_b++] = rank;
        }
        rank++;
    }
    prev_fps.swap(new_fps);
    lvl_met_m.n_rules = rank-1;
    prev_lvl_met = lvl_met_m;

    //TODO resize the stream with the merged rules
}

template<class gram_type>
void merge_two_grammars(gram_type& p_gram_a, gram_type& p_gram_b, std::vector<uint64_t>& fp_seeds) {

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

    gram_type p_gram_c;
    p_gram_c.rules.resize(std::max(p_gram_a.rules_metadata.size(), p_gram_b.rules_metadata.size()));
    p_gram_c.rules_metadata.push_back(p_gram_a.rules_metadata[0]);
    size_t shared_levels = std::min(p_gram_a.rules_metadata.size(), p_gram_b.rules_metadata.size());

    for(size_t i=0;i<shared_levels;i++) {

        merge_level(p_gram_a.rules[i], p_gram_a.rules_metadata[i+1], mt_map_a,
                    p_gram_b.rules[i], p_gram_b.rules_metadata[i+1], mt_map_b,
                    p_gram_c.rules[i], prev_lvl_met,
                    fp_seeds[i+1], lvl_fps, i);
        //free(p_gram_a.rules[i].stream);
        //free(p_gram_b.rules[i].stream);
        std::cout<<"level "<<i<<" "<<shared_levels<<std::endl;
    }
}

#endif //LCG_PARTIAL_GRAM_H
