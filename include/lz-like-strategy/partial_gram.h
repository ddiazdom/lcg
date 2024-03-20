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

    [[nodiscard]] inline size_t n_bits() const {
        return tot_symbols*sym_width;
    }
};

template<class ter_type>
struct partial_gram {

    std::vector<bitstream<size_t>> rules;
    std::vector<lvl_metadata_type> metadata;
    size_t lvl = 0;

    partial_gram(){
        lvl_metadata_type ter{};
        ter.sym_width = std::numeric_limits<ter_type>::digits;
        ter.tot_symbols = std::numeric_limits<ter_type>::max()+1;
        ter.n_rules = ter.tot_symbols;
        ter.terminals = true;
        metadata.push_back(ter);
        rules.resize(40);//strings up to 1TB are allowed
    }

    partial_gram& swap(partial_gram& other){
        rules.swap(other.rules);
        metadata.swap(other.rules_metadata);
        std::swap(lvl, other.lvl);
        return *this;
    }

    size_t serialize(std::ofstream &ofs){
        assert(lvl==rules.size());
        size_t written_bytes = serialize_elm(ofs, lvl);
        for(size_t i=0;i<lvl;i++){
            written_bytes += rules[i].serialize(ofs);
        }
        written_bytes+= serialize_plain_vector(ofs, metadata);
        return written_bytes;
    }

    void load(std::ifstream &ifs){
        load_elm(ifs, lvl);
        rules.resize(lvl);
        for(size_t i=0;i<lvl;i++){
            rules[i].load(ifs);
        }
        load_plain_vector(ifs, metadata);
    }

    size_t serialize_to_fd(int fd){
        assert(lvl==(metadata.size()-1));

        write(fd, &lvl, sizeof(size_t));
        size_t written_bytes = sizeof(size_t);

        for(size_t i=0;i<lvl;i++){
            write(fd, &rules[i].stream_size, sizeof(size_t));
            written_bytes += sizeof(size_t);
            write(fd, rules[i].stream, rules[i].stream_size*sizeof(size_t));
            written_bytes += rules[i].stream_size*sizeof(size_t);
        }

        write(fd, metadata.data(), metadata.size()*sizeof(lvl_metadata_type));
        written_bytes += metadata.size()*sizeof(lvl_metadata_type);
        return written_bytes;
    }

    size_t load_from_fd(int fd){
        read(fd, &lvl, sizeof(size_t));
        size_t read_bytes = sizeof(size_t);

        rules.resize(lvl);
        for(size_t i=0;i<lvl;i++){
            size_t tot_words;
            read(fd, &tot_words, sizeof(size_t));
            read_bytes+=sizeof(size_t);

            rules[i].reserve_in_words(tot_words);

            assert(rules[i].stream!= nullptr);
            read(fd, rules[i].stream, rules[i].stream_size*sizeof(size_t));
            read_bytes+=rules[i].stream_size*sizeof(size_t);
        }

        metadata.resize(lvl+1);
        read(fd, metadata.data(), metadata.size()*sizeof(lvl_metadata_type));
        read_bytes+=metadata.size()*sizeof(lvl_metadata_type);
        return read_bytes;
    }

    template<class sym_type>
    void append_new_lvl(sym_type* text, const lz_map::phrase_list_t &phrase_set, size_t tot_symbols,
                        std::vector<std::pair<uint32_t, uint64_t>>& perm){

        lvl_metadata_type lvl_met{};
        lvl_met.n_rules = phrase_set.size();

        //the extra bit is to mark the end of each rule
        lvl_met.sym_width = sym_width(metadata.back().n_rules)+1;
        lvl_met.tot_symbols = tot_symbols;
        lvl_met.terminals = false;
        rules[lvl].reserve_in_bits(lvl_met.n_bits());

        size_t acc_bits=0;
        for(size_t i=0;i<phrase_set.size();i++){
            size_t idx = perm[i].first;
            uint32_t source = phrase_set[idx].source/sizeof(sym_type);
            uint32_t len = (phrase_set[idx].len/sizeof(sym_type));
            uint32_t last = source + len-1;
            for(size_t j=source;j<last;j++){
                assert(text[j]>0);
                rules[lvl].write(acc_bits, acc_bits+lvl_met.sym_width-1, text[j]<<1UL);
                acc_bits+=lvl_met.sym_width;

                //TODO testing
                /*if(i<2){
                    std::cout<<int(text[j])<<" ";
                }*/
                //
            }

            /*if(i<2){
                std::cout<<int(text[last])<<" -> "<<i<<" fp:"<<perm[i].second<<std::endl;
            }*/

            rules[lvl].write(acc_bits, acc_bits+lvl_met.sym_width-1, ((text[last]<<1UL) | 1UL));
            acc_bits+=lvl_met.sym_width;
        }

        assert(acc_bits==lvl_met.n_bits());
        metadata.push_back(lvl_met);
        lvl++;
    }

    template<class sym_type>
    void add_compressed_string(sym_type* text, size_t size){

        lvl_metadata_type lvl_met{};
        lvl_met.n_rules = 1;
        //the extra bit is to mark the end of each rule
        lvl_met.sym_width = sym_width(metadata.back().n_rules)+1;
        lvl_met.tot_symbols = size/2;//we divide by 2 because each string is followed by a 0 (a byproduct of the compression)
        lvl_met.terminals = false;

        rules[lvl].reserve_in_bits(lvl_met.n_bits());

        //we skip the separator symbols (encoded as a 0)
        size_t last = size-2;
        size_t acc_bits=0;
        for(size_t j=0;j<last;j+=2){
            assert(text[j]>0);
            rules[lvl].write(acc_bits, acc_bits+lvl_met.sym_width-1, text[j]<<1UL);
            acc_bits+=lvl_met.sym_width;
        }
        rules[lvl].write(acc_bits, acc_bits+lvl_met.sym_width-1, ((text[last]<<1UL) | 1UL));
        acc_bits+=lvl_met.sym_width;

        assert(acc_bits==lvl_met.n_bits());
        metadata.push_back(lvl_met);
        lvl++;
    }

    void reset_grammar(){
        /*for(auto & lvl_rules : rules){
            memset(lvl_rules.stream, 0,  lvl_rules.stream_size*sizeof(size_t));
        }*/
        lvl=0;
        metadata.resize(1);
    }

    ~partial_gram(){
        for(auto & lvl_rules : rules){
            lvl_rules.destroy();
        }
        std::vector<lvl_metadata_type>().swap(metadata);
    }

    void print_stats(){
        for(size_t i=0;i<metadata.size();i++){
            std::cout<<"Level "<<i+1<<std::endl;
            std::cout<<"  Number of rules? "<<metadata[i].n_rules<<std::endl;
            std::cout<<"  Number of symbols? "<<metadata[i].tot_symbols<<std::endl;
            std::cout<<"  Is terminal? "<<metadata[i].terminals<<std::endl;
        }
    }
};

std::pair<uint64_t,size_t> get_rule_info(bitstream<size_t>& rule_stream, size_t pos, size_t width,
                                         std::vector<uint64_t>& prev_fps, std::vector<uint32_t>& mt_map, size_t fp_seed){
    size_t len=0, sym;
    uint64_t fingerprint;
    std::vector<uint64_t> tmp_phrase;

    do{
        sym = rule_stream.read(pos, pos+width-1);
        size_t tmp = sym>>1UL;
        assert(tmp<mt_map.size());
        tmp = mt_map[tmp];
        assert(tmp<prev_fps.size());
        tmp_phrase.push_back(prev_fps[tmp]);
        pos+=width;
        len++;
    }while(!(sym & 1UL));

    fingerprint = XXH64(tmp_phrase.data(), sizeof(uint64_t)*len, fp_seed);
    return {fingerprint, len};
}

int compare_rules(size_t pos_a, bitstream<size_t>& rule_a, uint8_t width_a,
                  size_t pos_b, bitstream<size_t>& rule_b, uint8_t width_b,
                  std::vector<uint32_t>& mt_map_a, std::vector<uint32_t>& mt_map_b){

    bool last_a=false, last_b, equal;
    size_t sym_a, sym_b;
    while(!last_a){

        assert((pos_a+width_a-1)<rule_a.capacity_in_bits());
        sym_a = rule_a.read(pos_a, pos_a+width_a-1);
        last_a = sym_a & 1UL;
        assert((sym_a>>1UL)>0);
        assert((sym_a>>1UL)<mt_map_a.size());
        sym_a = mt_map_a[sym_a>>1UL];

        assert((pos_b+width_b-1)<(rule_b.capacity_in_bits()));
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
    size_t min_bits = d_pos+(s_len*d_width);
    if(dest.capacity_in_bits()<=min_bits) {
        min_bits = INT_CEIL((min_bits * 120), 100);//20% increase in the stream size
        dest.reserve_in_bits(min_bits);
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
        assert((d_pos+d_width)<=dest.capacity_in_bits());
        dest.write(d_pos, d_pos+d_width-1, sym);
        d_pos+=d_width;
        s_pos+=s_width;
    }
    assert(last & 1UL);
}

void merge_level(bitstream<size_t> &stream_a, lvl_metadata_type &lvl_met_a, std::vector<uint32_t> &mt_map_a,
                 bitstream<size_t> &stream_b, lvl_metadata_type &lvl_met_b, std::vector<uint32_t> &mt_map_b,
                 bitstream<size_t> &stream_c, lvl_metadata_type &prev_lvl_met,
                 uint64_t &fp_seed, std::vector<uint64_t> &prev_fps) {

    assert(prev_lvl_met.n_rules>0);

    uint64_t fp_a, fp_b;
    size_t curr_rule_a=0, curr_rule_b=0;
    size_t curr_pos_a=0, curr_pos_b=0;
    size_t len_a, len_b;

    lvl_metadata_type lvl_met_c{};
    lvl_met_c.sym_width = sym_width(prev_lvl_met.n_rules)+1;
    lvl_met_c.terminals = lvl_met_a.terminals;

    uint8_t m_width = lvl_met_c.sym_width;
    size_t curr_pos_m=0;

    stream_c.reserve_in_bits(std::max(lvl_met_a.n_bits(), lvl_met_b.n_bits()));

    std::vector<uint8_t> merge_marks;
    std::vector<uint64_t> new_fps(1, 0);

    uint8_t a_width = lvl_met_a.sym_width;
    uint8_t b_width = lvl_met_b.sym_width;

    std::tie(fp_a, len_a) = get_rule_info(stream_a, curr_pos_a, a_width, prev_fps, mt_map_a, fp_seed);
    std::tie(fp_b, len_b) = get_rule_info(stream_b, curr_pos_b, b_width, prev_fps, mt_map_b, fp_seed);

    uint64_t prev_fp_a = fp_a;
    uint64_t prev_fp_b = fp_b;

    while(curr_rule_a<lvl_met_a.n_rules && curr_rule_b<lvl_met_b.n_rules){

        if(fp_a<fp_b){
            //write rule from A
            append_rule(curr_pos_a, a_width, len_a, stream_a, mt_map_a, curr_pos_m, m_width, stream_c);
            curr_rule_a++;

            new_fps.push_back(fp_a);
            lvl_met_c.tot_symbols+=len_a;
            merge_marks.push_back(1);
        } else if(fp_b<fp_a) {
            //write rule from B
            append_rule(curr_pos_b, b_width, len_b, stream_b, mt_map_b, curr_pos_m, m_width, stream_c);
            curr_rule_b++;

            new_fps.push_back(fp_b);
            lvl_met_c.tot_symbols+=len_b;
            merge_marks.push_back(2);
        } else {

            int eq_seq = compare_rules(curr_pos_a, stream_a, a_width,
                                       curr_pos_b, stream_b, b_width,
                                       mt_map_a, mt_map_b);

            //TODO asdasd
            if(eq_seq!=0){
                std::cout<<"Collision warning:  "<<fp_a<<" "<<fp_b<<" "<<curr_pos_a<<" "<<curr_pos_b<<std::endl;
                size_t pos_a = curr_pos_a;
                for(size_t i=0;i<len_a;i++){
                    size_t sym = stream_a.read(pos_a, pos_a+a_width-1);
                    sym>>=1UL;
                    pos_a+=a_width;
                    std::cout<<sym<<" ";
                }
                std::cout<<""<<std::endl;

                size_t pos_b = curr_pos_b;
                for(size_t i=0;i<len_b;i++){
                    size_t sym = stream_b.read(pos_b, pos_b+b_width-1);
                    sym>>=1UL;
                    pos_b+=b_width;
                    std::cout<<sym<<" ";
                }
                std::cout<<""<<std::endl;
            }
            //

            if(eq_seq==0){
                //write rule from A
                append_rule(curr_pos_a, a_width, len_a, stream_a, mt_map_a, curr_pos_m, m_width, stream_c);
                curr_rule_a++;
                curr_rule_b++;
                curr_pos_b +=len_b*b_width;

                new_fps.push_back(fp_a);
                lvl_met_c.tot_symbols+=len_a;
                merge_marks.push_back(3);
            } else if(eq_seq<0) { // collision: break ties by lex. rank
                //write rule from A
                append_rule(curr_pos_a, a_width, len_a, stream_a, mt_map_a, curr_pos_m, m_width, stream_c);
                curr_rule_a++;

                new_fps.push_back(fp_a);
                lvl_met_c.tot_symbols+=len_a;
                merge_marks.push_back(1);
            } else{
                //write rule from B
                append_rule(curr_pos_b, b_width, len_b, stream_b, mt_map_b, curr_pos_m, m_width, stream_c);
                curr_rule_b++;

                new_fps.push_back(fp_b);
                lvl_met_c.tot_symbols+=len_b;
                merge_marks.push_back(2);
            }
        }

        if((merge_marks.back() & 1) && curr_rule_a<lvl_met_a.n_rules){
            std::tie(fp_a, len_a) = get_rule_info(stream_a, curr_pos_a, a_width, prev_fps, mt_map_a, fp_seed);
            assert(fp_a>=prev_fp_a);
            prev_fp_a = fp_a;
        }

        if((merge_marks.back() & 2) && curr_rule_b<lvl_met_b.n_rules){
            std::tie(fp_b, len_b) = get_rule_info(stream_b, curr_pos_b, b_width, prev_fps, mt_map_b, fp_seed);
            assert(fp_b>=prev_fp_b);
            prev_fp_b = fp_b;
        }
    }

    while(curr_rule_a<lvl_met_a.n_rules){
        //write rule from A
        append_rule(curr_pos_a, a_width, len_a, stream_a, mt_map_a, curr_pos_m, m_width, stream_c);
        curr_rule_a++;

        new_fps.push_back(fp_a);
        lvl_met_c.tot_symbols+=len_a;
        merge_marks.push_back(1);

        if(curr_rule_a<lvl_met_a.n_rules){
            std::tie(fp_a, len_a) = get_rule_info(stream_a, curr_pos_a, a_width, prev_fps, mt_map_a, fp_seed);
            assert(fp_a>=prev_fp_a);
            prev_fp_a = fp_a;
        }
    }

    while(curr_rule_b<lvl_met_b.n_rules){
        //write rule from B
        append_rule(curr_pos_b, b_width, len_b, stream_b, mt_map_b, curr_pos_m, m_width, stream_c);
        curr_rule_b++;

        new_fps.push_back(fp_b);
        lvl_met_c.tot_symbols+=len_b;
        merge_marks.push_back(2);

        if(curr_rule_b<lvl_met_b.n_rules){
            std::tie(fp_b, len_b) = get_rule_info(stream_b, curr_pos_b, b_width, prev_fps, mt_map_b, fp_seed);
            assert(fp_b>=prev_fp_b);
            prev_fp_b = fp_b;
        }
    }

    //update mapping values
    //it is one-based as the metasymbols are also one-based because of the separator symbol
    size_t mt_sym_a=1, mt_sym_b=1, mg_mt_sym=1;
    mt_map_a.resize(lvl_met_a.n_rules+1);
    mt_map_b.resize(lvl_met_b.n_rules+1);
    mt_map_a[0] = 0;
    mt_map_b[0] = 0;
    for(unsigned char merge_mark : merge_marks){
        if(merge_mark==1){
            mt_map_a[mt_sym_a++] = mg_mt_sym;
        } else if(merge_mark==2){
            mt_map_b[mt_sym_b++] = mg_mt_sym;
        } else {
            mt_map_a[mt_sym_a++] = mg_mt_sym;
            mt_map_b[mt_sym_b++] = mg_mt_sym;
        }
        mg_mt_sym++;
    }

    prev_fps.swap(new_fps);
    lvl_met_c.n_rules = merge_marks.size();
    prev_lvl_met = lvl_met_c;
}

template<class gram_type>
void create_fake_level(gram_type& p_gram, size_t new_lvl, std::vector<uint64_t>& prev_fps,
                       uint64_t fp_seed, std::vector<uint32_t>& mt_map){

    assert(new_lvl>=(p_gram.lvl-1));
    lvl_metadata_type new_lvl_mt{};

    // new_level is to the previous level in the metadata because
    // the first element of the metadata vector has the terminal alphabet
    new_lvl_mt.sym_width = sym_width(p_gram.rules_metadata[new_lvl].n_rules)+1;
    new_lvl_mt.n_rules = p_gram.rules_metadata[new_lvl].n_rules;
    new_lvl_mt.tot_symbols = p_gram.rules_metadata[new_lvl].n_rules;
    new_lvl_mt.terminals = false;

    assert((new_lvl_mt.n_rules+1) == mt_map.size());

    p_gram.rules[new_lvl].stream.reserve_in_bits(new_lvl_mt.n_bits());

    std::vector<std::pair<uint64_t, uint64_t>> perm(new_lvl_mt.n_rules);
    for(size_t i=1;i<=new_lvl_mt.n_rules;i++){
        perm[i].first = i;
        uint64_t fp = prev_fps[mt_map[i]];
        perm[i].second = XXH64(&fp, sizeof(uint64_t), fp_seed);
    }

    std::sort(perm.begin(), perm.end(), [&](auto const& a, auto const &b) -> bool{
        if(a.second!=b.second){
            return a.second<b.second;
        }
        assert(prev_fps[mt_map[a.first]]!=prev_fps[mt_map[b.first]]);
        return prev_fps[mt_map[a.first]]<prev_fps[mt_map[b.first]];
    });

    size_t pos=0, width=new_lvl_mt.sym_width;
    for(auto & i : perm){
        p_gram.rules[new_lvl].write(pos, pos+width-1, (i.first<<1UL | 1));
        pos+=width;
    }
    assert(pos==new_lvl_mt.n_bits());
    p_gram.rules_metadata[new_lvl+1] = new_lvl_mt;
}

template<class gram_type, class sym_type>
void merge_two_grammars(gram_type& p_gram_a, gram_type& p_gram_b, std::vector<uint64_t>& fp_seeds) {

    //get the information about the terminal symbols
    lvl_metadata_type prev_lvl_met = p_gram_a.metadata[0];

    //in the first level, tot_symbols represent the alphabet of terminals
    std::vector<uint32_t> mt_map_a;
    mt_map_a.resize(prev_lvl_met.tot_symbols);

    std::vector<uint32_t> mt_map_b;
    mt_map_b.resize(prev_lvl_met.tot_symbols);

    std::vector<uint64_t> prev_fps;
    prev_fps.resize(prev_lvl_met.tot_symbols);

    for(size_t i=0;i<prev_fps.size();i++){
        prev_fps[i] = XXH64(&i, sizeof(sym_type), fp_seeds[0]);
        mt_map_a[i] = i;
        mt_map_b[i] = i;
    }

    gram_type p_gram_c;

    //we subtract one because the last level contains the compressed strings, and they are not merged
    size_t max_lvl = std::max(p_gram_a.lvl, p_gram_b.lvl)-1;
    size_t min_lvl = std::min(p_gram_a.lvl, p_gram_b.lvl)-1;

    //pointer to the grammar that still has levels to process
    gram_type *st_gram = nullptr;
    std::vector<uint32_t> *st_gram_map= nullptr;
    if(p_gram_a.lvl<p_gram_b.lvl){
        st_gram = &p_gram_a;
        st_gram_map = &mt_map_a;
    } else if(p_gram_b.lvl<p_gram_a.lvl) {
        st_gram = &p_gram_b;
        st_gram_map = &mt_map_b;
    }

    if(st_gram!=nullptr){
        st_gram->rules.resize(max_lvl+1);
        st_gram->metadata.resize(max_lvl+2);
        st_gram->rules[min_lvl].swap(st_gram->rules[max_lvl]);
        st_gram->lvl = st_gram->rules.size();
        //+1 because the first element is the metadata of the terminal symbols
        std::swap(st_gram->metadata[min_lvl+1], st_gram->metadata[max_lvl+1]);
    }

    size_t i=0;
    while(i<min_lvl) {

        //std::cout<<p_gram_a.rules_metadata[i].n_rules<<" "<<mt_map_a.size()<<std::endl;
        //std::cout<<p_gram_b.rules_metadata[i].n_rules<<" "<<mt_map_b.size()<<std::endl;

        merge_level(p_gram_a.rules[i], p_gram_a.metadata[i+1], mt_map_a,
                    p_gram_b.rules[i], p_gram_b.metadata[i+1], mt_map_b,
                    p_gram_c.rules[i], prev_lvl_met,
                    fp_seeds[i+1], prev_fps);

        p_gram_c.metadata.push_back(prev_lvl_met);

        std::cout<<"level "<<i<<" "<<min_lvl<<" "<<max_lvl<<std::endl;
        i++;

        /*if(i>=min_lvl && i<max_lvl){
            //TODO move the compressed string to the back
            create_fake_level(*st_gram, i, prev_fps, fp_seeds[i+1], *st_gram_map);
        }*/
    }

    //TODO checking
    std::cout<<"Merged grammar"<<std::endl;
    p_gram_c.lvl = p_gram_c.metadata.size()-1;
    p_gram_c.print_stats();
    //
}

#endif //LCG_PARTIAL_GRAM_H
