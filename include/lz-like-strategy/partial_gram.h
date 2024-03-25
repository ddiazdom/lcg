//
// Created by Diaz, Diego on 23.10.2023.
//

#ifndef LCG_PARTIAL_GRAM_H
#define LCG_PARTIAL_GRAM_H
#include "cds/bitstream.h"
#include "lz_like_map.h"

struct merge_data_t{

    std::vector<uint64_t> fps;
    std::vector<uint64_t> new_fps;
    std::vector<uint32_t> map_a;
    std::vector<uint32_t> map_b;
    size_t lvl_sigma=0;
    bitstream<size_t, true> buffer;
    std::vector<uint8_t> merge_marks;

    void initialize(size_t lvl_sigma_, size_t sym_bytes, uint64_t seed){
        lvl_sigma = lvl_sigma_;
        map_a.resize(lvl_sigma);
        map_b.resize(lvl_sigma);
        fps.resize(lvl_sigma);
        for(size_t i=0;i<fps.size();i++){
            fps[i] = XXH64(&i, sym_bytes, seed);
            map_a[i] = i;
            map_b[i] = i;
        }
    }

    [[nodiscard]] size_t space_usage() const {
        size_t bytes = fps.capacity()*sizeof(uint64_t);
        bytes += new_fps.capacity()*sizeof(uint64_t);
        bytes += map_a.capacity()*sizeof(uint32_t);
        bytes += map_b.capacity()*sizeof(uint32_t);
        bytes += merge_marks.capacity()*sizeof(uint8_t);
        return bytes;
    }

    ~merge_data_t(){
        destroy(fps);
        destroy(new_fps);
        destroy(map_a);
        destroy(map_b);
        buffer.destroy();
        destroy(merge_marks);
#ifdef __linux__
        malloc_trim(0);
#endif
    }
};

struct lvl_metadata_type{
    size_t n_rules=0;
    size_t tot_symbols=0;
    uint8_t sym_width=0;
    bool terminals=false;

    [[nodiscard]] inline size_t n_bits() const {
        return tot_symbols*sym_width;
    }
};

template<class ter_type, bool mmap_allocator=false>
struct partial_gram {

    typedef ter_type                            sym_type;
    typedef bitstream<size_t, mmap_allocator> stream_type;

    size_t text_size=0;
    ter_type max_tsym = std::numeric_limits<ter_type>::max();
    ter_type sep_tsym = 0;
    size_t lvl = 0;
    size_t longest_rule=0;
    std::vector<lvl_metadata_type> metadata;
    std::vector<stream_type> rules;

    partial_gram(){
        lvl_metadata_type ter{};
        ter.sym_width = std::numeric_limits<ter_type>::digits;
        ter.tot_symbols = std::numeric_limits<ter_type>::max()+1;
        ter.n_rules = ter.tot_symbols;
        ter.terminals = true;
        metadata.push_back(ter);
    }

    partial_gram& swap(partial_gram& other){
        rules.swap(other.rules);
        metadata.swap(other.rules_metadata);
        std::swap(lvl, other.lvl);
        return *this;
    }

    size_t serialize(std::ofstream &ofs){
        assert(lvl==rules.size());
        size_t written_bytes = serialize_elm(ofs, text_size);
        written_bytes += serialize_elm(ofs, max_tsym);
        written_bytes += serialize_elm(ofs, sep_tsym);
        written_bytes += serialize_elm(ofs, lvl);
        written_bytes += serialize_elm(ofs, longest_rule);

        written_bytes+= serialize_plain_vector(ofs, metadata);

        size_t n_words;
        for(size_t i=0;i<lvl;i++){
            n_words = rules[i].bits2words(metadata[i+1].n_bits());
            ofs.write((char *)rules[i].stream,  rules[i].words2bytes(n_words));
            written_bytes += rules[i].words2bytes(n_words);
        }
        return written_bytes;
    }

    void load_metadata(std::ifstream &ifs){
        load_elm(ifs, text_size);
        load_elm(ifs, max_tsym);
        load_elm(ifs, sep_tsym);
        load_elm(ifs, lvl);
        load_elm(ifs, longest_rule);
        load_plain_vector(ifs, metadata);
    }

    void load_next_rule_set(std::ifstream &ifs, size_t i, bitstream<size_t>& buffer){
        size_t n_words =  buffer.bits2words(metadata[i+1].n_bits());;
        buffer.reserve_in_words(n_words);
        assert(buffer.stream!= nullptr);
        ifs.read((char *)buffer.stream, long(buffer.words2bytes(n_words)));
    }

    void load(std::ifstream &ifs){
        load_metadata(ifs);
        rules.resize(lvl);
        for(size_t i=0;i<lvl;i++){
            rules[i].load(ifs);
        }
    }

    size_t serialize_to_fd(int fd){
        assert(lvl==(metadata.size()-1));

        write(fd, &text_size, sizeof(size_t));
        write(fd, &max_tsym, sizeof(size_t));
        write(fd, &sep_tsym, sizeof(size_t));
        write(fd, &lvl, sizeof(size_t));
        write(fd, &longest_rule, sizeof(size_t));
        size_t written_bytes = 5*sizeof(size_t);

        write(fd, metadata.data(), metadata.size()*sizeof(lvl_metadata_type));
        written_bytes += metadata.size()*sizeof(lvl_metadata_type);

        size_t n_words;
        for(size_t i=0;i<lvl;i++){
            n_words = rules[i].bits2words(metadata[i+1].n_bits());
            write(fd, rules[i].stream,  rules[i].words2bytes(n_words));
            written_bytes += rules[i].words2bytes(n_words);
        }
        return written_bytes;
    }

    size_t load_from_fd(int fd){
        read(fd, &text_size, sizeof(size_t));
        read(fd, &max_tsym, sizeof(size_t));
        read(fd, &sep_tsym, sizeof(size_t));
        read(fd, &lvl, sizeof(size_t));
        read(fd, &longest_rule, sizeof(size_t));
        size_t read_bytes = 5*sizeof(size_t);

        if(rules.size()>lvl){
            for(size_t i=lvl;i<rules.size();i++){
                rules[i].destroy();
            }
        }

        metadata.resize(lvl+1);
        read(fd, metadata.data(), metadata.size()*sizeof(lvl_metadata_type));
        read_bytes+=metadata.size()*sizeof(lvl_metadata_type);

        rules.resize(lvl);
        size_t n_words;
        for(size_t i=0;i<lvl;i++){
            n_words = rules[i].bits2words(metadata[i+1].n_bits());
            rules[i].reserve_in_words(n_words);
            assert(rules[i].stream!= nullptr);
            read(fd, rules[i].stream, rules[i].words2bytes(n_words));
            read_bytes+=rules[i].words2bytes(n_words);
        }
        return read_bytes;
    }

    template<class sym_type>
    void append_new_lvl(sym_type* text, const lz_like_map::phrase_list_t &phrase_set, size_t tot_symbols,
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
            if(len>longest_rule) longest_rule = len;
            uint32_t last = source + len-1;
            for(size_t j=source;j<last;j++){
                assert(text[j]>0);
                rules[lvl].write(acc_bits, acc_bits+lvl_met.sym_width-1, text[j]<<1UL);
                acc_bits+=lvl_met.sym_width;
            }

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
            assert(text[j]>0 && text[j+1]==0);
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
        lvl=0;
        metadata.resize(1);
    }

    [[nodiscard]] inline size_t txt_size() const {
        return text_size;
    }

    [[nodiscard]] inline size_t gram_size() const {
        size_t g_size=0;
        for(const auto & i : metadata){
            g_size+=i.tot_symbols;
        }
        return g_size;
    }

    [[nodiscard]] inline size_t tot_gram_symbols() const {
        size_t n_rules=0;
        for(const auto & i : metadata){
           n_rules+=i.n_rules;
        }
        return n_rules;
    }

    [[nodiscard]] inline size_t tot_strings() const {
        return metadata.back().tot_symbols;
    }

    [[nodiscard]] inline ter_type max_terminal_symbol() const {
        return max_tsym;
    }

    [[nodiscard]] inline ter_type separator_symbol() const {
        return sep_tsym;
    }

    ~partial_gram(){
        for(auto & lvl_rules : rules){
            lvl_rules.destroy();
        }
        std::vector<lvl_metadata_type>().swap(metadata);
    }

    size_t space_usage(){
        size_t n_bytes=0;
        for(auto & lvl_rules : rules){
            n_bytes +=lvl_rules.capacity_in_bytes();
        }
        n_bytes+=metadata.size()+sizeof(lvl_metadata_type);
        return n_bytes;
    }

    size_t gram_size_in_bytes(){
        size_t n_bits=0;
        for(auto & mt : metadata){
            n_bits +=mt.n_bits();
        }
        return (INT_CEIL(n_bits, 8));
    }

    void print_stats(){
        for(size_t i=0;i<metadata.size();i++){
            if(i==0){
                std::cout<<"Number of terminals: "<<metadata[i].tot_symbols<<std::endl;
            }else if(i<(metadata.size()-1)){
                std::cout<<"Level "<<i+1<<": number of rules: "<<metadata[i].n_rules<<"  number of symbols: "<<metadata[i].tot_symbols<<std::endl;
            }else{
                std::cout<<"Compressed string: number of strings: "<<metadata[i].tot_symbols<<std::endl;
            }
        }
    }
};

template<class stream_type>
inline void get_rule_info(stream_type& rule_stream, size_t& pos, size_t width,
                          std::vector<uint64_t>& prev_fps, std::vector<uint32_t>& mt_map, size_t fp_seed,
                          std::vector<uint64_t>& phrase_fp_buff, std::vector<uint64_t>& phrase_buff,
                          uint64_t& fingerprint, size_t& len){
    size_t sym;
    len=0;

    do{
        sym = rule_stream.read(pos, pos+width-1);
        size_t tmp = sym>>1UL;
        //assert(tmp<mt_map.size());
        tmp = mt_map[tmp];
        //assert(tmp<prev_fps.size());
        phrase_buff[len] = tmp;
        phrase_fp_buff[len] = prev_fps[tmp];
        pos+=width;
        len++;
    } while(!(sym & 1UL));

    fingerprint = XXH64(phrase_fp_buff.data(), sizeof(uint64_t)*len, fp_seed);
}

bool compare_rules(std::vector<uint64_t>& phrase_a, size_t len_a, std::vector<uint64_t>& phrase_b, size_t len_b){
    size_t n_comparisons = std::min(len_a, len_b);
    for(size_t i=0;i<n_comparisons;i++){
        if(phrase_a[i]!=phrase_b[i]){
            return phrase_a[i]<phrase_b[i];
        }
    }
    return len_a<len_b;
}

template<class stream_type>
void append_rule(std::vector<uint64_t>& s_phrase, size_t s_len, size_t& d_pos, size_t d_width, stream_type& dest){

    //check the appended rule fits the buffer of merged rules
    size_t min_bits = d_pos+(s_len*d_width);
    if(dest.capacity_in_bits()<=min_bits) {
        min_bits = INT_CEIL((min_bits * 120), 100);//20% increase in the stream size
        dest.reserve_in_bits(min_bits);
    }

    size_t last=s_len-1;
    for(size_t i=0;i<s_len;i++){
        dest.write(d_pos, d_pos+d_width-1, (s_phrase[i]<<1UL | (i==last)));
        d_pos+=d_width;
    }
}

template<class stream_type>
lvl_metadata_type merge_level(stream_type &stream_a, lvl_metadata_type &lvl_met_a,
                              stream_type &stream_b, lvl_metadata_type &lvl_met_b,
                              uint64_t &fp_seed, merge_data_t& mg_data, size_t longest_rule) {


    uint64_t fp_a, fp_b;
    size_t curr_rule_a=0, curr_rule_b=0;
    size_t curr_pos_a=0, curr_pos_b=0;
    size_t len_a, len_b;

    lvl_metadata_type lvl_met_c{};
    lvl_met_c.sym_width = sym_width(mg_data.lvl_sigma)+1;//we use the extra bit to mark the end of each phrase in the stream of rules
    lvl_met_c.terminals = lvl_met_a.terminals;

    uint8_t m_width = lvl_met_c.sym_width;
    size_t curr_pos_m=0;


    std::cout<<"Buffer: "<<report_space(mg_data.buffer.capacity_in_bytes())<<" A:"<<report_space(stream_a.capacity_in_bytes())<<" B:"<<report_space(stream_b.capacity_in_bytes())<<std::endl;

    mg_data.buffer.reserve_in_bits(std::max(lvl_met_a.n_bits(), lvl_met_b.n_bits()));

    mg_data.merge_marks.clear();
    mg_data.merge_marks.reserve(lvl_met_a.n_rules);

    mg_data.new_fps.clear();
    mg_data.new_fps.reserve(lvl_met_a.n_rules);
    mg_data.new_fps.push_back(0);//fake

    uint8_t a_width = lvl_met_a.sym_width;
    uint8_t b_width = lvl_met_b.sym_width;

    std::vector<uint64_t> phrase_fp_buff(longest_rule, 0);
    std::vector<uint64_t> phrase_a(longest_rule, 0);
    std::vector<uint64_t> phrase_b(longest_rule, 0);

    get_rule_info(stream_a, curr_pos_a, a_width, mg_data.fps, mg_data.map_a, fp_seed, phrase_fp_buff,
                  phrase_a, fp_a, len_a);
    assert(len_a<=longest_rule);
    get_rule_info(stream_b, curr_pos_b, b_width, mg_data.fps, mg_data.map_b, fp_seed, phrase_fp_buff,
                  phrase_b, fp_b, len_b);
    assert(len_b<=longest_rule);

    uint64_t prev_fp_a = fp_a;
    uint64_t prev_fp_b = fp_b;

    while(curr_rule_a<lvl_met_a.n_rules && curr_rule_b<lvl_met_b.n_rules){

        if(fp_a<fp_b){
            //write rule from A
            append_rule(phrase_a, len_a, curr_pos_m, m_width, mg_data.buffer);
            curr_rule_a++;

            mg_data.new_fps.push_back(fp_a);
            lvl_met_c.tot_symbols+=len_a;
            mg_data.merge_marks.push_back(1);
        } else if(fp_b<fp_a) {
            //write rule from B
            append_rule(phrase_b, len_b, curr_pos_m, m_width, mg_data.buffer);
            curr_rule_b++;

            mg_data.new_fps.push_back(fp_b);
            lvl_met_c.tot_symbols+=len_b;
            mg_data.merge_marks.push_back(2);
        } else {

            bool eq = len_a==len_b && (memcmp(phrase_a.data(), phrase_b.data(), len_a*sizeof(uint64_t))==0);

            if(eq){
                //write rule from A
                append_rule(phrase_a, len_a, curr_pos_m, m_width, mg_data.buffer);
                curr_rule_a++;
                curr_rule_b++;

                mg_data.new_fps.push_back(fp_a);
                lvl_met_c.tot_symbols+=len_a;
                mg_data.merge_marks.push_back(3);
            }else {

                //a collision occurred (extremely unlikely, but not impossible)
                std::cout << "Collision warning:  " << fp_a << " " << fp_b << " " << curr_pos_a << " " << curr_pos_b<< std::endl;
                // break ties by lex. rank
                bool a_is_lex_smaller = compare_rules(phrase_a, len_a, phrase_b, len_b);

                if(a_is_lex_smaller) {
                    //write rule from A
                    append_rule(phrase_a, len_a, curr_pos_m, m_width, mg_data.buffer);
                    curr_rule_a++;

                    mg_data.new_fps.push_back(fp_a);
                    lvl_met_c.tot_symbols += len_a;
                    mg_data.merge_marks.push_back(1);
                } else {
                    //write rule from B
                    append_rule(phrase_b, len_b, curr_pos_m, m_width, mg_data.buffer);
                    curr_rule_b++;

                    mg_data.new_fps.push_back(fp_b);
                    lvl_met_c.tot_symbols += len_b;
                    mg_data.merge_marks.push_back(2);
                }
            }
            //
        }

        if((mg_data.merge_marks.back() & 1) && curr_rule_a<lvl_met_a.n_rules){
            get_rule_info(stream_a, curr_pos_a, a_width, mg_data.fps, mg_data.map_a, fp_seed, phrase_fp_buff,
                          phrase_a, fp_a, len_a);
            assert(fp_a>=prev_fp_a);
            prev_fp_a = fp_a;
        }

        if((mg_data.merge_marks.back() & 2) && curr_rule_b<lvl_met_b.n_rules){
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

        mg_data.new_fps.push_back(fp_a);
        lvl_met_c.tot_symbols+=len_a;
        mg_data.merge_marks.push_back(1);

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

        mg_data.new_fps.push_back(fp_b);
        lvl_met_c.tot_symbols+=len_b;
        mg_data.merge_marks.push_back(2);

        if(curr_rule_b<lvl_met_b.n_rules){
            get_rule_info(stream_b, curr_pos_b, b_width, mg_data.fps, mg_data.map_b, fp_seed, phrase_fp_buff,
                          phrase_b, fp_b, len_b);
            assert(fp_b>=prev_fp_b);
            prev_fp_b = fp_b;
        }
    }

    //update mapping values
    //the new metasymbols are one-based to differentiate them from the separator symbol (0) is the parse
    size_t mt_sym_a=1, mt_sym_b=1, mg_mt_sym=1;
    mg_data.map_a.resize(lvl_met_a.n_rules+1);
    mg_data.map_b.resize(lvl_met_b.n_rules+1);
    mg_data.map_a[0] = 0;
    mg_data.map_b[0] = 0;
    for(unsigned char merge_mark : mg_data.merge_marks){
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

    mg_data.fps.swap(mg_data.new_fps);
    lvl_met_c.n_rules = mg_data.merge_marks.size();
    mg_data.lvl_sigma = lvl_met_c.n_rules;//store this information for the next round
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

template<class stream_type>
lvl_metadata_type concatenate_strings(stream_type &stream_a, lvl_metadata_type &lvl_met_a,
                                      stream_type &stream_b, lvl_metadata_type &lvl_met_b, merge_data_t& mg_data){

    std::cout<<"Buffer: "<<report_space(mg_data.buffer.capacity_in_bytes())<<" A:"<<report_space(stream_a.capacity_in_bytes())<<" B:"<<report_space(stream_b.capacity_in_bytes())<<std::endl;

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

template<class gram_type>
void merge_two_grammars(gram_type& p_gram_a, gram_type& p_gram_b, std::vector<uint64_t>& fp_seeds, merge_data_t& mg_data) {

    //in the first level, tot_symbols represent the alphabet of terminals
    mg_data.initialize(p_gram_a.metadata[0].tot_symbols, sizeof(typename gram_type::sym_type), fp_seeds[0]);

    //we subtract one because the last level contains the compressed strings,
    // and we do not merge but concatenate them
    size_t max_lvl = std::max(p_gram_a.lvl, p_gram_b.lvl)-1;
    size_t min_lvl = std::min(p_gram_a.lvl, p_gram_b.lvl)-1;

    //pointer to the grammar with the least number of levels
    gram_type *st_gram = nullptr;
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

    size_t longest_rule = std::max(p_gram_a.longest_rule, p_gram_b.longest_rule);

    lvl_metadata_type buffer_metadata;
    size_t i=0;
    while(i<max_lvl) {
        buffer_metadata = merge_level(p_gram_a.rules[i], p_gram_a.metadata[i+1],
                                      p_gram_b.rules[i], p_gram_b.metadata[i+1],
                                      fp_seeds[i+1], mg_data, longest_rule);
        i++;
        if(i>=min_lvl && i<max_lvl){
            assert(st_gram!=nullptr && st_gram!= nullptr);
            create_fake_level(*st_gram, i, mg_data.fps, fp_seeds[i+1], *st_gram_map);
        }

        p_gram_a.metadata[i] = buffer_metadata;
        //p_gram_a.rules[i-1].swap(mg_data.buffer);
        mg_data.buffer.copy(p_gram_a.metadata[i].n_bits(), p_gram_a.rules[i-1]);
        std::cout<<"Extra data is using"<<report_space(off_t(mg_data.space_usage()))<<std::endl;
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

#endif //LCG_PARTIAL_GRAM_H
