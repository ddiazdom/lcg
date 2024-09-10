//
// Created by Diaz, Diego on 23.10.2023.
//

#ifndef LCG_PARTIAL_GRAM_H
#define LCG_PARTIAL_GRAM_H

#include "cds/bitstream.h"
#include "lz-like-build-strat/lz_like_map.h"

struct lvl_metadata_type{
    size_t n_rules=0;
    size_t tot_symbols=0;
    uint8_t sym_width=0;
    bool terminals=false;

    [[nodiscard]] inline size_t n_bits() const {
        return tot_symbols*sym_width;
    }

    [[nodiscard]] inline size_t uint32_bytes() const {
        return tot_symbols*4;
    }
};

struct str_perm_val{
    size_t txt_id;
    size_t first_str;
    size_t last_str;
};

template<class ter_type, bool mmap_allocator=false>
struct partial_gram {

    typedef ter_type                            sym_type;
    typedef bitstream<size_t, mmap_allocator> stream_type;

    size_t text_size=0;
    size_t txt_id=0;
    ter_type max_tsym = std::numeric_limits<ter_type>::max();
    ter_type sep_tsym = 0;
    size_t lvl = 0;
    size_t longest_rule=0;
    size_t longest_str=0;
    uint64_t par_seed=0;
    std::vector<lvl_metadata_type> metadata;
    std::vector<stream_type> rules;

    partial_gram(): text_size(0),
                    txt_id(0),
                    max_tsym( std::numeric_limits<ter_type>::max()),
                    sep_tsym(0),
                    lvl(0),
                    longest_rule(0),
                    longest_str(0),
                    par_seed(0){
        lvl_metadata_type ter{};
        ter.sym_width = std::numeric_limits<ter_type>::digits;
        ter.tot_symbols = std::numeric_limits<ter_type>::max()+1;
        ter.n_rules = ter.tot_symbols;
        ter.terminals = true;
        metadata.push_back(ter);
    }

    partial_gram& swap(partial_gram& other){
        std::swap(text_size, other.text_size);
        std::swap(txt_id, other.txt_id);
        std::swap(max_tsym, other.max_tsym);
        std::swap(sep_tsym, other.sep_tsym);
        std::swap(lvl, other.lvl);
        std::swap(longest_rule, other.longest_rule);
        std::swap(longest_str, other.longest_str);
        std::swap(par_seed, other.par_seed);
        metadata.swap(other.rules_metadata);
        rules.swap(other.rules);
        return *this;
    }

    size_t serialize(std::ofstream &ofs){
        assert(lvl==rules.size());
        size_t written_bytes = serialize_elm(ofs, text_size);
        written_bytes += serialize_elm(ofs, txt_id);
        written_bytes += serialize_elm(ofs, max_tsym);
        written_bytes += serialize_elm(ofs, sep_tsym);
        written_bytes += serialize_elm(ofs, lvl);
        written_bytes += serialize_elm(ofs, longest_rule);
        written_bytes += serialize_elm(ofs, longest_str);
        written_bytes += serialize_elm(ofs, par_seed);

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
        load_elm(ifs, txt_id);
        load_elm(ifs, max_tsym);
        load_elm(ifs, sep_tsym);
        load_elm(ifs, lvl);
        load_elm(ifs, longest_rule);
        load_elm(ifs, longest_str);
        load_elm(ifs, par_seed);
        load_plain_vector(ifs, metadata);
    }

    size_t bytes_metadata(){
        size_t bytes=0;
        bytes+=sizeof(text_size);
        bytes+=sizeof(txt_id);
        bytes+=sizeof(max_tsym);
        bytes+=sizeof(sep_tsym);
        bytes+=sizeof(lvl);
        bytes+=sizeof(longest_rule);
        bytes+=sizeof(longest_str);
        bytes+=sizeof(par_seed);

        bytes+=metadata.size()*sizeof(lvl_metadata_type);
        bytes+=sizeof(size_t);
        return bytes;
    }

    void load_concat_string(std::ifstream &ifs, bitstream<size_t>& buffer){
        size_t n_words=0;
        for(size_t i=1;i<metadata.size()-1;i++){
            n_words +=  buffer.bits2words(metadata[i].n_bits());
        }
        size_t bytes = buffer.words2bytes(n_words);
        bytes+=ifs.tellg();
        ifs.seekg(long(bytes));

        buffer.reserve_in_words(n_words);
        assert(buffer.stream!= nullptr);
        ifs.read((char *)buffer.stream, long(buffer.words2bytes(n_words)));
    }

    /*void reset_file_to_first_level(std::ifstream &ifs){
        size_t n_words =  buffer.bits2words(metadata[i+1].n_bits());;
        buffer.reserve_in_words(n_words);
        assert(buffer.stream!= nullptr);
        ifs.read((char *)buffer.stream, long(buffer.words2bytes(n_words)));
    }*/

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
        write(fd, &txt_id, sizeof(size_t));
        write(fd, &max_tsym, sizeof(size_t));
        write(fd, &sep_tsym, sizeof(size_t));
        write(fd, &lvl, sizeof(size_t));
        write(fd, &longest_rule, sizeof(size_t));
        write(fd, &longest_str, sizeof(size_t));
        write(fd, &par_seed, sizeof(uint64_t));
        size_t written_bytes = 8*sizeof(size_t);

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
        read(fd, &txt_id, sizeof(size_t));
        read(fd, &max_tsym, sizeof(size_t));
        read(fd, &sep_tsym, sizeof(size_t));
        read(fd, &lvl, sizeof(size_t));
        read(fd, &longest_rule, sizeof(size_t));
        read(fd, &longest_str, sizeof(size_t));
        read(fd, &par_seed, sizeof(uint64_t));
        size_t read_bytes = 8*sizeof(size_t);

        for(size_t i=lvl;i<rules.size();i++){
            rules[i].destroy();
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
    off_t append_new_lvl(sym_type* text, const typename lz_like_map<sym_type>::phrase_list_t &phrase_set, size_t tot_symbols,
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
            uint32_t source = phrase_set[idx].source;
            uint32_t len = phrase_set[idx].len;
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
        return rules[lvl-1].capacity_in_bytes();
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

    //string_ranges encodes the original order of the strings that were merged in this partial grammar
    //the pair is (txt_chunk_id, n_str):
    void reorder_strings(std::vector<std::pair<size_t, size_t>>& string_ranges){

        stream_type permuted_strings;
        permuted_strings.reserve_in_bits(metadata.back().n_bits());

        std::vector<str_perm_val> sorted_str_ranges(string_ranges.size());
        size_t acc_strings=0;
        for(size_t i=0;i<string_ranges.size();i++){
            sorted_str_ranges[i].txt_id = string_ranges[i].first;
            sorted_str_ranges[i].first_str = acc_strings;
            sorted_str_ranges[i].last_str = acc_strings+string_ranges[i].second-1;
            acc_strings+=string_ranges[i].second;
        }

        std::sort(sorted_str_ranges.begin(), sorted_str_ranges.end(), [](auto const& a, auto const&b){
            return a.txt_id<b.txt_id;
        });

        size_t pos =0;
        size_t width = metadata.back().sym_width;
        size_t str_sym;
        for(auto const& str_range : sorted_str_ranges){
            for(size_t str=str_range.first_str;str<=str_range.last_str;str++){
                str_sym = rules[lvl-1].read(str*width, (str+1)*width-1);
                permuted_strings.write(pos*width, (pos+1)*width-1, str_sym);
                pos++;
            }
        }
        assert(pos==metadata.back().tot_symbols);
        rules[lvl-1].swap(permuted_strings);
    }

    void reset_grammar(){
        text_size=0;
        txt_id=0;
        sep_tsym=0;
        lvl=0;
        longest_rule=0;
        longest_str=0;
        metadata.resize(1);
    }

    void destroy_gram(){
        for(auto & lvl_rules : rules){
            lvl_rules.destroy();
        }
        destroy(metadata);
#ifdef __linux__
        malloc_trim(0);
#endif
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

    [[nodiscard]] inline size_t gram_uint32_bytes() const {
        size_t bytes=0;
        for(const auto & i : metadata){
            bytes+=i.uint32_bytes();
        }
        return bytes;
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

    [[nodiscard]] inline size_t longest_string() const {
        return longest_str;
    }

    [[nodiscard]] inline ter_type max_terminal_symbol() const {
        return max_tsym;
    }

    [[nodiscard]] inline ter_type separator_symbol() const {
        return sep_tsym;
    }

    ~partial_gram(){
        destroy_gram();
    }

    [[nodiscard]] size_t gram_size_in_bytes() const {
        size_t n_bits=0;
        for(auto & mt : metadata){
            n_bits +=mt.n_bits();
        }
        n_bits += metadata.size()*sizeof(lvl_metadata_type);
        return (INT_CEIL(n_bits, 8));
    }

    /*void print_stats(){
        for(size_t i=0;i<metadata.size();i++){
            if(i==0){
                std::cout<<"Number of terminals: "<<metadata[i].tot_symbols<<std::endl;
            }else if(i<(metadata.size()-1)){
                std::cout<<"Level "<<i+1<<": number of rules: "<<metadata[i].n_rules<<"  number of symbols: "<<metadata[i].tot_symbols<<std::endl;
            }else{
                std::cout<<"Compressed string: number of strings: "<<metadata[i].tot_symbols<<std::endl;
            }
        }
    }*/

    void stats(size_t pad) const {

        size_t g = gram_size();
        size_t s = tot_strings();
        size_t n = text_size;
        size_t g_syms = tot_gram_symbols();
        uint8_t r_bits = sym_width(g_syms);
        size_t n_ter = max_terminal_symbol()+1;
        size_t n_nter = g_syms - n_ter;
        auto g_bytes = INT_CEIL(g * r_bits, 8);

        float comp_ratio = float(n) / float(g_bytes);

        std::string pad_string(pad, ' ');
        std::cout << pad_string << "Parsing seed:                 " << par_seed<<std::endl;
        std::cout << pad_string << "Number of compressed symbols: " << n << " (" << report_space((off_t) n) << ")"<< std::endl;
        std::cout << pad_string << "Number of compressed strings: " << s <<std::endl;
        std::cout << pad_string << "Separator symbol:             " << (int) sep_tsym << std::endl;
        std::cout << pad_string << "Number of terminals:          " << n_ter << std::endl;
        std::cout << pad_string << "Number of nonterminals:       " << n_nter<< std::endl;
        std::cout << pad_string << "Grammar size:                 " << g <<std::endl;
        std::cout << pad_string << "Space usage partial grammar:  " << report_space((off_t) gram_size_in_bytes())<< std::endl;
        std::cout << pad_string << "Approx. compression ratio:    " << comp_ratio << std::endl;
    }

    void breakdown(size_t pad) {
        std::cout<<"Breakdown of the partial gram"<<std::endl;
        stats(pad);
        std::string pad_string(pad, ' ');
        size_t n_rules;
        std::cout << pad_string << "Grammar rules per level" << std::endl;
        for (size_t i = 0; i < metadata.size() - 1; i++) {
            n_rules = metadata[i].n_rules;
            if (n_rules > 0) {
                std::cout << pad_string << "  Level " << (i + 1) << ": number of rules: " << n_rules << std::endl;
            }
        }
    }
};

uint64_t get_par_seed_par_gram(std::string& p_gram_file){
    partial_gram<uint8_t> p_gram;
    std::ifstream ifs(p_gram_file);
    p_gram.load_metadata(ifs);
    return p_gram.par_seed;
}

size_t vbyte_size(partial_gram<uint8_t>& p_gram){

    size_t pos, n_bits;
    uint8_t width;
    uint32_t sym;
    size_t n_bytes=0;

    for(size_t i=0;i<p_gram.lvl;i++) {

        pos = 0;
        width = p_gram.metadata[i+1].sym_width;
        n_bits = p_gram.metadata[i+1].n_bits();
        bool last;

        while(pos<n_bits){
            sym = p_gram.rules[i].read(pos, pos+width-1);
            last = sym & 1UL;
            sym = sym>>1;
            n_bytes+= vbyte_len(sym);
            pos+=width;

            if(last){
                n_bytes++;
            }
        }
        assert(pos==n_bits);
    }
    return n_bytes;
}

void get_breakdown(std::string& p_gram_file){
    partial_gram<uint8_t> p_gram;
    std::ifstream ifs(p_gram_file);
    p_gram.load_metadata(ifs);
    p_gram.breakdown(2);
}

template<class stream_type>
inline void get_rule_info(stream_type& rule_stream, size_t& pos, size_t width,
                          std::vector<uint64_t>& prev_fps, std::vector<uint32_t>& mt_map, size_t fp_seed,
                          std::vector<uint64_t>& fp_seq, std::vector<uint32_t>& phrase,
                          uint64_t& fp, size_t& len){
    size_t sym;
    len=0;

    do{
        sym = rule_stream.read(pos, pos+width-1);
        phrase[len] = mt_map[sym>>1UL];
        pos+=width;
        len++;
    } while(!(sym & 1UL));

    for(size_t i=0;i<len;i++){
        fp_seq[i] = prev_fps[phrase[i]];
    }

    fp = XXH64(fp_seq.data(), sizeof(uint64_t)*len, fp_seed);
}

bool rules_lex_comp(std::vector<uint32_t>& phrase_a, size_t len_a, std::vector<uint32_t>& phrase_b, size_t len_b){
    size_t n_comparisons = std::min(len_a, len_b);
    for(size_t i=0;i<n_comparisons;i++){
        if(phrase_a[i]!=phrase_b[i]){
            return phrase_a[i]<phrase_b[i];
        }
    }
    return len_a<len_b;
}

template<class stream_type>
void append_rule(std::vector<uint32_t>& s_phrase, size_t s_len, size_t& d_pos, size_t d_width, stream_type& dest){

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

template<class comp_gram_type>
void partial2complete_gram(comp_gram_type& new_gram, std::string& p_gram_file, uint64_t par_seed_){

    partial_gram<uint8_t> p_gram;
    std::ifstream ifs(p_gram_file, std::ios::binary);
    p_gram.load_metadata(ifs);

    new_gram.n = p_gram.txt_size();
    new_gram.r = p_gram.tot_gram_symbols();
    new_gram.g = p_gram.gram_size();
    new_gram.s = p_gram.tot_strings();
    new_gram.c = p_gram.tot_strings();
    new_gram.par_seed = par_seed_;

    new_gram.max_tsym = size_t(p_gram.max_terminal_symbol());
    new_gram.sep_tsym = p_gram.separator_symbol();

    new_gram.lvl_rules.reserve(p_gram.metadata.size()-1);
    size_t acc=new_gram.max_tsym+1, tmp;
    //-1 because the last element of lvl_rules is the compressed string
    for(size_t i=1;i<p_gram.metadata.size()-1;i++){
        tmp = p_gram.metadata[i].n_rules;
        new_gram.lvl_rules.push_back(acc);
        acc+=tmp;
    }
    new_gram.lvl_rules.push_back(acc);

    off_t r_bits = sym_width(new_gram.r);
    new_gram.r_bits = r_bits;
    new_gram.rule_stream.reserve_in_bits(new_gram.r_bits*new_gram.g);

    new_gram.rl_ptr.set_width(sym_width(new_gram.g));
    new_gram.rl_ptr.resize(new_gram.r-new_gram.max_tsym);

    size_t bit_pos =0;
    new_gram.terminals.resize(new_gram.max_tsym+1);
    for(size_t ter=0;ter<=new_gram.max_tsym;ter++){
        new_gram.rule_stream.write(bit_pos, bit_pos+r_bits-1, ter);
        bit_pos+=r_bits;
        new_gram.terminals[ter] = ter;
    }

    bool last;
    size_t rule=0, pos, width, n_bits, rule_start_ptr=bit_pos, sym, min_sym=0, max_sym=new_gram.max_tsym;
    bitstream<size_t> rules_buffer;

    for(size_t i=0;i<p_gram.lvl;i++) {

        p_gram.load_next_rule_set(ifs, i, rules_buffer);
        pos = 0;
        width = p_gram.metadata[i+1].sym_width;
        n_bits = p_gram.metadata[i+1].n_bits();

        while(pos<n_bits){
            sym = rules_buffer.read(pos, pos+width-1);
            //std::cout<<(sym>>1)<<" ";
            last = sym & 1UL;
            sym = min_sym+(sym>>1);//sym is one-based
            //new_gram.rules.write(j, sym);
            new_gram.rule_stream.write(bit_pos, bit_pos+r_bits-1, sym);
            bit_pos+=r_bits;
            pos+=width;

            if(last){
                new_gram.rl_ptr.write(rule, rule_start_ptr/r_bits);
                rule_start_ptr = bit_pos;
                rule++;
                //std::cout<<""<<std::endl;
            }
        }
        min_sym = max_sym;
        max_sym+=p_gram.metadata[i+1].n_rules;
        assert((max_sym-new_gram.max_tsym)==rule);
        assert(pos==n_bits);
    }

    rules_buffer.destroy();
    assert((bit_pos/r_bits)==new_gram.g);
    assert(rule==(new_gram.r-(new_gram.max_tsym+1)));
    new_gram.rl_ptr.write(rule, bit_pos/r_bits);

    size_t offset = new_gram.g-new_gram.c;
    new_gram.str_boundaries.resize(new_gram.s+1);
    for(size_t str=0;str<=new_gram.s;str++){
        new_gram.str_boundaries[str] = offset+str;
    }
    assert(new_gram.str_boundaries[0]==offset);
    assert(new_gram.str_boundaries.back()==bit_pos/r_bits);
}

template<class com_gram_type>
void complete2partial_gram(com_gram_type& gram, std::string& new_p_gram_file){

}

#endif //LCG_PARTIAL_GRAM_H
