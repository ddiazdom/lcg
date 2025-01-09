//
// Created by Diaz, Diego on 3.10.2024.
//

#ifndef LCG_PLAIN_GRAM_H
#define LCG_PLAIN_GRAM_H
#include "phrase_set.h"

//a string subset is a subset of consecutive strings compressed in a text chunk.
// Text chunks are small (a couple of hundreds MBs), so they should not contain
// that many strings. We need this struct to sort the compressed strings in input
// order after the compression algorithm finishes
struct string_subset{
    uint64_t txt_id:48;
    uint8_t lvl;
    uint64_t offset:38; //a grammar can have up to 274,877 million of strings
    uint32_t n_strings:26;//a string subset (i.e., those in chunk) can have up to 67 million of strings

    string_subset(size_t _txt_id, uint8_t _lvl, uint32_t _offset, uint32_t _n_strings): txt_id(_txt_id),
                                                                                        lvl(_lvl),
                                                                                        offset(_offset),
                                                                                        n_strings(_n_strings){}
    string_subset(): txt_id(0),
                     lvl(0),
                     offset(0),
                     n_strings(0){}
};

struct plain_gram {

    uint64_t text_size=0;
    size_t n_levels=0;
    uint8_t s_sym='\n';

    std::vector<uint32_t> comp_string;
    std::vector<string_subset> str_orders;

    std::vector<size_t> fps_len;
    std::vector<uint64_t *> fps;

    phrase_set<uint8_t> ter_dict;
    std::vector<phrase_set<uint32_t>> nt_dicts;

    explicit plain_gram(size_t lvl_cap, uint8_t sep_sym, float load_factor=0.85){
        assert(lvl_cap>2);
        ter_dict.set_load_factor(load_factor);
        fps.resize(lvl_cap);
        fps_len.resize(lvl_cap);
        nt_dicts.resize(lvl_cap-2);
        for(auto &nt_lvl: nt_dicts){
            nt_lvl.set_load_factor(load_factor);
        }
        s_sym = sep_sym;

        size_t alpha_size = std::numeric_limits<uint8_t>::max()+1;
        fps[0] = mem<uint64_t>::allocate(alpha_size);
        fps_len[0] = alpha_size;
        for(size_t i=0;i<alpha_size;i++){
            fps[0][i] = XXH3_64bits(&i, sizeof(uint8_t));
            assert(fps[0][i]!=0);
        }
        fps[0][sep_sym]=0;//small hack

        //small hack as the metasymbols are one-based
        for(size_t i=1;i<fps.size();i++){
            fps[i] = mem<uint64_t>::allocate(1);
            fps[i][0]=0;
            fps_len[i]=1;
        }
    }

    plain_gram(plain_gram&& other) noexcept {
        fps.swap(other.fps);
        fps_len.swap(other.fps_len);
        ter_dict.swap(other.ter_dict);
        nt_dicts.swap(other.nt_dicts);
        comp_string.swap(other.comp_string);
        str_orders.swap(other.str_orders);
        std::swap(n_levels, other.n_levels);
        std::swap(s_sym, other.s_sym);
    }

    plain_gram(const plain_gram& other) noexcept {
        copy(other);
    }

    plain_gram& operator=(const plain_gram& other){
        copy(other);
        return *this;
    }

    void copy(const plain_gram& other){

        for(auto& fps_set : fps){
            if(fps_set!= nullptr){
                mem<uint64_t>::deallocate(fps_set);
            }
        }

        fps.resize(other.fps.size());
        for(size_t i=0;i<fps.size();i++){
            fps[i] = mem<uint64_t>::allocate(other.fps_len[i]);
            memcpy(fps[i], other.fps[i], other.fps_len[i]*sizeof(uint64_t));
        }

        fps_len = other.fps_len;
        ter_dict = other.ter_dict;
        nt_dicts = other.nt_dicts;
        comp_string = other.comp_string;
        str_orders = other.str_orders;
        n_levels = other.n_levels;
        s_sym = other.s_sym;
    }

    void get_gram_levels() {
        n_levels = !ter_dict.empty();
        size_t l=0;
        while(!nt_dicts[l].empty()){
            n_levels++;
            l++;
        }
    }

    void swap(plain_gram& other){
        fps.swap(other.fps);
        fps_len.swap(other.fps_len);
        ter_dict.swap(other.ter_dict);
        nt_dicts.swap(other.nt_dicts);
        comp_string.swap(other.comp_string);
        str_orders.swap(other.str_orders);
        std::swap(n_levels, other.n_levels);
        std::swap(s_sym, other.s_sym);
    }

    [[nodiscard]] inline bool empty() const {
        return ter_dict.empty();
    }

    size_t mem_usage(){
        size_t bytes = ter_dict.mem_usage();
        for(auto const& dict: nt_dicts){
            bytes+=dict.mem_usage();
        }

        for(auto const& f_len : fps_len){
            bytes+=f_len*sizeof(uint64_t);
            bytes+=sizeof(uint64_t);//the length
        }

        bytes+=comp_string.capacity()*sizeof(uint32_t);
        bytes+=str_orders.capacity()*sizeof(string_subset);

        return bytes;
    }

    size_t eff_mem_usage(){
        size_t bytes = ter_dict.eff_mem_usage();
        for(auto const& dict: nt_dicts){
            bytes+=dict.eff_mem_usage();
        }

        for(auto const& f_len : fps_len){
            bytes+=f_len*sizeof(uint64_t);
            bytes+=sizeof(uint64_t);//the length
        }

        bytes+=comp_string.size()*sizeof(uint32_t);
        bytes+=str_orders.size()*sizeof(string_subset);

        return bytes;
    }

    size_t av_bytes(){
        size_t bytes = ter_dict.buff_bytes_available();
        for(auto const& dict: nt_dicts){
            bytes+=dict.buff_bytes_available();
        }
        return bytes;
    }

    [[nodiscard]] inline uint8_t sep_sym() const {
        return s_sym;
    }

    [[nodiscard]] inline size_t lvl_cap() const {
        return fps.size();
    }

    inline void update_fps(size_t round){
        if(round>0){
            nt_dicts[round-1].update_fps(fps[round], fps_len[round], fps[round+1], fps_len[round+1]);
        }else{
            ter_dict.update_fps(fps[round], fps_len[round], fps[round+1], fps_len[round+1]);
        }
    }

    inline void update_fps() {
        size_t round=0;
        update_fps(round++);
        while(!nt_dicts[round-1].empty()){
            update_fps(round++);
        }
    }

    inline void update_fps_with_sink(size_t round, const plain_gram& sink){
        if(round>0){
            nt_dicts[round-1].update_fps_with_sink(sink.fps[round], sink.fps_len[round], sink.alphabet(round),
                                                   fps[round], fps_len[round],
                                                   fps[round+1], fps_len[round+1]);
        }else{
            ter_dict.update_fps_with_sink(sink.fps[round], sink.fps_len[round], sink.alphabet(round),
                                          fps[round], fps_len[round],
                                          fps[round+1], fps_len[round+1]);
        }
    }

    [[nodiscard]] inline size_t alphabet(size_t round) const {
        switch (round) {
            case 0:
                return std::numeric_limits<uint8_t>::max();
            case 1:
                return ter_dict.size();
            default:
                return nt_dicts[round-2].size();
        }
    }

    [[nodiscard]] inline size_t vbyte_usage() const {
        size_t bytes= ter_dict.vbyte_usage();
        for(auto const& nt_dict : nt_dicts){
            bytes+=nt_dict.vbyte_usage();
        }

        for(unsigned int i : comp_string){
            bytes+= vbyte_len(i);
        }
        bytes+=str_orders.capacity()*sizeof(string_subset);

        for(auto const& f_len : fps_len){
            bytes+=f_len*5;
            bytes+=sizeof(uint64_t);//the length
        }
        return bytes;
    }

    void clear_fps(){
        for(size_t i=1;i<fps.size();i++){
            fps[i] = mem<uint64_t>::reallocate(fps[i], 1);
            fps[i][0]=0;
            fps_len[i]=1;
        }
    }

    void destroy_tables(){
        ter_dict.destroy_table();
        for(auto & nt_dict : nt_dicts){
            nt_dict.destroy_table();
        }
    }

    [[nodiscard]] inline size_t gram_alphabet() const {
        size_t tot_symbols = std::numeric_limits<uint8_t>::max()+1;//alphabet of terminals

        //alphabet of nonterminals (different from the grammar's start symbol)
        tot_symbols+=ter_dict.size();
        for(auto const& nt_dict : nt_dicts){
            tot_symbols+=nt_dict.size();
        }
        //

        return tot_symbols+1;//the +1 counts for the grammar's start symbol, where comp_string is the right-hand side
    }

    [[nodiscard]] inline size_t gram_size() const {
        size_t size = std::numeric_limits<uint8_t>::max()+1;
        size+=ter_dict.tot_symbols();
        for(auto const& nt_dict : nt_dicts){
            size+=nt_dict.tot_symbols();
        }
        size+=comp_string.size();
        return size;
    }

    [[nodiscard]] inline uint8_t tot_rounds() const {
        uint8_t max_lvl = 0;
        for(auto const& subset : str_orders){
            if(subset.lvl>max_lvl) max_lvl = subset.lvl;
        }
        return max_lvl;
    }

    void clear(){
        ter_dict.clear();
        for(auto &dict: nt_dicts){
            dict.clear();
        }

        for(size_t i=1;i<fps.size();i++){
            fps[i] = mem<uint64_t>::reallocate(fps[i], 1);
            //small hack as the metasymbols are one-based
            fps[i][0]=0;
            fps_len[i]=1;
        }
    }

    void reorder_strings(){
        std::vector<uint32_t> new_string(comp_string.size(), 0);
        std::sort(str_orders.begin(), str_orders.end(), [](auto const& a, auto const& b) -> bool{
            return a.txt_id < b.txt_id;
        });

        size_t i=0;
        for(auto & subset : str_orders){
            memcpy(&new_string[i], &comp_string[subset.offset], subset.n_strings*sizeof(uint32_t));
            subset.offset = i;
            i+=subset.n_strings;
        }

        assert(i==new_string.size());
        new_string.swap(comp_string);
    }

    [[nodiscard]] inline size_t tot_strings() const {
        return comp_string.size();
    }

    [[nodiscard]] inline size_t txt_size() const {
        return text_size;
    }

    [[nodiscard]] inline size_t phrases_round(size_t round) const {
        if(round>0){
            return nt_dicts[round-1].size();
        }else{
            return ter_dict.size();
        }
    }

    std::string get_stats(size_t pad=0){
        std::string pad_str(pad, ' ');

        std::string stats=pad_str+"Level 1, number of phrases: "+std::to_string(ter_dict.size())+",  number of symbols: "+std::to_string(ter_dict.tot_symbols());
        size_t round=2;
        for(auto const& lvl_set: nt_dicts){
            if(lvl_set.empty()) break;
            stats+="\n"+pad_str+"Level "+std::to_string(round)+", number of phrases: "+std::to_string(lvl_set.size())+", number of symbols: "+std::to_string(lvl_set.tot_symbols());
            round++;
        }
        stats+="\n"+pad_str+"Tot. strings: "+std::to_string(comp_string.size());
        stats+="\n"+pad_str+"Space space: real:"+report_space((off_t)mem_usage())+" eff:"+report_space((off_t)eff_mem_usage());
        return stats;
    }

    ~plain_gram() {

        /*if(!ter_dict.empty()){
            std::cout<<" dict ter "<<ter_dict.size()<<", "<<ter_dict.load_factor()<<std::endl;
            ter_dict.psl_dist();
        }
        size_t i=0;
        for(auto &nt_dict : nt_dicts){
            if(!nt_dict.empty()){
                std::cout<<i<<" "<<nt_dict.size()<<", "<<nt_dict.load_factor()<<std::endl;
                nt_dict.psl_dist();
            }
            i++;
        }*/

        for(auto fps_set : fps){
            if(fps_set!= nullptr){
                mem<uint64_t>::deallocate(fps_set);
            }
        }
    }

    size_t serialize(std::ofstream &ofs){
        size_t written_bytes = 0;
        written_bytes+= serialize_elm(ofs, text_size);
        written_bytes+= serialize_elm(ofs, n_levels);
        written_bytes+= serialize_elm(ofs, s_sym);
        //assert(fps.size()==fps_len.size());
        //written_bytes+= serialize_plain_vector(ofs, fps_len);
        //for(size_t i=0;i<fps_len.size();i++){
        //    written_bytes+= serialize_raw_vector(ofs, fps[i], fps_len[i]);
        //}
        written_bytes+= ter_dict.serialize(ofs);
        size_t n_dicts = nt_dicts.size();
        written_bytes+= serialize_elm(ofs, n_dicts);
        for(auto & nt_dict : nt_dicts){
            written_bytes+= nt_dict.serialize(ofs);
        }
        written_bytes+= serialize_plain_vector(ofs, comp_string);
        written_bytes+= serialize_plain_vector(ofs, str_orders);

        return written_bytes;
    }

    void load(std::ifstream &ifs){
        load_elm(ifs, text_size);
        load_elm(ifs, n_levels);
        load_elm(ifs, s_sym);
        //load_plain_vector(ifs, fps_len);
        //for(size_t i=0;i<fps_len.size();i++){
        //    load_raw_vector(ifs, fps[i], fps_len[i]);
        //}
        ter_dict.load(ifs);
        size_t n_dicts;
        load_elm(ifs, n_dicts);
        nt_dicts.resize(n_dicts);
        for(auto & nt_dict : nt_dicts){
            nt_dict.load(ifs);
        }
        load_plain_vector(ifs, comp_string);
        load_plain_vector(ifs, str_orders);
    }
};

struct plain_gram_buffer{

    uint64_t text_size=0;
    size_t n_levels=0;
    size_t active_level = 0;
    uint8_t s_sym='\n';

    std::ofstream ofs;
    std::ifstream ifs;

    std::vector<uint32_t> comp_string;
    std::vector<string_subset> str_orders;

    phrase_set<uint8_t> ter_dict;
    phrase_set<uint32_t> nt_dict;

    plain_gram_buffer(std::string& i_file, std::string& o_file){

        ifs = std::ifstream(i_file, std::ios::binary);
        load_elm(ifs, text_size);
        load_elm(ifs, n_levels);
        load_elm(ifs, s_sym);
        assert(n_levels>0);

        ter_dict.load(ifs);
        //TODO check o_file does not exist
        ofs = std::ofstream(o_file, std::ios::binary);

        std::streamoff header_bytes = sizeof(text_size)+sizeof(n_levels)+sizeof(s_sym);
        ofs.seekp(header_bytes);
    }

    void next_level(){
        if(active_level==0){
            ter_dict.serialize(ofs);
            ter_dict.destroy();
        }else{
            nt_dict.serialize(ofs);
            nt_dict.destroy();
        }

        active_level++;
        if(active_level<n_levels){
            nt_dict.load(ifs);
        }
    }

    void load_comp_string(){
        load_plain_vector(ifs, comp_string);
        load_plain_vector(ifs, str_orders);
    }

    void close(){
        serialize_plain_vector(ofs, comp_string);
        serialize_plain_vector(ofs, str_orders);

        //store the header information at the beginning of the file.
        //we are assuming the grammar was modified
        ofs.seekp(0);
        serialize_elm(ofs, text_size);
        serialize_elm(ofs, n_levels);
        serialize_elm(ofs, s_sym);
    }
};
#endif //LCG_PLAIN_GRAM_H
