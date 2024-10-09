//
// Created by Diaz, Diego on 17.9.2023.
//

#ifndef LCG_LZ_LIKE_LC_PARSING_H
#define LCG_LZ_LIKE_LC_PARSING_H

#include <random>
#include <fcntl.h>
#include <cstring>
#include "xxhash.h"

#include "cds/macros.h"
#include "cds/ts_queue.h"
#include "cds/ts_priority_queue.h"
#include "cds/utils.h"
#include "cds/cdt_common.hpp"
#include "cds/mmap_allocator.h"
#include "cds/vbyte_encoding.h"

#include "text_handler.h"
#include "phrase_set.h"
#include "grammar.h"
#include "merge_grams.h"

struct parsing_opts{
    size_t n_threads{};
    unsigned long n_chunks{};
    off_t chunk_size{};
    off_t page_cache_limit{};
    size_t sep_sym{};
    uint64_t orig_seed=0;//this is the seed to create random seeds for p_seeds
    std::vector<uint64_t> p_seeds;
    float i_frac{};
};

struct parsing_state{
    size_t chunk_size;
    uint8_t sep_sym;
    int fd_r;
    size_t f_size=0;
    off_t rem_bytes=0;
    off_t read_bytes=0;
    off_t r_page_cache_bytes=0;
    off_t page_cache_limit=0;
    size_t chunk_id=0;
    size_t proc_syms=0;
    size_t sink_gram_mem_usage=0;
    float max_frac;

    parsing_state(std::string& i_file, uint8_t s_sym, size_t c_size, off_t p_cache_lim, float _max_frac): chunk_size(c_size),
                                                                                                          sep_sym(s_sym),
                                                                                                          page_cache_limit(p_cache_lim),
                                                                                                          max_frac(_max_frac){
        fd_r = open(i_file.c_str(), O_RDONLY);
        f_size = file_size(i_file);
#ifdef __linux__
        posix_fadvise(fd_r, 0, f_size, POSIX_FADV_SEQUENTIAL);
#endif
        rem_bytes = (off_t)f_size;
    }

    ~parsing_state(){
#ifdef __linux__
        posix_fadvise(fd_r, 0, f_size, POSIX_FADV_DONTNEED);
#endif
        close(fd_r);
    }

    void flush_page_cache(){
#ifdef __linux__
        std::cout<<"removing from page cache "<<r_page_cache_bytes<<" "<<read_bytes<<std::endl;
        posix_fadvise(fd_r, read_bytes-r_page_cache_bytes, r_page_cache_bytes, POSIX_FADV_DONTNEED);
        r_page_cache_bytes=0;
#endif
    }
};

struct phrase_overflow{
    uint32_t position;
    uint32_t length;
    uint32_t metasymbol;

    [[nodiscard]] inline uint32_t right_end() const {
        return position+length-1;
    }
};

/*struct inv_perm_elm{
    uint32_t orig_mt;
    uint64_t fp;
};*/

/*template<class sym_type, bool p_round>
void create_meta_sym(text_chunk& chunk, uint64_t pf_seed, const phrase_set<sym_type>& dict,
                      buff_vector<uint64_t>& prev_fps){

    sym_type* text;
    off_t txt_size;
    if constexpr (p_round){
        text = chunk.text;
        txt_size = chunk.text_bytes;
    }else{
        text = chunk.parse;
        txt_size = chunk.parse_size;
    }

    auto avb_addr = chunk.get_free_mem_area();
    buff_vector<uint64_t> new_fps(avb_addr.first, avb_addr.second);
    new_fps.resize(dict.size()+1);
    chunk.add_used_bytes((off_t)new_fps.static_buff_usage());

    size_t source, end, tot_symbols=0;
    std::vector<uint64_t> fp_sequence;
    new_fps[0] = 0;
    for(size_t i=0, mt_sym=1;i<dict.size();i++, mt_sym++) {
        source = phrase_set[i].source;
        end = source + phrase_set[i].len;
        for(size_t j=source;j<end;j++){
            assert(text[j]>0);
            assert(text[j]<prev_fps.size());
            fp_sequence.push_back(prev_fps[text[j]]);
        }
        //inv_mt_perm[i].fp = XXH3_64bits_withSeed(fp_sequence.data(), fp_sequence.size()*sizeof(uint64_t), pf_seed);
        new_fps[mt_sym] = XXH3_64bits(fp_sequence.data(), fp_sequence.size()*sizeof(uint64_t));
        tot_symbols+=fp_sequence.size();
        fp_sequence.clear();
    }

    new_fps.swap(prev_fps);
    new_fps.destroy();
    //TODO handle this later
    //chunk.p_gram.template append_new_lvl<sym_type>(text, phrase_set, tot_symbols);
}*/

void mul_thread_ter_collapse(plain_gram& sink_gram, std::vector<text_chunk>& chunks){

    phrase_set<uint8_t>& sink_set = sink_gram.ter_dict;
    size_t size_before = sink_set.size();

    auto ter_worker = [&](plain_gram& gram, size_t& new_syms) {

        uint64_t* o_map = gram.fps[1];
        uint64_t o_map_len = gram.fps_len[1];
        phrase_set<uint8_t>& coll_set = gram.ter_dict;

        assert(o_map_len==(coll_set.size()+1));
        coll_set.destroy_table();
        new_syms=0;

        auto it = coll_set.begin();
        auto it_end = coll_set.end();

        uint32_t mt;
        size_t phrase=1;
        while(it!=it_end){
            auto phr = *it;
            if(sink_set.find(phr.phrase, phr.len, mt)){
                o_map[phrase] = mt;
            }else{
                o_map[phrase] = it.pos() | 0x8000000000000000UL;
                new_syms+= phr.len+(2*sizeof(uint32_t));
            }
            phrase++;
            ++it;
        }
        assert(phrase==o_map_len);
    };

    std::vector<size_t> nb(chunks.size(), 0);
    std::vector<std::thread> threads;
    threads.reserve(chunks.size()-1);
    for(size_t i=0;i<chunks.size()-1;i++){
        threads.emplace_back(ter_worker, std::ref(chunks[i].gram), std::ref(nb[i]));
    }
    ter_worker(chunks.back().gram, nb.back());

    for (auto &thread: threads) {
        thread.join();
    }

    size_t tot_new_sym=0;
    for(auto b : nb){
        tot_new_sym+=b;
    }
    sink_set.set_stream_capacity(sink_set.stream_len()+tot_new_sym);

    uint32_t len;
    size_t pos;
    for(auto & chunk : chunks) {
        const uint8_t *stream = chunk.gram.ter_dict.phr_stream();
        for(size_t j=0;j<chunk.gram.fps_len[1];j++){
            if(chunk.gram.fps[1][j] & 0x8000000000000000UL){
                pos = chunk.gram.fps[1][j] & 0x7FFFFFFFFFFFFFFF;
                memcpy(&len, &stream[pos], sizeof(uint32_t));
                pos+=sizeof(uint32_t);
                chunk.gram.fps[1][j] = sink_set.insert(&stream[pos], len);
            }
        }
        chunk.gram.ter_dict.clear();//the table is fully destroyed
    }

    if(!sink_set.empty()){
        std::cout<<"  Growth factor 1 "<<float(size_before)/float(sink_gram.ter_dict.size())<<std::endl;
    }
}

void sin_thread_ter_collapse(plain_gram& sink_gram, std::vector<text_chunk>& chunks) {

    phrase_set<uint8_t>& sink_set = sink_gram.ter_dict;
    size_t size_before = sink_set.size();

    for(auto & chunk : chunks) {
        uint64_t* o_map = chunk.gram.fps[1];
        uint64_t o_map_len = chunk.gram.fps_len[1];
        phrase_set<uint8_t>& coll_set = chunk.gram.ter_dict;

        assert(o_map_len==(coll_set.size()+1));
        coll_set.destroy_table();

        auto it = coll_set.begin();
        auto it_end = coll_set.end();

        uint32_t phrase = 1;
        while(it!=it_end){
            auto phr = *it;
            o_map[phrase++] = sink_set.insert(phr.phrase, phr.len);
            ++it;
        }
        assert(phrase==o_map_len);
        coll_set.clear();
    }

    if(!sink_set.empty()){
        std::cout<<"  Growth factor 1 "<<float(size_before)/float(sink_set.size())<<std::endl;
    }
}

void mul_thread_nt_collapse(plain_gram& sink_gram, std::vector<text_chunk>& chunks, size_t round, uint32_t prev_alpha_sink) {

    size_t size_before = sink_gram.nt_dicts[round-1].size();
    phrase_set<uint32_t>& sink_set = sink_gram.nt_dicts[round-1];

    auto nt_worker = [&](plain_gram& gram, size_t& new_syms) {

        const uint64_t* i_map = gram.fps[round];
        const uint64_t i_map_len = gram.fps_len[round];
        uint64_t* o_map = gram.fps[round+1];
        uint64_t o_map_len = gram.fps_len[round+1];
        phrase_set<uint32_t>& coll_set = gram.nt_dicts[round-1];
        new_syms = 0;

        assert(o_map_len==(coll_set.size()+1));
        coll_set.destroy_table();

        auto it = coll_set.begin();
        auto it_end = coll_set.end();

        uint32_t mt;
        size_t phrase=1;
        uint64_t flag = 0x8000000000000000UL;

        while(it!=it_end){
            auto phr = *it;
            for(size_t j=0;j<phr.len;j++){
                if(phr.phrase[j]>prev_alpha_sink){
                    size_t rank = phr.phrase[j]-prev_alpha_sink;
                    assert(rank<i_map_len);
                    phr.phrase[j] = i_map[rank]+1;
                }
                assert(phr.phrase[j]>0);
            }

            if(sink_set.find(phr.phrase, phr.len, mt)){
                o_map[phrase] = mt;
            } else {
                o_map[phrase] = it.pos() | flag;
                new_syms += phr.len+2;//the +2 considers the length and the metasymbol
            }
            phrase++;
            ++it;
        }
        assert(phrase==o_map_len);
    };

    std::vector<std::thread> threads;
    threads.reserve(chunks.size()-1);
    std::vector<size_t> nb(chunks.size(), 0);

    for(size_t i=0;i<chunks.size()-1;i++) {
        threads.emplace_back(nt_worker, std::ref(chunks[i].gram), std::ref(nb[i]));
    }
    nt_worker(chunks.back().gram, nb.back());

    for (auto &thread: threads) {
        thread.join();
    }

    size_t tot_new_sym=0;
    for(auto b : nb){
        tot_new_sym+=b;
    }
    sink_set.set_stream_capacity(sink_set.stream_len()+tot_new_sym);

    uint32_t len;
    size_t pos;
    for(auto & chunk : chunks) {
        const uint32_t *stream = chunk.gram.nt_dicts[round-1].phr_stream();
        uint64_t* o_map = chunk.gram.fps[round+1];
        uint64_t o_map_len = chunk.gram.fps_len[round+1];
        for(size_t j=0;j<o_map_len;j++){
            if(o_map[j] & 0x8000000000000000UL){
                pos = o_map[j] & 0x7FFFFFFFFFFFFFFF;
                memcpy(&len, &stream[pos], sizeof(uint32_t));
                pos++;
                o_map[j] = sink_set.insert(&stream[pos], len);
            }
        }
        chunk.gram.nt_dicts[round-1].clear();
        chunk.gram.fps[round]= mem<uint64_t>::reallocate(chunk.gram.fps[round], 1);
        chunk.gram.fps_len[round] = 1;
    }
}

void sin_thread_nt_collapse(plain_gram& sink_gram, std::vector<text_chunk>& chunks, size_t round, uint32_t prev_alpha_sink){

    assert(round>=1);
    size_t size_before = sink_gram.nt_dicts[round-1].size();
    phrase_set<uint32_t>& sink_set = sink_gram.nt_dicts[round-1];

    for(auto & chunk : chunks){

        const uint64_t* i_map = chunk.gram.fps[round];
        const uint64_t i_map_len = chunk.gram.fps_len[round];
        uint64_t* o_map = chunk.gram.fps[round+1];
        uint64_t o_map_len = chunk.gram.fps_len[round+1];
        phrase_set<uint32_t>& coll_set = chunk.gram.nt_dicts[round-1];

        assert(o_map_len==(coll_set.size()+1));
        coll_set.destroy_table();

        auto it = coll_set.begin();
        auto it_end = coll_set.end();

        uint32_t phrase = 1;
        while(it!=it_end) {
            auto phr = *it;
            for(size_t j=0;j<phr.len;j++){
                if(phr.phrase[j]>prev_alpha_sink){
                    size_t rank = phr.phrase[j]-prev_alpha_sink;
                    assert(rank<i_map_len);
                    phr.phrase[j] = i_map[rank]+1;
                }
                assert(phr.phrase[j]>0);
            }
            o_map[phrase++] = sink_set.insert(phr.phrase, phr.len);
            ++it;
        }
        assert(phrase==o_map_len);
        coll_set.clear();
        chunk.gram.fps[round] = mem<uint64_t>::reallocate(chunk.gram.fps[round], 1);
        chunk.gram.fps_len[round] = 1;
    }

    if(!sink_set.empty()){
        std::cout<<"  Growth factor "<<round+1<<" "<<float(size_before)/float(sink_set.size())<<std::endl;
    }
}

void collapse_grams(plain_gram& sink_gram, std::vector<text_chunk>& chunks) {

    //compute the length of each rule set
    std::vector<uint32_t> prev_lvl_alpha(sink_gram.fps.size(), 0);
    prev_lvl_alpha[1] = sink_gram.ter_dict.size();
    for(size_t i=2;i<prev_lvl_alpha.size();i++){
        prev_lvl_alpha[i] = sink_gram.nt_dicts[i-2].size();
    }

    //just swap if the sink grammar is empty
    if(sink_gram.empty()){
        sink_gram.swap(chunks[0].gram);
    }

    size_t tot_strings = 0;
    for(auto & text_chunk : chunks){
        tot_strings+=text_chunk.gram.comp_string.size();
        text_chunk.gram.get_gram_levels();
    }

    if(chunks.size()>1){
        mul_thread_ter_collapse(sink_gram, chunks);
    }else{
        sin_thread_ter_collapse(sink_gram, chunks);
    }

    for(size_t round=1;round<sink_gram.nt_dicts.size();round++){
        if(round<=4 && chunks.size()>1){
            mul_thread_nt_collapse(sink_gram, chunks, round, prev_lvl_alpha[round]);
        } else {
            sin_thread_nt_collapse(sink_gram, chunks, round, prev_lvl_alpha[round]);
        }
    }

    size_t pos = chunks[0].gram.comp_string.size();
    sink_gram.comp_string.resize(sink_gram.comp_string.size()+tot_strings);
    for(auto & chunk : chunks){
        size_t n = chunk.gram.n_levels;
        for(size_t j=0;j<chunk.gram.comp_string.size();j++){
            //TODO fix this
            //assert(text_chunks[i].comp_string[j]>0 && text_chunks[i].comp_string[j]<text_chunks[i].fps_len[n]);
            //text_chunks[0].comp_string[pos++] = text_chunks[i].fps[n][text_chunks[i].comp_string[j]]+1;
        }
        chunk.gram.comp_string.clear();
    }

    size_t round=0;
    sink_gram.update_fps(round++);
    while(!sink_gram.nt_dicts[round-1].empty()){
        sink_gram.update_fps(round++);
    }
}

void finish_byte_parse(text_chunk& chunk, off_t &parse_distance, std::vector<phrase_overflow>& phr_with_ovf){

    uint32_t mt_sym=1;
    off_t txt_size = chunk.text_bytes;

    off_t ovf_idx=0, next_ovf=-1;
    if(!phr_with_ovf.empty()){
        next_ovf = phr_with_ovf[ovf_idx].right_end();
    }

    off_t pos = txt_size-1;
    chunk.parse = (uint32_t *)(chunk.text+parse_distance);
    chunk.parse--;

    while(pos>=0){
        if(pos==next_ovf){
            pos -= phr_with_ovf[ovf_idx].length;
            *chunk.parse = phr_with_ovf[ovf_idx].metasymbol;
            next_ovf=-1;
            if(ovf_idx<off_t(phr_with_ovf.size()-1)){
                next_ovf = phr_with_ovf[++ovf_idx].right_end();
            }
        }  else {
            pos-=vbyte_decoder<uint32_t>::read_right2left(&chunk.text[pos], mt_sym);
            *chunk.parse = mt_sym;
        }
        chunk.parse--;
        while(pos>next_ovf && chunk.text[pos]==0) pos--;
    }

    chunk.parse++;
    assert(pos==-1);
    assert(chunk.parse[chunk.parse_size-1]==0);
    assert((uintptr_t) chunk.text<= (uintptr_t)chunk.parse && (uintptr_t)chunk.parse<= (uintptr_t)&chunk.text[txt_size-1]);
    assert((uintptr_t) chunk.text<= (uintptr_t)(chunk.parse+chunk.parse_size) &&
           (uintptr_t) (chunk.parse+chunk.parse_size) <= (uintptr_t)(chunk.text+chunk.buffer_bytes));
    assert(parse_distance==(off_t)chunk.dist(reinterpret_cast<uint8_t*>(chunk.parse+chunk.parse_size)));
}

template<bool query_sink>
void byte_par_r2l(text_chunk& chunk, off_t& n_strings, size_t sep_sym);

template<>
void byte_par_r2l<true>(text_chunk& chunk, off_t& n_strings, size_t sep_sym) {

    size_t next_av_in_sink = chunk.sink_gram.ter_dict.size();
    uint64_t * fps = chunk.gram.fps[chunk.round];
    uint8_t * text = chunk.text;
    uint64_t hash;
    off_t lb, rb = chunk.text_bytes-1, i=chunk.text_bytes-2, byte_offset, parse_size;
    uint8_t v_len;
    assert(text[i+1]==sep_sym && text[i]>text[i+1]);

    text[rb] = 128;//the vbyte code of the metasymbol mt=0 representing a separator in the next round of parsing
    n_strings=1;
    parse_size=4;//we will count in bytes of sizeof(uint32_t)
    off_t mbo[2]={4};

    std::vector<phrase_overflow> phr_with_ovf;
    bool r_cmp = true, l_cmp, new_str;
    uint32_t mt_sym, phrase_len;
    uint8_t mid_sym = text[i];
    while(--i>0 && text[i]==mid_sym);

    while(i>=0){
        l_cmp = fps[text[i]]>fps[mid_sym];
        if(l_cmp && !r_cmp){
            lb = i+1;
            new_str = mid_sym==sep_sym;
            lb += new_str;
            phrase_len=rb-lb;

            byte_offset = rb + parse_size;
            mbo[byte_offset > mbo[1]] = byte_offset;

            hash = XXH3_64bits(&text[lb], phrase_len);
            bool found = chunk.sink_gram.ter_dict.find(&text[lb], phrase_len, mt_sym, hash);
            if(!found){
                mt_sym = next_av_in_sink + chunk.gram.ter_dict.insert(&text[lb], phrase_len, hash);
            }
            mt_sym++;

            v_len = vbyte_len(mt_sym);
            if(__builtin_expect(v_len>phrase_len, 0)){
                //metasymbol does not fit its phrase
                phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
            } else {
                vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
                memset(&text[lb+v_len], 0, phrase_len-v_len);
            }

            parse_size+=4;
            rb = i+1;
            if(new_str){
                assert(text[rb]==sep_sym);
                text[rb] = 128;//vbyte code for 0 (the separator symbol in the next levels)
                n_strings++;
                parse_size+=4;
            }
        }
        r_cmp = l_cmp;
        mid_sym = text[i];
        while(--i>0 && text[i]==mid_sym);
    }

    lb = 0;
    phrase_len=rb-lb;
    byte_offset = rb + parse_size;
    mbo[byte_offset > mbo[1]] = byte_offset;

    hash = XXH3_64bits(&text[lb], phrase_len);
    bool found = chunk.sink_gram.ter_dict.find(&text[lb], phrase_len, mt_sym, hash);
    if(!found){
        mt_sym = next_av_in_sink + chunk.gram.ter_dict.insert(&text[lb], phrase_len, hash);
    }
    mt_sym++;

    v_len = vbyte_len(mt_sym);
    if(__builtin_expect(v_len>phrase_len, 0)){
        phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
    }else {
        vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
        memset(&text[lb+v_len], 0, phrase_len-v_len);
    }
    parse_size+=4;
    chunk.parse_size = parse_size/4;//divide by 4=sizeof(uint32_t) to avoid the <<

    //computes how many bytes we require for the parse
    off_t max_byte_offset = INT_CEIL(mbo[1], sizeof(uint32_t))*sizeof(uint32_t);
    chunk.increase_capacity(max_byte_offset);

    chunk.gram.update_fps_with_sink(chunk.round, chunk.sink_gram);
    finish_byte_parse(chunk, max_byte_offset, phr_with_ovf);

    //std::cout<<"\nThe alphabet is "<<chunk.gram.alphabet(chunk.round)<<" and the size is "<<chunk.gram.ter_dict.size()<<" "<<chunk.parse_size<<std::endl;
}

template<>
void byte_par_r2l<false>(text_chunk& chunk, off_t& n_strings, size_t sep_sym) {

    uint64_t * fps = chunk.gram.fps[chunk.round];
    uint8_t * text = chunk.text;
    off_t lb, rb = chunk.text_bytes-1, i=chunk.text_bytes-2, byte_offset, parse_size;
    uint8_t v_len;
    assert(text[i+1]==sep_sym && text[i]>text[i+1]);

    text[rb] = 128;//the vbyte code of the metasymbol mt=0 representing a separator in the next round of parsing
    n_strings=1;
    parse_size=4;//we will count in bytes of sizeof(uint32_t)
    off_t mbo[2]={4};

    std::vector<phrase_overflow> phr_with_ovf;
    bool r_cmp = true, l_cmp, new_str;
    uint32_t mt_sym, phrase_len;
    uint8_t mid_sym = text[i];
    while(--i>0 && text[i]==mid_sym);

    while(i>=0){
        l_cmp = fps[text[i]]>fps[mid_sym];
        if(l_cmp && !r_cmp){
            lb = i+1;
            new_str = mid_sym==sep_sym;
            lb += new_str;
            phrase_len=rb-lb;

            byte_offset = rb + parse_size;
            mbo[byte_offset > mbo[1]] = byte_offset;
            mt_sym = chunk.gram.ter_dict.insert(&text[lb], phrase_len)+1;

            v_len = vbyte_len(mt_sym);
            if(__builtin_expect(v_len>phrase_len, 0)){
                //metasymbol does not fit its phrase
                phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
            } else {
                vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
                memset(&text[lb+v_len], 0, phrase_len-v_len);
            }

            parse_size+=4;
            rb = i+1;
            if(new_str){
                assert(text[rb]==sep_sym);
                text[rb] = 128;//vbyte code for 0 (the separator symbol in the next levels)
                n_strings++;
                parse_size+=4;
            }
        }
        r_cmp = l_cmp;
        mid_sym = text[i];
        while(--i>0 && text[i]==mid_sym);
    }

    lb = 0;
    phrase_len=rb-lb;
    byte_offset = rb + parse_size;
    mbo[byte_offset > mbo[1]] = byte_offset;
    mt_sym = chunk.gram.ter_dict.insert(&text[lb], phrase_len) + 1;

    v_len = vbyte_len(mt_sym);
    if(__builtin_expect(v_len>phrase_len, 0)){
        phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
    }else {
        vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
        memset(&text[lb+v_len], 0, phrase_len-v_len);
    }
    parse_size+=4;
    chunk.parse_size = parse_size/4;//divide by 4=sizeof(uint32_t) to avoid the <<

    //computes how many bytes we require for the parse
    off_t max_byte_offset = INT_CEIL(mbo[1], sizeof(uint32_t))*sizeof(uint32_t);
    chunk.increase_capacity(max_byte_offset);

    chunk.gram.update_fps(chunk.round);
    finish_byte_parse(chunk, max_byte_offset, phr_with_ovf);
}

template<bool query_sink>
void int_par_l2r(text_chunk& chunk);

template<>
void int_par_l2r<true>(text_chunk& chunk){

    assert(chunk.round>0);
    const uint64_t* fps[2] = {chunk.gram.fps[chunk.round], chunk.sink_gram.fps[chunk.round]};

    uint64_t hash;
    uint32_t next_av_mt_in_sink = chunk.sink_gram.nt_dicts[chunk.round-1].size();
    uint32_t alpha_sink = chunk.sink_gram.alphabet(chunk.round);
    uint64_t sym_offset[2] = {alpha_sink, 0};

    uint32_t *text = chunk.parse;
    uint32_t mt_sym, sep_sym=0, txt_size = chunk.parse_size;
    uint32_t left_sym, middle_sym;
    uint64_t left_fp, middle_fp, right_fp;
    off_t i=0, parse_size = 0, phrase_len, lb, rb;
    bool new_str=false, in_sink;

    lb = 0;
    left_sym = text[i];
    in_sink = left_sym<=alpha_sink;
    left_fp = fps[in_sink][left_sym-sym_offset[in_sink]];

    while(++i<txt_size && text[i]==left_sym);
    assert(i<txt_size);

    middle_sym = text[i];
    in_sink = middle_sym<=alpha_sink;
    middle_fp = fps[in_sink][middle_sym-sym_offset[in_sink]];
    rb=i;
    while(++i<txt_size && text[i]==middle_sym);

    while(i<txt_size) {
        in_sink = text[i]<=alpha_sink;
        right_fp = fps[in_sink][text[i]-sym_offset[in_sink]];

        if(left_fp>middle_fp && middle_fp<right_fp){//local minimum
            phrase_len = rb-lb;

            hash = XXH3_64bits(&text[lb], phrase_len*sizeof(uint32_t));
            bool found = chunk.sink_gram.nt_dicts[chunk.round-1].find(&text[lb], phrase_len, mt_sym, hash);
            if(!found){
                mt_sym = next_av_mt_in_sink + chunk.gram.nt_dicts[chunk.round-1].insert(&text[lb], phrase_len, hash);
            }

            text[parse_size]=sep_sym;
            parse_size+=new_str;
            text[parse_size++] = mt_sym+1;
            new_str = text[rb]==sep_sym;
            lb = rb+new_str;
        }
        left_fp = middle_fp;
        middle_fp = right_fp;
        middle_sym = text[i];
        rb = i;
        while(++i<txt_size && text[i]==middle_sym);
    }
    assert(rb==(txt_size-1) && text[rb]==sep_sym);

    phrase_len = rb-lb;
    hash = XXH3_64bits(&text[lb], phrase_len*sizeof(uint32_t));
    bool found = chunk.sink_gram.nt_dicts[chunk.round-1].find(&text[lb], phrase_len, mt_sym, hash);
    if(!found){
        mt_sym = next_av_mt_in_sink+chunk.gram.nt_dicts[chunk.round-1].insert(&text[lb], phrase_len, hash);
    }

    text[parse_size]=sep_sym;
    parse_size+=new_str;
    text[parse_size++] = mt_sym+1;
    text[parse_size++] = sep_sym;

    chunk.gram.update_fps_with_sink(chunk.round, chunk.sink_gram);
    chunk.parse_size = parse_size;
}

template<>
void int_par_l2r<false>(text_chunk& chunk){

    uint64_t *fps = chunk.gram.fps[chunk.round];
    uint32_t *text = chunk.parse;

    uint32_t mt_sym, sep_sym=0, txt_size = chunk.parse_size;
    uint32_t left_sym, middle_sym;
    uint64_t left_fp, middle_fp, right_fp;
    off_t i=0, parse_size = 0, phrase_len, lb, rb;
    bool new_str=false;

    lb = 0;
    left_sym = text[i];
    left_fp = fps[left_sym];
    while(++i<txt_size && text[i]==left_sym);
    assert(i<txt_size);

    middle_sym = text[i];
    middle_fp = fps[middle_sym];
    rb=i;
    while(++i<txt_size && text[i]==middle_sym);

    while(i<txt_size) {
        right_fp = fps[text[i]];

        if(left_fp>middle_fp && middle_fp<right_fp){//local minimum
            phrase_len = rb-lb;
            mt_sym = chunk.gram.nt_dicts[chunk.round-1].insert(&text[lb], phrase_len);

            text[parse_size]=sep_sym;
            parse_size+=new_str;
            text[parse_size++] = mt_sym+1;
            new_str = text[rb]==sep_sym;
            lb = rb+new_str;
        }
        left_fp = middle_fp;
        middle_fp = right_fp;
        middle_sym = text[i];
        rb = i;
        while(++i<txt_size && text[i]==middle_sym);
    }
    assert(rb==(txt_size-1) && text[rb]==sep_sym);

    phrase_len = rb-lb;
    mt_sym = chunk.gram.nt_dicts[chunk.round-1].insert(&text[lb], phrase_len);

    text[parse_size]=sep_sym;
    parse_size+=new_str;
    text[parse_size++] = mt_sym+1;
    text[parse_size++] = sep_sym;

    chunk.gram.update_fps(chunk.round);
    chunk.parse_size = parse_size;
}

template<bool query_sink>
void compress_text_chunk(text_chunk& chunk){

    off_t n_strings=0;
    size_t sep_sym = chunk.sep_sym;
    chunk.round = 0;

    //auto start = std::chrono::steady_clock::now();
    byte_par_r2l<query_sink>(chunk, n_strings, sep_sym);
    //auto end = std::chrono::steady_clock::now();
    //report_time(start, end , 2);

    off_t size_limit = n_strings*2;
    chunk.round++;

    while(chunk.parse_size!=size_limit){
        assert(chunk.parse_size>=size_limit);
        //start = std::chrono::steady_clock::now();
        int_par_l2r<query_sink>(chunk);
        //end = std::chrono::steady_clock::now();
        //report_time(start, end , 2);
        chunk.round++;
    }

    //start = std::chrono::steady_clock::now();
    size_t pos = chunk.gram.comp_string.size();
    chunk.gram.comp_string.resize(pos+(chunk.parse_size/2));
    for(off_t i=0;i<chunk.parse_size;i+=2){
        assert(i==0 || chunk.parse[i-1]==0);
        chunk.gram.comp_string[pos++] = chunk.parse[i];
    }
    //chunk.p_gram.add_compressed_string(chunk.parse, chunk.parse_size);
    //end = std::chrono::steady_clock::now();
    //report_time(start, end , 2);
}


template<bool query_sink>
void fill_chunk_grammars(std::vector<text_chunk>& text_chunks, parsing_state& p_state){

    ts_queue<size_t> buffers_to_process;
    ts_queue<size_t> buffers_to_reuse;
    std::atomic<size_t> parser_finished{0};

    auto compressor_worker = [&]() {
        size_t buff_id;
        bool res;

        while (true) {
            res = buffers_to_process.pop(buff_id);
            if (!res){
                parser_finished.fetch_add(1, std::memory_order_acq_rel);
                break;
            }
            assert(text_chunks[buff_id].text_bytes > 0);

            text_chunks[buff_id].t_start = std::chrono::steady_clock::now();
            compress_text_chunk<query_sink>(text_chunks[buff_id]);
            memset(text_chunks[buff_id].text, 0, text_chunks[buff_id].buffer_bytes);
            text_chunks[buff_id].t_end = std::chrono::steady_clock::now();
            buffers_to_reuse.push(buff_id);
        }
    };

    std::vector<std::thread> threads;
    for (size_t i = 0; i < text_chunks.size(); i++) {
        threads.emplace_back(compressor_worker);
    }

    auto tmp_ck_size = off_t(INT_CEIL(p_state.chunk_size, sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type));
    size_t buff_id = 0;
    std::vector<size_t> byte_counts(text_chunks.size(), 0);
    size_t acc_bytes = p_state.sink_gram_mem_usage;
    float input_frac = 0;//we set it to 0 to enter the loop in case the sink grammar already exceeds the fraction threshold

    while(buff_id < text_chunks.size() && p_state.rem_bytes > 0 && input_frac<p_state.max_frac) {

        tmp_ck_size = std::min(tmp_ck_size, p_state.rem_bytes);
        text_chunks[buff_id].text_bytes = tmp_ck_size;
        text_chunks[buff_id].sep_sym = (text_chunk::size_type) p_state.sep_sym;
        text_chunks[buff_id].increase_capacity((tmp_ck_size*115)/100);
        text_chunks[buff_id].id = p_state.chunk_id++;

        read_chunk_from_file(p_state.fd_r, p_state.rem_bytes, p_state.read_bytes, text_chunks[buff_id]);
        buffers_to_process.push(buff_id);

#ifdef __linux__
        p_state.r_page_cache_bytes+=text_chunks[buff_id].e_bytes;
        if(p_state.r_page_cache_bytes>p_state.page_cache_limit){
            p_state.flush_page_cache();
        }
#endif
        buff_id++;
    }

    while (p_state.rem_bytes > 0 && input_frac<p_state.max_frac) {

        buffers_to_reuse.pop(buff_id);
        p_state.proc_syms+=text_chunks[buff_id].text_bytes;

        size_t new_byte_count = text_chunks[buff_id].gram.mem_usage();
        acc_bytes -= byte_counts[buff_id];
        acc_bytes += new_byte_count;
        byte_counts[buff_id] = new_byte_count;
        input_frac = float(acc_bytes)/float(p_state.f_size);

        std::cout<<"Parsed_input "<<report_space((off_t)p_state.proc_syms)<<", malloc_peak "<<report_space((off_t)malloc_count_peak())<<", byte_usage "<<report_space((off_t)text_chunks[buff_id].gram.mem_usage())<<", input_fraction "<<input_frac<<" "<<std::endl;
        //malloc_count_reset_peak();

        text_chunks[buff_id].text_bytes = tmp_ck_size;
        text_chunks[buff_id].id = p_state.chunk_id++;
        read_chunk_from_file(p_state.fd_r, p_state.rem_bytes, p_state.read_bytes, text_chunks[buff_id]);
        buffers_to_process.push(buff_id);

#ifdef __linux__
        p_state.r_page_cache_bytes+=text_chunks[buff_id].e_bytes;
        if(p_state.r_page_cache_bytes>p_state.page_cache_limit){
            p_state.flush_page_cache();
        }
#endif
    }

    //wait for the queue to be empty and close the input file
    while (!buffers_to_process.empty());
    buffers_to_process.done();

    //wait for all the parsers to finish
    while(parser_finished.load(std::memory_order_acquire)!=text_chunks.size());

    while(!buffers_to_reuse.empty()){
        buffers_to_reuse.pop(buff_id);
        p_state.proc_syms+=text_chunks[buff_id].text_bytes;

        size_t new_byte_count = text_chunks[buff_id].gram.mem_usage();
        acc_bytes -= byte_counts[buff_id];
        acc_bytes += new_byte_count;
        byte_counts[buff_id] = new_byte_count;
        input_frac = float(acc_bytes)/float(p_state.f_size);

        std::cout<<"Parsed_input "<<report_space((off_t)p_state.proc_syms)<<", malloc_peak "<<report_space((off_t)malloc_count_peak())<<", byte_usage "<<report_space((off_t)text_chunks[buff_id].gram.mem_usage())<<", input_fraction "<<input_frac<<std::endl;
        //malloc_count_reset_peak();

#ifdef __linux__
        p_state.r_page_cache_bytes+=text_chunks[buff_id].e_bytes;
        if(p_state.r_page_cache_bytes>p_state.page_cache_limit){
            p_state.flush_page_cache();
        }
#endif
    }
    buffers_to_reuse.done();

    for (auto &thread: threads) {
        thread.join();
    }
}

void build_partial_grammars(parsing_opts& p_opts, std::string& text_file) {

    plain_gram sink_gram(40, p_opts.sep_sym);
    parsing_state par_state(text_file, p_opts.sep_sym, p_opts.chunk_size, p_opts.page_cache_limit, p_opts.i_frac);

    //std::vector<text_chunk> text_chunks(p_opts.n_chunks, text_chunk(par_state.sink_gram));
    std::vector<text_chunk> text_chunks;
    text_chunks.reserve(p_opts.n_threads);
    for(size_t i=0;i<p_opts.n_threads;i++){
        text_chunks.emplace_back(sink_gram);
    }

    fill_chunk_grammars<false>(text_chunks, par_state);
    collapse_grams(sink_gram, text_chunks);
    par_state.sink_gram_mem_usage=sink_gram.eff_mem_usage();
    std::cout<<"Sink grammar space: tot:"<<report_space((off_t)sink_gram.mem_usage())<<" eff:"<<report_space((off_t)sink_gram.eff_mem_usage())<<std::endl;

    while(par_state.rem_bytes>0){
        fill_chunk_grammars<true>(text_chunks, par_state);
        collapse_grams(sink_gram, text_chunks);
        par_state.sink_gram_mem_usage=sink_gram.eff_mem_usage();
        std::cout<<"Sink grammar space: tot:"<<report_space((off_t)sink_gram.mem_usage())<<" eff:"<<report_space((off_t)sink_gram.eff_mem_usage())<<std::endl;
    }
    sink_gram.print_stats();
}

template<class sym_type>
void merge_partial_grammars(std::string& ct_p_grams_file, std::string& mg_p_gram_file,
                            std::vector<uint64_t>& p_seeds, size_t n_threads) {

    using p_gram_type = partial_gram<sym_type>;

    std::vector<std::pair<p_gram_type, std::vector<std::pair<size_t, size_t>>>> initial_grams(n_threads);

    ts_queue<size_t> gram_to_merge_queue;
    ts_queue<size_t> av_buff_queue;

    //std::string p_grams_file = tmp_ws.get_file("concatenated_grams");
    int fd_r = open(ct_p_grams_file.c_str(), O_RDONLY);
    size_t tot_bytes = file_size(ct_p_grams_file);
    size_t rem_bytes =  tot_bytes;
    size_t read_bytes;
    size_t i=0;

    while(i<n_threads && rem_bytes>0){

        read_bytes = initial_grams[i].first.load_from_fd(fd_r);
        //store the range of strings this partial gram covers
        initial_grams[i].second.push_back({initial_grams[i].first.txt_id,
                                           initial_grams[i].first.tot_strings()});
        rem_bytes-=read_bytes;

        //this is just to indicate that the first n_thread buffers are ready to be
        // loaded with n_thread grammars, which we will merge the initial grammars
        av_buff_queue.push(i);
        //

        i++;
    }

    n_threads = i;
    std::vector<p_gram_type> grams_to_merge(n_threads);
    size_t prog_bytes=0;

    auto gram_read_worker = [&](){

        if(rem_bytes==0){
            gram_to_merge_queue.done();
        }else{
            size_t buff_id;
            while(rem_bytes > 0){
                av_buff_queue.pop(buff_id);
                read_bytes = grams_to_merge[buff_id].load_from_fd(fd_r);
                gram_to_merge_queue.push(buff_id);
                rem_bytes-=read_bytes;
                prog_bytes+=read_bytes;
                //std::cout<<"Processed data: "<<double(prog_bytes)/double(tot_bytes)*100<<std::endl;
                std::cout<<"  Merged data: "<<double(prog_bytes)/double(tot_bytes)*100<<"% with peak "<<report_space((off_t)malloc_count_peak())<<std::endl;
            }
            while(av_buff_queue.size()<n_threads);
            gram_to_merge_queue.done();
        }

        //no longer needed
        for(auto &gram : grams_to_merge){
            gram.destroy_gram();
        }

        if(n_threads>1) {
            //merge the partial grams of the different threads into one
            size_t n_ranges=0;
            for (size_t i = 1; i < n_threads; i++)  n_ranges+=initial_grams[i].second.size();
            initial_grams[i].second.reserve(n_ranges);

            //this step is incorrect as it assumes all the initial grams are the same size,
            // but it gives us an estimate of the remaining %
            size_t rem_prog_bytes = INT_CEIL((tot_bytes-prog_bytes), n_threads);

            for (size_t i = 1; i < n_threads; i++) {
                merge_two_partial_grammars_in_memory<p_gram_type>(initial_grams[0].first, initial_grams[i].first, p_seeds);
                //store the range of strings this partial gram covers
                initial_grams[0].second.insert(initial_grams[0].second.end(), initial_grams[i].second.begin(), initial_grams[i].second.end());
                initial_grams[i].first.destroy_gram();
                destroy(initial_grams[i].second);

                prog_bytes+=rem_prog_bytes;
                std::cout<<"  Merged data: "<<double(prog_bytes)/double(tot_bytes)*100<<"% with peak "<<report_space((off_t)malloc_count_peak())<<std::endl;
                malloc_count_reset_peak();
            }
            initial_grams[0].second.shrink_to_fit();
#ifdef __linux__
             malloc_trim(0);
#endif
             initial_grams[0].first.reorder_strings(initial_grams[0].second);
        }

        store_to_file(mg_p_gram_file, initial_grams[0].first);
        initial_grams[0].first.destroy_gram();
#ifdef __linux__
        malloc_trim(0);
#endif
    };

    auto gram_merge_worker = [&](size_t idx){

        size_t buff_id;
        bool res;

        while (true) {
            res = gram_to_merge_queue.pop(buff_id);
            if (!res) break;
            merge_two_partial_grammars_in_memory<p_gram_type>(initial_grams[idx].first, grams_to_merge[buff_id], p_seeds);
            if(n_threads>1){//this step is to later reorder the strings in the final grammar
                initial_grams[idx].second.push_back({grams_to_merge[buff_id].txt_id,
                                                     grams_to_merge[buff_id].tot_strings()});
            }
            av_buff_queue.push(buff_id);
        }
    };

    std::vector<std::thread> threads;
    threads.emplace_back(gram_read_worker);
    for(size_t j=0;j<n_threads;j++){
        threads.emplace_back(gram_merge_worker, j);
    }

    for (auto &thread: threads) {
        thread.join();
    }
}

void lc_parsing_algo(std::string& i_file, std::string& o_file, size_t n_threads,
                     size_t n_chunks, size_t chunk_size, size_t par_seed, bool par_gram, float i_frac) {

    off_t f_size = file_size(i_file);
    parsing_opts p_opts;
    p_opts.chunk_size = chunk_size==0 ? off_t(ceil(0.005 * double(f_size))) : (off_t)chunk_size;
    //p_opts.chunk_size = std::min<off_t>(p_opts.chunk_size, std::numeric_limits<uint32_t>::max());//the chunks cannot exceed the 4GB by design
    p_opts.chunk_size = std::min<off_t>(p_opts.chunk_size, 1024*1024*200);//the chunks doest not exceed the 200MB by design

    size_t tot_chunks = INT_CEIL(f_size, p_opts.chunk_size);
    n_threads = std::min(n_threads, tot_chunks);

    //TODO remove the control over the number of chunks
    p_opts.n_chunks = n_chunks==0? n_threads+1 : n_chunks;
    p_opts.n_chunks = std::min<unsigned long>(p_opts.n_chunks, tot_chunks);

    p_opts.n_threads = n_threads;

    p_opts.page_cache_limit = 1024*1024*1024;
    p_opts.sep_sym = '\n';
    p_opts.orig_seed = par_seed;
    p_opts.i_frac = i_frac;

    std::mt19937 gen(par_seed); // Standard mersenne_twister_engine seeded with a fixed value
    std::uniform_int_distribution<uint64_t> distrib(1, std::numeric_limits<uint64_t>::max());
    p_opts.p_seeds.resize(32);
    for(size_t i=0;i<32;i++){
        //p_opts.p_seeds[i] = distrib(gen);
        p_opts.p_seeds[i] = 0;
    }

    std::cout<<"  Settings"<<std::endl;
    //std::cout<<"    Parsing mode              : short strings (<= 4 GBs)"<<std::endl;
    std::cout<<"    Parsing threads           : "<<p_opts.n_threads<<std::endl;
    //std::cout<<"    Parsing seed              : "<<p_opts.orig_seed<<std::endl;
    std::cout<<"    Active text chunks in RAM : "<<p_opts.n_chunks<<std::endl;
    std::cout<<"    Size of each chunk        : "<<report_space(p_opts.chunk_size)<<std::endl;
    std::cout<<"    Chunks' approx. mem usage : "<<report_space(off_t(((p_opts.chunk_size*115)/100)*p_opts.n_chunks))<<"\n"<<std::endl;

    build_partial_grammars(p_opts, i_file);

    //merge_many_grams_in_serial(ct_p_grams_file, p_opts.n_threads);
    /*std::string mg_p_gram_file = par_gram? o_file : tmp_ws.get_file("merged_p_grams");
    merge_partial_grammars<sym_type>(ct_p_grams_file, mg_p_gram_file, p_opts.p_seeds, p_opts.n_threads);

    if(!par_gram){
        gram_type final_grammar;
        partial2complete_gram(final_grammar, mg_p_gram_file, par_seed);
        store_to_file(o_file, final_grammar);
    }*/
    exit(0);
}
#endif //LCG_LZ_LIKE_LC_PARSING_H
