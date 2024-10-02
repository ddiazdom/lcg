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

void mul_thread_ter_collapse(std::vector<text_chunk>& chunks){

    auto ter_worker = [&](phrase_set<uint8_t>& sink, phrase_set<uint8_t>& coll_set, uint64_t*& o_map, size_t& o_map_len) {
        assert(o_map_len==(coll_set.size()+1));
        coll_set.destroy_table();
        auto it = coll_set.begin();
        auto it_end = coll_set.end();
        bool found;
        uint32_t mt;
        size_t phrase_id=1;
        uint64_t flag = 0x8000000000000000UL;
        while(it!=it_end){
            auto phr = *it;
            found = sink.find(phr.phrase, phr.len, mt);
            if(found){
                o_map[phrase_id] = mt;
            }else{
                o_map[phrase_id] = it.pos() | flag;
            }
            phrase_id++;
            ++it;
        }
        assert(phrase_id==o_map_len);
    };

    std::vector<std::thread> threads;
    std::vector<size_t> found_per_thread(chunks.size(), 0);
    for (size_t i = 1; i < chunks.size(); i++) {
        threads.emplace_back(ter_worker,
                             std::ref(chunks[0].ter_dict), std::ref(chunks[i].ter_dict),
                             std::ref(chunks[i].fps[1]), std::ref(chunks[i].fps_len[1]));
    }

    for (auto &thread: threads) {
        thread.join();
    }

    uint32_t len;
    size_t pos;
    size_t size_before = chunks[0].ter_dict.size();
    for (size_t i = 1; i < chunks.size(); i++) {
        if(!chunks[i].ter_dict.empty()){
            const uint8_t *stream = chunks[i].ter_dict.phr_stream();
            for(size_t j=0;j<chunks[i].fps_len[1];j++){
                if(chunks[i].fps[1][j] & 0x8000000000000000UL){
                    pos = chunks[i].fps[1][j] & 0x7FFFFFFFFFFFFFFF;
                    memcpy(&len, &stream[pos], sizeof(uint32_t));
                    pos+=sizeof(uint32_t);
                    chunks[i].fps[1][j] = chunks[0].ter_dict.insert(&stream[pos], len);
                }
            }
            chunks[i].ter_dict.clear();
        }
    }

    if(!chunks[0].ter_dict.empty()){
        std::cout<<"  Growth factor 1 "<<float(size_before)/float(chunks[0].ter_dict.size())<<std::endl;
    }
}

void mul_thread_nt_collapse(std::vector<text_chunk>& chunks, size_t round){

    auto nt_worker = [&](
            phrase_set<uint32_t>& sink_set, phrase_set<uint32_t>& coll_set,
            uint64_t*& i_map, size_t& i_map_len,
            uint64_t*& o_map, size_t& o_map_len) {

        assert(o_map_len==(coll_set.size()+1));
        coll_set.destroy_table();
        auto it = coll_set.begin();
        auto it_end = coll_set.end();
        bool found;
        uint32_t mt;
        size_t phrase_id=1;
        uint64_t flag = 0x8000000000000000UL;

        while(it!=it_end){
            auto phr = *it;
            for(size_t j=0;j<phr.len;j++){
                assert(phr.phrase[j]>0 && phr.phrase[j]<i_map_len);
                phr.phrase[j] = i_map[phr.phrase[j]]+1;
            }
            found = sink_set.find(phr.phrase, phr.len, mt);
            if(found){
                o_map[phrase_id] = mt;
            }else{
                o_map[phrase_id] = it.pos() | flag;
            }
            phrase_id++;
            ++it;
        }
        assert(phrase_id==o_map_len);
    };

    std::vector<std::thread> threads;
    std::vector<size_t> found_per_thread(chunks.size(), 0);
    for (size_t i = 1; i < chunks.size(); i++) {
        threads.emplace_back(nt_worker,
                             std::ref(chunks[0].nt_dicts[round-1]), std::ref(chunks[i].nt_dicts[round-1]),
                             std::ref(chunks[i].fps[round]), std::ref(chunks[i].fps_len[round]),
                             std::ref(chunks[i].fps[round+1]), std::ref(chunks[i].fps_len[round+1]));
    }

    for (auto &thread: threads) {
        thread.join();
    }


    uint32_t len;
    size_t pos;
    size_t size_before = chunks[0].nt_dicts[round-1].size();
    for (size_t i = 1; i < chunks.size(); i++) {
        if(!chunks[i].nt_dicts[round-1].empty()){

            const uint32_t *stream = chunks[i].nt_dicts[round-1].phr_stream();
            for(size_t j=0;j<chunks[i].fps_len[round+1];j++){
                if(chunks[i].fps[round+1][j] & 0x8000000000000000UL){
                    pos = chunks[i].fps[round+1][j] & 0x7FFFFFFFFFFFFFFF;
                    memcpy(&len, &stream[pos], sizeof(uint32_t));
                    pos++;
                    chunks[i].fps[round+1][j] = chunks[0].nt_dicts[round-1].insert(&stream[pos], len);
                }
            }

            chunks[i].nt_dicts[round-1].clear();
            chunks[i].fps[round]= mem<uint64_t>::reallocate(chunks[i].fps[round], 1);
            chunks[i].fps_len[round] = 1;
        }
    }

    if(!chunks[0].nt_dicts[round-1].empty()){
        std::cout<<"  Growth factor "<<round+1<<" "<<float(size_before)/float(chunks[0].nt_dicts[round-1].size())<<std::endl;
    }
}

void sin_thread_nt_collapse(std::vector<text_chunk>& chunks, size_t round){
    size_t size_before = chunks[0].nt_dicts[round-1].size();
    for(size_t i=1;i<chunks.size();i++){
        if(!chunks[i].nt_dicts[round-1].empty()){
            chunks[0].nt_dicts[round-1].absorb_set(chunks[i].nt_dicts[round-1],
                                                   chunks[i].fps[round], chunks[i].fps_len[round],
                                                   chunks[i].fps[round+1], chunks[i].fps_len[round+1]);
            chunks[i].nt_dicts[round-1].clear();
            chunks[i].fps[round]= mem<uint64_t>::reallocate(chunks[i].fps[round], 1);
            chunks[i].fps_len[round] = 1;
        }
    }

    if(!chunks[0].nt_dicts[round-1].empty()){
        std::cout<<"  Growth factor "<<round+1<<" "<<float(size_before)/float(chunks[0].nt_dicts[round-1].size())<<std::endl;
    }
}

void collapse_grams(std::vector<text_chunk>& text_chunks) {

    size_t tot_strings = 0;
    for(auto & text_chunk : text_chunks){
        tot_strings+=text_chunk.comp_string.size();
        text_chunk.get_gram_levels();
    }

    mul_thread_ter_collapse(text_chunks);
    for(size_t round=1;round<=text_chunks[0].nt_dicts.size();round++){
        if(round<=4 && text_chunks.size()>1){
            mul_thread_nt_collapse(text_chunks, round);
        }else{
            sin_thread_nt_collapse(text_chunks, round);
        }
    }

    size_t pos = text_chunks[0].comp_string.size();
    text_chunks[0].comp_string.resize(tot_strings);
    for(size_t i=1;i<text_chunks.size();i++){
        size_t n = text_chunks[i].n_levels;
        for(size_t j=0;j<text_chunks[i].comp_string.size();j++){
            //TODO fix this
            //assert(text_chunks[i].comp_string[j]>0 && text_chunks[i].comp_string[j]<text_chunks[i].fps_len[n]);
            //text_chunks[0].comp_string[pos++] = text_chunks[i].fps[n][text_chunks[i].comp_string[j]]+1;
        }
    }

    std::cout<<"Final stats: "<<std::endl;
    std::cout<<"Level 1\t Phrases: "<<text_chunks[0].ter_dict.size()<<" Size: "<<text_chunks[0].ter_dict.tot_symbols()<<std::endl;

    text_chunks[0].ter_dict.update_fps(text_chunks[0].fps[0], text_chunks[0].fps_len[0], text_chunks[0].fps[1], text_chunks[0].fps_len[1]);
    size_t round=1;
    while(!text_chunks[0].nt_dicts[round-1].empty()){
        text_chunks[0].nt_dicts[round-1].update_fps(text_chunks[0].fps[round], text_chunks[0].fps_len[round],
                                                    text_chunks[0].fps[round+1], text_chunks[0].fps_len[round+1]);
        std::cout<<"Level "<<round+1<<"\t Phrases: "<<text_chunks[0].nt_dicts[round-1].size()<<" Size: "<<text_chunks[0].nt_dicts[round-1].tot_symbols()<<std::endl;
        round++;
    }
    std::cout<<"Tot. strings: "<<text_chunks[0].comp_string.size()<<std::endl;

    /*std::cout<<"Merging the sets"<<std::endl;
    std::vector<uint32_t> sym_map;
    for(size_t i=1;i<text_chunks.size();i++){
        size_t size_before = text_chunks[0].ter_dict.size();
        size_t size_after = text_chunks[0].ter_dict.absorb_set(text_chunks[i].ter_dict, sym_map);
        std::cout<<"  Growth factor 1 "<<float(size_before)/float(size_after)<<std::endl;
        size_t l=0;
        while(!text_chunks[i].nt_dicts[l].empty()){
            size_before = text_chunks[0].nt_dicts[l].size();
            size_after = text_chunks[0].nt_dicts[l].absorb_set(text_chunks[i].nt_dicts[l], sym_map);
            std::cout<<"  Growth factor "<<l+2<<" "<<float(size_before)/float(size_after)<<std::endl;
            l++;
        }

        size_t pos = text_chunks[0].comp_string.size();
        text_chunks[0].comp_string.resize(pos+text_chunks[i].comp_string.size());
        for(size_t j=0;j<text_chunks[i].comp_string.size();j++){
            text_chunks[0].comp_string[pos++] = sym_map[text_chunks[i].comp_string[j]]+1;
        }
        destroy(sym_map);
    }

    text_chunks[0].ter_dict.update_fps(text_chunks[0].fps[0], text_chunks[0].fps_len[0], text_chunks[0].fps[1], text_chunks[0].fps_len[1]);
    std::cout<<"Final stats: "<<std::endl;
    std::cout<<"Level 1\t Phrases: "<<text_chunks[0].ter_dict.size()<<" Size: "<<text_chunks[0].ter_dict.tot_symbols()<<std::endl;
    size_t l=0, lvl=1;
    while(!text_chunks[0].nt_dicts[l].empty()){

        text_chunks[0].nt_dicts[l].update_fps(text_chunks[0].fps[lvl], text_chunks[0].fps_len[lvl],
                                              text_chunks[0].fps[lvl+1], text_chunks[0].fps_len[lvl+1]);

        std::cout<<"Level "<<l+2<<"\t Phrases: "<<text_chunks[0].nt_dicts[l].size()<<" Size: "<<text_chunks[0].nt_dicts[l].tot_symbols()<<std::endl;
        l++;
        lvl++;
    }
    std::cout<<"Tot. strings: "<<text_chunks[0].comp_string.size()<<std::endl;*/
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

void byte_par_r2l(text_chunk& chunk, off_t& n_strings, size_t sep_sym) {

    uint8_t * text = chunk.text;
    off_t lb, rb = chunk.text_bytes-1, i=chunk.text_bytes-2, byte_offset, parse_size;
    uint8_t v_len;
    assert(text[i+1]==sep_sym && text[i]>text[i+1]);

    text[rb] = 128;//the vbyte code of the metasymbol mt=0 representing a separator in the next round of parsing
    n_strings=1;
    parse_size=4;//we will count in bytes of sizeof(uint32_t)
    off_t mbo[2]={4};

    //given a phrase T[a..b-1], byte_offset tells the number of bytes (4 per mt) used by the metasymbols after T[b-1].
    //This value defines the new size for text, now with the parse
    //max_byte_offset tells the maximum offset

    std::vector<phrase_overflow> phr_with_ovf;
    bool r_cmp = true, l_cmp, new_str;
    uint32_t mt_sym, phrase_len;
    uint8_t mid_sym = text[i];
    while(--i>0 && text[i]==mid_sym);

    while(i>=0){
        l_cmp = chunk.fps[chunk.round][text[i]]>chunk.fps[chunk.round][mid_sym];
        if(l_cmp && !r_cmp){
            lb = i+1;
            new_str = mid_sym==sep_sym;
            lb += new_str;
            phrase_len=rb-lb;

            byte_offset = rb + parse_size;
            mbo[byte_offset > mbo[1]] = byte_offset;

            mt_sym = chunk.ter_dict.insert(&text[lb], phrase_len) + 1;
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

    mt_sym = chunk.ter_dict.insert(&text[lb], phrase_len)+1;
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

    chunk.ter_dict.update_fps(chunk.fps[chunk.round], chunk.fps_len[chunk.round],
                              chunk.fps[chunk.round+1], chunk.fps_len[chunk.round+1]);
    finish_byte_parse(chunk, max_byte_offset, phr_with_ovf);
}

// this method parses the text and store the parse in the text itself.
// It only works for parsing rounds other than the first one because the length of symbol each
// cell is the same as the length of cell where we store the metasymbols, so there is no overflow
void int_par_l2r(text_chunk& chunk){

    uint32_t *text = chunk.parse;
    uint32_t mt_sym, sep_sym=0, txt_size = chunk.parse_size;
    uint32_t left_sym, middle_sym;
    uint64_t left_hash, middle_hash, right_hash;
    off_t i=0, parse_size = 0, phrase_len, lb, rb;
    bool new_str=false;

    lb = 0;
    left_sym = text[i];
    left_hash = chunk.fps[chunk.round][left_sym];
    while(++i<txt_size && text[i]==left_sym);
    assert(i<txt_size);

    middle_sym = text[i];
    middle_hash = chunk.fps[chunk.round][middle_sym];
    rb=i;
    while(++i<txt_size && text[i]==middle_sym);

    while(i<txt_size) {
        right_hash = chunk.fps[chunk.round][text[i]];
        if(left_hash>middle_hash && middle_hash<right_hash){//local minimum

            phrase_len = rb-lb;
            mt_sym = chunk.nt_dicts[chunk.round-1].insert(&text[lb], phrase_len);

            text[parse_size]=sep_sym;
            parse_size+=new_str;
            text[parse_size++] = mt_sym+1;
            new_str = text[rb]==sep_sym;
            lb = rb+new_str;
        }

        left_hash = middle_hash;
        middle_hash = right_hash;
        middle_sym = text[i];
        rb = i;
        while(++i<txt_size && text[i]==middle_sym);
    }
    assert(rb==(txt_size-1) && text[rb]==sep_sym);

    phrase_len = rb-lb;
    mt_sym = chunk.nt_dicts[chunk.round-1].insert(&text[lb], phrase_len);

    text[parse_size]=sep_sym;
    parse_size+=new_str;
    text[parse_size++] = mt_sym+1;
    text[parse_size++] = sep_sym;
    chunk.nt_dicts[chunk.round-1].update_fps(chunk.fps[chunk.round], chunk.fps_len[chunk.round],
                                             chunk.fps[chunk.round+1], chunk.fps_len[chunk.round+1]);
    chunk.parse_size = parse_size;
}

void compress_text_chunk(text_chunk& chunk){

    off_t n_strings=0;
    size_t sep_sym = chunk.sep_sym;
    chunk.round = 0;

    //auto start = std::chrono::steady_clock::now();
    byte_par_r2l(chunk, n_strings, sep_sym);
    //auto end = std::chrono::steady_clock::now();
    //report_time(start, end , 2);

    off_t size_limit = n_strings*2;
    chunk.round++;

    while(chunk.parse_size!=size_limit){
        assert(chunk.parse_size>=size_limit);
        //start = std::chrono::steady_clock::now();
        int_par_l2r(chunk);
        //end = std::chrono::steady_clock::now();
        //report_time(start, end , 2);
        chunk.round++;
    }

    //start = std::chrono::steady_clock::now();
    size_t pos = chunk.comp_string.size();
    chunk.comp_string.resize(pos+(chunk.parse_size/2));
    for(off_t i=0;i<chunk.parse_size;i+=2){
        assert(i==0 || chunk.parse[i-1]==0);
        chunk.comp_string[pos++] = chunk.parse[i];
    }
    //chunk.p_gram.add_compressed_string(chunk.parse, chunk.parse_size);
    //end = std::chrono::steady_clock::now();
    //report_time(start, end , 2);
}

void build_partial_grammars(parsing_opts& p_opts, std::string& text_file, std::string& ct_p_grams_file) {

    ts_queue<size_t> buffers_to_process;
    ts_queue<size_t> buffers_to_reuse;

    std::atomic<size_t> parser_finished{0};
    std::vector<text_chunk> text_chunks(p_opts.n_chunks);

    auto io_worker = [&]() -> void {

        int fd_r = open(text_file.c_str(), O_RDONLY);

        size_t f_size = file_size(text_file);
#ifdef __linux__
        off_t r_page_cache_bytes = 0, w_page_cache_bytes = 0;
        posix_fadvise(fd_r, 0, f_size, POSIX_FADV_SEQUENTIAL);
#endif

        off_t rem_bytes = (off_t) f_size, r_acc_bytes = 0;
        size_t chunk_id = 0;

        auto tmp_ck_size = off_t(INT_CEIL(p_opts.chunk_size, sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type));

        while (chunk_id < p_opts.n_chunks && rem_bytes > 0) {

            tmp_ck_size = std::min(tmp_ck_size, rem_bytes);

            text_chunks[chunk_id].text_bytes = tmp_ck_size;
            text_chunks[chunk_id].sep_sym = (text_chunk::size_type) p_opts.sep_sym;
            text_chunks[chunk_id].increase_capacity((tmp_ck_size*115)/100);
            text_chunks[chunk_id].id = chunk_id;

            read_chunk_from_file(fd_r, rem_bytes, r_acc_bytes, text_chunks[chunk_id]);
            buffers_to_process.push(chunk_id);
#ifdef __linux__
            r_page_cache_bytes+=text_chunks[chunk_id].e_bytes;
            if(r_page_cache_bytes>p_opts.page_cache_limit){
                std::cout<<"removing from page cache "<<r_page_cache_bytes<<" "<<r_acc_bytes<<std::endl;
                posix_fadvise(fd_r, r_acc_bytes-r_page_cache_bytes, r_page_cache_bytes, POSIX_FADV_DONTNEED);
                r_page_cache_bytes=0;
            }
#endif
            chunk_id++;
       }

       size_t buff_id;
       size_t proc_syms=0;
       while (rem_bytes > 0) {
           buffers_to_reuse.pop(buff_id);
           proc_syms+=text_chunks[buff_id].text_bytes;
           std::cout<<"Parsed input "<<report_space((off_t)proc_syms)<<" with peak "<<report_space((off_t)malloc_count_peak())<<" and byte usage "<<report_space((off_t)text_chunks[buff_id].mem_usage())<<" and eff_byte usage "<<report_space((off_t)text_chunks[buff_id].eff_mem_usage())<<" buff available "<<report_space((off_t)text_chunks[buff_id].dict_buff_av())<<std::endl;
           malloc_count_reset_peak();
           text_chunks[buff_id].text_bytes = tmp_ck_size;
           text_chunks[buff_id].id = chunk_id++;

           read_chunk_from_file(fd_r, rem_bytes, r_acc_bytes, text_chunks[buff_id]);
           buffers_to_process.push(buff_id);
#ifdef __linux__
           r_page_cache_bytes+=text_chunks[buff_id].e_bytes;
           if(r_page_cache_bytes>p_opts.page_cache_limit){
               std::cout<<"removing from page cache "<<r_page_cache_bytes<<" "<<r_acc_bytes<<std::endl;
               posix_fadvise(fd_r, r_acc_bytes-r_page_cache_bytes, r_page_cache_bytes, POSIX_FADV_DONTNEED);
               r_page_cache_bytes=0;
           }
#endif
       }

       //wait for the queue to be empty and close the input file
       while (!buffers_to_process.empty());
       buffers_to_process.done();
#ifdef __linux__
       posix_fadvise(fd_r, 0, f_size, POSIX_FADV_DONTNEED);
#endif
       close(fd_r);

       //wait for all the parsers to finish
       while(parser_finished.load(std::memory_order_acquire)!=p_opts.n_threads);

       //store the remaining grammars in the temporary file
       while(!buffers_to_reuse.empty()){
           buffers_to_reuse.pop(buff_id);
           proc_syms+=text_chunks[buff_id].text_bytes;
           std::cout<<"Parsed input "<<report_space((off_t)proc_syms)<<" with peak "<<report_space((off_t)malloc_count_peak())<<" and byte usage "<<report_space((off_t)text_chunks[buff_id].mem_usage())<<" and eff_byte usage "<<report_space((off_t)text_chunks[buff_id].eff_mem_usage())<<" buff available "<<report_space((off_t)text_chunks[buff_id].dict_buff_av())<<std::endl;
           malloc_count_reset_peak();
       }
       buffers_to_reuse.done();
       collapse_grams(text_chunks);
       std::cout<<"\nParsed finished "<<std::endl;
    };

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
            compress_text_chunk(text_chunks[buff_id]);
            memset(text_chunks[buff_id].text, 0, text_chunks[buff_id].buffer_bytes);
            text_chunks[buff_id].t_end = std::chrono::steady_clock::now();
            buffers_to_reuse.push(buff_id);
        }
    };

    std::vector<std::thread> threads;
    threads.emplace_back(io_worker);

    for (size_t i = 0; i < p_opts.n_threads; i++) {
        threads.emplace_back(compressor_worker);
    }

    for (auto &thread: threads) {
        thread.join();
    }

    for(auto &chunk : text_chunks){
        //free(chunk.text);
        mem<uint8_t>::deallocate(chunk.text);
        chunk.text=nullptr;
    }
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
                     size_t n_chunks, size_t chunk_size, size_t par_seed, bool par_gram) {

    off_t f_size = file_size(i_file);
    parsing_opts p_opts;
    p_opts.chunk_size = chunk_size==0 ? off_t(ceil(0.005 * double(f_size))) : (off_t)chunk_size;
    //p_opts.chunk_size = std::min<off_t>(p_opts.chunk_size, std::numeric_limits<uint32_t>::max());//the chunks cannot exceed the 4GB by design
    p_opts.chunk_size = std::min<off_t>(p_opts.chunk_size, 1024*1024*200);//the chunks doest not exceed the 200MB by design

    size_t tot_chunks = INT_CEIL(f_size, p_opts.chunk_size);
    n_threads = std::min(n_threads, tot_chunks);

    p_opts.n_chunks = n_chunks==0? n_threads+1 : n_chunks;
    p_opts.n_chunks = std::min<unsigned long>(p_opts.n_chunks, tot_chunks);

    p_opts.n_threads = n_threads;

    p_opts.page_cache_limit = 1024*1024*1024;
    p_opts.sep_sym = '\n';
    p_opts.orig_seed = par_seed;

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

    build_partial_grammars(p_opts, i_file, o_file);

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

/*
    off_t byte_p_round_bck(uint8_t*& text, off_t txt_size, off_t& buffer_size, uint32_t*& parse,
                           off_t& n_strings, size_t sep_sym, uint64_t fp_seed,
                           std::vector<uint64_t>& prev_fps, partial_gram<uint8_t>& p_gram){

        size_t prev_sym, curr_sym, next_sym;
        uint64_t prev_hash, curr_hash, next_hash;

        text_chunk::size_type mt_sym;
        off_t txt_pos = txt_size-1, parse_size = 0, lb, rb;
        lz_like_map map(text);
        uint32_t phrase_len;

        off_t parse_distance=0, offset;
        std::vector<phrase_overflow> phr_with_ovf;
        uint8_t v_len;

        bool inserted;
        n_strings = 0;

        while(txt_pos>0){

            assert(text[txt_pos]==sep_sym);

            offset = txt_pos+1+(parse_size*4);
            parse_distance = offset>parse_distance ? offset : parse_distance ;

            text[txt_pos] = 128;//vbyte code for 0 (the separator symbol in the next levels)
            parse_size++;
            n_strings++;

            rb = txt_pos;
            prev_sym = text[--txt_pos];

            curr_sym = text[--txt_pos];
            while(curr_sym==prev_sym && txt_pos>0) curr_sym = text[--txt_pos];

            if(curr_sym==sep_sym){
                lb = txt_pos+1;
                phrase_len = rb-lb;

                offset = rb+(parse_size*4);
                parse_distance = offset>parse_distance ? offset : parse_distance ;
                //std::cout<<lb<<" "<<phrase_len<<" "<<map.size()<<" "<<parse_size<<std::endl;
                mt_sym = map.insert(lb, phrase_len, inserted)+1;

                v_len = vbyte_len(mt_sym);
                if(__builtin_expect(v_len>phrase_len, 0)){
                    phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
                }else if(!inserted){
                    vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
                    memset(&text[lb+v_len], 0, phrase_len-v_len);
                    //test
                    //uint32_t tmp;
                    //vbyte_decoder<uint32_t>::read_right2left(&text[lb+v_len-1], tmp);
                    //assert(tmp==mt_sym);
                    //breaks.emplace_back(lb+v_len-1, parse_size);
                    //
                }
                parse_size++;
                continue;
            }

            next_sym = text[--txt_pos];
            while(next_sym==curr_sym && txt_pos>0) next_sym = text[--txt_pos];

            prev_hash = prev_fps[prev_sym];
            curr_hash = prev_fps[curr_sym];

            bool local_minimum;

            while(next_sym!=sep_sym && txt_pos>0){

                next_hash = prev_fps[next_sym];
                local_minimum = prev_hash>curr_hash && curr_hash<next_hash;
                if(local_minimum){
                    lb = txt_pos+1;
                    phrase_len = rb-lb;

                    offset = rb+(parse_size*4);
                    parse_distance = offset>parse_distance ? offset : parse_distance ;

                    //std::cout<<lb<<" "<<phrase_len<<" "<<map.size()<<" "<<parse_size<<std::endl;
                    mt_sym = map.insert(lb, phrase_len, inserted)+1;
                    v_len = vbyte_len(mt_sym);
                    if(__builtin_expect(v_len>phrase_len, 0)){
                        phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
                    }else if(!inserted){
                        vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
                        //test
                        //uint32_t tmp;
                        //vbyte_decoder<uint32_t>::read_right2left(&text[lb+v_len-1], tmp);
                        //assert(tmp==mt_sym);
                        //breaks.emplace_back(lb+v_len-1, parse_size);
                        //
                        memset(&text[lb+v_len], 0, phrase_len-v_len);
                    }
                    parse_size++;
                    rb = lb;
                }

                prev_hash = curr_hash;
                curr_hash = next_hash;

                curr_sym = next_sym;
                next_sym = text[--txt_pos];
                while(next_sym==curr_sym && txt_pos>0) next_sym = text[--txt_pos];
            }

            if(txt_pos==0){
                next_hash = prev_fps[next_sym];
                local_minimum = prev_hash>curr_hash && curr_hash<next_hash;
                if(local_minimum){
                    lb = txt_pos+1;
                    phrase_len = rb-lb;

                    offset = rb+(parse_size*4);
                    parse_distance = offset>parse_distance ? offset : parse_distance ;

                    //std::cout<<lb<<" "<<phrase_len<<" "<<map.size()<<" "<<parse_size<<std::endl;
                    mt_sym = map.insert(lb, phrase_len, inserted)+1;
                    v_len = vbyte_len(mt_sym);
                    if(__builtin_expect(v_len>phrase_len, 0)){
                        phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
                    }else if(!inserted){
                        vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
                        memset(&text[lb+v_len], 0, phrase_len-v_len);
                        //test
                        //uint32_t tmp;
                        //vbyte_decoder<uint32_t>::read_right2left(&text[lb+v_len-1], tmp);
                        //assert(tmp==mt_sym);
                        //breaks.emplace_back(lb+v_len-1, parse_size);
                        //
                    }
                    parse_size++;
                    rb = lb;
                }
                lb = 0;
            }else{
                lb=txt_pos+1;
            }

            phrase_len = rb-lb;

            offset = rb+(parse_size*4);
            parse_distance = offset>parse_distance ? offset : parse_distance ;

            //std::cout<<lb<<" "<<phrase_len<<" "<<map.size()<<" "<<parse_size<<std::endl;
            mt_sym = map.insert(lb, phrase_len, inserted)+1;
            v_len = vbyte_len(mt_sym);
            if(__builtin_expect(v_len>phrase_len, 0)){
                phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
            }else if(!inserted){
                vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
                memset(&text[lb+v_len], 0, phrase_len-v_len);
                //test
                //uint32_t tmp;
                //vbyte_decoder<uint32_t>::read_right2left(&text[lb+v_len-1], tmp);
                //assert(tmp==mt_sym);
                //breaks.emplace_back(lb+v_len-1, parse_size);
                //
            }

            parse_size++;

            //str_len = prev_sep_sym-txt_pos-1;
            //longest_str = str_len>longest_str ? str_len : longest_str;
            //prev_sep_sym = txt_pos;
        }

        map.shrink_to_fit();
        map.destroy_table();

        std::vector<uint32_t> perm;
        create_meta_sym<uint8_t, true>(perm, fp_seed, map.phrase_set, text, txt_size, prev_fps, p_gram);
        finish_byte_parse(text, txt_size, buffer_size, parse_size, map, parse_distance, perm, parse, phr_with_ovf);
        return parse_size;
    }
    */

/*
// this method parses the text and store the parse in the text itself.
// It only works for parsing rounds other than the first one because the length of symbol each
// cell is the same as the length of cell where we store the metasymbols, so there is no overflow
off_t int_p_round_fwd(uint32_t* text, off_t txt_size, off_t& n_strings, uint64_t fp_seed,
                      std::vector<uint64_t>& prev_fps, partial_gram<uint8_t>& p_gram) {

    uint32_t mt_sym, sep_sym=0;
    size_t prev_sym, curr_sym, next_sym, dummy_sym=std::numeric_limits<text_chunk::size_type>::max();
    off_t txt_pos = 0, parse_size = 0, phrase_len, lb, rb;
    lz_like_map<uint32_t> map(text);

    bool inserted;
    n_strings = 0;
    off_t sym_bytes = sizeof(uint32_t);

    while(txt_pos<txt_size) {

        lb = txt_pos;
        prev_sym = text[txt_pos++];

        curr_sym = text[txt_pos++];
        while(curr_sym==prev_sym) curr_sym = text[txt_pos++];
        rb = txt_pos-1;

        if(curr_sym==sep_sym){
            phrase_len = rb-lb;
            mt_sym = map.insert(lb, phrase_len, inserted);
            assert(text[rb]==sep_sym);
            if(!inserted){
                //we can not replace the first phrase occurrence as we use it as source for the dictionary
                assert(text[lb]!=dummy_sym);
                text[lb] = mt_sym+1;//store the metasymbol in the first phrase position
                memset(&text[lb+1], (int)dummy_sym, sym_bytes*(phrase_len-1));//pad the rest of the phrase with dummy symbols
            }
            text[rb] = sep_sym;
            parse_size+=2;//+1 for the separator symbol
            n_strings++;
            continue;
        }

        next_sym = text[txt_pos++];
        while(next_sym==curr_sym) next_sym = text[txt_pos++];

        while(next_sym!=sep_sym){
            uint32_t tmp_p = prev_sym;
            uint32_t tmp_c = curr_sym;
            uint32_t tmp_n = next_sym;

            if(tmp_p>tmp_c && tmp_c<tmp_n){//local minimum
                phrase_len = rb-lb;
                mt_sym = map.insert(lb, phrase_len, inserted);
                //we can not replace the first phrase occurrence as we use it as source for the dictionary
                if(!inserted){
                    assert(text[lb]!=dummy_sym);
                    text[lb] = mt_sym+1;
                    memset(&text[lb+1], (int)dummy_sym, sym_bytes*(phrase_len-1));
                }
                parse_size++;
                lb = rb;
            }
            rb = txt_pos-1;
            prev_sym = curr_sym;
            curr_sym = next_sym;
            next_sym = text[txt_pos++];
            while(next_sym==curr_sym) next_sym = text[txt_pos++];
        }

        phrase_len = txt_pos-1-lb;
        assert(text[lb+phrase_len]==sep_sym);
        mt_sym = map.insert(lb, phrase_len, inserted);
        //we can not replace the first phrase occurrence as we use it as source for the dictionary
        if(!inserted){
            assert(text[lb]!=dummy_sym);
            text[lb] = mt_sym+1;//store the metasymbol in the first phrase position
            memset(&text[lb+1], (int)dummy_sym, sym_bytes*(phrase_len-1));//pad the rest of the phrase with dummy symbols
        }
        text[lb+phrase_len] = sep_sym;
        parse_size+=2;//+1 for the separator symbol
        n_strings++;
    }

    map.shrink_to_fit();
    map.destroy_table();

    assert(map.phrase_set.size()<dummy_sym);
    std::vector<uint32_t> mt_perm;
    create_meta_sym<uint32_t, false>(mt_perm, fp_seed, map.phrase_set, text, txt_size, prev_fps, p_gram);

    // create the parse in place
    map.insert_dummy_entry({uint32_t(txt_size), 0, false, false});
    size_t tot_phrases = map.phrase_set.size()-1;//do not count the dummy
    mt_sym = 0, lb = 0;
    off_t i=0, k=0;

    while(mt_sym<tot_phrases) {
        assert(i==lb);
        text[k++] = mt_perm[mt_sym+1];
        i+= map.phrase_set[mt_sym].len;//move out of the phrase boundary

        mt_sym++;
        lb = map.phrase_set[mt_sym].source;//position for the next phrase

        while(i<lb){//process the text area between consecutive phrases
            assert(text[i]<mt_perm.size());
            text[k++] = mt_perm[text[i]];
            i++;
            while(text[i]==dummy_sym && i<lb) i++;
        }
    }

    assert(k==parse_size);
    return parse_size;
}*/
#endif //LCG_LZ_LIKE_LC_PARSING_H
