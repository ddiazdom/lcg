//
// Created by Diaz, Diego on 9.10.2024.
//

#ifndef LCG_COLLAPSE_TEMP_GRAM_H
#define LCG_COLLAPSE_TEMP_GRAM_H

#include "plain_gram.h"
#include "text_handler.h"

#ifdef DEBUG_MODE
    #define COLL_REPORT(lvl, output_set) \
        do{\
            if(!sink_set.empty()){\
                std::cout<<"  Growth_factor level_"<<lvl<<": "<<float(size_before)/float(output_set.size())<<std::endl;\
            }\
        }while(0);
#else
    #define COLL_REPORT(lvl, output_set) do{}while(0);
#endif

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
    COLL_REPORT(1, sink_set)
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
    COLL_REPORT(1, sink_set)
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
    }
    COLL_REPORT((round+1), sink_set)
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
    }
    COLL_REPORT((round+1), sink_set)
}

void add_compressed_strings(plain_gram& sink_gram, plain_gram& gram, std::vector<uint32_t>& prev_alphas_sink){
    size_t pos = sink_gram.comp_string.size();
    sink_gram.comp_string.resize(sink_gram.comp_string.size()+gram.comp_string.size());
    sink_gram.str_orders.reserve(sink_gram.str_orders.size()+gram.str_orders.size());

    size_t n=0;
    for(auto &str_data : gram.str_orders){
        uint64_t *i_map = gram.fps[str_data.lvl];
        uint64_t i_map_len = gram.fps_len[str_data.lvl];
        str_data.offset = pos;
        for(size_t i=0;i<str_data.n_strings;i++){
            if(gram.comp_string[n]>prev_alphas_sink[str_data.lvl]){
                size_t rank = gram.comp_string[n]-prev_alphas_sink[str_data.lvl];
                assert(rank<i_map_len);
                gram.comp_string[n] = i_map[rank]+1;
            }
            assert(gram.comp_string[n]>0);
            sink_gram.comp_string[pos++] = gram.comp_string[n++];
        }
        sink_gram.str_orders.push_back(str_data);
    }
    assert(pos==sink_gram.comp_string.size());
    assert(gram.comp_string.size()==n);

    gram.comp_string.clear();
    gram.str_orders.clear();
}

void collapse_grams(plain_gram& sink_gram, std::vector<text_chunk>& chunks) {

    //compute the length of each rule set
    std::vector<uint32_t> prev_lvl_alpha(sink_gram.fps.size(), 0);
    prev_lvl_alpha[1] = sink_gram.ter_dict.size();
    for(size_t i=2;i<prev_lvl_alpha.size();i++){
        prev_lvl_alpha[i] = sink_gram.nt_dicts[i-2].size();
    }

    //just swap grammars if the sink grammar is empty
    if(sink_gram.empty()){
        sink_gram.swap(chunks[0].gram);
    }

    //destroy the tables to free some working memory
    for(auto & text_chunk : chunks){
        text_chunk.gram.destroy_tables();
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

    //add the compressed strings
    for(auto &chunk : chunks){
        add_compressed_strings(sink_gram, chunk.gram, prev_lvl_alpha);
    }

    //reset the fingerprints. I can't do it before because the symbols
    // in the compressed strings have different levels
    for(auto &chunk : chunks){
        chunk.gram.clear_fps();
    }

    //TODO remove later
    std::cout<<"The new sink grammar "<<std::endl;
    sink_gram.print_stats();
    //
}
#endif //LCG_COLLAPSE_TEMP_GRAM_H
