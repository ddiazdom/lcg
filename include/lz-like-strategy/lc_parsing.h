//
// Created by Diaz, Diego on 17.9.2023.
//

#ifndef LCG_LZ_LIKE_LC_PARSING_H
#define LCG_LZ_LIKE_LC_PARSING_H

#include <random>
#include <fcntl.h>
#include <cstring>

#include "cds/macros.h"
#include "cds/ts_queue.h"
#include "cds/ts_priority_queue.h"
#include "cds/utils.h"
#include "cds/cdt_common.hpp"

#include "text_handler.h"
#include "lz_like_map.h"
#include "grammar.h"
#include "cds/mmap_allocator.h"

namespace lz_like_strat {

    struct parsing_opts{
        size_t n_threads{};
        unsigned long n_chunks{};
        off_t chunk_size{};
        off_t page_cache_limit{};
        size_t sep_sym{};
        std::vector<uint64_t> p_seeds;
    };

    struct heap_node{
        size_t chunk_id;
        size_t buff_idx;
    };
    struct heap_node_compare{
        bool operator()(heap_node& a, heap_node& b){
            return a.chunk_id>b.chunk_id;
        }
    };
    typedef ts_priority_queue<heap_node, std::vector<heap_node>, heap_node_compare> min_heap_type;

    template<class sym_type, bool p_round>
    off_t create_meta_sym(std::vector<uint32_t>& mt_perm, uint64_t pf_seed,
                         const lz_like_map::phrase_list_t & phrase_set,
                         sym_type * text, off_t txt_size,
                         std::vector<uint64_t>& prev_fps,
                         partial_gram<uint8_t>& p_gram){

        size_t source, end, tot_symbols=0;
        std::vector<std::pair<uint32_t, uint64_t>> perm(phrase_set.size());
        std::vector<uint64_t> fp_sequence;

        for(size_t i=0;i<phrase_set.size();i++) {
            perm[i].first = i;
            source = phrase_set[i].source/sizeof(sym_type);
            end = source + (phrase_set[i].len/sizeof(sym_type));
            for(size_t j=source;j<end;j++){
                assert(text[j]>0);
                assert(text[j]<prev_fps.size());
                fp_sequence.push_back(prev_fps[text[j]]);
            }
            perm[i].second = XXH64(fp_sequence.data(), fp_sequence.size()*sizeof(uint64_t), pf_seed);
            tot_symbols+=fp_sequence.size();
            fp_sequence.clear();
        }

        //sort the phrases according their hash values
        std::sort(perm.begin(), perm.end(), [&](auto const& a, auto const &b) -> bool{
            return a.second<b.second;
        });

        size_t prev_hash = perm[0].second;
        off_t prev_pos=0, tot_phrases=off_t(perm.size());
        off_t n_cols =0;
        for(off_t i=1; i<tot_phrases; i++){
            if(prev_hash!=perm[i].second){
                if((i-prev_pos)>1){

                    n_cols +=(i-prev_pos);
                    std::cout<<"Warning: we have "<<(i-prev_pos)<<" colliding phrases"<<std::endl;
                    //TODO testing
                    for(off_t k=prev_pos;k<i;k++) {
                        std::cout<<"sorted pos: "<<k<<" orig pos: "<<perm[k].first<<" len "<<phrase_set[perm[k].first].len/sizeof(sym_type)<< " source " << phrase_set[perm[k].first].source/sizeof(sym_type)<< " -> ";
                        off_t s = phrase_set[perm[k].first].source/sizeof(sym_type);
                        off_t len = phrase_set[perm[k].first].len/sizeof(sym_type);
                        assert(s<txt_size);

                        fp_sequence.clear();
                        for (off_t u = s, l=0; l < len; u++, l++) {
                            fp_sequence.push_back(prev_fps[text[u]]);
                        }
                        uint64_t test_hash = XXH64(fp_sequence.data(), fp_sequence.size()*sizeof(uint64_t), pf_seed);
                        assert(perm[k].second==test_hash);
                    }
                    //

                    //sort the range [prev_pos..i-1]
                    std::sort(perm.begin()+prev_pos, perm.begin()+i, [&](auto a, auto b) -> bool{

                        lz_like_map::phrase_t phrase_a = phrase_set[a.first];
                        auto data_a = (sym_type *)&text[phrase_a.source/sizeof(sym_type)];
                        off_t len_a = phrase_a.len/sizeof(sym_type);

                        lz_like_map::phrase_t phrase_b = phrase_set[b.first];
                        auto data_b = (sym_type *)&text[phrase_b.source/sizeof(sym_type)];
                        off_t len_b = phrase_b.len/sizeof(sym_type);

                        size_t len = std::min(len_a, len_b);
                        size_t j=0;

                        while(j<len && data_a[j]==data_b[j]) j++;

                        if constexpr (p_round){
                            return prev_fps[data_a[j]] < prev_fps[data_b[j]];
                        }else{
                            return data_a[j]<data_b[j];
                        }
                    });
                }
                prev_pos = i;
                prev_hash = perm[i].second;
            }
        }

        prev_fps.resize(perm.size()+1);
        mt_perm.resize(perm.size()+1);
        prev_fps[0] = 0;
        mt_perm[0] = 0;
        for(size_t i=0, mt_sym=1;i<perm.size();i++, mt_sym++){
            size_t perm_mt_sym =  perm[i].first+1;
            assert(perm_mt_sym<mt_perm.size());
            mt_perm[perm_mt_sym] = mt_sym;
            prev_fps[mt_sym] = perm[i].second;
        }
        prev_fps.shrink_to_fit();
        p_gram.template append_new_lvl<sym_type>(text, phrase_set, tot_symbols, perm);
        return n_cols;
    }

    template<class sym_type, bool p_round>
    off_t parsing_round(sym_type* text, off_t txt_size, text_chunk::size_type* parse,
                        off_t& n_strings, size_t sep_sym, uint64_t fp_seed,
                        std::vector<uint64_t>& prev_fps, partial_gram<uint8_t>& p_gram){

        size_t prev_sym, curr_sym, next_sym, dummy_sym=std::numeric_limits<text_chunk::size_type>::max();
        uint64_t prev_hash=0, curr_hash=0, next_hash=0;

        text_chunk::size_type mt_sym;
        off_t txt_pos = 0, parse_size = 0, phrase_len, lb, rb;
        lz_like_map map((uint8_t *)text);

        bool inserted;
        n_strings = 0;
        off_t sym_bytes = sizeof(sym_type);

        while(txt_pos<txt_size){

            lb = txt_pos;
            prev_sym = text[txt_pos++];

            curr_sym = text[txt_pos++];
            while(txt_pos<txt_size && curr_sym==prev_sym) curr_sym = text[txt_pos++];
            rb = txt_pos-1;

            if(curr_sym==sep_sym){
                phrase_len = rb-lb;
                mt_sym = map.insert(lb*sym_bytes, phrase_len*sym_bytes, inserted);
                assert(text[rb]==sep_sym);

                if constexpr (p_round){
                    parse[parse_size++] = mt_sym+1;
                    parse[parse_size++] = 0;
                }else{//replace the phase with the metasymbol in place
                    assert(text[lb]!=dummy_sym);
                    if(!inserted){
                        text[lb] = mt_sym+1;//store the metabmol in the first phrase position
                        memset(&text[lb+1], dummy_sym, sym_bytes*(phrase_len-1));//pad the rest of the phrase with dummy symbols
                    }
                    text[rb] = 0;
                    parse_size+=2;//+1 for the separator symbol
                }
                n_strings++;
                continue;
            }

            next_sym = text[txt_pos++];
            while(txt_pos<txt_size && next_sym==curr_sym) next_sym = text[txt_pos++];

            if constexpr (p_round) {
                prev_hash = prev_fps[prev_sym];
                curr_hash = prev_fps[curr_sym];
            }

            bool local_minimum;

            while(txt_pos<txt_size && next_sym!=sep_sym){

                if constexpr (p_round){
                    next_hash = prev_fps[next_sym];
                    local_minimum = prev_hash>curr_hash && curr_hash<next_hash;
                }else{
                    local_minimum = prev_sym>curr_sym && curr_sym<next_sym;
                }

                if(local_minimum){
                    phrase_len = rb-lb;
                    mt_sym = map.insert(lb*sym_bytes, phrase_len*sym_bytes, inserted);
                    if constexpr (p_round){
                        parse[parse_size++] = mt_sym+1;
                    }else{//store the metasymbol in place
                        assert(text[lb]!=dummy_sym);
                        if(!inserted){
                            text[lb] = mt_sym+1;
                            memset(&text[lb+1], dummy_sym, sym_bytes*(phrase_len-1));
                        }
                        parse_size++;
                    }
                    lb = rb;
                }

                rb = txt_pos-1;
                if constexpr (p_round){
                    prev_hash = curr_hash;
                    curr_hash = next_hash;
                } else {
                    prev_sym = curr_sym;
                }

                curr_sym = next_sym;
                next_sym = text[txt_pos++];
                while(txt_pos<txt_size && next_sym==curr_sym) next_sym = text[txt_pos++];
            }

            phrase_len = txt_pos-1-lb;
            assert(text[lb+phrase_len]==sep_sym);
            mt_sym = map.insert(lb*sym_bytes, phrase_len*sym_bytes, inserted);

            if constexpr (p_round){
                parse[parse_size++] = mt_sym+1;
                parse[parse_size++] = 0;
            } else{//replace the phase with the metasymbol in place
                assert(text[lb]!=dummy_sym);
                if(!inserted){
                    text[lb] = mt_sym+1;//store the metabmol in the first phrase position
                    memset(&text[lb+1], dummy_sym, sym_bytes*(phrase_len-1));//pad the rest of the phrase with dummy symbols
                }
                text[lb+phrase_len] = 0;
                parse_size+=2;//+1 for the separator symbol
            }
            n_strings++;
        }

        map.shrink_to_fit();
        map.destroy_table();

        assert(map.phrase_set.size()<dummy_sym);

        std::vector<uint32_t> perm;
        create_meta_sym<sym_type, p_round>(perm, fp_seed, map.phrase_set, text, txt_size, prev_fps, p_gram);

        if constexpr (p_round){
            for(off_t i=0;i<parse_size;i++){
                assert(parse[i]<perm.size());
                parse[i] = perm[parse[i]];
            }
        }else{
            // replace the parts of the text that
            // we use as sources for the phases
            mt_sym =0;
            while(mt_sym<map.phrase_set.size()) {
                lb = map.phrase_set[mt_sym].source/sizeof(sym_type);
                phrase_len = map.phrase_set[mt_sym].len/sizeof(sym_type);
                text[lb] = mt_sym+1;//store the metasymol in the first phrase position
                memset(&text[lb+1], dummy_sym, sym_bytes*(phrase_len-1));//pad the rest of the phrase with dummy symbols
                mt_sym++;
            }

            off_t i=0, k=0;
            while(i<txt_size){
                assert(text[i]<perm.size());
                text[k++] = perm[text[i]];
                i++;
                while(i<txt_size && text[i]==dummy_sym) i++;
            }
            assert(k==parse_size);
        }
        return parse_size;
    }

    template<class sym_type>
    void compress_text_chunk(text_chunk& chunk, std::vector<uint64_t>& fp_seeds){

        size_t alpha_size = std::numeric_limits<sym_type>::max()+1;
        std::vector<uint64_t> prev_fps(alpha_size);

        for(size_t i=0;i<alpha_size;i++){
            prev_fps[i] = XXH64(&i, sizeof(sym_type), fp_seeds[0]);
        }

        off_t n_strings=0;
        size_t p_round=0;
        size_t sep_sym = chunk.sep_sym;

        chunk.p_gram.max_tsym = std::numeric_limits<sym_type>::max();
        chunk.p_gram.sep_tsym = chunk.sep_sym;
        chunk.p_gram.text_size = chunk.text_bytes/sizeof(sym_type);

        off_t parse_size = parsing_round<sym_type, true>(chunk.text, chunk.text_bytes/sizeof(sym_type), chunk.parse,
                                                         n_strings, sep_sym, fp_seeds[p_round+1], prev_fps, chunk.p_gram);
        sep_sym = 0;
        off_t size_limit = n_strings*2;
        auto * new_parse = (text_chunk::size_type *)chunk.text;
        p_round++;

        while(parse_size!=size_limit){
            assert(parse_size>=size_limit);
            parse_size = parsing_round<text_chunk::size_type, false>(chunk.parse, parse_size, new_parse, n_strings,
                                                                     sep_sym, fp_seeds[p_round+1], prev_fps, chunk.p_gram);
            //std::swap(chunk.parse, new_parse);
            p_round++;
        }

        chunk.p_gram.add_compressed_string(chunk.parse, parse_size);
    }

    template<class sym_type>
    void build_grammars(parsing_opts& p_opts, std::string& i_file, tmp_workspace& tmp_ws){
        ts_queue<size_t> chunks_to_read;
        min_heap_type chunks_to_merge;

        std::atomic<size_t> parser_finished{0};
        std::vector<text_chunk> text_chunks(p_opts.n_chunks);
        struct stat st{};
        if (stat(i_file.c_str(), &st) != 0) return;

        auto text_read_worker = [&]() -> void {

            int fd_r = open(i_file.c_str(), O_RDONLY);

#ifdef __linux__
            off_t r_page_cache_bytes = 0, w_page_cache_bytes = 0;
            posix_fadvise(fd_r, 0, st.st_size, POSIX_FADV_SEQUENTIAL);
#endif

            off_t rem_bytes = st.st_size, r_acc_bytes = 0;
            size_t chunk_id = 0;

            auto tmp_ck_size = off_t(INT_CEIL(p_opts.chunk_size, sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type));

            while (chunk_id < p_opts.n_chunks && rem_bytes > 0) {

                tmp_ck_size = std::min(tmp_ck_size, rem_bytes);

                text_chunks[chunk_id].text_bytes = tmp_ck_size;
                text_chunks[chunk_id].sep_sym = (text_chunk::size_type) p_opts.sep_sym;

                //the parse size is (text_len/2)*(sizeof(size_type)/sizeof(sym_type)),
                // where ``text_len'' is the number of input symbols that fits the buffer
                off_t parse_bytes = INT_CEIL((tmp_ck_size/sizeof(sym_type)), 2)*(sizeof(text_chunk::size_type)/sizeof(sym_type));

                text_chunks[chunk_id].buffer_bytes = off_t(tmp_ck_size + parse_bytes);
                //text_chunks[chunk_id].buffer = (text_chunk::size_type *) malloc(text_chunks[chunk_id].buffer_bytes);
                text_chunks[chunk_id].buffer = (text_chunk::size_type *) mmap_allocate(text_chunks[chunk_id].buffer_bytes);
                text_chunks[chunk_id].id = chunk_id;

                read_chunk_from_file<sym_type>(fd_r, rem_bytes, r_acc_bytes, text_chunks[chunk_id]);

                //next aligned position within the buffer
                size_t parse_start =  INT_CEIL(text_chunks[chunk_id].text_bytes, sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type);
                text_chunks[chunk_id].parse = (text_chunk::size_type *) &text_chunks[chunk_id].text[parse_start/sizeof(sym_type)];

                //this is enough for 4GB strings
                text_chunks[chunk_id].p_gram.rules.resize(32);

                chunks_to_read.push(chunk_id);
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

            //data to write
            std::string concat_grams = tmp_ws.get_file("concatenated_grams");
            heap_node h_node{};
            size_t next_chunk=0;
            int fd_w = open(concat_grams.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
#ifdef __linux__
            off_t w_acc_bytes=0;
#endif
            //

            size_t proc_syms=0;
            size_t acc_strings=0, n_strings;
            while (rem_bytes > 0) {

                while(true) {
                    chunks_to_merge.pop(h_node);
                    if(h_node.chunk_id==next_chunk){

                        n_strings = text_chunks[h_node.buff_idx].p_gram.tot_strings();
                        text_chunks[h_node.buff_idx].p_gram.first_string = acc_strings;
                        text_chunks[h_node.buff_idx].p_gram.last_string = acc_strings+n_strings-1;
                        acc_strings+=n_strings;

                        size_t g_bytes = text_chunks[h_node.buff_idx].p_gram.serialize_to_fd(fd_w);
                        std::cout<<report_space(text_chunks[h_node.buff_idx].text_bytes)<<" compressed to "<<report_space((off_t)g_bytes)<<std::endl;
#ifdef __linux__
                        w_page_cache_bytes += g_bytes;
                        if(w_page_cache_bytes>p_opts.page_cache_limit){
                            std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                            posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
                            w_page_cache_bytes=0;
                        }
#endif
                        text_chunks[h_node.buff_idx].p_gram.reset_grammar();
                        next_chunk++;
                        break;
                    }else{
                        chunks_to_merge.push(h_node);
                    }
                }

                proc_syms+=text_chunks[h_node.buff_idx].text_bytes;

                //std::cout<<"\n  Processed input "<<report_space((off_t)proc_syms)<<"    "<<std::flush;
                std::cout<<"  Processed input "<<report_space((off_t)proc_syms)<<"    "<<rem_bytes<<std::endl;

                text_chunks[h_node.buff_idx].text_bytes = tmp_ck_size;
                text_chunks[h_node.buff_idx].id = chunk_id++;

                read_chunk_from_file<sym_type>(fd_r, rem_bytes, r_acc_bytes, text_chunks[h_node.buff_idx]);

                //next aligned position
                size_t parse_start =  INT_CEIL(text_chunks[h_node.buff_idx].text_bytes, sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type);
                text_chunks[h_node.buff_idx].parse = (text_chunk::size_type *) &text_chunks[h_node.buff_idx].text[parse_start/sizeof(sym_type)];

                chunks_to_read.push(h_node.buff_idx);
#ifdef __linux__
                r_page_cache_bytes+=text_chunks[h_node.buff_idx].e_bytes;
                if(r_page_cache_bytes>p_opts.page_cache_limit){
                    std::cout<<"removing from page cache "<<r_page_cache_bytes<<" "<<r_acc_bytes<<std::endl;
                    posix_fadvise(fd_r, r_acc_bytes-r_page_cache_bytes, r_page_cache_bytes, POSIX_FADV_DONTNEED);
                    r_page_cache_bytes=0;
                }
#endif
            }

            //wait for the stack to be empty and close the input file
            while (!chunks_to_read.empty());
            chunks_to_read.done();
#ifdef __linux__
            posix_fadvise(fd_r, 0, st.st_size, POSIX_FADV_DONTNEED);
#endif
            close(fd_r);

            //wait for all the parsers to finish
            while(parser_finished.load(std::memory_order_acquire)!=p_opts.n_threads);

            //store the remaining grammars in the temporary file
            while(next_chunk<chunk_id) {
                while(true){
                    chunks_to_merge.pop(h_node);

                    if(h_node.chunk_id==next_chunk) {

                        n_strings = text_chunks[h_node.buff_idx].p_gram.tot_strings();
                        text_chunks[h_node.buff_idx].p_gram.first_string = acc_strings;
                        text_chunks[h_node.buff_idx].p_gram.last_string = acc_strings+n_strings-1;
                        acc_strings+=n_strings;

                        size_t g_bytes = text_chunks[h_node.buff_idx].p_gram.serialize_to_fd(fd_w);
                        std::cout<<report_space(text_chunks[h_node.buff_idx].text_bytes)<<" compressed to "<<report_space((off_t)g_bytes)<<std::endl;
#ifdef __linux__
                        w_page_cache_bytes += g_bytes;
                        if(w_page_cache_bytes>p_opts.page_cache_limit){
                            std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                            posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
                            w_page_cache_bytes=0;
                        }
#endif
                        next_chunk++;
                        break;
                    } else {
                        chunks_to_merge.push(h_node);
                    }
                }
                proc_syms+=text_chunks[h_node.buff_idx].text_bytes;
                //std::cout<<"\n  Processed input "<<report_space((off_t)proc_syms)<<"     "<<std::flush;
                std::cout<<"  Processed input "<<report_space((off_t)proc_syms)<<"     "<<std::endl;
            }
            chunks_to_merge.done();

#ifdef __linux__
            if(w_page_cache_bytes>0){
                std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
            }
#endif
            close(fd_w);
        };

        auto parser_worker = [&]() {
            size_t buff_id;
            bool res;

            while (true) {
                res = chunks_to_read.pop(buff_id);
                assert(text_chunks[buff_id].text_bytes > 0);
                if (!res){
                    parser_finished.fetch_add(1, std::memory_order_acq_rel);
                    break;
                }

                compress_text_chunk<sym_type>(text_chunks[buff_id], p_opts.p_seeds);
                memset(text_chunks[buff_id].buffer, 0, text_chunks[buff_id].buffer_bytes);

                chunks_to_merge.push({text_chunks[buff_id].id, buff_id});
            }
        };

        std::vector<std::thread> threads;
        threads.emplace_back(text_read_worker);

        for (size_t i = 0; i < p_opts.n_threads; i++) {
            threads.emplace_back(parser_worker);
        }

        for (auto &thread: threads) {
            thread.join();
        }

        for(auto &chunk : text_chunks){
            mmap_deallocate(chunk.buffer, chunk.buffer_bytes);
            chunk.buffer=nullptr;
        }
    }

    template<class sym_type>
    void merge_grammars(parsing_opts& p_opts, std::string& o_file, tmp_workspace& tmp_ws){

        using p_gram_type = partial_gram<sym_type, true>;

        size_t n_threads = p_opts.n_threads;

        std::vector<std::pair<p_gram_type,
                    std::vector<std::pair<size_t, size_t>>>
                    > initial_grams(n_threads);

        ts_queue<size_t> gram_to_merge_queue;
        ts_queue<std::pair<size_t, off_t>> av_buff_queue;

        std::string p_grams_file = tmp_ws.get_file("concatenated_grams");
        int fd_r = open(p_grams_file.c_str(), O_RDONLY);
        size_t tot_bytes = file_size(p_grams_file);
        size_t rem_bytes =  tot_bytes;
        size_t read_bytes;
        size_t i=0;

        while(i<n_threads && rem_bytes>0){
            read_bytes = initial_grams[i].first.load_from_fd(fd_r);
            //store the range of strings this partial gram covers
            initial_grams[i].second.push_back({initial_grams[i].first.first_string,
                                               initial_grams[i].first.last_string});
            rem_bytes-=read_bytes;
            av_buff_queue.push({i, 0});
            i++;
        }
        n_threads = i;
        std::vector<p_gram_type> grams_to_merge(n_threads);

        auto gram_read_worker = [&](){

            if(rem_bytes==0){
                gram_to_merge_queue.done();
            }else{
                std::pair<size_t , off_t> proc_input;
                while(rem_bytes > 0){
                    std::cout<<"Processed data: "<<double(rem_bytes)/double(tot_bytes)*100<<std::endl;
                    av_buff_queue.pop(proc_input);
                    read_bytes = grams_to_merge[proc_input.first].load_from_fd(fd_r);
                    gram_to_merge_queue.push(proc_input.first);
                    rem_bytes-=read_bytes;
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

                for (size_t i = 1; i < n_threads; i++) {
                    merge_two_grammars<p_gram_type>(initial_grams[0].first, initial_grams[i].first, p_opts.p_seeds);
                    //store the range of strings this partial gram covers
                    initial_grams[0].second.insert(initial_grams[0].second.end(), initial_grams[i].second.begin(), initial_grams[i].second.end());
                    initial_grams[i].first.destroy_gram();
                }
                initial_grams[0].second.shrink_to_fit();
#ifdef __linux__
                malloc_trim(0);
#endif
                initial_grams[0].first.reorder_strings(initial_grams[0].second);
            }

            std::string mg_p_gram_file = tmp_ws.get_file("merged_partial_grams");
            store_to_file(mg_p_gram_file, initial_grams[0].first);

            initial_grams[0].first.destroy_gram();

#ifdef __linux__
            malloc_trim(0);
#endif

            lc_gram_t final_grammar(mg_p_gram_file, p_opts.p_seeds);
            store_to_file(o_file, final_grammar);
        };

        auto gram_merge_worker = [&](size_t idx){

            size_t buff_id;
            bool res;

            while (true) {
                res = gram_to_merge_queue.pop(buff_id);
                if (!res) break;
                //report_mem_peak();
                //std::cout<<grams_to_merge[buff_id].first_string<<" --- "<<grams_to_merge[buff_id].last_string<<std::endl;
                merge_two_grammars<p_gram_type>(initial_grams[idx].first, grams_to_merge[buff_id], p_opts.p_seeds);

                if(n_threads>1){
                    initial_grams[idx].second.push_back({grams_to_merge[buff_id].first_string,
                                                         grams_to_merge[buff_id].last_string});
                }
                av_buff_queue.push({buff_id, initial_grams[idx].first.text_size});
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

    template<class sym_type>
    void lc_parsing_algo(std::string& i_file, std::string& p_file, std::string& o_file,
                         tmp_workspace& tmp_ws, size_t n_threads,
                         size_t n_chunks, size_t chunk_size) {

        parsing_opts p_opts;
        p_opts.n_threads = n_threads;
        p_opts.n_chunks = n_chunks==0? n_threads*2 : n_chunks;
        p_opts.chunk_size = chunk_size==0 ? off_t(ceil(0.025 * double(file_size(i_file)))) : (off_t)chunk_size;
        p_opts.chunk_size = std::min<off_t>(p_opts.chunk_size, std::numeric_limits<uint32_t>::max());//the chunks cannot exceed the 4GB by design
        //p_opts.chunk_size = std::min<off_t>(1020*1024*100, file_size(i_file));
        p_opts.page_cache_limit = 1024*1024*1024;
        p_opts.sep_sym = 10;

        std::string seed_source;
        if(!p_file.empty()) {
            load_pl_vector(p_file, p_opts.p_seeds);
            assert(!p_opts.p_seeds.empty());
            seed_source = p_file;
        } else {
            std::random_device rd;  // Will be used to obtain a seed for the random number engine
            //std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
            std::mt19937 gen(5489U); // Standard mersenne_twister_engine seeded with a fixed value
            std::uniform_int_distribution<uint64_t> distrib(1, std::numeric_limits<uint64_t>::max());
            p_opts.p_seeds.resize(32);
            for(size_t i=0;i<32;i++){
                p_opts.p_seeds[i] = distrib(gen);
            }
            //TODO create the seeds
            //store_pl_vector("fixed_hash_functions", p_opts.p_functions);
            seed_source = "randomly chosen";
        }

        std::cout<<"  Settings"<<std::endl;
        std::cout<<"    Parsing seeds             : "<<seed_source<<std::endl;
        std::cout<<"    Parsing threads           : "<<p_opts.n_threads<<std::endl;
        std::cout<<"    Active text chunks in RAM : "<<p_opts.n_chunks<<std::endl;
        std::cout<<"    Size of each chunk        : "<<report_space(p_opts.chunk_size)<<std::endl;
        std::cout<<"    Chunks' approx. mem usage : "<<report_space(off_t(p_opts.chunk_size*p_opts.n_chunks*3))<<"\n"<<std::endl;

        build_grammars<sym_type>(p_opts, i_file, tmp_ws);
        merge_grammars<sym_type>(p_opts, o_file, tmp_ws);
    }
}

#endif //LCG_LZ_LIKE_LC_PARSING_H
