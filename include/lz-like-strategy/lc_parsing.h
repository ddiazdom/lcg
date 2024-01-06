//
// Created by Diaz, Diego on 17.9.2023.
//

#ifndef LCG_LZ_LIKE_LC_PARSING_H
#define LCG_LZ_LIKE_LC_PARSING_H

#include <random>
#include <fcntl.h>

#include "cds/macros.h"
#include "cds/ts_queue.h"
#include "cds/ts_priority_queue.h"
#include "cds/utils.h"
#include "cds/cdt_common.hpp"

#include "../external/xxHash-dev/xxhash.h"

#include "text_handler.h"
#include "lz_map.h"

namespace lzstrat {

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
                         const lz_map::phrase_list_t & phrase_set,
                         sym_type * text, off_t txt_size,
                         std::vector<uint64_t>& hash_values_vector,
                         partial_gram<uint8_t>& p_gram){

        size_t source, end, tot_symbols=0;
        std::vector<std::pair<uint32_t, uint64_t>> perm(phrase_set.size());
        std::vector<uint64_t> tmp_phrase;

        for(size_t i=0;i<phrase_set.size();i++) {
            perm[i].first = i;
            source = phrase_set[i].source/sizeof(sym_type);
            end = source + (phrase_set[i].len/sizeof(sym_type));
            for(size_t j=source;j<end;j++){
                assert(text[j]>0);
                if constexpr (p_round){
                    uint64_t fp_sym = XXH64(&text[j], sizeof(sym_type), pf_seed);
                    tmp_phrase.push_back(fp_sym);
                }else{
                    assert((text[j]-1)<hash_values_vector.size());
                    tmp_phrase.push_back(hash_values_vector[text[j]-1]);
                }
            }
            perm[i].second = XXH64(tmp_phrase.data(), tmp_phrase.size()*sizeof(uint64_t), pf_seed);
            tot_symbols+=tmp_phrase.size();
            tmp_phrase.clear();
        }

        //sort the phrases according their hash values
        std::sort(perm.begin(), perm.end(), [&](auto a, auto b) -> bool{
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

                        tmp_phrase.clear();
                        for (off_t u = s, l=0; l < len; u++, l++) {
                            uint64_t fp;
                            if constexpr (p_round){
                                fp = XXH64(&text[u], sizeof(sym_type), pf_seed);
                            }else{
                                fp = hash_values_vector[text[u]-1];
                            }
                            tmp_phrase.push_back(fp);
                        }
                        uint64_t test_hash = XXH64(tmp_phrase.data(), tmp_phrase.size()*sizeof(uint64_t), pf_seed);
                        assert(perm[k].second==test_hash);
                    }
                    //

                    //sort the range [prev_pos..i-1]
                    std::sort(perm.begin()+prev_pos, perm.begin()+i, [&](auto a, auto b) -> bool{

                        lz_map::phrase_t phrase_a = phrase_set[a.first];
                        auto data_a = (sym_type *)&text[phrase_a.source/sizeof(sym_type)];
                        off_t len_a = phrase_a.len/sizeof(sym_type);

                        lz_map::phrase_t phrase_b = phrase_set[b.first];
                        auto data_b = (sym_type *)&text[phrase_b.source/sizeof(sym_type)];
                        off_t len_b = phrase_b.len/sizeof(sym_type);

                        size_t len = std::min(len_a, len_b);
                        size_t j=0;

                        while(j<len && data_a[j]==data_b[j]) j++;

                        if constexpr (p_round){
                            return XXH64(&data_a[j], sizeof(sym_type), pf_seed) < XXH64(&data_b[j], sizeof(sym_type), pf_seed);
                        }else{
                            return data_a[j]<data_b[j];
                        }
                    });
                }
                prev_pos = i;
                prev_hash = perm[i].second;
            }
        }

        hash_values_vector.resize(perm.size());
        for(size_t i=0;i<hash_values_vector.size();i++){
            hash_values_vector[i] = perm[i].second;
        }
        hash_values_vector.shrink_to_fit();

        mt_perm.resize(perm.size());
        for(size_t i=0, rank=0;i<mt_perm.size();i++, rank++){
            assert(perm[i].first<mt_perm.size());
            mt_perm[perm[i].first] = rank;
        }
        p_gram.template append_new_lvl<sym_type>(text, phrase_set, tot_symbols, perm);
        return n_cols;
    }

    template<class sym_type, bool p_round>
    off_t parsing_round(sym_type* text, off_t txt_size, text_chunk::size_type* parse,
                        off_t& n_strings, size_t sep_sym, uint64_t fp_seed,
                        std::vector<uint64_t>& prev_hash_values, partial_gram<uint8_t>& p_gram){

        size_t prev_sym, curr_sym, next_sym;
        text_chunk::size_type mt_sym;
        off_t txt_pos = 0, parse_pos = 0, phrase_len, lb, rb;
        lz_map map((uint8_t *)text);

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
                n_strings++;
                mt_sym = map.insert((txt_pos-2)*sym_bytes, sym_bytes, inserted);
                //assert(text[txt_pos-2]!=0);

                parse[parse_pos++] = mt_sym+1;
                parse[parse_pos++] = 0;
                continue;
            }

            next_sym = text[txt_pos++];
            while(txt_pos<txt_size && next_sym==curr_sym) next_sym = text[txt_pos++];

            while(txt_pos<txt_size && next_sym!=sep_sym){

                if(prev_sym>curr_sym && curr_sym<next_sym){
                    phrase_len = rb-lb;
                    //assert(text[lb+phrase_len-1]>0);
                    //assert(text[lb]>0);
                    mt_sym = map.insert(lb*sym_bytes, phrase_len*sym_bytes, inserted);
                    parse[parse_pos++] = mt_sym+1;
                    lb = rb;
                }

                rb = txt_pos-1;
                prev_sym = curr_sym;
                curr_sym = next_sym;
                next_sym = text[txt_pos++];
                while(txt_pos<txt_size && next_sym==curr_sym) next_sym = text[txt_pos++];
            }

            phrase_len = txt_pos-1-lb;
            assert(text[lb+phrase_len]==sep_sym);
            //assert(text[lb+phrase_len-1]>0);
            //assert(text[lb]>0);
            //std::cout<<lb<<" "<<txt_pos<<" "<<(lb+phrase_len)<<" "<<txt_size<<std::endl;
            //for(off_t k=lb, u=0;u<phrase_len;u++, k++){
            //    std::cout<<text[k]<<" ";
            //}
            //std::cout<<" -> "<<int(text[lb+phrase_len])<<std::endl;
            //assert((lb+phrase_len)<=txt_size);

            mt_sym = map.insert(lb*sym_bytes, phrase_len*sym_bytes, inserted);
            parse[parse_pos++] = mt_sym+1;
            parse[parse_pos++] = 0;
            n_strings++;
        }

        map.shrink_to_fit();
        map.destroy_table();

        std::vector<uint32_t> perm;
        create_meta_sym<sym_type, p_round>(perm, fp_seed, map.phrase_set, text, txt_size, prev_hash_values, p_gram);
        for(off_t i=0;i<parse_pos;i++){
            if(parse[i]==0) continue;
            assert(parse[i]<=perm.size());
            parse[i] = perm[parse[i]-1]+1;
        }
        return parse_pos;
    }

    template<class sym_type>
    void compress_text_chunk(text_chunk& chunk, std::vector<uint64_t>& fp_seeds){

        std::vector<uint64_t> prev_hash_values;
        off_t n_strings=0;
        size_t p_round=0;
        size_t sep_sym = chunk.sep_sym;

        off_t parse_size = parsing_round<sym_type, true>(chunk.text, chunk.text_bytes/sizeof(sym_type), chunk.parse, n_strings, sep_sym, fp_seeds[p_round++], prev_hash_values, chunk.p_gram);
        sep_sym = 0;
        off_t size_limit = n_strings*2;
        auto * new_parse = (text_chunk::size_type *)chunk.text;

        while(parse_size!=size_limit){
            assert(parse_size>=size_limit);
            parse_size = parsing_round<text_chunk::size_type, false>(chunk.parse, parse_size, new_parse, n_strings, sep_sym, fp_seeds[p_round++], prev_hash_values, chunk.p_gram);
            std::swap(chunk.parse, new_parse);
        }

        chunk.p_gram.add_compressed_string(chunk.text, parse_size);
    }

    template<class sym_type>
    void lc_parsing_algo(std::string& i_file, std::string& p_file, std::string& o_file,
                         tmp_workspace& tmp_ws, size_t n_threads,
                         size_t n_chunks, size_t chunk_size) {

        parsing_opts p_opts;
        p_opts.n_threads = n_threads;
        p_opts.n_chunks = n_chunks==0? n_threads*2 : n_chunks;
        p_opts.chunk_size = chunk_size==0 ? off_t(ceil(0.025 * double(file_size(i_file)))) : (off_t)chunk_size; //std::min<off_t>(1020*1024*100, file_size(i_file));
        p_opts.page_cache_limit = 1024*1024*1024;
        p_opts.sep_sym = 10;

        std::string seed_source;
        if(!p_file.empty()) {
            load_pl_vector(p_file, p_opts.p_seeds);
            assert(!p_opts.p_seeds.empty());
            seed_source = p_file;
        } else {
            std::random_device rd;  // Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
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


        ts_queue<size_t> chunks_to_read;
        min_heap_type chunks_to_merge;

        std::atomic<size_t> parser_finished{0};
        std::vector<text_chunk> text_chunks(p_opts.n_chunks);
        struct stat st{};
        if (stat(i_file.c_str(), &st) != 0) return;

        auto read_worker = [&]() -> void {

            int fd_r = open(i_file.c_str(), O_RDONLY);

#ifdef __linux__
            off_t r_page_cache_bytes = 0, w_page_cache_bytes = 0;
            posix_fadvise(r_fd, 0, st.st_size, POSIX_FADV_SEQUENTIAL);
#endif

            off_t rem_bytes = st.st_size, r_acc_bytes = 0;
            size_t chunk_id = 0;

            auto tmp_ck_size = off_t((p_opts.chunk_size/sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type));

            while (chunk_id < p_opts.n_chunks && rem_bytes > 0) {

                tmp_ck_size = std::min(tmp_ck_size, rem_bytes);

                text_chunks[chunk_id].text_bytes = tmp_ck_size;
                text_chunks[chunk_id].sep_sym = (text_chunk::size_type) p_opts.sep_sym;

                //the parse size is (text_len/2)*(sizeof(size_type)/sizeof(sym_type)),
                // where ``text_len'' is the number of input symbols that fits the buffer
                off_t parse_bytes = INT_CEIL((tmp_ck_size/sizeof(sym_type)), 2)*(sizeof(text_chunk::size_type)/sizeof(sym_type));

                text_chunks[chunk_id].buffer_bytes = off_t(tmp_ck_size + parse_bytes);
                text_chunks[chunk_id].buffer = (text_chunk::size_type *) malloc(text_chunks[chunk_id].buffer_bytes);
                text_chunks[chunk_id].id = chunk_id;

                read_chunk_from_file<sym_type>(fd_r, rem_bytes, r_acc_bytes, text_chunks[chunk_id]);

                //next aligned position within the buffer
                size_t parse_start =  INT_CEIL(text_chunks[chunk_id].text_bytes, sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type);
                text_chunks[chunk_id].parse = (text_chunk::size_type *) &text_chunks[chunk_id].text[parse_start/sizeof(sym_type)];

                chunks_to_read.push(chunk_id);
#ifdef __linux__
                r_page_cache_bytes+=text_chunks[chunk_id].eff_buff_bytes();
                if(r_page_cache_bytes>p_opts.page_cache_limit){
                    std::cout<<"removing from page cache "<<r_page_cache_bytes<<" "<<r_acc_bytes<<std::endl;
                    posix_fadvise(fd_r, acc_bytes-r_page_cache_bytes, r_page_cache_bytes, POSIX_FADV_DONTNEED);
                    r_page_cache_bytes=0;
                }
#endif
                chunk_id++;
            }

            //data to write
            std::string concat_grams = tmp_ws.get_file("concatenated_grams");
            heap_node node{};
            size_t next_chunk=0;
            int fd_w = open(concat_grams.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
            off_t w_acc_bytes=0;
            //

            size_t proc_syms=0;
            while (rem_bytes > 0) {

                while(true) {
                    chunks_to_merge.pop(node);
                    if(node.chunk_id==next_chunk){
                        size_t g_bytes = text_chunks[node.buff_idx].p_gram.serialize_to_fd(fd_w);
                        std::cout<<report_space(text_chunks[node.buff_idx].text_bytes)<<" compressed to "<<report_space((off_t)g_bytes)<<std::endl;
#ifdef __linux__
                        w_page_cache_bytes += text_chunks[node.buff_idx].gram_bytes();
                        if(w_page_cache_bytes>opts.page_cache_limit){
                            std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                            posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
                            w_page_cache_bytes=0;
                        }
#endif
                        text_chunks[node.buff_idx].p_gram.reset_grammar();
                        next_chunk++;
                        break;
                    }else{
                        chunks_to_merge.push(node);
                    }
                }

                proc_syms+=text_chunks[node.buff_idx].text_bytes;

                std::cout<<"\n  Processed input "<<report_space((off_t)proc_syms)<<"    "<<std::flush;

                text_chunks[node.buff_idx].text_bytes = tmp_ck_size;
                text_chunks[node.buff_idx].id = chunk_id++;

                read_chunk_from_file<sym_type>(fd_r, rem_bytes, r_acc_bytes, text_chunks[node.buff_idx]);

                //next aligned position
                size_t parse_start =  INT_CEIL(text_chunks[node.buff_idx].text_bytes, sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type);
                text_chunks[node.buff_idx].parse = (text_chunk::size_type *) &text_chunks[node.buff_idx].text[parse_start/sizeof(sym_type)];

                chunks_to_read.push(node.buff_idx);
#ifdef __linux__
                r_page_cache_bytes+=text_chunks[node.buff_idx].eff_buff_bytes();
                if(r_page_cache_bytes>p_opts.page_cache_limit){
                    std::cout<<"removing from page cache "<<r_page_cache_bytes<<" "<<r_acc_bytes<<std::endl;
                    posix_fadvise(fd_r, acc_bytes-r_page_cache_bytes, r_page_cache_bytes, POSIX_FADV_DONTNEED);
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
                    chunks_to_merge.pop(node);

                    if(node.chunk_id==next_chunk) {
                        //TODO store grammar to disk
                        size_t g_bytes = text_chunks[node.buff_idx].p_gram.serialize_to_fd(fd_w);
                        std::cout<<report_space(text_chunks[node.buff_idx].text_bytes)<<" compressed to "<<report_space((off_t)g_bytes)<<std::endl;
#ifdef __linux__
                        w_page_cache_bytes += text_chunks[node.buff_idx].gram_bytes();
                        if(w_page_cache_bytes>opts.page_cache_limit){
                            std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                            posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
                            w_page_cache_bytes=0;
                        }
#endif
                        next_chunk++;
                        break;
                    } else {
                        chunks_to_merge.push(node);
                    }
                }
                proc_syms+=text_chunks[node.buff_idx].text_bytes;
                std::cout<<"\n  Processed input "<<report_space((off_t)proc_syms)<<"     "<<std::flush;
            }
            chunks_to_merge.done();
            std::cout<<""<<std::endl;

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
                //chunks_to_reuse.push(buff_id);
                //chunks_to_merge.push(text_chunks[buff_id].id);
            }
        };

        std::vector<std::thread> threads;
        threads.emplace_back(read_worker);

        for (size_t i = 0; i < p_opts.n_threads; i++) {
            threads.emplace_back(parser_worker);
        }

        for (auto &thread: threads) {
            thread.join();
        }
    }
}

#endif //LCG_LZ_LIKE_LC_PARSING_H
