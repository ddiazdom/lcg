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
#include "cds/vbyte_encoding.h"

namespace lz_like_strat {

    struct parsing_opts{
        size_t n_threads{};
        unsigned long n_chunks{};
        off_t chunk_size{};
        off_t page_cache_limit{};
        size_t sep_sym{};
        std::vector<uint64_t> p_seeds;
    };

    /*struct heap_node{
        size_t chunk_id;
        size_t buff_idx;
    };
    struct heap_node_compare{
        bool operator()(heap_node& a, heap_node& b){
            return a.chunk_id>b.chunk_id;
        }
    };
    typedef ts_priority_queue<heap_node, std::vector<heap_node>, heap_node_compare> min_heap_type;*/

    template<class sym_type, bool p_round>
    off_t create_meta_sym(std::vector<uint32_t>& mt_perm, uint64_t pf_seed,
                         const typename lz_like_map<sym_type>::phrase_list_t& phrase_set,
                         sym_type * text, off_t txt_size,
                         std::vector<uint64_t>& prev_fps,
                         partial_gram<uint8_t>& p_gram){

        size_t source, end, tot_symbols=0;
        std::vector<std::pair<uint32_t, uint64_t>> perm(phrase_set.size());
        std::vector<uint64_t> fp_sequence;

        for(size_t i=0;i<phrase_set.size();i++) {
            perm[i].first = i;
            source = phrase_set[i].source;
            end = source + phrase_set[i].len;
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
                        std::cout<<"sorted pos: "<<k<<" orig pos: "<<perm[k].first<<" len "<<phrase_set[perm[k].first].len<< " source " << phrase_set[perm[k].first].source<< " -> ";
                        off_t s = phrase_set[perm[k].first].source;
                        off_t len = phrase_set[perm[k].first].len;
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

                        typename lz_like_map<sym_type>::phrase_t phrase_a = phrase_set[a.first];
                        sym_type * data_a = &text[phrase_a.source];
                        off_t len_a = phrase_a.len;

                        typename lz_like_map<sym_type>::phrase_t phrase_b = phrase_set[b.first];
                        sym_type * data_b = &text[phrase_b.source];
                        off_t len_b = phrase_b.len;

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

    // this method parses the text and store the parse in the text itself.
    // It only works for parsing rounds other than the first one because the length of symbol each
    // cell is the same as the length of cell where we store the metasymbols, so there is no overflow
    off_t inplace_parsing_round(uint32_t* text, off_t txt_size, off_t& n_strings, uint64_t fp_seed,
                                std::vector<uint64_t>& prev_fps, partial_gram<uint8_t>& p_gram){

        uint32_t mt_sym, sep_sym =0;
        size_t prev_sym, curr_sym, next_sym, dummy_sym=std::numeric_limits<text_chunk::size_type>::max();
        off_t txt_pos = 0, parse_size = 0, phrase_len, lb, rb;
        lz_like_map<uint32_t> map(text);

        bool inserted;
        n_strings = 0;
        off_t sym_bytes = sizeof(uint32_t);

        while(txt_pos<txt_size){

            lb = txt_pos;
            prev_sym = text[txt_pos++];

            curr_sym = text[txt_pos++];
            while(txt_pos<txt_size && curr_sym==prev_sym) curr_sym = text[txt_pos++];
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
            while(txt_pos<txt_size && next_sym==curr_sym) next_sym = text[txt_pos++];
            bool local_minimum;

            while(txt_pos<txt_size && next_sym!=sep_sym){

                local_minimum = prev_sym>curr_sym && curr_sym<next_sym;

                if(local_minimum){
                    phrase_len = rb-lb;
                    mt_sym = map.insert(lb, phrase_len, inserted);
                    if(!inserted){
                        //we can not replace the first phrase occurrence as we use it as source for the dictionary
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
                while(txt_pos<txt_size && next_sym==curr_sym) next_sym = text[txt_pos++];
            }

            phrase_len = txt_pos-1-lb;
            assert(text[lb+phrase_len]==sep_sym);
            mt_sym = map.insert(lb, phrase_len, inserted);
            if(!inserted){
                //we can not replace the first phrase occurrence as we use it as source for the dictionary
                assert(text[lb]!=dummy_sym);
                text[lb] = mt_sym+1;//store the metabmol in the first phrase position
                memset(&text[lb+1], (int)dummy_sym, sym_bytes*(phrase_len-1));//pad the rest of the phrase with dummy symbols
            }
            text[lb+phrase_len] = sep_sym;
            parse_size+=2;//+1 for the separator symbol
            n_strings++;
        }

        map.shrink_to_fit();
        map.destroy_table();

        assert(map.phrase_set.size()<dummy_sym);

        std::vector<uint32_t> perm;
        create_meta_sym<uint32_t, false>(perm, fp_seed, map.phrase_set, text, txt_size, prev_fps, p_gram);

        // create the parse in place
        map.insert_dummy_entry({uint32_t(txt_size), 0});
        size_t tot_phrases = map.phrase_set.size()-1;//do not count the dummy
        mt_sym = 0, lb = 0;
        off_t i=0, k=0;

        while(mt_sym<tot_phrases) {
            assert(i==lb);
            text[k++] = perm[mt_sym+1];
            i+= map.phrase_set[mt_sym].len;//move out of the phrase boundary

            mt_sym++;
            lb = map.phrase_set[mt_sym].source;//position for the next phrase

            while(i<lb){//process the text area between consecutive phrases
                assert(text[i]<perm.size());
                text[k++] = perm[text[i]];
                i++;
                while(i<lb && text[i]==dummy_sym) i++;
            }
        }
        assert(k==parse_size);
        return parse_size;
    }

    // this method parses the text and store the parse in the text itself.
    // It only works for parsing rounds other than the first one because the length of symbol each
    // cell is the same as the length of cell where we store the metasymbols, so there is no overflow
    off_t inplace_first_par_round(uint8_t* text, off_t txt_size, off_t& n_strings, uint64_t fp_seed,
                                  std::vector<uint64_t>& prev_fps, partial_gram<uint8_t>& p_gram){

        uint32_t mt_sym, sep_sym =0;
        size_t prev_sym, curr_sym, next_sym, dummy_sym=std::numeric_limits<text_chunk::size_type>::max();
        off_t txt_pos = 0, parse_size = 0, phrase_len, lb, rb;
        lz_like_map<uint8_t> map(text);

        bool inserted;
        n_strings = 0;
        off_t sym_bytes = sizeof(uint32_t);

        while(txt_pos<txt_size){

            lb = txt_pos;
            prev_sym = text[txt_pos++];

            curr_sym = text[txt_pos++];
            while(txt_pos<txt_size && curr_sym==prev_sym) curr_sym = text[txt_pos++];
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
            while(txt_pos<txt_size && next_sym==curr_sym) next_sym = text[txt_pos++];
            bool local_minimum;

            while(txt_pos<txt_size && next_sym!=sep_sym){

                local_minimum = prev_sym>curr_sym && curr_sym<next_sym;

                if(local_minimum){
                    phrase_len = rb-lb;
                    mt_sym = map.insert(lb, phrase_len, inserted);
                    if(!inserted){
                        //we can not replace the first phrase occurrence as we use it as source for the dictionary
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
                while(txt_pos<txt_size && next_sym==curr_sym) next_sym = text[txt_pos++];
            }

            phrase_len = txt_pos-1-lb;
            assert(text[lb+phrase_len]==sep_sym);
            mt_sym = map.insert(lb, phrase_len, inserted);
            if(!inserted){
                //we can not replace the first phrase occurrence as we use it as source for the dictionary
                assert(text[lb]!=dummy_sym);
                text[lb] = mt_sym+1;//store the metabmol in the first phrase position
                memset(&text[lb+1], (int)dummy_sym, sym_bytes*(phrase_len-1));//pad the rest of the phrase with dummy symbols
            }
            text[lb+phrase_len] = sep_sym;
            parse_size+=2;//+1 for the separator symbol
            n_strings++;
        }

        map.shrink_to_fit();
        map.destroy_table();

        assert(map.phrase_set.size()<dummy_sym);

        std::vector<uint32_t> perm;
        create_meta_sym<uint8_t, false>(perm, fp_seed, map.phrase_set, text, txt_size, prev_fps, p_gram);

        // create the parse in place
        map.insert_dummy_entry({uint32_t(txt_size), 0});
        size_t tot_phrases = map.phrase_set.size()-1;//do not count the dummy
        mt_sym = 0, lb = 0;
        off_t i=0, k=0;

        while(mt_sym<tot_phrases) {
            assert(i==lb);
            text[k++] = perm[mt_sym+1];
            i+= map.phrase_set[mt_sym].len;//move out of the phrase boundary

            mt_sym++;
            lb = map.phrase_set[mt_sym].source;//position for the next phrase

            while(i<lb){//process the text area between consecutive phrases
                assert(text[i]<perm.size());
                text[k++] = perm[text[i]];
                i++;
                while(i<lb && text[i]==dummy_sym) i++;
            }
        }
        assert(k==parse_size);
        return parse_size;
    }

    off_t first_parsing_round(uint8_t* text, off_t txt_size, text_chunk::size_type* parse,
                              off_t& n_strings, size_t sep_sym, uint64_t fp_seed,
                              std::vector<uint64_t>& prev_fps, partial_gram<uint8_t>& p_gram){

        size_t prev_sym, curr_sym, next_sym;
        uint64_t prev_hash, curr_hash, next_hash;

        text_chunk::size_type mt_sym;
        off_t txt_pos = 0, parse_size = 0, phrase_len, lb, rb;
        lz_like_map map(text);

        bool inserted;
        n_strings = 0;

        auto start = std::chrono::steady_clock::now();
        while(txt_pos<txt_size){

            lb = txt_pos;
            prev_sym = text[txt_pos++];

            curr_sym = text[txt_pos++];
            while(curr_sym==prev_sym && txt_pos<txt_size) curr_sym = text[txt_pos++];
            rb = txt_pos-1;

            if(curr_sym==sep_sym){
                phrase_len = rb-lb;
                mt_sym = map.insert(lb, phrase_len, inserted);
                //assert(text[rb]==sep_sym);

                parse[parse_size++] = mt_sym+1;
                parse[parse_size++] = 0;
                n_strings++;
                continue;
            }

            next_sym = text[txt_pos++];
            while(next_sym==curr_sym && txt_pos<txt_size) next_sym = text[txt_pos++];

            prev_hash = prev_fps[prev_sym];
            curr_hash = prev_fps[curr_sym];

            bool local_minimum;

            while(next_sym!=sep_sym && txt_pos<txt_size){

                next_hash = prev_fps[next_sym];
                local_minimum = prev_hash>curr_hash && curr_hash<next_hash;

                if(local_minimum){
                    phrase_len = rb-lb;
                    mt_sym = map.insert(lb, phrase_len, inserted);

                    parse[parse_size++] = mt_sym+1;
                    lb = rb;
                }

                rb = txt_pos-1;
                prev_hash = curr_hash;
                curr_hash = next_hash;

                curr_sym = next_sym;
                next_sym = text[txt_pos++];
                while(next_sym==curr_sym && txt_pos<txt_size) next_sym = text[txt_pos++];
            }

            phrase_len = txt_pos-1-lb;
            assert(text[lb+phrase_len]==sep_sym);
            mt_sym = map.insert(lb, phrase_len, inserted);

            parse[parse_size++] = mt_sym+1;
            parse[parse_size++] = 0;
            n_strings++;
        }

        map.shrink_to_fit();
        map.destroy_table();

        std::vector<uint32_t> perm;
        create_meta_sym<uint8_t, true>(perm, fp_seed, map.phrase_set, text, txt_size, prev_fps, p_gram);

        for(off_t i=0;i<parse_size;i++){
            assert(parse[i]<perm.size());
            parse[i] = perm[parse[i]];
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
        chunk.p_gram.txt_id = chunk.id;

        off_t parse_size = first_parsing_round(chunk.text, chunk.text_bytes/sizeof(sym_type), chunk.parse,
                                               n_strings, sep_sym, fp_seeds[p_round+1], prev_fps, chunk.p_gram);
        off_t size_limit = n_strings*2;
        p_round++;

        while(parse_size!=size_limit){
            assert(parse_size>=size_limit);
            parse_size = inplace_parsing_round(chunk.parse, parse_size, n_strings, fp_seeds[p_round+1], prev_fps, chunk.p_gram);
            p_round++;
        }
        chunk.p_gram.add_compressed_string(chunk.parse, parse_size);
    }

    template<class sym_type>
    void build_partial_grammars(parsing_opts& p_opts, std::string& text_file, std::string& ct_p_grams_file) {

        ts_queue<size_t> buffers_to_process;
        ts_queue<size_t> buffers_to_reuse;

        std::atomic<size_t> parser_finished{0};
        std::vector<text_chunk> text_chunks(p_opts.n_chunks);

        struct stat st{};
        if (stat(text_file.c_str(), &st) != 0) return;

        auto io_worker = [&]() -> void {

            int fd_r = open(text_file.c_str(), O_RDONLY);

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

                //this value is enough for 4GB text chunks
                text_chunks[chunk_id].p_gram.rules.resize(32);

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

            //data to write
            size_t buff_id;
            int fd_w = open(ct_p_grams_file.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
#ifdef __linux__
            off_t w_acc_bytes=0;
#endif
            //

            size_t proc_syms=0;
            while (rem_bytes > 0) {
                buffers_to_reuse.pop(buff_id);
                size_t g_bytes = text_chunks[buff_id].p_gram.serialize_to_fd(fd_w);
                std::cout<<report_space(text_chunks[buff_id].text_bytes)<<" compressed to "<<report_space((off_t)g_bytes)<<std::endl;
#ifdef __linux__
               w_page_cache_bytes += g_bytes;
                   if(w_page_cache_bytes>p_opts.page_cache_limit){
                       std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                       posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
                       w_page_cache_bytes=0;
                   }
#endif
                proc_syms+=text_chunks[buff_id].text_bytes;
                text_chunks[buff_id].p_gram.reset_grammar();

                //std::cout<<"\n  Processed input "<<report_space((off_t)proc_syms)<<"    "<<std::flush;
                std::cout<<"  Processed input "<<report_space((off_t)proc_syms)<<"    "<<rem_bytes<<std::endl;

                text_chunks[buff_id].text_bytes = tmp_ck_size;
                text_chunks[buff_id].id = chunk_id++;

                read_chunk_from_file<sym_type>(fd_r, rem_bytes, r_acc_bytes, text_chunks[buff_id]);

                //next aligned position
                size_t parse_start =  INT_CEIL(text_chunks[buff_id].text_bytes, sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type);
                text_chunks[buff_id].parse = (text_chunk::size_type *) &text_chunks[buff_id].text[parse_start/sizeof(sym_type)];

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
            posix_fadvise(fd_r, 0, st.st_size, POSIX_FADV_DONTNEED);
#endif
            close(fd_r);

            //wait for all the parsers to finish
            while(parser_finished.load(std::memory_order_acquire)!=p_opts.n_threads);

            //store the remaining grammars in the temporary file
            while(!buffers_to_reuse.empty()){

                buffers_to_reuse.pop(buff_id);
                size_t g_bytes = text_chunks[buff_id].p_gram.serialize_to_fd(fd_w);
                std::cout<<report_space(text_chunks[buff_id].text_bytes)<<" compressed to "<<report_space((off_t)g_bytes)<<std::endl;
#ifdef __linux__
                w_page_cache_bytes += g_bytes;
                if(w_page_cache_bytes>p_opts.page_cache_limit){
                    std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                    posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
                    w_page_cache_bytes=0;
                }
#endif
                proc_syms+=text_chunks[buff_id].text_bytes;
                //std::cout<<"\n  Processed input "<<report_space((off_t)proc_syms)<<"     "<<std::flush;
                std::cout<<"  Processed input "<<report_space((off_t)proc_syms)<<"     "<<std::endl;
            }
            buffers_to_reuse.done();

#ifdef __linux__
            if(w_page_cache_bytes>0){
                std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
            }
#endif
            close(fd_w);
        };

        auto compressor_worker = [&]() {
            size_t buff_id;
            bool res;

            while (true) {
                res = buffers_to_process.pop(buff_id);
                assert(text_chunks[buff_id].text_bytes > 0);

                if (!res){
                    parser_finished.fetch_add(1, std::memory_order_acq_rel);
                    break;
                }

                compress_text_chunk<sym_type>(text_chunks[buff_id], p_opts.p_seeds);
                memset(text_chunks[buff_id].buffer, 0, text_chunks[buff_id].buffer_bytes);

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
            mmap_deallocate(chunk.buffer, chunk.buffer_bytes);
            chunk.buffer=nullptr;
        }
    }

    template<class sym_type>
    void merge_partial_grammars(std::string& ct_p_grams_file, std::string& mg_p_gram_file,
                                std::vector<uint64_t>& p_seeds, size_t n_threads) {

        using p_gram_type = partial_gram<sym_type, true>;


        std::vector<std::pair<p_gram_type,
                    std::vector<std::pair<size_t, size_t>>>
                    > initial_grams(n_threads);

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
                    std::cout<<"Processed data: "<<double(prog_bytes)/double(tot_bytes)*100<<std::endl;
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
                    merge_two_grammars<p_gram_type>(initial_grams[0].first, initial_grams[i].first, p_seeds);
                    //store the range of strings this partial gram covers
                    initial_grams[0].second.insert(initial_grams[0].second.end(), initial_grams[i].second.begin(), initial_grams[i].second.end());
                    initial_grams[i].first.destroy_gram();
                    destroy(initial_grams[i].second);

                    prog_bytes+=rem_prog_bytes;
                    std::cout<<"Processed data: "<<double(prog_bytes)/double(tot_bytes)*100<<std::endl;
                }
                initial_grams[0].second.shrink_to_fit();
#ifdef __linux__
                malloc_trim(0);
#endif
                initial_grams[0].first.reorder_strings(initial_grams[0].second);
            }

            store_to_file(mg_p_gram_file, initial_grams[0].first);
            initial_grams[0].first.destroy_gram();

            std::cout<<"Processed data: "<<100<<std::endl;

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
                merge_two_grammars<p_gram_type>(initial_grams[idx].first, grams_to_merge[buff_id], p_seeds);

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

        std::string ct_p_grams_file = tmp_ws.get_file("concat_p_grams");
        build_partial_grammars<sym_type>(p_opts, i_file, ct_p_grams_file);

        std::string mg_p_gram_file = tmp_ws.get_file("merged_p_grams");
        merge_partial_grammars<sym_type>(ct_p_grams_file, mg_p_gram_file, p_opts.p_seeds, p_opts.n_threads);

        lc_gram_t final_grammar(mg_p_gram_file, p_opts.p_seeds);
        store_to_file(o_file, final_grammar);
    }
}

#endif //LCG_LZ_LIKE_LC_PARSING_H
