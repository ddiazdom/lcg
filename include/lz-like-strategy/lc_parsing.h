//
// Created by Diaz, Diego on 17.9.2023.
//

#ifndef LCG_LZ_LIKE_LC_PARSING_H
#define LCG_LZ_LIKE_LC_PARSING_H

#include "lz-like-strategy/text_handler.h"
#include "cds/utils.h"
#include "lz_map.h"

namespace lzstrat {

    struct parsing_opts{
        size_t n_threads{};
        size_t n_chunks{};
        off_t chunk_size{};
        off_t page_cache_limit{};
        size_t sep_sym{};
    };

    template<class sym_type>
    off_t parsing_round(sym_type* text, off_t txt_size, text_chunk::size_type* parse, off_t& n_strings, size_t sep_sym){

        size_t prev_sym, curr_sym, next_sym;
        text_chunk::size_type mt_sym;

        off_t txt_pos = 0;
        off_t parse_pos = 0;
        off_t phrase_len, lb, rb;
        lz_map map((uint8_t *)text);
        bool inserted;
        n_strings = 0;

        while(txt_pos<txt_size){

            lb = txt_pos;
            prev_sym = text[txt_pos++];

            curr_sym = text[txt_pos++];
            while(txt_pos<txt_size && curr_sym==prev_sym) curr_sym = text[txt_pos++];
            rb = txt_pos-1;

            if(curr_sym==sep_sym){
                n_strings++;
                mt_sym = map.insert(txt_pos-2, sizeof(sym_type), inserted);
                parse[parse_pos++] = mt_sym;
                parse[parse_pos++] = 0;
                continue;
            }

            next_sym = text[txt_pos++];
            while(txt_pos<txt_size && next_sym==curr_sym) next_sym = text[txt_pos++];

            while(next_sym!=sep_sym){

                if(prev_sym>curr_sym && curr_sym<next_sym){
                    phrase_len = rb-lb;
                    mt_sym = map.insert(lb, phrase_len*sizeof(sym_type), inserted);
                    parse[parse_pos++] = mt_sym;
                    lb = rb;
                }

                rb = txt_pos-1;
                prev_sym = curr_sym;
                curr_sym = next_sym;
                next_sym = text[txt_pos++];
                while(txt_pos<txt_size && next_sym==curr_sym) next_sym = text[txt_pos++];
            }

            phrase_len = txt_pos-1-lb;
            mt_sym = map.insert(lb, phrase_len*sizeof(sym_type), inserted);
            parse[parse_pos++] = mt_sym;
            parse[parse_pos++] = 0;

            n_strings++;
        }
        std::cout<<map.table_mem_usage()<<" "<<map.size()<<" "<<map.capacity()<<std::endl;
        return parse_pos;
    }

    template<class sym_type>
    void compress_text_chunk(text_chunk& chunk, tmp_workspace& tmp_ws){

        off_t n_strings=0;
        size_t sep_sym = chunk.sep_sym;

        auto *text = (sym_type *)chunk.data_start;
        off_t text_size = chunk.data_bytes/sizeof(sym_type);
        text_chunk::size_type * parse = chunk.buffer;
        std::cout<<"Parsing round 1"<<std::endl;
        off_t parse_size = parsing_round<sym_type>(text, text_size, parse, n_strings, sep_sym);
        sep_sym = 0;
        off_t size_limit = n_strings*2;

        std::cout<<parse_size<<" "<<size_limit<<" "<<text_size<<" "<<n_strings<<std::endl;
        while(parse_size!=size_limit){
            assert(parse_size>=size_limit);
            parse_size = parsing_round<text_chunk::size_type>(parse, parse_size, parse, n_strings, sep_sym);
            std::cout<<parse_size<<" "<<size_limit<<" "<<text_size<<" "<<n_strings<<std::endl;
        }
    }

    template<class sym_type>
    void lc_parsing_algo(std::string& i_file, std::vector<hashing>& phf, std::string& o_file,
                         tmp_workspace& tmp_ws, parsing_opts& p_opts) {

        ts_queue<size_t> in_strings;
        ts_queue<size_t> out_strings;

        std::vector<text_chunk> text_chunks(p_opts.n_chunks, phf);
        struct stat st{};
        if (stat(i_file.c_str(), &st) != 0) return;

        auto read_worker = [&]() -> void {

            int fd = open(i_file.c_str(), O_RDONLY);

#ifdef __linux__
            posix_fadvise(fd, 0, st.st_size, POSIX_FADV_SEQUENTIAL);
#endif
            off_t rem_bytes = st.st_size, acc_bytes = 0;
            size_t chunk_id = 0;
            off_t tmp_ck_size;

            while (chunk_id < p_opts.n_chunks && rem_bytes > 0) {

                tmp_ck_size = std::min(p_opts.chunk_size, rem_bytes);

                text_chunks[chunk_id].data_bytes = tmp_ck_size;
                text_chunks[chunk_id].sep_sym = (text_chunk::size_type) p_opts.sep_sym;
                text_chunks[chunk_id].buffer_bytes = off_t(tmp_ck_size*sizeof(text_chunk::size_type));
                text_chunks[chunk_id].buffer = (text_chunk::size_type *) malloc(text_chunks[chunk_id].buffer_bytes);
                text_chunks[chunk_id].id = chunk_id;

                read_chunk_from_file(fd, rem_bytes, acc_bytes, text_chunks[chunk_id]);
                in_strings.push(chunk_id);
                chunk_id++;

#ifdef __linux__
                page_cache_bytes+=text_chunks[chunk_id].eff_buff_bytes();
                if(page_cache_bytes>page_cache_limit){
                    std::cout<<"removing from page cache "<<page_cache_bytes<<" "<<acc_bytes<<std::endl;
                    posix_fadvise(fd, acc_bytes-page_cache_bytes, page_cache_bytes, POSIX_FADV_DONTNEED);
                    page_cache_bytes=0;
                }
#endif
            }

            size_t buff_idx;
            while (rem_bytes > 0) {
                out_strings.pop(buff_idx);
                text_chunks[buff_idx].id = chunk_id++;
                read_chunk_from_file(fd, rem_bytes, acc_bytes, text_chunks[buff_idx]);
                in_strings.push(buff_idx);
#ifdef __linux__
                page_cache_bytes+=text_chunks[chunk_id].eff_buff_bytes();
                if(page_cache_bytes>page_cache_limit){
                    std::cout<<"removing from page cache "<<page_cache_bytes<<" "<<acc_bytes<<std::endl;
                    posix_fadvise(fd, acc_bytes-page_cache_bytes, page_cache_bytes, POSIX_FADV_DONTNEED);
                    page_cache_bytes=0;
                }
#endif
            }
            while (!in_strings.empty());

            while (!out_strings.empty()) {
                out_strings.pop(buff_idx);
            }

            in_strings.done();
            out_strings.done();
#ifdef __linux__
            posix_fadvise(fd, 0, st.st_size, POSIX_FADV_DONTNEED);
#endif
            close(fd);
        };

        auto parser_worker = [&]() {

            size_t buff_id;
            bool res;

            while (true) {
                res = in_strings.pop(buff_id);
                assert(text_chunks[buff_id].data_bytes > 0);

                if (!res) break;
                compress_text_chunk<sym_type>(text_chunks[buff_id], tmp_ws);
                out_strings.push(buff_id);
            }
        };

        auto gram_merge_worker = [&]() {

        };

        std::vector<std::thread> threads;
        threads.emplace_back(read_worker);
        for (size_t i = 0; i < p_opts.n_threads; i++) {
            threads.emplace_back(parser_worker);
        }

        for (auto &thread: threads) {
            thread.join();
        }
    };
}

#endif //LCG_LZ_LIKE_LC_PARSING_H
