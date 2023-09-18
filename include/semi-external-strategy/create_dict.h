//
// Created by Diaz, Diego on 26.4.2023.
//

#ifndef SIMPLE_EXAMPLE_CREATE_DICTIONARY_H
#define SIMPLE_EXAMPLE_CREATE_DICTIONARY_H

#include "local_minima_parser.hpp"
#include "text_handler.h"
#include "cds/ts_queue.h"
#include <unistd.h>

template<class parser_t, class map_t>
struct create_dict{

    typedef typename parser_t::text_chunk text_chunk;

    void operator()(std::string& input_file, map_t& map,
                    std::vector<off_t>& str_ptrs, parsing_opts& opts) {


        ts_queue<size_t> in_strings;
        ts_queue<size_t> out_strings;

        std::vector<text_chunk> text_chunks(opts.active_chunks);
        struct stat st{};
        if(stat(input_file.c_str(), &st) != 0)  return;

        auto read_worker = [&]() -> void {

            int fd = open(input_file.c_str(), O_RDONLY);

#ifdef __linux__
            posix_fadvise(fd, 0, st.st_size, POSIX_FADV_SEQUENTIAL);
#endif

            off_t rem_bytes = st.st_size, acc_bytes=0;
            off_t * chunk_str_ptr = str_ptrs.data();

            size_t chunk_id=0;
            off_t tmp_ck_size;

            while(chunk_id<opts.active_chunks && rem_bytes>0){

                tmp_ck_size = std::min(opts.chunk_size, rem_bytes);
                text_chunks[chunk_id].bytes = tmp_ck_size;
                //text_chunks[chunk_id].buff_cap = tmp_ck_size/sizeof(sym_type);
                text_chunks[chunk_id].buffer = (uint8_t *)malloc(tmp_ck_size);
                text_chunks[chunk_id].id = chunk_id;

                text_chunks[chunk_id].template read_chunk_from_file<parser_t>(fd, rem_bytes, acc_bytes, chunk_str_ptr);
                in_strings.push(chunk_id);
                chunk_id++;
            }

            size_t buff_idx;
            while(rem_bytes>0){
                out_strings.pop(buff_idx);
                text_chunks[buff_idx].id = chunk_id++;
                text_chunks[buff_idx].template read_chunk_from_file<parser_t>(fd, rem_bytes, acc_bytes, chunk_str_ptr);
                in_strings.push(buff_idx);
            }
            while(!in_strings.empty());

            while(!out_strings.empty()){
                out_strings.pop(buff_idx);
            }

            in_strings.done();
            out_strings.done();

#ifdef __linux__
            posix_fadvise(fd, 0, st.st_size, POSIX_FADV_DONTNEED);
#endif
            close(fd);
        };

        auto lms_worker = [&]() {

            map.register_thread();

            size_t buff_id;
            bool res;

            while(true){
                res = in_strings.pop(buff_id);
                assert(text_chunks[buff_id].bytes>0);

                if(!res) break;
                parser_t::scan(text_chunks[buff_id], map);
                //lms_breaks_left2right<text_chunk, false, map_t>(text_chunks[buff_id], map);
                out_strings.push(buff_id);
            }
        };

        std::vector<std::thread> threads;
        threads.emplace_back(read_worker);
        for(size_t i=0;i<opts.n_threads;i++){
            threads.emplace_back(lms_worker);
        }

        for(auto & thread : threads){
            thread.join();
        }

        std::vector<text_chunk>().swap(text_chunks);
        std::cout<<"Map has "<<map.size()<<" phrases "<<std::endl;
    }
};
#endif //SIMPLE_EXAMPLE_CREATE_DICTIONARY_H
