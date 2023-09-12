//
// Created by Diaz, Diego on 26.4.2023.
//

#ifndef SIMPLE_EXAMPLE_CREATE_DICT_FULL_SCAN_H
#define SIMPLE_EXAMPLE_CREATE_DICT_FULL_SCAN_H

#include "text_handler.h"
#include <unistd.h>

template<class parser_t, class map_t>
struct create_dict_full_scan {

    typedef typename parser_t::text_chunk_t text_chunk_t;

    size_t operator()(std::string& input_file, map_t& map, parsing_opts& opts) {

        //using chunk_type = text_chunk<sym_type, sym_type>;
        std::vector<text_chunk_t> text_chunks(opts.active_chunks);
        struct stat st{};
        if(stat(input_file.c_str(), &st) != 0)  return 0;

        std::atomic_long read_chunks{-1};
        std::atomic_bool read_finished{false};
        std::atomic<size_t> max_mt_syms_global{0};

        auto read_worker = [&]() -> void {

            int fd = open(input_file.c_str(), O_RDONLY);

#ifdef __linux__
            off_t page_cache_bytes=0;
            posix_fadvise(fd, 0, st.st_size, POSIX_FADV_SEQUENTIAL);
#endif
            off_t rem_bytes = st.st_size, acc_bytes=0;
            off_t * chunk_str_ptr = opts.str_ptr->data();

            size_t chunk_id=0, circular_id=0;
            off_t tmp_ck_size;

            while(chunk_id<opts.active_chunks && rem_bytes>0){
                tmp_ck_size = std::min(opts.chunk_size, rem_bytes);
                text_chunks[chunk_id].bytes = tmp_ck_size;
                text_chunks[chunk_id].buffer = (uint8_t *)malloc(tmp_ck_size);
                text_chunks[chunk_id].sym_perm = &opts.sym_perm;

                text_chunks[chunk_id].id = chunk_id;
                text_chunks[chunk_id].template read_chunk_from_file<parser_t>(fd, rem_bytes, acc_bytes, chunk_str_ptr);

#ifdef __linux__
                page_cache_bytes+=text_chunks[chunk_id].eff_buff_bytes();
                if(page_cache_bytes>opts.page_cache_limit){
                    std::cout<<"removing from page cache "<<page_cache_bytes<<" "<<acc_bytes<<std::endl;
                    posix_fadvise(fd, acc_bytes-page_cache_bytes, page_cache_bytes, POSIX_FADV_DONTNEED);
                    page_cache_bytes=0;
                }
#endif
                read_chunks.store((long)chunk_id++, std::memory_order_release);
            }

            long long n_scans;
            while(rem_bytes>0){
                n_scans = text_chunks[circular_id].scans.load(std::memory_order_acquire);
                if(n_scans<(long)opts.n_threads) continue;//wait for all the lms workers to mark the chunk
                assert(n_scans==(long)opts.n_threads);

                text_chunks[circular_id].id = chunk_id;
                text_chunks[circular_id].template read_chunk_from_file<parser_t>(fd, rem_bytes, acc_bytes, chunk_str_ptr);
#ifdef __linux__
                page_cache_bytes+=text_chunks[circular_id].eff_buff_bytes();
                if(page_cache_bytes>opts.page_cache_limit){
                    std::cout<<"removing from page cache "<<page_cache_bytes<<" "<<acc_bytes<<" "<<text_chunks[circular_id].eff_buff_bytes()<<std::endl;
                    posix_fadvise(fd, acc_bytes-page_cache_bytes, page_cache_bytes, POSIX_FADV_DONTNEED);
                    page_cache_bytes=0;
                }
#endif
                text_chunks[circular_id].scans.store(0, std::memory_order_release);
                read_chunks.store((long)chunk_id++, std::memory_order_release);

                circular_id = (circular_id+1) % opts.active_chunks;
            }

            //while for all the chunks to be read
            size_t eff_active_chunks = std::min(opts.active_chunks, chunk_id);
            while(true){
                size_t comp=0;
                for(size_t i=0;i<eff_active_chunks;i++){
                    comp +=text_chunks[i].scans.load(std::memory_order_acquire)==(long)opts.n_threads;
                }
                if(comp==eff_active_chunks) break;
            }
            read_finished.store(true, std::memory_order_release);

#ifdef __linux__
            if(page_cache_bytes!=0){
                std::cout<<"removing from page cache "<<page_cache_bytes<<" "<<acc_bytes<<std::endl;
                posix_fadvise(fd, acc_bytes-page_cache_bytes, page_cache_bytes, POSIX_FADV_DONTNEED);
            }
#endif
            close(fd);
        };

        auto lms_worker = [&]() {

            map.register_thread();
            size_t circular_id=0;

            //chunk_id is the current chunk this thread is reading
            //r_chunk is the last chunk that was successfully read
            long chunk_id=0, r_chunks;
            size_t max_mt_syms=0;
            while(true){
                r_chunks = read_chunks.load(std::memory_order_acquire);
                if(chunk_id>r_chunks){
                    if(read_finished.load(std::memory_order_acquire)){
                        break;
                    }
                    continue; //wait until the next chunk is available
                }

                assert(text_chunks[circular_id].id == (size_t)chunk_id);
                size_t mt_syms = parser_t::template get_phrases<map_t>(text_chunks[circular_id], map);

                if(mt_syms>max_mt_syms){
                    max_mt_syms = mt_syms;
                }

                text_chunks[circular_id].scans.fetch_add(1, std::memory_order_release);
                circular_id = (circular_id+1)%opts.active_chunks;
                chunk_id++;
            }

            if(max_mt_syms>max_mt_syms_global.load(std::memory_order::memory_order_acquire)){
                max_mt_syms_global.store(max_mt_syms, std::memory_order_release);
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
        std::vector<text_chunk_t>().swap(text_chunks);

        return max_mt_syms_global.load(std::memory_order_acquire);
    }
};
#endif //SIMPLE_EXAMPLE_CREATE_DICT_FULL_SCAN_H
