//
// Created by Diaz, Diego on 26.4.2023.
//

#ifndef SE_STRAT_PARSE_TEXT_H
#define SE_STRAT_PARSE_TEXT_H

#include "parsing_opts.h"
#include "cds/ts_priority_queue.h"
#include "cds/ts_queue.h"
#include "text_handler.h"

struct heap_node{
    size_t chunk_id;
    size_t chunk_idx;
};

struct heap_node_compare{
    bool operator()(heap_node& a, heap_node& b){
        return a.chunk_id>b.chunk_id;
    }
};

typedef ts_priority_queue<heap_node, std::vector<heap_node>, heap_node_compare> min_heap_type;

template<class parser_t, class map_t, class out_encoder, class text_chunk_t>
struct parse_text {

    //typedef typename parser_t::text_chunk_t text_chunk_t;
    typedef typename out_encoder::sym_type out_sym_type;

    uint64_t operator()(std::string& input_file, std::string& output_file, map_t& map, parsing_opts& opts){

        ts_queue<size_t> chunks_to_read;//chunks containing unprocessed data
        min_heap_type chunks_to_write;//chunks containing a parsed sequence
        std::atomic_uint64_t n_syms{0};

        std::vector<text_chunk_t> text_chunks(opts.active_chunks);
        struct stat st{};
        if(stat(input_file.c_str(), &st) != 0)  return 0;

        auto io_worker = [&]() -> void {

            //data to read
            int fd_r = open(input_file.c_str(), O_RDONLY); //file descriptor to the input file
            off_t rem_bytes = st.st_size, r_acc_bytes=0; //r_acc_bytes = number of bytes read so far
            off_t* chunk_str_ptr = opts.str_ptr->data(); //pointer to the first string start
            size_t chunk_id=0; //an id for a chunk that was read but not processed
            off_t tmp_ck_size;
            //

#ifdef __linux__
            off_t r_page_cache_bytes=0, w_page_cache_bytes=0;
            posix_fadvise(fd_r, 0, st.st_size, POSIX_FADV_SEQUENTIAL);
#endif

            //data to write
            heap_node node{};
            size_t next_chunk=0;
            int fd_w = open(output_file.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
            off_t w_acc_bytes=0;
            //

            while(chunk_id<opts.active_chunks && rem_bytes>0){

                tmp_ck_size = std::min(opts.chunk_size, rem_bytes);
                text_chunks[chunk_id].bytes_cap = tmp_ck_size;
                text_chunks[chunk_id].buffer = (uint8_t *)malloc(tmp_ck_size);
                text_chunks[chunk_id].vbyte_sym_perm = &opts.vbyte_sym_perm;
                text_chunks[chunk_id].sep_sym = opts.sep_sym;
                text_chunks[chunk_id].p_seed = opts.p_seed;
                memset(text_chunks[chunk_id].buffer, 0,  tmp_ck_size);

                text_chunks[chunk_id].id = chunk_id;
                text_chunks[chunk_id].acc_sec_bytes = 0;
                text_chunks[chunk_id].template read_chunk_from_file<parser_t>(fd_r, rem_bytes, r_acc_bytes, chunk_str_ptr);

                if constexpr (std::is_same<out_encoder, vbyte_decoder<out_sym_type>>::value){
                    text_chunks[chunk_id].sec_bytes = opts.max_mt_sym_in_buff*vbyte_len(map.size()-1);
                } else{
                    text_chunks[chunk_id].sec_bytes = opts.max_mt_sym_in_buff*sizeof(out_sym_type);
                }
                text_chunks[chunk_id].sec_buffer = (uint8_t *)malloc(text_chunks[chunk_id].sec_bytes);

#ifdef __linux__
                r_page_cache_bytes+=text_chunks[chunk_id].eff_buff_bytes();
                if(r_page_cache_bytes>opts.page_cache_limit) {
                    std::cout<<"+ removing from page cache "<<r_page_cache_bytes<<" "<<r_acc_bytes<<std::endl;
                    posix_fadvise(fd_r, r_acc_bytes-r_page_cache_bytes, r_page_cache_bytes, POSIX_FADV_DONTNEED);
                    r_page_cache_bytes=0;
                }
#endif
                chunks_to_read.push(chunk_id++);
            }

            size_t buff_idx;
            while(rem_bytes>0) {

                while(true){
                    chunks_to_write.pop(node);
                    assert(text_chunks[node.chunk_idx].bytes>0);
                    if(node.chunk_id==next_chunk){
                        buff_idx = node.chunk_idx;
                        text_chunks[buff_idx].write_chunk_to_file(fd_w, w_acc_bytes);
#ifdef __linux__
                        w_page_cache_bytes += text_chunks[buff_idx].eff_mt_buff_bytes();
                        if(w_page_cache_bytes>opts.page_cache_limit){
                            std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                            posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
                            w_page_cache_bytes=0;
                        }
#endif
                        next_chunk++;
                        break;
                    }else{
                        chunks_to_write.push(node);
                    }
                }

                text_chunks[buff_idx].id = chunk_id++;
                text_chunks[buff_idx].acc_sec_bytes = 0;
                text_chunks[buff_idx].template read_chunk_from_file<parser_t>(fd_r, rem_bytes, r_acc_bytes, chunk_str_ptr);

                //TODO realloc the secondary buffer
                /*if(text_chunks[buff_idx].bytes!=text_chunks[buff_idx].sec_bytes){
                    text_chunks[buff_idx].sec_bytes = text_chunks[buff_idx].bytes*2;
                    text_chunks[buff_idx].sec_buffer = (uint8_t*)realloc(text_chunks[buff_idx].sec_buffer, text_chunks[buff_idx].sec_bytes);
                }*/

#ifdef __linux__
                r_page_cache_bytes+=text_chunks[buff_idx].eff_buff_bytes();
                if(r_page_cache_bytes>opts.page_cache_limit){
                    std::cout<<"+ removing from page cache "<<r_page_cache_bytes<<" "<<r_acc_bytes<<" "<<text_chunks[buff_idx].eff_buff_bytes()<<std::endl;
                    posix_fadvise(fd_r, r_acc_bytes-r_page_cache_bytes, r_page_cache_bytes, POSIX_FADV_DONTNEED);
                    r_page_cache_bytes=0;
                }
#endif
                chunks_to_read.push(buff_idx);
            }

            while(!chunks_to_read.empty());
            chunks_to_read.done();

#ifdef __linux__
            if(r_page_cache_bytes!=0){
                std::cout<<"+ removing from page cache "<<r_page_cache_bytes<<" "<<r_acc_bytes<<std::endl;
                posix_fadvise(fd_r, r_acc_bytes-r_page_cache_bytes, r_page_cache_bytes, POSIX_FADV_DONTNEED);
            }
#endif
            close(fd_r);

            while(next_chunk<chunk_id) {
                while(true){
                    chunks_to_write.pop(node);

                    assert(text_chunks[node.chunk_idx].bytes>0);
                    if(node.chunk_id==next_chunk){
                        buff_idx = node.chunk_idx;
                        text_chunks[buff_idx].write_chunk_to_file(fd_w, w_acc_bytes);
#ifdef __linux__
                        w_page_cache_bytes += text_chunks[buff_idx].eff_mt_buff_bytes();
                        if(w_page_cache_bytes>opts.page_cache_limit){
                            std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                            posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
                            w_page_cache_bytes=0;
                        }
#endif
                        next_chunk++;
                        break;
                    } else {
                        chunks_to_write.push(node);
                    }
                }
            }
            chunks_to_write.done();

#ifdef __linux__
            if(w_page_cache_bytes!=0){
                std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
            }
#endif
            close(fd_w);
        };

        auto lms_worker = [&]() {
            size_t buff_id, n_mt_syms=0;
            bool res;

            while(true) {
                res = chunks_to_read.pop(buff_id);
                if (!res) break;
                assert(text_chunks[buff_id].bytes > 0);
                off_t n_mt = parser_t::template parse_text<map_t, out_encoder>(text_chunks[buff_id], map);

                n_mt_syms +=n_mt;
                chunks_to_write.push({text_chunks[buff_id].id, buff_id});
            }

            n_syms.fetch_add(n_mt_syms, std::memory_order_relaxed);
        };

        std::vector<std::thread> threads;
        threads.emplace_back(io_worker);
        for(size_t i=0;i<opts.n_threads;i++){
            threads.emplace_back(lms_worker);
        }

        for(auto & thread : threads){
            thread.join();
        }
        std::vector<text_chunk_t>().swap(text_chunks);
        return n_syms.load(std::memory_order_acquire);
    }
};
#endif //SE_STRAT_PARSE_TEXT_H
