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
        unsigned long n_chunks{};
        off_t chunk_size{};
        off_t page_cache_limit{};
        size_t sep_sym{};
    };

    struct heap_node_compare{
        bool operator()(const size_t& a, const size_t& b){
            return a>b;
        }
    };
    typedef ts_priority_queue<size_t, std::vector<size_t>, heap_node_compare> min_heap_type;

    template<class sym_type, bool p_round>
    void create_met(std::vector<uint32_t>& mt_perm, hashing& pf,
                    const lz_map::phrase_list_t & phrase_set,
                    sym_type * text, std::vector<uint64_t>& hash_values_vector){

        std::vector<uint64_t> tmp_phrase;
        size_t source, end;
        std::vector<std::pair<uint32_t, uint64_t>> perm(phrase_set.size());

        for(size_t i=0;i<phrase_set.size();i++){

            perm[i].first = i;

            source = phrase_set[i].source/sizeof(sym_type);
            end = source + (phrase_set[i].len/sizeof(sym_type));

            //TODO testing
            /*if(phrase_set[i].source<1000){
                std::cout<<phrase_set[i].len<<" "<<phrase_set[i].source<<" -> ";
                for(size_t j=source;j<end;j++){
                    std::cout<<text[j]<<" ";
                }
                std::cout<<""<<std::endl;
            }*/
            //std::cout<<source<<" "<<end<<" "<<i<<" "<<perm.size()<<std::endl;
            //

            for(size_t j=source;j<end;j++){
                if constexpr (p_round){
                    tmp_phrase.push_back(pf.symbol_hash(text[j]));
                }else{
                    assert(text[j]>0);
                    tmp_phrase.push_back(hash_values_vector[text[j]-1]);
                }
            }
            perm[i].second = pf.string_hash((char *)tmp_phrase.data(), tmp_phrase.size()*sizeof(uint64_t), 64);
            tmp_phrase.clear();
        }
        std::vector<uint64_t>().swap(tmp_phrase);

        //sort the phrases according their hash values
        std::sort(perm.begin(), perm.end(), [&](auto a, auto b) -> bool{
            return a.second<b.second;
        });

        size_t prev_hash = perm[0].second;
        off_t prev_pos=0, tot_phrases=off_t(perm.size());
        for(off_t i=0; i<tot_phrases; i++){

            if(prev_hash!=perm[i].second){
                if((i-prev_pos)>1){
                    std::cout<<"Warning: we have "<<(i-prev_pos)<<" colliding phrases"<<std::endl;

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
                            return pf.symbol_hash(data_a[j])<pf.symbol_hash(data_b[j]);
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
        for(size_t i=0;i<perm.size();i++){
            hash_values_vector[i] = perm[i].second;
        }
        hash_values_vector.shrink_to_fit();

        mt_perm.resize(perm.size());
        for(size_t i=0, rank=0;i<perm.size();i++, rank++){
            mt_perm[perm[i].first] = rank;
        }
    }

    template<class sym_type, bool p_round>
    off_t parsing_round(sym_type* text, off_t txt_size, text_chunk::size_type* parse,
                        off_t& n_strings, size_t sep_sym, hashing& hf,
                        std::vector<uint64_t>& prev_hash_values){

        size_t prev_sym, curr_sym, next_sym;
        text_chunk::size_type mt_sym;
        off_t txt_pos = 0;
        off_t parse_pos = 0;
        off_t phrase_len, lb, rb;
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

            while(next_sym!=sep_sym){

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
            //assert(text[lb+phrase_len-1]>0);
            //assert(text[lb]>0);
            mt_sym = map.insert(lb*sym_bytes, phrase_len*sym_bytes, inserted);
            parse[parse_pos++] = mt_sym+1;
            parse[parse_pos++] = 0;

            //std::cout<<n_strings<<" "<<txt_pos<<" "<<parse_pos<<" "<<txt_pos/(n_strings+1)<<std::endl;
            n_strings++;
        }

        //map.shrink_to_fit();
        //map.destroy_table();

        std::vector<uint32_t> perm;

        /*size_t j = 0;
        while(true){
            std::cout<<text[j]<<" ";
            if((p_round && text[j]==10) || text[j]==0){
                break;
            }
            j++;
        }
        std::cout<<""<<std::endl;*/

        create_met<sym_type, p_round>(perm, hf, map.phrase_set, text, prev_hash_values);

        for(off_t i=0;i<parse_pos;i++){
            if(parse[i]==0) continue;
            parse[i] = perm[parse[i]-1]+1;
        }

        return parse_pos;
    }

    template<class sym_type>
    void compress_text_chunk(text_chunk& chunk, tmp_workspace& tmp_ws, std::vector<hashing>& pf){

        std::vector<uint64_t> prev_hash_values;
        off_t n_strings=0;
        size_t sep_sym = chunk.sep_sym;
        off_t parse_size = parsing_round<sym_type, true>(chunk.text, chunk.text_bytes/sizeof(sym_type), chunk.parse, n_strings, sep_sym, pf[0], prev_hash_values);
        sep_sym = 0;
        off_t size_limit = n_strings*2;
        //std::cout<<"psize: "<<parse_size<<", size_limit: "<<size_limit<<". text_size: "<<chunk.text_bytes/sizeof(sym_type)<<", n_strings: "<<n_strings<<std::endl;
        size_t i=1;

        auto * new_parse = (text_chunk::size_type *)chunk.text;

        while(parse_size!=size_limit){
            assert(parse_size>=size_limit);
            parse_size = parsing_round<text_chunk::size_type, false>(chunk.parse, parse_size, new_parse, n_strings, sep_sym, pf[i++], prev_hash_values);
            std::swap(chunk.parse, new_parse);
            //std::cout<<"psize: "<<parse_size<<", size_limit: "<<size_limit<<". text_size: "<<chunk.text_bytes/sizeof(sym_type)<<", n_strings: "<<n_strings<<std::endl;
        }
    }

    template<class sym_type>
    void lc_parsing_algo(std::string& i_file, std::vector<hashing>& phf, std::string& o_file,
                         tmp_workspace& tmp_ws, parsing_opts& p_opts) {

        std::cout<<"  Settings"<<std::endl;
        std::cout<<"    Parsing threads           : "<<p_opts.n_threads<<std::endl;
        std::cout<<"    Active text chunks in RAM : "<<p_opts.n_chunks<<std::endl;
        std::cout<<"    Size of each chunk        : "<<report_space(p_opts.chunk_size)<<std::endl;
        std::cout<<"    Chunks' approx. mem usage : "<<report_space(off_t(p_opts.chunk_size*p_opts.n_chunks*3))<<"\n"<<std::endl;

        ts_queue<size_t> chunks_to_read;
        ts_queue<size_t> chunks_to_reuse;
        //min_heap_type chunks_to_merge;

        std::atomic<size_t> parser_finished{0};
        std::vector<text_chunk> text_chunks(p_opts.n_chunks);
        struct stat st{};
        if (stat(i_file.c_str(), &st) != 0) return;

        auto read_worker = [&]() -> void {

            int fd = open(i_file.c_str(), O_RDONLY);

#ifdef __linux__
            posix_fadvise(fd, 0, st.st_size, POSIX_FADV_SEQUENTIAL);
#endif
            off_t rem_bytes = st.st_size, acc_bytes = 0;
            size_t chunk_id = 0;

            auto tmp_ck_size = off_t((p_opts.chunk_size/sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type));

            while (chunk_id < p_opts.n_chunks && rem_bytes > 0) {

                tmp_ck_size = std::min(tmp_ck_size, rem_bytes);

                text_chunks[chunk_id].phf = &phf;
                text_chunks[chunk_id].text_bytes = tmp_ck_size;
                text_chunks[chunk_id].sep_sym = (text_chunk::size_type) p_opts.sep_sym;

                //the parse size is (text_len/2)*(sizeof(size_type)/sizeof(sym_type)),
                // where ``text_len'' is the number of input symbols that fits the buffer
                off_t parse_bytes = INT_CEIL((tmp_ck_size/sizeof(sym_type)), 2)*(sizeof(text_chunk::size_type)/sizeof(sym_type));

                text_chunks[chunk_id].buffer_bytes = off_t(tmp_ck_size + parse_bytes);
                text_chunks[chunk_id].buffer = (text_chunk::size_type *) malloc(text_chunks[chunk_id].buffer_bytes);
                text_chunks[chunk_id].id = chunk_id;

                read_chunk_from_file(fd, rem_bytes, acc_bytes, text_chunks[chunk_id]);

                //next aligned position within the buffer
                text_chunks[chunk_id].text = (sym_type *) text_chunks[chunk_id].buffer;
                size_t parse_start =  (text_chunks[chunk_id].text_bytes/sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type);
                text_chunks[chunk_id].parse = (text_chunk::size_type *) &text_chunks[chunk_id].text[parse_start/sizeof(sym_type)];

                chunks_to_read.push(chunk_id);
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

            size_t proc_syms=0;
            size_t buff_idx;

            while (rem_bytes > 0) {
                chunks_to_reuse.pop(buff_idx);

                proc_syms+=text_chunks[buff_idx].text_bytes;
                std::cout<<"\rProcessed symbols "<<report_space((off_t)proc_syms)<<"    "<<std::flush;

                text_chunks[buff_idx].id = chunk_id++;
                read_chunk_from_file(fd, rem_bytes, acc_bytes, text_chunks[buff_idx]);

                text_chunks[buff_idx].text = (sym_type *) text_chunks[buff_idx].buffer;
                //next aligned position
                size_t parse_start =  (text_chunks[buff_idx].text_bytes/sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type);
                text_chunks[buff_idx].parse = (text_chunk::size_type *) &text_chunks[buff_idx].text[parse_start/sizeof(sym_type)];

                chunks_to_read.push(buff_idx);
#ifdef __linux__
                page_cache_bytes+=text_chunks[chunk_id].eff_buff_bytes();
                if(page_cache_bytes>page_cache_limit){
                    std::cout<<"removing from page cache "<<page_cache_bytes<<" "<<acc_bytes<<std::endl;
                    posix_fadvise(fd, acc_bytes-page_cache_bytes, page_cache_bytes, POSIX_FADV_DONTNEED);
                    page_cache_bytes=0;
                }
#endif
            }
            while (!chunks_to_read.empty());
            chunks_to_read.done();

            //wait for all the parsers to finish
            while(parser_finished.load(std::memory_order_acquire)!=p_opts.n_threads);

            while (!chunks_to_reuse.empty()) {
                chunks_to_reuse.pop(buff_idx);
                proc_syms+=text_chunks[buff_idx].text_bytes;
                std::cout<<"\r"<<"Processed symbols "<<report_space((off_t)proc_syms)<<"     "<<std::flush;
            }
            chunks_to_reuse.done();
            std::cout<<""<<std::endl;
#ifdef __linux__
            posix_fadvise(fd, 0, st.st_size, POSIX_FADV_DONTNEED);
#endif
            close(fd);
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

                compress_text_chunk<sym_type>(text_chunks[buff_id], tmp_ws, phf);
                chunks_to_reuse.push(buff_id);
                //chunks_to_merge.push(text_chunks[buff_id].id);
            }
        };

        /*auto gram_merge_worker = [&]() {

            //partial grammars to marge
            size_t chunk_id;
            size_t next_chunk=0;

            while(true){
                chunks_to_merge.pop(chunk_id);
                assert(text_chunks[chunk_id].text_bytes>0);
                if(chunk_id==next_chunk){
                    //TODO merge the grammars
                    next_chunk++;
                    break;
                }else{
                    chunks_to_merge.push(chunk_id);
                }
            }

            //wait for all the chunks to be processed
            while(!chunks_to_reuse.empty());

            while(!chunks_to_merge.empty()) {
                while(true){
                    chunks_to_merge.pop(chunk_id);
                    assert(text_chunks[chunk_id].text_bytes>0);
                    if(chunk_id==next_chunk){
                        //TODO merge grammars
                        next_chunk++;
                        break;
                    } else {
                        chunks_to_merge.push(chunk_id);
                    }
                }
            }
            chunks_to_merge.done();
        };*/

        std::vector<std::thread> threads;
        threads.emplace_back(read_worker);
        //threads.emplace_back(gram_merge_worker);

        for (size_t i = 0; i < p_opts.n_threads; i++) {
            threads.emplace_back(parser_worker);
        }

        for (auto &thread: threads) {
            thread.join();
        }
    };
}

#endif //LCG_LZ_LIKE_LC_PARSING_H
