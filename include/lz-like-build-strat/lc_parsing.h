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
#include "lz_like_map.h"
#include "grammar.h"
#include "merge_grams.h"

namespace lz_like_strat {

    struct parsing_opts{
        size_t n_threads{};
        unsigned long n_chunks{};
        off_t chunk_size{};
        off_t page_cache_limit{};
        size_t sep_sym{};
        uint64_t orig_seed;//this is the seed to create random seeds for p_seeds
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

    struct inv_perm_elm{
        uint32_t orig_mt;
        uint64_t fp;
    };

    template<class sym_type, bool p_round>
    off_t create_meta_sym(text_chunk& chunk, buff_vector<uint32_t>& mt_perm, uint64_t pf_seed,
                         const typename lz_like_map<sym_type>::phrase_list_t& phrase_set,
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
        buff_vector<inv_perm_elm> inv_mt_perm(avb_addr.first, avb_addr.second);
        inv_mt_perm.resize(phrase_set.size());
        chunk.add_used_bytes((off_t)inv_mt_perm.static_buff_usage());


        size_t source, end, tot_symbols=0;
        std::vector<uint64_t> fp_sequence;
        for(size_t i=0;i<phrase_set.size();i++) {
            inv_mt_perm[i].orig_mt = i;
            source = phrase_set[i].source;
            end = source + phrase_set[i].len;
            for(size_t j=source;j<end;j++){
                assert(text[j]>0);
                assert(text[j]<prev_fps.size());
                fp_sequence.push_back(prev_fps[text[j]]);
            }
            inv_mt_perm[i].fp = XXH64(fp_sequence.data(), fp_sequence.size()*sizeof(uint64_t), pf_seed);
            tot_symbols+=fp_sequence.size();
            fp_sequence.clear();
        }

        //sort the phrases according their hash values
        std::sort(inv_mt_perm.get_data(), inv_mt_perm.get_data()+inv_mt_perm.size(), [&](auto const &a, auto const& b) -> bool{
            return a.fp<b.fp;
        });

        size_t prev_hash = inv_mt_perm[0].fp;
        off_t prev_pos=0, tot_phrases=off_t(inv_mt_perm.size());
        off_t n_cols =0;
        for(off_t i=1; i<tot_phrases; i++){

            if(prev_hash!=inv_mt_perm[i].fp){
                if((i-prev_pos)>1){
                    n_cols +=(i-prev_pos);
                    std::cout<<"Warning: we have "<<(i-prev_pos)<<" colliding phrases"<<std::endl;
                    //TODO testing
                    for(off_t k=prev_pos;k<i;k++) {
                        std::cout<<"\tsorted pos:"<<k<<" orig pos:"<<inv_mt_perm[k].orig_mt<<", "<<phrase_set[inv_mt_perm[k].orig_mt].to_string()<<", phrase: ";
                        off_t s = phrase_set[inv_mt_perm[k].orig_mt].source;
                        off_t len = phrase_set[inv_mt_perm[k].orig_mt].len;
                        assert(s<txt_size);

                        fp_sequence.clear();
                        for (off_t u = s, l=0; l < len; u++, l++) {
                            std::cout<<text[u]<<" ";
                            fp_sequence.push_back(prev_fps[text[u]]);
                        }
                        std::cout<<""<<std::endl;
                        uint64_t test_hash = XXH64(fp_sequence.data(), fp_sequence.size()*sizeof(uint64_t), pf_seed);
                        assert(inv_mt_perm[k].fp==test_hash);
                    }
                    std::cout<<""<<std::endl;
                    //

                    //sort the range [prev_pos..i-1]
                    std::sort(inv_mt_perm.get_data()+prev_pos, inv_mt_perm.get_data()+i, [&](auto const &a, auto const&b) -> bool{

                        typename lz_like_map<sym_type>::phrase_t phrase_a = phrase_set[a.orig_mt];
                        sym_type * data_a = &text[phrase_a.source];
                        off_t len_a = phrase_a.len;

                        typename lz_like_map<sym_type>::phrase_t phrase_b = phrase_set[b.orig_mt];
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
                prev_hash = inv_mt_perm[i].fp;
            }
        }

        avb_addr = chunk.get_free_mem_area();
        buff_vector<uint64_t> new_fps(avb_addr.first, avb_addr.second);
        new_fps.resize(inv_mt_perm.size() + 1);
        //size_t bytes_in_ram = inv_mt_perm.eff_mem_usage() + mt_perm.eff_mem_usage() + phrase_set.eff_mem_usage() + prev_fps.eff_mem_usage() + chunk.p_gram.buff_bytes() + new_fps.eff_mem_usage() + chunk.buffer_bytes;
        //std::cout <<"\n\tI have approximately " << report_space((off_t) bytes_in_ram) << " in main memory versus " << report_space((off_t) malloc_count_current()) << "/" << report_space((off_t) malloc_count_peak()) << std::endl;
        //std::cout<<"inv_perm="<<inv_mt_perm.called_malloc()<<", perm="<<mt_perm.called_malloc()<<", phrases="<<phrase_set.called_malloc()<<" prev_fps="<<prev_fps.called_malloc()<<" av_bytes="<<report_space((off_t)chunk.av_bytes())<<"/"<< report_space((off_t) chunk.buffer_bytes)<<", mem_peak:"<< report_space((off_t)malloc_count_peak())<<std::endl;
        //std::cout<<report_space((off_t)inv_mt_perm.mem_usage())<<", "<<report_space((off_t)mt_perm.mem_usage())<<", "<<report_space((off_t)phrase_set.mem_usage())<<", "<<report_space((off_t) prev_fps.mem_usage())<<std::endl;
        //malloc_count_reset_peak();
        //std::cout<<prev_fps.eff_mem_usage()<<" "<<new_fps.eff_mem_usage()<<std::endl;

        new_fps[0] = 0;
        mt_perm[0] = 0;
        for (size_t i = 0, mt_sym = 1; i < inv_mt_perm.size(); i++, mt_sym++) {
            size_t perm_mt_sym = inv_mt_perm[i].orig_mt + 1;
            assert(perm_mt_sym < mt_perm.size());
            mt_perm[perm_mt_sym] = mt_sym;
            new_fps[mt_sym] = inv_mt_perm[i].fp;
        }
        new_fps.swap(prev_fps);
        new_fps.destroy();
        //std::cout<<prev_fps.eff_mem_usage()<<" "<<new_fps.eff_mem_usage()<<std::endl;

        //std::cout <<"\t\tDeberia ser igual o menor " << report_space((off_t) malloc_count_current()) << "/" << report_space((off_t) malloc_count_peak())<<std::endl;
        chunk.p_gram.template append_new_lvl<sym_type>(text, phrase_set, tot_symbols, inv_mt_perm);
        //std::cout <<"\t\tDeberia haber crecido " << report_space((off_t) malloc_count_current()) << "/" << report_space((off_t) malloc_count_peak())<<std::endl;
        //std::cout<<report_space(off_t(inv_mt_perm.size())*off_t(sizeof(std::pair<uint32_t, uint64_t>)))<<" "<<report_space(off_t(prev_fps.size())*8)<<" "<<report_space(off_t(mt_perm.size())*4)<<" "<<report_space(off_t(phrase_set.size())*8)<<" "<<report_space(lvl_bytes)<<std::endl;
        return n_cols;
    }

    void finish_byte_parse(text_chunk& chunk, lz_like_map<uint8_t>& dict,
                           off_t &parse_distance, buff_vector<uint32_t>& mt_perm,
                           std::vector<phrase_overflow>& phr_with_ovf){

        uint32_t mt_sym=1;
        uint8_t v_len;
        off_t txt_size = chunk.text_bytes;

        for(auto const& phrase : dict.phrase_set){
            v_len = vbyte_len(mt_sym);
            if(v_len<=phrase.len){
                vbyte_decoder<uint32_t>::write_right2left(&chunk.text[phrase.source], mt_sym, v_len);
                memset(&chunk.text[phrase.source+v_len], 0, phrase.len-v_len);
            }
            mt_sym++;
        }

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
                *chunk.parse = mt_perm[phr_with_ovf[ovf_idx].metasymbol];
                next_ovf=-1;
                if(ovf_idx<off_t(phr_with_ovf.size()-1)){
                    next_ovf = phr_with_ovf[++ovf_idx].right_end();
                }
            }  else {
                pos-=vbyte_decoder<uint32_t>::read_right2left(&chunk.text[pos], mt_sym);
                assert(mt_sym<mt_perm.size());
                *chunk.parse = mt_perm[mt_sym];
            }
            chunk.parse--;
            while(pos>next_ovf && chunk.text[pos]==0) pos--;
        }
        assert(pos==-1);
        chunk.parse++;
        assert(chunk.parse[chunk.parse_size-1]==0);
        assert((uintptr_t) chunk.text<= (uintptr_t)chunk.parse && (uintptr_t)chunk.parse<= (uintptr_t)&chunk.text[txt_size-1]);
        assert((uintptr_t) chunk.text<= (uintptr_t)(chunk.parse+chunk.parse_size) &&
               (uintptr_t) (chunk.parse+chunk.parse_size) <= (uintptr_t)(chunk.text+chunk.buffer_bytes));
        assert(parse_distance==(off_t)chunk.dist(reinterpret_cast<uint8_t*>(chunk.parse+chunk.parse_size)));
    }

    //this function scans a text with a byte alphabet from right to left and parses it according its randomized LMS phrases.
    // we perform the parsing in place, meaning that we store the parse in directly in the input text reinterpreted as an unit32_t sequence.
    // this strategy is space-efficient but challenging to implement because the metasymbol mt we assign to a phrase F could use more bytes
    // if n_bytes(F)>4 and 4<=n_bytes(mt)> n_bytes(F). We encode the metasymbols using vbyte encoding and those codes not fitting their phrases are treated as
    // with overflow and handled separately
    //Parameters:
    // text: text to be parsed
    // txt_size: number of symbols in the text
    // buffer_size: number of bytes for the buffer containing the text
    // parse: uin32_t pointer that will point to the resulting parsing
    // n_strings: number of strings in text. It will be filled during the parsing
    // sep_sym: byte symbol in text representing the sentinel separating the strings
    // fp_seed: integer that we use as seed to create fingerprints for the phrases resulted from the parse
    // prev_fps: the fingerprints for the alphabet symbols that we use to compute the parsing breaks
    // p_gram: partial gram that will store the set of phrases resulted from the parsing
    void byte_par_r2l(text_chunk& chunk, off_t& n_strings, size_t sep_sym, uint64_t fp_seed, buff_vector<uint64_t>& prev_fps) {

        uint8_t * text = chunk.text;
        off_t lb, rb = chunk.text_bytes-1, i=chunk.text_bytes-2, max_byte_offset, byte_offset, parse_size;
        uint8_t v_len;
        lz_like_map dict(text);
        assert(text[i+1]==sep_sym && text[i]>text[i+1]);

        text[rb] = 128;//the vbyte code of the metasymbol mt=0 representing a separator in the next round of parsing
        n_strings=1;
        parse_size=1;
        max_byte_offset=4;

        //given a phrase T[a..b-1], byte_offset tells the number of bytes (4 per mt) used by the metasymbols after T[b-1].
        //This value defines the new size for text, now with the parse
        //max_byte_offset tells the maximum offset

        std::vector<phrase_overflow> phr_with_ovf;
        bool r_cmp = true, l_cmp, inserted, new_str;
        uint32_t mt_sym, phrase_len;
        uint8_t mid_sym = text[i];
        while(--i>0 && text[i]==mid_sym);

        while(i>=0){
            l_cmp = prev_fps[text[i]]>prev_fps[mid_sym];
            if(l_cmp && !r_cmp){
                lb = i+1;
                new_str = mid_sym==sep_sym;
                lb += new_str;
                phrase_len=rb-lb;

                byte_offset = rb + (parse_size << 2);
                max_byte_offset = byte_offset > max_byte_offset ? byte_offset : max_byte_offset;

                mt_sym = dict.insert(lb, phrase_len, inserted) + 1;

                v_len = vbyte_len(mt_sym);
                if(__builtin_expect(v_len>phrase_len, 0)){
                    //metasymbol does not fit its phrase
                    phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
                }else if(!inserted){
                    vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
                    memset(&text[lb+v_len], 0, phrase_len-v_len);
                }

                parse_size++;
                rb = i+1;
                if(new_str){
                    assert(text[rb]==sep_sym);
                    text[rb] = 128;//vbyte code for 0 (the separator symbol in the next levels)
                    n_strings++;
                    parse_size++;
                }
            }

            r_cmp = l_cmp;
            mid_sym = text[i];
            while(--i>0 && text[i]==mid_sym);
        }

        lb = 0;
        phrase_len=rb-lb;
        byte_offset = rb + (parse_size << 2);
        max_byte_offset = byte_offset > max_byte_offset ? byte_offset : max_byte_offset;

        mt_sym = dict.insert(lb, phrase_len, inserted)+1;
        v_len = vbyte_len(mt_sym);
        if(__builtin_expect(v_len>phrase_len, 0)){
            phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
        }else if(!inserted){
            vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
            memset(&text[lb+v_len], 0, phrase_len-v_len);
        }
        parse_size++;

        dict.shrink_to_fit();
        dict.destroy_table();
        chunk.parse_size = parse_size;

        //computes how many bytes we require for the parse
        max_byte_offset = INT_CEIL(max_byte_offset, sizeof(uint32_t))*sizeof(uint32_t);
        chunk.increase_capacity(max_byte_offset);
        chunk.update_used_bytes(max_byte_offset);//update the amount of used bytes

        auto avb_addr = chunk.get_free_mem_area();
        //buff_vector<uint32_t> mt_perm(avb_addr.first, avb_addr.second);
        buff_vector<uint32_t> mt_perm;
        mt_perm.resize(dict.size()+1);
        chunk.add_used_bytes((off_t)mt_perm.static_buff_usage());

        create_meta_sym<uint8_t, true>(chunk, mt_perm, fp_seed, dict.phrase_set, prev_fps);
        finish_byte_parse(chunk, dict, max_byte_offset, mt_perm, phr_with_ovf);
        chunk.update_used_bytes(max_byte_offset);

        avb_addr = chunk.get_free_mem_area();
        prev_fps.move_buffer(avb_addr.first, avb_addr.second);
        chunk.add_used_bytes((off_t)prev_fps.static_buff_usage());
        //std::cout<<"main memory usage "<<report_space((off_t)malloc_count_current())<<"/"<<report_space((off_t)malloc_count_peak())<<std::endl;
    }

    // this method parses the text and store the parse in the text itself.
    // It only works for parsing rounds other than the first one because the length of symbol each
    // cell is the same as the length of cell where we store the metasymbols, so there is no overflow
    void int_par_l2r(text_chunk& chunk, off_t& n_strings,
                      uint64_t fp_seed, buff_vector<uint64_t>& prev_fps){

        //std::cout<<"We have "<<report_space((off_t)chunk.av_bytes())<<" for the round"<<std::endl;

        uint32_t *text = chunk.parse;
        uint32_t mt_sym, sep_sym=0, txt_size = chunk.parse_size;
        size_t left_sym, middle_sym, dummy_sym=std::numeric_limits<text_chunk::size_type>::max();
        off_t i=0, parse_size = 0, phrase_len, lb, rb;

        auto avb_addr = chunk.get_free_mem_area();
        lz_like_map<uint32_t> dict(text, avb_addr.first, avb_addr.second);

        bool inserted, new_str;
        n_strings = 0;
        off_t sym_bytes = sizeof(uint32_t);

        lb = 0;
        left_sym = text[i];
        while(++i<txt_size && text[i]==left_sym);
        assert(i<txt_size);

        middle_sym = text[i];
        rb=i;
        while(++i<txt_size && text[i]==middle_sym);

        while(i<txt_size) {
            if(left_sym>middle_sym && middle_sym<text[i]){//local minimum

                phrase_len = rb-lb;
                mt_sym = dict.insert(lb, phrase_len, inserted);

                //we can not replace the first phrase occurrence as we use it as source for the dictionary
                if(!inserted){
                    assert(text[lb]!=dummy_sym);
                    text[lb] = mt_sym+1;
                    memset(&text[lb+1], (int)dummy_sym, sym_bytes*(phrase_len-1));
                }

                new_str = text[rb]==sep_sym;
                parse_size+=1+new_str;
                n_strings+=new_str;

                lb = rb+new_str;
            }

            left_sym = middle_sym;
            middle_sym = text[i];
            rb = i;
            while(++i<txt_size && text[i]==middle_sym);
        }
        assert(rb==(txt_size-1) && text[rb]==sep_sym);

        phrase_len = rb-lb;
        mt_sym = dict.insert(lb, phrase_len, inserted);
        //we can not replace the first phrase occurrence as we use it as source for the dictionary
        if(!inserted){
            assert(text[lb]!=dummy_sym);
            text[lb] = mt_sym+1;//store the metasymbol in the first phrase position
            memset(&text[lb+1], (int)dummy_sym, sym_bytes*(phrase_len-1));//pad the rest of the phrase with dummy symbols
        }
        parse_size+=2;//+1 for the separator symbol
        n_strings++;

        dict.shrink_to_fit();
        dict.destroy_table();

        //std::cout<<"we called malloc? "<<dict.called_malloc()<<" "<<dict.phrases_buff_usage()<<" "<<chunk.boundary()<<" "<<dict.size()<<std::endl;
        assert(dict.phrase_set.size()<dummy_sym);

        chunk.add_used_bytes(off_t(dict.phrase_set.static_buff_usage()));

        avb_addr = chunk.get_free_mem_area();
        buff_vector<uint32_t> mt_perm(avb_addr.first, avb_addr.second);
        mt_perm.resize(dict.size()+1);
        chunk.add_used_bytes((off_t)mt_perm.static_buff_usage());

        create_meta_sym<uint32_t, false>(chunk, mt_perm, fp_seed, dict.phrase_set, prev_fps);

        // create the parse in place
        size_t last_phrase = dict.phrase_set.size()-1;
        mt_sym = 0, lb = 0;
        off_t k=0;
        i=0;

        while(mt_sym<last_phrase) {
            assert(i==lb);
            text[k++] = mt_perm[mt_sym+1];
            i+= dict.phrase_set[mt_sym].len;//move out of the phrase boundary

            mt_sym++;
            lb = dict.phrase_set[mt_sym].source;//position for the next phrase

            while(i<lb){//process the text area between consecutive phrases
                assert(text[i]<mt_perm.size());
                text[k++] = mt_perm[text[i]];
                while(text[++i]==dummy_sym && i<lb);
            }
        }
        //process the rightmost phrase and the parsed text that follows
        assert(mt_sym==last_phrase && i==dict.phrase_set[mt_sym].source);
        text[k++] = mt_perm[mt_sym+1];
        i+= dict.phrase_set[mt_sym].len;
        lb = txt_size;
        while(i<lb){
            assert(text[i]<mt_perm.size());
            text[k++] = mt_perm[text[i]];
            while(text[++i]==dummy_sym && i<lb);
        }

        assert(k==parse_size);
        chunk.parse_size = parse_size;
        chunk.update_used_bytes((off_t)chunk.dist(reinterpret_cast<uint8_t*>(chunk.parse+chunk.parse_size)));

        avb_addr = chunk.get_free_mem_area();
        prev_fps.move_buffer(avb_addr.first, avb_addr.second);
        chunk.add_used_bytes((off_t)prev_fps.static_buff_usage());
        //std::cout<<"\tmain memory usage "<<report_space((off_t)malloc_count_current())<<"/"<<report_space((off_t)malloc_count_peak())<<std::endl;
    }

    template<class sym_type>
    void compress_text_chunk(text_chunk& chunk, std::vector<uint64_t>& fp_seeds){

        //std::cout<<"\nThis is our peak "<<report_space((off_t)malloc_count_peak())<<std::endl;
        size_t alpha_size = std::numeric_limits<sym_type>::max()+1;
        buff_vector<uint64_t> prev_fps;
        prev_fps.resize(alpha_size);

        for(size_t i=0;i<alpha_size;i++){
            prev_fps[i] = XXH64(&i, sizeof(sym_type), fp_seeds[0]);
            assert(prev_fps[i]!=0);
        }
        prev_fps[chunk.sep_sym]=0;

        off_t n_strings=0;
        size_t p_round=0;
        size_t sep_sym = chunk.sep_sym;

        chunk.p_gram.max_tsym = std::numeric_limits<sym_type>::max();
        chunk.p_gram.sep_tsym = chunk.sep_sym;
        chunk.p_gram.text_size = chunk.text_bytes/sizeof(sym_type);
        chunk.p_gram.txt_id = chunk.id;

        //auto start = std::chrono::steady_clock::now();
        byte_par_r2l(chunk, n_strings, sep_sym, fp_seeds[p_round + 1], prev_fps);
        //auto end = std::chrono::steady_clock::now();
        //report_time(start, end , 2);

        off_t size_limit = n_strings*2;
        p_round++;

        while(chunk.parse_size!=size_limit){
            assert(chunk.parse_size>=size_limit);

            //start = std::chrono::steady_clock::now();
            int_par_l2r(chunk, n_strings, fp_seeds[p_round+1], prev_fps);
            //end = std::chrono::steady_clock::now();
            //report_time(start, end , 2);

            p_round++;
        }
        //start = std::chrono::steady_clock::now();
        chunk.p_gram.add_compressed_string(chunk.parse, chunk.parse_size);
        //end = std::chrono::steady_clock::now();
        //report_time(start, end , 2);
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
                text_chunks[chunk_id].increase_capacity((tmp_ck_size*115)/100);
                text_chunks[chunk_id].id = chunk_id;

                read_chunk_from_file(fd_r, rem_bytes, r_acc_bytes, text_chunks[chunk_id]);

                //this value is enough for 4GB text chunks
                text_chunks[chunk_id].p_gram.rules.resize(32);
                text_chunks[chunk_id].p_gram.par_seed = p_opts.orig_seed;

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
                text_chunks[buff_id].p_gram.serialize_to_fd(fd_w);
                //std::cout<<"\r  "<<report_space(text_chunks[buff_id].text_bytes)<<" compressed to "<<report_space((off_t)g_bytes)<<" (speed: "<<report_speed(text_chunks[buff_id].text_bytes, text_chunks[buff_id].t_start, text_chunks[buff_id].t_end)<<", ratio: )"<<std::endl;
#ifdef __linux__
               w_page_cache_bytes += g_bytes;
                   if(w_page_cache_bytes>p_opts.page_cache_limit){
                       std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                       posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
                       w_page_cache_bytes=0;
                   }
#endif
                proc_syms+=text_chunks[buff_id].text_bytes;

                std::cout<<"\n  Parsed input "<<report_space((off_t)proc_syms)<<" with peak "<<report_space((off_t)malloc_count_peak())<<" and grammar size "<<report_space((off_t)text_chunks[buff_id].p_gram.gram_size_in_bytes())<<"/ in int_words:"<<report_space((off_t)text_chunks[buff_id].p_gram.gram_uint32_bytes())<<" / buff_bytes grammar: "<< report_space((off_t)text_chunks[buff_id].p_gram.buff_bytes())<<std::flush;
                malloc_count_reset_peak();

                text_chunks[buff_id].p_gram.reset_grammar();

                //std::cout<<"\r  Processed input "<<report_space((off_t)proc_syms)<<"/"<<report_space(rem_bytes)<<std::endl;

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
            posix_fadvise(fd_r, 0, st.st_size, POSIX_FADV_DONTNEED);
#endif
            close(fd_r);

            //wait for all the parsers to finish
            while(parser_finished.load(std::memory_order_acquire)!=p_opts.n_threads);

            //store the remaining grammars in the temporary file
            while(!buffers_to_reuse.empty()){

                buffers_to_reuse.pop(buff_id);
                text_chunks[buff_id].p_gram.serialize_to_fd(fd_w);
                //std::cout<<report_space(text_chunks[buff_id].text_bytes)<<" compressed to "<<report_space((off_t)g_bytes)<<std::endl;
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
                std::cout<<"\n  Parsed input "<<report_space((off_t)proc_syms)<<" with peak "<<report_space((off_t)malloc_count_peak())<<std::flush;
                malloc_count_reset_peak();
            }
            buffers_to_reuse.done();

#ifdef __linux__
            if(w_page_cache_bytes>0){
                std::cout<<"- removing from page cache "<<w_page_cache_bytes<<" "<<w_acc_bytes<<std::endl;
                posix_fadvise(fd_w, w_acc_bytes-w_page_cache_bytes, w_page_cache_bytes, POSIX_FADV_DONTNEED);
            }
#endif
            close(fd_w);
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
                compress_text_chunk<sym_type>(text_chunks[buff_id], p_opts.p_seeds);
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

            //std::cout<<"Processed data: "<<100<<std::endl;

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

    template<class sym_type, class gram_type>
    void lc_parsing_algo(std::string& i_file, std::string& o_file,
                         tmp_workspace& tmp_ws, size_t n_threads,
                         size_t n_chunks, size_t chunk_size, size_t par_seed, bool par_gram) {

        off_t f_size = file_size(i_file);
        parsing_opts p_opts;
        p_opts.chunk_size = chunk_size==0 ? off_t(ceil(0.025 * double(f_size))) : (off_t)chunk_size;
        p_opts.chunk_size = std::min<off_t>(p_opts.chunk_size, std::numeric_limits<uint32_t>::max());//the chunks cannot exceed the 4GB by design

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
            p_opts.p_seeds[i] = distrib(gen);
        }

        std::cout<<"  Settings"<<std::endl;
        std::cout<<"    Parsing mode              : short strings (<= 4 GBs)"<<std::endl;
        std::cout<<"    Parsing threads           : "<<p_opts.n_threads<<std::endl;
        std::cout<<"    Parsing seed              : "<<p_opts.orig_seed<<std::endl;
        std::cout<<"    Active text chunks in RAM : "<<p_opts.n_chunks<<std::endl;
        std::cout<<"    Size of each chunk        : "<<report_space(p_opts.chunk_size)<<std::endl;
        std::cout<<"    Chunks' approx. mem usage : "<<report_space(off_t(((p_opts.chunk_size*115)/100)*p_opts.n_chunks))<<"\n"<<std::endl;

        std::string ct_p_grams_file = tmp_ws.get_file("concat_p_grams");
        build_partial_grammars<sym_type>(p_opts, i_file, ct_p_grams_file);

        std::string mg_p_gram_file = par_gram? o_file : tmp_ws.get_file("merged_p_grams");
        merge_partial_grammars<sym_type>(ct_p_grams_file, mg_p_gram_file, p_opts.p_seeds, p_opts.n_threads);

        if(!par_gram){
            gram_type final_grammar;
            partial2complete_gram(final_grammar, mg_p_gram_file, par_seed);
            store_to_file(o_file, final_grammar);
        }
    }
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
