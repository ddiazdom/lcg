//
// Created by Diaz, Diego on 31.3.2023.
//

#ifndef PARALLEL_PARSING_PARALLEL_PARSER_HPP
#define PARALLEL_PARSING_PARALLEL_PARSER_HPP

#include "text_handler.h"
#include "cds/utils.h"
#include "cds/int_decoders.h"
#include "cds/int_array.h"
#include "cds/rank_support.h"
#include "cds/vbyte_encoding.h"
#include "grammar.h"
#include <functional>

#ifdef __linux__
#include <malloc.h>
#endif

template<class decoder_type>
struct lms_parsing {

    typedef decoder_type decoder_t;
    typedef typename decoder_t::sym_type sym_type;
    typedef text_chunk<decoder_t> text_chunk_t;

    //the scan is performed in the range [ps, end_ps-1] is the leftmost byte out of the range
    template<bool parse,
             class encoder_t=decoder_t,
             class map_type>
    static size_t scan(text_chunk_t& text_chunk, off_t ps, off_t end_ps, map_type& map) {

        off_t len, lb, rb;
        typename decoder_t::sym_type prev_sym, curr_sym, next_sym;
        uint8_t *lb_ptr, *rb_ptr;
        size_t n_phrases=0;

        auto *ptr = &text_chunk.buffer[ps];
        auto const *l_boundary = ptr-1;
        auto const *r_boundary = ptr+(end_ps-ps);

        //leftmost byte of the rightmost valid symbol
        off_t lb_rm_s = end_ps - decoder_t::backward(&text_chunk.buffer[end_ps], l_boundary, 1);

        lb = ps;
        lb_ptr = ptr;

        if constexpr (std::is_same<decoder_t, plain_decoder<sym_type>>::value) {
            decoder_t::read_forward(ptr, prev_sym);
        }else{
            text_chunk.read_forward(ptr, prev_sym);
        }

        assert(ps<=lb_rm_s);

        if constexpr (std::is_same<decoder_t, plain_decoder<sym_type>>::value){
            ps += decoder_t::mov_to_next_diff_sym(ptr, r_boundary, curr_sym);
        }else{
            ps += text_chunk.mov_to_next_diff_sym(ptr, r_boundary, curr_sym);
        }

        rb = ps;
        rb_ptr = ptr;

        while(ps<lb_rm_s) {

            //read the next symbol and the next that differs
            if constexpr (std::is_same<decoder_t, plain_decoder<sym_type>>::value){
                ps += decoder_t::mov_to_next_diff_sym(ptr, r_boundary, next_sym);
            }else{
                ps += text_chunk.mov_to_next_diff_sym(ptr, r_boundary, next_sym);
            }

            if(prev_sym>curr_sym && curr_sym<next_sym){

                len = rb - lb;
                if constexpr (parse){
                    typename map_type::value_type val=0;
                    map.find(lb_ptr, len, val);
                    typename encoder_t::sym_type mt_sym = val;
                    text_chunk.acc_sec_bytes+= encoder_t::write_forward(&text_chunk.sec_buffer[text_chunk.acc_sec_bytes], mt_sym);
                    assert(text_chunk.acc_sec_bytes<=text_chunk.sec_bytes);
                }else{
                    map.value_add_submap(lb_ptr, len, 1);
                }
                n_phrases++;

                lb = rb;
                lb_ptr = rb_ptr;
            }

            prev_sym = curr_sym;
            curr_sym = next_sym;
            rb = ps;
            rb_ptr = ptr;
        }
        assert(ps==rb);
        len = end_ps-lb;
        assert(len>0);

        if constexpr (parse){
            typename map_type::value_type mt_sym;
            map.find(lb_ptr, len, mt_sym);
            typename encoder_t::sym_type sym = mt_sym;
            text_chunk.acc_sec_bytes+= encoder_t::write_forward(&text_chunk.sec_buffer[text_chunk.acc_sec_bytes], sym);
            assert(text_chunk.acc_sec_bytes<=text_chunk.sec_bytes);
        }else{
            map.value_add_submap(lb_ptr, len, 1);
        }
        n_phrases++;
        return n_phrases;
    }

    //returns the rightmost byte of the rightmost symbol within text_chunk
    // that is a local minimium (i.e., a break). It also returns the number
    // of bytes of the break
    static off_t rm_break(text_chunk_t& text_chunk){

        off_t lb = text_chunk.rm_str();//start of the rightmost string in the chunk
        off_t ps = text_chunk.bytes-1;

        auto *ptr = &text_chunk.buffer[ps];
        auto const *l_boundary = (&text_chunk.buffer[std::max<off_t>(lb, 0)])-1;

        //position of ps within the file
        off_t global_ps = text_chunk.n_bytes_before+text_chunk.bytes-1;

        //get the rightmost byte of the rightmost valid symbol within the buffer
        ps -= decoder_t::mov_to_prev_valid_byte(ptr, global_ps);

        //this means that the bytes of the rm string do not contain an entire symbol
        if(ptr==l_boundary) return -1;

        sym_type right_sym, curr_sym, left_sym;

        if constexpr (std::is_same<decoder_t, plain_decoder<sym_type>>::value){
            ps -= decoder_t::read_backwards(ptr, l_boundary, right_sym);//read the rightmost symbol from right to left
            ps -= decoder_t::mov_to_prev_diff_sym(ptr, l_boundary, curr_sym);
        }else{
            ps -= text_chunk.read_backwards(ptr, l_boundary, right_sym);//read the rightmost symbol from right to left
            ps -= text_chunk.mov_to_prev_diff_sym(ptr, l_boundary, curr_sym);
        }

        while(ps>lb && ps>0) {

            if constexpr (std::is_same<decoder_t, plain_decoder<sym_type>>::value){
                ps -= decoder_t::mov_to_prev_diff_sym(ptr, l_boundary, left_sym);
            }else{
                ps -= text_chunk.mov_to_prev_diff_sym(ptr, l_boundary, left_sym);
            }

            if(left_sym>curr_sym && curr_sym<right_sym){
                return ps+decoder_t::forward(ptr, 1)-1;
            }
            right_sym = curr_sym;
            curr_sym = left_sym;
        }
        return -1;
    }

    template<class map_type>
    static size_t get_phrases(text_chunk_t& text_chunk, map_type& map){

        off_t ps, end_ps;
        ps = text_chunk.str_buff_start(-1);//leftmost byte of str
        end_ps = text_chunk.str_buff_start(0);
        size_t n_mt_syms = 0;

        // ps > end_ps means the first symbol in the buffer is the start of str
        if (ps < end_ps){
            //scan the suffix of a string
            n_mt_syms+=scan<false>(text_chunk, ps, end_ps, map);
        }

        for(off_t str=0;str<text_chunk.n_str;str++) {
            ps = text_chunk.str_buff_start(str);//leftmost byte of str
            end_ps = text_chunk.str_buff_start(str + 1);
            n_mt_syms+=scan<false>(text_chunk, ps, end_ps, map);
        }
        return n_mt_syms;
    }

    template<class map_type, class encoder_t>
    static size_t parse_text(text_chunk_t& text_chunk, map_type& map){
        off_t ps, end_ps;
        ps = text_chunk.str_buff_start(-1);//leftmost byte of str
        end_ps = text_chunk.str_buff_start(0);
        size_t n_mt_syms=0;

        // ps > end_ps means the first symbol in the buffer is the start of str
        if (ps < end_ps){
            //scan the suffix of a string
            n_mt_syms+=scan<true, encoder_t>(text_chunk, ps, end_ps, map);
        }

        for(off_t str=0; str<text_chunk.n_str; str++) {
            ps = text_chunk.str_buff_start(str);//leftmost byte of str
            end_ps = text_chunk.str_buff_start(str + 1);

            //update the start of the string
            text_chunk.str_ptr[str] = text_chunk.acc_sec_bytes;
            n_mt_syms+=scan<true, encoder_t>(text_chunk, ps, end_ps, map);
        }
        return n_mt_syms;
    }

    static off_t overlap(const text_chunk_t& chunk, off_t ps){
        return ps+1;
    }

    template<class map_t>
    static void p_set_stats(map_t& map, parsing_opts& p_opts) {
        p_opts.n_sym = 0, p_opts.max_sym=0, p_opts.tot_phrases=0;
        typename decoder_t::sym_type  sym=0;
        off_t len, bytes;
        for(auto const& pair : map){
            len = pair.first.second;
            auto * ptr = (uint8_t*)pair.first.first;
            while(len>0){
                bytes = decoder_t::read_forward(ptr, sym);
                ptr += bytes;
                len -= bytes;
                p_opts.n_sym++;

                if constexpr (std::is_same<decoder_t, vbyte_decoder<sym_type>>::value){
                    sym = p_opts.sym_perm[sym];
                }

                if(sym>p_opts.max_sym){
                    p_opts.max_sym = sym;
                }
            }
            assert(len==0);
            p_opts.tot_phrases+=pair.second;
        }
        p_opts.p_alph_bytes = INT_CEIL(sym_width(map.size()), 8);
    }

    template<class map_t>
    static void p_set_stats_and_vbyte_codes(map_t& map, parsing_opts& p_opts){

        p_opts.n_sym = 0, p_opts.max_sym=0, p_opts.tot_phrases=0;
        typename decoder_t::sym_type  sym=0;

        off_t len, bytes;
        std::vector<std::pair<size_t, size_t>> freqs;
        size_t rank = 0;

        freqs.reserve(map.size());
        for(auto const& pair : map){

            len = pair.first.second;
            auto * ptr = (uint8_t*)pair.first.first;

            while(len>0){

                bytes = decoder_t::read_forward(ptr, sym);
                ptr += bytes;
                len -= bytes;
                p_opts.n_sym++;

                if constexpr (std::is_same<decoder_t, vbyte_decoder<sym_type>>::value){
                    sym = p_opts.sym_perm[sym];
                }

                if(sym>p_opts.max_sym){
                    p_opts.max_sym = sym;
                }
            }
            assert(len==0);
            size_t freq = pair.second;
            p_opts.tot_phrases+=freq;
            freqs.emplace_back(rank++, freq);
        }

        std::sort(freqs.begin(), freqs.end(), [&](auto const& a, auto const& b){
            return a.second>b.second;
        });

        p_opts.vbyte_size=0;
        for(size_t i=0;i<freqs.size();i++){
            p_opts.vbyte_size+= vbyte_len(i)*freqs[i].second;
        }

        p_opts.p_alph_bytes = INT_CEIL(sym_width(map.size()), 8);
        p_opts.vbyte_comp_ratio = float(p_opts.vbyte_size)/float(round_to_power_of_two(p_opts.p_alph_bytes)*p_opts.tot_phrases);
        p_opts.parse_compressible &= p_opts.vbyte_comp_ratio <= p_opts.vbyte_threshold;

        if(p_opts.parse_compressible){
            std::vector<size_t> perm(freqs.size());
            for(size_t i=0;i<perm.size();i++){
                perm[freqs[i].first] = i;
            }
            size_t i=0;
            for(auto pair : map){
                pair.second = perm[i++];
            }
        }else{
            p_opts.vbyte_comp_ratio = 0;
        }
    }

    template<class map_t>
    static void assign_met(map_t &map, parsing_opts& p_opts, lc_gram_buffer_t& gram_buffer, tmp_workspace& ws){

        map.store_tables(ws.get_file("tmp_tables"));
        map.destroy_tables();

        //get statistics about the parse
        p_opts.parse_compressible = map.size()<=p_opts.vbyte_alphabet_threshold;
        if(p_opts.parse_compressible){
            p_set_stats_and_vbyte_codes<map_t>(map, p_opts);
        }else{
            p_opts.vbyte_size=0;
            p_opts.vbyte_comp_ratio = 0;
            p_set_stats(map, p_opts);
        }

        typename decoder_t::sym_type  sym=0;
        size_t len, bytes, k=0;
        int_array<size_t> parsing_set(p_opts.n_sym, sym_width(p_opts.max_sym));
        std::vector<rand_order> r_order(map.size());

        for(auto const& pair : map){
            len = pair.first.second;
            auto * ptr = (uint8_t*)pair.first.first;

            r_order[k].str_ptr = parsing_set.size();
            r_order[k].orig_order = k;
            r_order[k].str_len = 0;

            while(len>0){
                bytes = decoder_t::read_forward(ptr, sym);
                ptr+=bytes;
                len-=bytes;

                if constexpr (std::is_same<decoder_t, vbyte_decoder<sym_type>>::value){
                    sym = p_opts.sym_perm[sym];
                }
                parsing_set.push_back(sym);
                r_order[k].str_len++;
            }
            k++;
        }
        assert(p_opts.n_sym==parsing_set.size());
        map.store_buffers(ws.get_file("tmp_buffers"));
        map.destroy_buffers();

        gram_buffer.create_lc_rules(r_order, parsing_set);

        int_array<size_t> ranks(map.size(), 0, sym_width(map.size()));
        size_t rank=0;
        for(auto & i : r_order){
            ranks[i.orig_order] = rank++;
        }
        std::vector<rand_order>().swap(r_order);

        map.load_buffers(ws.get_file("tmp_buffers"));
        size_t i=0;
        if(p_opts.parse_compressible){
            std::vector<uint32_t> new_sym_perm(map.size(), 0);
            for(auto pair : map){
                new_sym_perm[pair.second] = ranks[i++];
            }
            p_opts.new_sym_perm.swap(new_sym_perm);
        }else{
            for(auto pair : map){
                pair.second = ranks[i++];
            }
        }
        ranks.erase();
#ifdef __linux__
        malloc_trim(0);
#endif
        map.load_tables(ws.get_file("tmp_tables"));
    }
};
#endif //PARALLEL_PARSING_PARALLEL_PARSER_HPP
