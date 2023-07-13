//
// Created by Diaz, Diego on 15.12.2022.
//

#ifndef GRLBWT_PARSING_STRATEGIES_H
#define GRLBWT_PARSING_STRATEGIES_H

#include <thread>
#include "common.h"

template<class stream_t, bool first_round>
struct lms_parsing{

    typedef stream_t                       stream_type;
    typedef typename stream_type::sym_type sym_type;

    inline void operator()(stream_t& ifs,
                           size_t f_string, size_t l_string, size_t max_symbol,
                           std::function<void(string_t&)>&& process_phrase,
                           std::function<std::pair<long, long>(size_t&)>&& init_str,
                           std::function<void(size_t& sym)>&& forward_comp_str) const {

        sym_type curr_sym, prev_sym;
        string_t phrase(2, sym_width(max_symbol));
        size_t end_ps, start_ps;
        uint8_t type;

        for(size_t str=l_string+1;str-->f_string;) {

            auto range = init_str(str);
            assert(range.first<=range.second);

            if(range.first<range.second) { //if this is not true, it means the string was fully compressed
                start_ps = range.first;
                end_ps = range.second;

                if(first_round){
                    //assert(ifs.read(end_ps)=='\n');
                    end_ps--;
                }

                prev_sym = ifs.read(end_ps);

                phrase.push_back(prev_sym);
                type = 0;

                for (size_t i = end_ps; i-- > start_ps;) {

                    curr_sym = ifs.read(i);

                    //if constexpr (!first_round){
                    //    rep = (rep << 1UL) | (curr_sym & 1UL);
                    //    curr_sym >>=1UL;
                    //}

                    if (curr_sym != prev_sym) {
                        type = (type<<1UL) | (curr_sym < prev_sym);
                        //if ((type & 3U) == 2 && (rep & 3U)==3U) {//LMS suffix
                        if ((type & 3U) == 2) {//LMS suffix
                            //process the previous phrase
                            assert(!phrase.empty());
                            process_phrase(phrase);

                            //create the new phrase
                            phrase.clear();
                        }
                    } else {
                        type = (type<<1UL) | (type & 1UL);
                    }

                    phrase.push_back(curr_sym);
                    prev_sym = curr_sym;
                }

                assert(!phrase.empty());
                process_phrase(phrase);
                phrase.clear();
            } else {
                assert(range.first==range.second);
                size_t tmp = ifs.read(range.first);
                forward_comp_str(tmp);
            }
        }
    }
};

template<class parser_type,
         template<class, class> class hash_functor,
         template<class, class, class> class parse_functor>
struct mt_parse_strat_t {//multi thread strategy

    typedef typename parser_type::stream_type istream_t;

    struct thread_worker_data_t{

        size_t               start_str{};
        size_t               end_str{};
        istream_t            ifs;
        std::string          o_file;
        phrase_map_t&        map;
        buffered_map_t       inner_map;
        std::vector<long>&   str_ptr;
        size_t               max_symbol{};
        size_t               rb{};

        thread_worker_data_t(size_t start_, size_t end_, std::string& i_file, std::string& o_file_,
                             phrase_map_t& map_, std::vector<long>& str_ptr_,
                             const size_t &hb_size, void *hb_addr,
                             size_t max_symbol_): start_str(start_),
                                                  end_str(end_),
                                                  ifs(i_file, BUFFER_SIZE),
                                                  o_file(o_file_),
                                                  map(map_),
                                                  inner_map(hb_size, o_file + "_phrases", 0.7, hb_addr, map.description_bits()),
                                                  str_ptr(str_ptr_),
                                                  max_symbol(max_symbol_),
                                                  rb(str_ptr[end_str+1]-1){};
        thread_worker_data_t()=default;

        inline std::pair<long, long> str2range(size_t str) const {
            assert(start_str<=str && str<=end_str);
            size_t bg = str_ptr[str];
            size_t end = str ==end_str ? rb : str_ptr[str+1]-1;
            assert(bg<=end);
            return {bg, end};
        }
    };

    std::string                       i_file;
    std::string                       o_file;
    phrase_map_t                      map;
    parsing_info&                     p_info;
    std::vector<thread_worker_data_t> threads_data;
    char *                            buff_addr=nullptr;

    mt_parse_strat_t(std::string &i_file_, std::string& o_file_,
                     parsing_info& p_info_, size_t hbuff_size,
                     size_t n_threads) : i_file(i_file_),
                                         o_file(o_file_),
                                         map(0.8, sym_width(INT_CEIL(p_info_.longest_str*sym_width(p_info_.max_symbol),8)*8)),
                                         p_info(p_info_) {

        std::vector<std::pair<size_t, size_t>> thread_ranges;
        size_t str_per_thread = INT_CEIL(p_info.str_ptrs.size()-1, n_threads);
        n_threads = INT_CEIL((p_info.str_ptrs.size()-1), str_per_thread);

        for(size_t i=0;i<n_threads;i++){
            thread_ranges.emplace_back(str_per_thread*i, std::min(str_per_thread*(i+1)-1, p_info.str_ptrs.size()-2));
        }
        threads_data.reserve(thread_ranges.size());

        // each thread has a hast table with a buffer of at least 8MB
        hbuff_size = std::max<size_t>(hbuff_size, BUFFER_SIZE*n_threads);

        // each thread has a buffer with the same size, and with an integral number of size_t words
        size_t hb_bytes = INT_CEIL(INT_CEIL(hbuff_size, n_threads), sizeof(size_t)) * sizeof(size_t);
        hbuff_size = hb_bytes*n_threads;

        buff_addr = (char *)malloc(hbuff_size);
        size_t k = 0;
        for (auto &range: thread_ranges) {
            std::string tmp_o_file = o_file.substr(0, o_file.size() - 5);
            tmp_o_file.append("_range_" + std::to_string(range.first) + "_" + std::to_string(range.second));
            threads_data.emplace_back(range.first, range.second, i_file, tmp_o_file, map, p_info.str_ptrs, hb_bytes,
                                      buff_addr+(hb_bytes*k), p_info.max_symbol);
            k++;
        }
        assert(!threads_data.empty());
    };

    mt_parse_strat_t()=default;

    std::pair<size_t, size_t> get_phrases() {

        std::vector<std::thread> threads(threads_data.size());
        hash_functor<thread_worker_data_t, parser_type> hf;

        for (size_t i = 0; i < threads_data.size(); i++) {
            threads[i] = std::thread(hf, std::ref(threads_data[i]));
        }

        for (size_t i = 0; i < threads_data.size(); i++) {
            threads[i].join();
        }

        for (size_t i = 0; i < threads_data.size(); i++) {
            threads_data[i].inner_map.flush();
        }

        free(buff_addr);
#ifdef __linux__
        malloc_trim(0);
#endif
        auto join_res = join_thread_phrases();
        return join_res;
    }

    std::pair<size_t, size_t> join_thread_phrases() {

        size_t dic_bits=0, freq, max_freq=0;
        std::string file;

        if(p_info.p_round>0){
            map.resize_table(prev_power_of_two(p_info.lms_phrases));
        }

        for(auto const& thread : threads_data) {

            file = thread.inner_map.dump_file();
            size_t tot_bytes = std::filesystem::file_size(file);
            if(tot_bytes==0) continue;

            size_t d_bits = thread.inner_map.description_bits();
            size_t value_bits = thread.inner_map.value_bits();
            size_t longest_key = thread.inner_map.longest_key();//in bits

            size_t buffer_size = std::max<size_t>(INT_CEIL(longest_key, 8), std::min<size_t>(tot_bytes, BUFFER_SIZE));
            buffer_size = next_power_of_two(buffer_size);
            i_file_stream<size_t> data_disk_buffer(file, buffer_size);

            //there is a region of the file that does not contain data
            // * * * | * 0 0 0 0 0 0 0 <- tot_bits
            //             | <- if next_bit falls in this region, the loop is still valid, but there is no more data
            size_t tot_bits = (tot_bytes*8)-8;
            size_t key_bits;
            size_t next_bit = 0;
            auto * key= (uint8_t*) malloc(INT_CEIL(longest_key, bitstream<ht_buff_t>::word_bits)*sizeof(ht_buff_t));

            while(next_bit<tot_bits) {

                key_bits = data_disk_buffer.read_bits(next_bit, next_bit+d_bits-1);

                assert(key_bits>0 && key_bits<=longest_key);

                key[INT_CEIL(key_bits, 8)-1] = 0;
                next_bit+=d_bits;
                data_disk_buffer.read_bit_chunk(key, next_bit, next_bit+key_bits-1);

                next_bit+=key_bits;
                freq = data_disk_buffer.read_bits(next_bit, next_bit+value_bits-1);
                next_bit+=value_bits;

                auto res = map.increment_value(key, key_bits, freq);
                if(res==freq) dic_bits+=key_bits;
                freq = res;
                if(freq>max_freq) max_freq = freq;
            }

            data_disk_buffer.close(true);
            free(key);
        }
        map.shrink_databuff();
        return {dic_bits/ sym_width(p_info.lms_phrases), max_freq};
    }

    template<class o_sym_type>
    size_t parse_text() {

        std::vector<std::thread> threads(threads_data.size());
        parse_functor<thread_worker_data_t, parser_type, o_file_stream<o_sym_type>> pf;

        for(size_t i=0;i<threads_data.size();i++){
            threads[i] = std::thread(pf, std::ref(threads_data[i]));
        }
        for(size_t i=0;i<threads_data.size();i++) {
            threads[i].join();
        }

        //invert the chunks
        o_file_stream<o_sym_type> of(o_file, BUFFER_SIZE, std::ios::out);
        for(auto const& thread: threads_data){
            i_file_stream<o_sym_type> inv_chunk(thread.o_file, BUFFER_SIZE);
            for(size_t i = inv_chunk.size();i-->0;){
                of.push_back(inv_chunk.read(i));
            }
            inv_chunk.close(true);
        }
        of.close();
        size_t psize = of.size();

        //join the phrases in one single file
        //size_t psize = join_parse_chunks();

        //update string pointers
        long acc=0, prev, str_len=0;
        for(size_t i=0;i<threads_data.size();i++){
            prev = p_info.str_ptrs[threads_data[i].start_str];
            for(size_t j=threads_data[i].start_str; j<=threads_data[i].end_str;j++){
                acc += (prev-p_info.str_ptrs[j]);
                prev = p_info.str_ptrs[j];

                p_info.str_ptrs[j] = acc;
                if(j>0 && (p_info.str_ptrs[j]-p_info.str_ptrs[j-1])>str_len){
                    str_len = p_info.str_ptrs[j]-p_info.str_ptrs[j-1];
                }
            }
            acc +=prev+1;
        }

        p_info.str_ptrs.back() = (long)psize;
        if((p_info.str_ptrs.back() - p_info.str_ptrs[p_info.str_ptrs.size()-2])>str_len){
            str_len = (p_info.str_ptrs.back() - p_info.str_ptrs[p_info.str_ptrs.size()-2]);
        }
        p_info.longest_str = str_len;
        return psize;
    }
};

template<class parser_type,
         template<class, class> class hash_functor,
         template<class, class, class> class parse_functor>
struct st_parse_strat_t {//parse data for single thread

    typedef parser_type                    parser_t;
    typedef typename parser_t::stream_type istream_t;

    istream_t           ifs;
    std::string         o_file;
    std::string         tmp_o_file;

    phrase_map_t&       inner_map;
    parsing_info&       p_info;
    std::vector<long>&  str_ptr;
    size_t              max_symbol{};
    size_t              start_str;
    size_t              end_str;
    phrase_map_t        map;

    st_parse_strat_t(std::string &i_file_, std::string& o_file_,
                     parsing_info& p_info_): ifs(i_file_, BUFFER_SIZE),
                                             o_file(o_file_),
                                             inner_map(map),
                                             p_info(p_info_),
                                             str_ptr(p_info.str_ptrs),
                                             max_symbol(p_info.max_symbol),
                                             start_str(0),
                                             end_str(str_ptr.size()-2),
                                             map(0.8, sym_width(INT_CEIL(p_info_.longest_str*sym_width(max_symbol),8)*8)){
        tmp_o_file = o_file.substr(0, o_file.size() - 4);
        tmp_o_file.append("_inv");

        if(p_info.p_round>1){
            size_t n_buckets = prev_power_of_two(p_info.lms_phrases);
            map.resize_table(n_buckets);
        }
    }

    [[nodiscard]] inline std::pair<long, long> str2range(size_t str) const {
        assert(start_str<=str && str<=end_str);
        size_t bg = str_ptr[str];
        size_t end = str_ptr[str+1]-1;
        assert(bg<=end);
        return {bg, end};
    }

    std::pair<size_t, size_t> get_phrases() {

        hash_functor<st_parse_strat_t, parser_type>()(*this);

        map.shrink_databuff();
        key_wrapper key_w{sym_width(max_symbol), map.description_bits(), map.get_data()};
        size_t n_syms=0, max_freq=0, freq;

        for(auto const &ptr : map){
            n_syms += key_w.size(ptr);
            freq = 0;
            map.get_value_from(ptr, freq);
            assert(freq>0);
            if(freq>max_freq) max_freq = freq;
        }
        return {n_syms, max_freq};
    }

    template<class o_sym_type>
    size_t parse_text() {

        size_t parse_size = parse_functor<st_parse_strat_t, parser_type, o_file_stream<o_sym_type>>()(*this);

        //update string pointers
        long acc=0, str_len=0;
        size_t prev = p_info.str_ptrs[start_str];
        for(size_t j=start_str; j<=end_str; j++){
            acc += (prev-p_info.str_ptrs[j]);
            prev = p_info.str_ptrs[j];
            p_info.str_ptrs[j] = acc;
            if(j>0 && (p_info.str_ptrs[j]-p_info.str_ptrs[j-1])>str_len){
                str_len = p_info.str_ptrs[j]-p_info.str_ptrs[j-1];
            }
        }
        p_info.str_ptrs.back() = (long)parse_size;
        if((p_info.str_ptrs.back() - p_info.str_ptrs[p_info.str_ptrs.size()-2])>str_len){
            str_len = (p_info.str_ptrs.back() - p_info.str_ptrs[p_info.str_ptrs.size()-2]);
        }
        p_info.longest_str = str_len;

        i_file_stream<o_sym_type> inv_parse(o_file, BUFFER_SIZE);
        o_file_stream<o_sym_type> final_parse(tmp_o_file, BUFFER_SIZE, std::ios::out);
        for(size_t i=inv_parse.size();i-->0;){
            final_parse.push_back(inv_parse.read(i));
        }
        inv_parse.close(true);
        final_parse.close();
        rename(tmp_o_file.c_str(), o_file.c_str());

        return parse_size;
    }
};

typedef i_file_stream<uint8_t>                        uint8t_i_stream;
typedef i_file_stream<uint16_t>                       uint16t_i_stream;
typedef i_file_stream<uint32_t>                       uint32t_i_stream;
typedef i_file_stream<uint64_t>                       uint64t_i_stream;

typedef lms_parsing<uint8t_i_stream, false>   uint8t_parser_t;
typedef lms_parsing<uint16t_i_stream, false>  uint16t_parser_t;
typedef lms_parsing<uint32t_i_stream, false>  uint32t_parser_t;
typedef lms_parsing<uint64t_i_stream, false>  uint64t_parser_t;

#endif //GRLBWT_PARSING_STRATEGIES_H
