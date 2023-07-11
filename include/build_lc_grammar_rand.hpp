//
// Created by Diaz, Diego on 26.1.2023.
//

#ifndef GRLBWT_EXACT_PAR_PHASE_H
#define GRLBWT_EXACT_PAR_PHASE_H

#include "common.h"
#include "utils.h"
#include <random>

struct parsing_info{
    size_t lms_phrases=0; //number of LMS phrases in the current parsing round
    size_t p_round=0; //parsing round
    size_t max_symbol=0;
    size_t min_symbol=0;
    size_t longest_str=0; //longest string in the parsing
    std::vector<long> str_ptrs;
};

//TODO this is for the locally-consistent parsing
template<bool identity=false>
struct perm_type{
    std::vector<size_t> perm;
    size_t min_sym;
    size_t max_sym;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    inline perm_type(size_t min_sym_, size_t max_sym_) : min_sym(min_sym_),
                                                         max_sym(max_sym_){
        size_t len = max_sym-min_sym+1;
        perm.resize(len);
        for(size_t i=0;i<perm.size();i++) perm[i] = i;
        if constexpr(!identity){
            std::shuffle(perm.begin(), perm.end(), std::default_random_engine(seed));
        }
    }

    inline size_t operator()(size_t sym) const {
        assert(min_sym<=sym && sym<=max_sym);
        return min_sym + perm[sym-min_sym];
    }

    inline void new_perm(){
        if constexpr (!identity){
            std::shuffle(perm.begin(), perm.end(), std::default_random_engine(seed));
        }
    }

    inline size_t serialize(std::ostream &out){
        size_t written_bytes = serialize_plain_vector(out, perm);
        written_bytes += serialize_elm(out, min_sym);
        written_bytes += serialize_elm(out, max_sym);
        return  written_bytes;
    }

    void load(std::istream &in){
        load_plain_vector(in, perm);
        load_elm(in, min_sym);
        load_elm(in, max_sym);
    }
};

template<class sym_type,
         bool first_round>
struct lc_parser_t{

    typedef i_file_stream<sym_type>             stream_t;
    typedef perm_type<true>                     perm_t;

    size_t min_sym;
    size_t max_sym;
    phrase_map_t map;
    perm_t perm;
    stream_t ifs;
    std::vector<long>& str_boundaries;
    tmp_workspace& ws;

    lc_parser_t(std::string& i_file, std::vector<long>& str_boundaries_,
                size_t min_sym_, size_t max_sym_,
                tmp_workspace& ws_) : min_sym(min_sym_),
                                      max_sym(max_sym_),
                                      perm(min_sym, max_sym),
                                      ifs(i_file, BUFFER_SIZE),
                                      str_boundaries(str_boundaries_),
                                      ws(ws_) {
    }

    lc_parser_t(std::string& i_file, std::vector<long>& str_boundaries_,
                size_t min_sym_, size_t max_sym_, perm_t& ext_perm,
                tmp_workspace& ws_) : min_sym(min_sym_),
                                      max_sym(max_sym_),
                                      perm(ext_perm),
                                      ifs(i_file, BUFFER_SIZE),
                                      str_boundaries(str_boundaries_),
                                      ws(ws_){
    }

    struct hash_functor{

        phrase_map_t in_map;
        std::vector<long>& str_boundaries;
        size_t first_str;
        size_t last_str;

        explicit hash_functor(std::vector<long>& str_boundaries_): str_boundaries(str_boundaries_),
                                                                   first_str(0),
                                                                   last_str(str_boundaries.size()-2){};

        void process_phrase(string_t& phrase){
            phrase.mask_tail();
            in_map.increment_value(phrase.data(), phrase.n_bits(), 1);
        }

        void forward_comp_str(__attribute__((unused)) size_t& sym){}

        std::pair<long, long> init_str (size_t& str) {
            if constexpr (first_round){
                return {str_boundaries[str], str_boundaries[str+1]-2};
            }else{
                return {str_boundaries[str], str_boundaries[str+1]-1};
            }
        }

        [[nodiscard]] inline size_t get_first_str() const {
            return first_str;
        }

        [[nodiscard]] inline size_t get_last_str() const {
            return last_str;
        }

        inline void finish_scan(){}

        std::pair<size_t, size_t> get_stats(size_t& m_sym){
            key_wrapper key_w{sym_width(m_sym), in_map.description_bits(), in_map.get_data()};
            size_t tot_syms=0;
            for (auto const &ptr: in_map) {
                tot_syms+=key_w.size(ptr);
            }
            return {tot_syms, in_map.size()};
        }

        void store_to_file(std::string& map_file){
            std::ofstream ofs(map_file, std::ios::binary);
            in_map.serialize(ofs);
        }
    };

    template<class o_sym_type>
    struct parse_functor{

        typedef o_file_stream<o_sym_type> o_stream_t;
        phrase_map_t& map;
        std::vector<long>& str_boundaries;
        o_stream_t ofs;
        size_t first_str;
        size_t last_str;

        parse_functor(phrase_map_t& map_, std::vector<long>& str_boundaries_,
                      std::string& o_file): map(map_),
                                            str_boundaries(str_boundaries_),
                                            ofs(o_file, BUFFER_SIZE, std::ios::out),
                                            first_str(0),
                                            last_str(str_boundaries.size()-2){};

        void process_phrase(string_t& phrase) {
            phrase.mask_tail();
            size_t sym = 0;
            auto res = map.key2value(phrase.data(), phrase.n_bits(), sym);

            //TODO testing
            /*for(size_t i=0;i<phrase.size();i++){
                std::cout<<phrase[i]<<" ";
            }
            std::cout<<""<<std::endl;*/
            //

            assert(res);
            ofs.push_back(sym);
        };

        std::pair<long, long> init_str(size_t& str)  {
            std::pair<long, long> range;
            if constexpr (first_round){
                range = {str_boundaries[str], str_boundaries[str+1]-2};
            }else{
                range = {str_boundaries[str], str_boundaries[str+1]-1};
            }

            if((str+1)<=last_str){
                str_boundaries[str+1] = ofs.size()-1;
            }

            return range;
        };

        inline void finish_scan(){
            str_boundaries[first_str] = ofs.size()-1;
        }

        [[nodiscard]] inline size_t get_first_str() const {
            return first_str;
        }

        [[nodiscard]] inline size_t get_last_str() const {
            return last_str;
        }

        void forward_comp_str(size_t& sym) {
            ofs.push_back(sym);
        };

        [[nodiscard]] inline size_t parse_size() const {
            return ofs.size();
        }

        ~parse_functor(){
            ofs.close(false);
        }
    };

    template<class functor>
    inline void lc_scan(functor& func, size_t max_symbol) {

        sym_type curr_sym, prev_sym;
        string_t phrase(2, sym_width(max_symbol));
        size_t end_ps, start_ps;
        uint8_t type;

        size_t first_str=func.get_first_str();
        size_t end_str=func.get_last_str();

        for(size_t str=end_str+1;str-->first_str;) {

            auto range = func.init_str(str);
            assert(range.first<=range.second);

            if(range.first<range.second) { //if this is not true, it means the string was fully compressed
                start_ps = range.first;
                end_ps = range.second;

                prev_sym = ifs.read(end_ps);

                phrase.push_back(prev_sym);
                type = 0;

                for (size_t i = end_ps; i-- > start_ps;) {

                    curr_sym = ifs.read(i);
                    if (curr_sym != prev_sym) {
                        type = (type<<1UL) | (perm(curr_sym) < perm(prev_sym));
                        if ((type & 3U) == 2) {//LMS suffix
                            //process the previous phrase
                            assert(!phrase.empty());
                            func.process_phrase(phrase);

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
                func.process_phrase(phrase);
                phrase.clear();
            } else {
                assert(range.first==range.second);
                size_t tmp = ifs.read(range.first);
                func.forward_comp_str(tmp);
            }
        }
        func.finish_scan();
    };

    std::pair<size_t, size_t> partition_text(size_t n_it){
        std::vector<std::pair<size_t, size_t>> parse_size;
        parse_size.reserve(n_it);
        for(size_t i=0;i<n_it;i++){
            hash_functor hf(str_boundaries);
            lc_scan(hf, max_sym);
            auto res = hf.get_stats(max_sym);
            //std::cout<<"it "<<(i+1)<<" : "<<res.first<<" "<<res.second<<std::endl;
            parse_size.push_back(res);
            std::string map_file = ws.get_file("map_"+std::to_string(i));
            store_to_file(map_file, hf.in_map);
            std::string perm_file = ws.get_file("perm_"+std::to_string(i));
            store_to_file(perm_file, perm);
            perm.new_perm();
        }

        //get the iteration with the best performance
        size_t arg_min=0;
        size_t min=std::numeric_limits<size_t>::max();
        for(size_t i=0;i<n_it;i++){
            if(parse_size[i].second<min){
                min = parse_size[i].second;
                arg_min = i;
            }
        }
        //

        std::string best_par_set = ws.get_file("map_"+std::to_string(arg_min));
        load_from_file(best_par_set, map);
        //std::cout<<arg_min<<" "<<map.size()<<" "<<std::endl;
        std::string best_perm = ws.get_file("perm_"+std::to_string(arg_min));
        load_from_file(best_perm, perm);

        return parse_size[arg_min];
    }

    phrase_map_t & get_map(){
        return map;
    }

    size_t produce_next_string(std::string& o_file){
        size_t new_max_sym = max_sym+map.size()-1;
        size_t bps = sym_width(new_max_sym);
        size_t p_size;
        if(bps<=8){
            parse_functor<uint8_t> pf(map, str_boundaries, o_file);
            lc_scan(pf, max_sym);
            p_size = pf.parse_size();
        }else if(bps<=16){
            parse_functor<uint16_t> pf(map, str_boundaries, o_file);
            lc_scan(pf, max_sym);
            p_size = pf.parse_size();
        }else if(bps<=32){
            parse_functor<uint32_t> pf(map, str_boundaries, o_file);
            lc_scan(pf, max_sym);
            p_size = pf.parse_size();
        }else{
            parse_functor<uint64_t> pf(map, str_boundaries, o_file);
            lc_scan(pf, max_sym);
            p_size = pf.parse_size();
        }
        return p_size;
    }
};
//

template<typename parse_data_t,
         typename parser_t>
struct hash_functor{

    void operator()(parse_data_t& data) {

        auto hash_phrase = [&](string_t& phrase) -> void {
            phrase.mask_tail();
            data.inner_map.increment_value(phrase.data(), phrase.n_bits(), 1);
        };

        auto init_str = [&](size_t& str) -> std::pair<long, long>{
            auto range = data.str2range(str);
            return range;
        };

        auto forward_comp_str = [&](__attribute__((unused)) size_t& sym) -> void{ };

        parser_t()(data.ifs, data.start_str, data.end_str, data.max_symbol, hash_phrase, init_str, forward_comp_str);
    };
};

template<class parse_data_t,
         class parser_t,
         class o_stream_type>
struct parse_functor{

    size_t operator()(parse_data_t& data) {

        o_stream_type ofs(data.o_file, BUFFER_SIZE, std::ios::out);

        auto phrase2symbol = [&](string_t& phrase) -> void {
            phrase.mask_tail();
            size_t sym = 0;
            auto res = data.map.key2value(phrase.data(), phrase.n_bits(), sym);
            assert(res);
            ofs.push_back(sym);
        };

        auto init_str = [&](size_t& str) -> std::pair<long, long>{
            auto range = data.str2range(str);
            if((str+1)<=data.end_str){
                data.str_ptr[str+1] = ofs.size()-1;
            }
            return range;
        };

        auto forward_comp_str = [&](size_t& sym) -> void{
            ofs.push_back(sym);
        };

        parser_t()(data.ifs, data.start_str, data.end_str, data.max_symbol, phrase2symbol, init_str, forward_comp_str);

        data.str_ptr[data.start_str] = ofs.size()-1;

        data.ifs.close();
        ofs.close();
        return ofs.size();
    };
};

template<class sym_type>
size_t build_lc_grammar(std::string &i_file, std::string & o_file, size_t n_threads, tmp_workspace &ws);

template<class parse_strategy>
size_t par_round(parse_strategy &p_strat, parsing_info &p_info, tmp_workspace &ws);

/***
 *
 * @param i_file : input text file
 * @param n_threads : number of working threads
 * @param hbuff_size : buffer size for the hashing step
 */
template<class sym_type>
void gram_algo(std::string &i_file, std::string& o_file, tmp_workspace & tmp_ws, size_t n_threads){
    build_lc_grammar<sym_type>(i_file, o_file, n_threads, tmp_ws);
    std::cout<<"The resulting grammar was stored in "<<o_file<<std::endl;
}
#endif //GRLBWT_EXACT_PAR_PHASE_H
