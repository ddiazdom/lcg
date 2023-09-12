//
// Created by Diaz, Diego on 26.1.2023.
//

#ifndef GRLBWT_EXACT_PAR_PHASE_H
#define GRLBWT_EXACT_PAR_PHASE_H

#include "common.h"
#include "utils.h"
#include <random>
#include "hashing.h"

struct rand_order{
    size_t str_ptr;
    size_t str_len;
    size_t hash;
    size_t orig_order;
};

struct parsing_info{
    size_t par_phrases=0; //number of parsing phrases in the current round
    size_t par_symbols=0; //number of symbols in the parsing set of the current round
    size_t p_round=0; //parsing round
    std::vector<hashing> p_functions;
    size_t max_symbol=0;
    size_t min_symbol=0;
    size_t longest_str=0; //longest string in the parsing
    std::vector<long> str_ptrs;
};

/*
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
};*/

template<class sym_type,
         bool first_round>
struct lc_parser_t{

    typedef i_file_stream<sym_type>             stream_t;
    //typedef perm_type<false>                    perm_t;
    hashing                                  par_function;

    size_t min_sym;
    size_t max_sym;
    phrase_map_t map;
    //perm_t perm;
    stream_t ifs;
    std::vector<long>& str_boundaries;
    tmp_workspace& ws;

    lc_parser_t(std::string& i_file, parsing_info& p_info, tmp_workspace& ws_) : min_sym(p_info.min_symbol),
                                                                                 max_sym(p_info.max_symbol),
                                                                                 ifs(i_file, BUFFER_SIZE),
                                                                                 str_boundaries(p_info.str_ptrs),
                                                                                 ws(ws_) {
        if(p_info.p_round<p_info.p_functions.size()){
            par_function = p_info.p_functions[p_info.p_round];
        }
    }

    /*lc_parser_t(std::string& i_file, std::vector<long>& str_boundaries_,
                size_t min_sym_, size_t max_sym_, hashing& par_function_,
                tmp_workspace& ws_) : par_function(par_function_),
                                      min_sym(min_sym_),
                                      max_sym(max_sym_),
                                      ifs(i_file, BUFFER_SIZE),
                                      str_boundaries(str_boundaries_),
                                      ws(ws_){}*/

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

        void close(){
            ofs.close();
        }

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

                        if constexpr (first_round){
                            type = (type<<1UL) | (par_function.symbol_hash(curr_sym) < par_function.symbol_hash(prev_sym));
                        }else{
                            type = (type<<1UL) | (curr_sym < prev_sym);
                        }

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

        std::pair<size_t, size_t> parse_res;
        std::string map_file = ws.get_file("map_file");
        {
            hash_functor hf(str_boundaries);
            lc_scan(hf, max_sym);
            parse_res = hf.get_stats(max_sym);
            store_to_file(map_file, hf.in_map);
        }
        load_from_file(map_file, map);

        /*size_t min=std::numeric_limits<size_t>::max(), arg_min=0;
        for(size_t i=0;i<n_it;i++){
            hash_functor hf(str_boundaries);
            lc_scan(hf, max_sym);
            auto res = hf.get_stats(max_sym);
            std::cout<<"try "<<(i+1)<<" : "<<res.first<<" "<<res.second<<std::endl;
            if(res.first<min){
                std::string map_file = ws.get_file("map_"+std::to_string(i));
                store_to_file(map_file, hf.in_map);
                std::string perm_file = ws.get_file("perm_"+std::to_string(i));
                store_to_file(perm_file, perm);
                arg_min = i;
                min = res.first;
                parse_res = res;
            }
            perm.new_perm();
        }
        std::string best_par_set = ws.get_file("map_"+std::to_string(arg_min));
        load_from_file(best_par_set, map);
        std::cout<<(arg_min+1)<<" "<<map.size()<<" "<<std::endl;
        std::string best_perm = ws.get_file("perm_"+std::to_string(arg_min));
        load_from_file(best_perm, perm);*/
        return parse_res;
    }

    phrase_map_t & get_map(){
        return map;
    }

    hashing& get_par_function(){
        return par_function;
    }

    void update_str_positions(size_t p_size){
        //update string pointers
        long acc=0, prev;
        prev = str_boundaries[0];
        for(size_t j=0; j<str_boundaries.size()-1;j++){
            acc += (prev-str_boundaries[j]);
            prev = str_boundaries[j];
            str_boundaries[j] = acc;
        }
        str_boundaries.back() = (long)p_size;
    }

    template<typename o_sym_type>
    size_t produce_next_string_int(std::string& o_file){
        parse_functor<o_sym_type> pf(map, str_boundaries, o_file);
        lc_scan(pf, max_sym);
        size_t p_size = pf.parse_size();
        update_str_positions(p_size);
        pf.close();

        std::string tmp_o_file = ws.get_file("tmp_parse_file");
        i_file_stream<o_sym_type> inv_parse(o_file, BUFFER_SIZE);
        o_file_stream<o_sym_type> final_parse(tmp_o_file, BUFFER_SIZE, std::ios::out);
        for(size_t i=inv_parse.size();i-->0;){
            final_parse.push_back(inv_parse.read(i));
        }
        inv_parse.close(true);
        final_parse.close();
        rename(tmp_o_file.c_str(), o_file.c_str());

        return p_size;
    }

    size_t produce_next_string(std::string& o_file){
        size_t new_max_sym = max_sym+map.size()-1;
        size_t bps = sym_width(new_max_sym);
        if(bps<=8){
            return produce_next_string_int<uint8_t>(o_file);
        }else if(bps<=16){
            return produce_next_string_int<uint16_t>(o_file);
        }else if(bps<=32){
            return produce_next_string_int<uint32_t>(o_file);
        }else{
            return produce_next_string_int<uint64_t>(o_file);
        }
    }
};
//

/*
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
};*/

template<class sym_type>
size_t build_lc_grammar(std::string &i_file, std::string& pf_file, std::string & o_file, size_t n_tries, size_t n_threads, tmp_workspace &ws);

#endif //GRLBWT_EXACT_PAR_PHASE_H
