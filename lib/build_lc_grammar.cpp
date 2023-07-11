//
// Created by Diaz, Diego on 26.1.2023.
//

#include "malloc_count.h"
#include "build_lc_grammar.hpp"
#include "parsing_strategies.h"
#include "grammar.h"

/*
    template<class vector_type, class sa_type>
    size_t process_dictionary_int(dictionary &dict, parsing_info &p_info, tmp_workspace &ws) {

        sa_type sa;
        if constexpr (std::is_same<sa_type, vector_t>::value) {
            size_t width = sym_width(dict.dict.size())+2;
            sa.set_width(width);
            sa.resize(dict.dict.size());
            sa.initialize(0, sa.size());
        }else{
            sa.resize(dict.dict.size());
            memset(sa.data(), 0, sa.size()*sizeof(typename sa_type::value_type));
        }

        std::cout <<"    Sorting the dictionary and constructing the preliminary BWT" << std::flush;
        auto start = std::chrono::steady_clock::now();
        suffix_induction<vector_type, sa_type>(dict, sa);

        phrase_map_t new_phrases_ht;
        bv_t phr_marks(dict.dict.size(), false);
        size_t n_phrases = produce_pre_bwt<sa_type>(dict, sa, new_phrases_ht, phr_marks, p_info, ws);
        auto end = std::chrono::steady_clock::now();
        report_time(start, end, 2);

        //malloc_count_print_status();
        //malloc_count_reset_peak();

        std::cout << "    Compressing the dictionary" << std::flush;
        start = std::chrono::steady_clock::now();
        produce_grammar<sa_type>(dict, sa, new_phrases_ht, phr_marks, p_info, ws);
        end = std::chrono::steady_clock::now();
        report_time(start, end, 35);

#ifdef __linux__
        malloc_trim(0);
#endif
        return n_phrases;
    }

    size_t process_dictionary(dictionary &dict, parsing_info &p_info, tmp_workspace &ws) {

        uint8_t width = sym_width(dict.dict.size()) + 1;
        size_t n_phrases;

        if (width <= 8) {
            using uint8_vector_t = std::vector<uint8_t, mallocator<uint8_t>>;
            n_phrases = process_dictionary_int<uint8_vector_t, uint8_vector_t>(dict, p_info, ws);
        } else if (width <= 16) {
            using uint16_vector_t = std::vector<uint16_t, mallocator<uint16_t>>;
            n_phrases = process_dictionary_int<uint16_vector_t, uint16_vector_t>(dict, p_info, ws);
        } else if (width <= 32) {
            using uint32_vector_t = std::vector<uint32_t, mallocator<uint32_t>>;
            n_phrases = process_dictionary_int<uint32_vector_t, uint32_vector_t>(dict, p_info, ws);
        } else {
            using uint64_vector_t = std::vector<uint64_t, mallocator<uint64_t>>;
            n_phrases = process_dictionary_int<uint64_vector_t, vector_t>(dict, p_info, ws);
        }
        return n_phrases;
    }
*/

void check_plain_grammar(lc_gram_t& gram, std::string& uncomp_file) {

    std::cout<<"Checking the grammar produces the exact input string"<<std::endl;
    std::cout<<"  This step is optional and for debugging purposes"<<std::endl;
    std::cout<<"  Terminals:              "<<gram.n_terminals()<<std::endl;
    std::cout<<"  Number of nonterminals: "<<gram.n_nonterminals()<<std::endl;
    std::cout<<"  Compressed string:      "<<gram.comp_str_size()<<std::endl;


    i_file_stream<uint8_t> if_stream(uncomp_file, BUFFER_SIZE);

    size_t start_symbol = gram.start_symbol();
    auto res = gram.nt2phrase(start_symbol);

    std::stack<size_t> stack;

    size_t f = res.first;
    size_t l = res.second;
    size_t idx=0;

    std::string decompression;
    size_t str=0;

    for(size_t i=f; i <= l; i++) {

        stack.emplace(gram.pos2symbol(i));
        assert(stack.size()<=if_stream.size());

        while(!stack.empty()){

            auto curr_sym = stack.top() ;
            stack.pop();

            if(gram.is_terminal(curr_sym)){
                decompression.push_back((char)curr_sym);
            }else{
                auto res2 = gram.nt2phrase(curr_sym);
                if(gram.is_rl_sym(curr_sym)){
                    assert(res2.second-res2.first+1==2);
                    size_t len = gram.pos2symbol(res2.second);
                    for(size_t j=0;j<len;j++){
                        stack.emplace(gram.pos2symbol(res2.first));
                    }
                }else{
                    for(size_t j=res2.second+1; j-->res2.first;){
                        stack.emplace(gram.pos2symbol(j));
                    }
                }
            }
        }

        for(char sym : decompression){
            if(sym!=(char)if_stream.read(idx)){
                std::cout<<(int)sym<<" "<<if_stream.read(idx)<<" "<<str<<" "<<gram.str_boundaries.size()-1<<std::endl;
            }
            assert(sym==(char)if_stream.read(idx));
            idx++;
        }
        if(gram.str_boundaries[str+1]==(i+1)){
            idx++;
            str++;
        }
        decompression.clear();
    }
    std::cout<<"\tGrammar is correct!!"<<std::endl;
}

void run_length_compress(lc_gram_buffer_t &gram_buff) {

    std::cout<<"Run-length compressing the grammar"<<std::endl;

    i_file_stream<size_t> rules(gram_buff.rules_file, BUFFER_SIZE);

    std::vector<size_t> rl_rules;
    rl_rules.reserve(gram_buff.g);

    size_t new_id = gram_buff.r;
    size_t tmp_sym;
    phrase_map_t ht;
    string_t pair(2,sdsl::bits::hi(rules.size())+1);

    size_t i=gram_buff.max_tsym+2;
    size_t run_len=1;
    bool prev_flag = rules.read(i-1) & 1UL;
    size_t prev_sym = rules.read(i-1) >> 1UL;
    size_t sym;
    bool flag;
    assert(prev_flag);

    //put the terminals
    for(size_t j=0;j<=gram_buff.max_tsym;j++){
        rl_rules.push_back(rules.read(j));
        assert(rules.read(j) & 1UL);
    }

    while(i<rules.size()) {

        sym = rules.read(i)>>1UL;
        flag =  rules.read(i) & 1UL;

        if(sym!=prev_sym || flag){

            if(run_len>1){
                pair.write(0, prev_sym);
                pair.write(1, run_len);
                auto res = ht.insert(pair.data(), pair.n_bits(), 0);

                if(res.second){
                    tmp_sym = new_id++;
                    ht.insert_value_at(res.first, tmp_sym);
                }else{
                    tmp_sym = 0;
                    ht.get_value_from(res.first, tmp_sym);
                }
            } else{
                tmp_sym = prev_sym;
            }
            //n_rules+=prev_flag;
            tmp_sym = (tmp_sym<<1UL | prev_flag);
            rl_rules.push_back(tmp_sym);

            prev_flag = flag;
            prev_sym = sym;
            run_len=0;
        }
        run_len++;
        i++;
    }

    if(run_len>1){
        pair.write(0, prev_sym);
        pair.write(1, run_len);
        auto res = ht.insert(pair.data(), pair.n_bits(), 0);
        if(res.second){
            tmp_sym = new_id;
            ht.insert_value_at(res.first, tmp_sym);
        }else{
            tmp_sym = 0;
            ht.get_value_from(res.first, tmp_sym);
        }
    } else{
        tmp_sym = prev_sym;
    }
    tmp_sym = (tmp_sym<<1UL | prev_flag);
    rl_rules.push_back(tmp_sym);

    //
    const bitstream<phrase_map_t::buff_t>& stream = ht.get_data();
    key_wrapper key_w{pair.width(), ht.description_bits(), stream};
    for(auto const& phrase : ht){
        rl_rules.push_back((key_w.read(phrase, 0)<<1UL) | 1UL);
        rl_rules.push_back(key_w.read(phrase, 1)<<1UL);
    }

    gram_buff.rl_rules.first = gram_buff.r;
    gram_buff.rl_rules.second = ht.size();
    gram_buff.r += ht.size();
    gram_buff.g  = rl_rules.size();

    std::cout<<"    Stats:"<<std::endl;
    std::cout<<"      Grammar size before:        "<<rules.size()<<std::endl;
    std::cout<<"      Grammar size after:         "<<rl_rules.size()<<std::endl;
    std::cout<<"      Number of new nonterminals: "<<ht.size()<<std::endl;
    std::cout<<"      Compression ratio:          "<<float(rl_rules.size())/float(rules.size())<<std::endl;

    basic_store_vector_to_file(gram_buff.rules_file, rl_rules);
    rules.close();
}

void create_lc_rules(lc_gram_buffer_t& gram, phrase_map_t& map, size_t txt_alphabet){
    key_wrapper key_w{sym_width(txt_alphabet), map.description_bits(), map.get_data()};
    size_t sym, tot_syms=0;
    bool first;
    for (auto const &ptr: map) {
        first = true;
        for (size_t i = key_w.size(ptr); i-- > 0;) {
            sym = key_w.read(ptr, i);
            assert(sym>0);
            sym = (sym<<1UL) | first;
            gram.rules_buffer.push_back(sym);
            first = false;
        }
        tot_syms+=key_w.size(ptr);
    }
    gram.r += map.size();
    gram.g += tot_syms;
    gram.lvl_rules.push_back(map.size());
    gram.lvl_size.push_back(tot_syms);
}

template<class sym_type>
void insert_comp_string(lc_gram_buffer_t& gram, std::string& input_file){

    i_file_stream<sym_type> ifs(input_file, BUFFER_SIZE);

    std::vector<size_t> rules_buffer;
    basic_load_vector_from_file(gram.rules_file, rules_buffer);
    assert(rules_buffer.size()==gram.g);

    rules_buffer.reserve(rules_buffer.size()+ifs.size());
    for(size_t i=0;i<ifs.size();i++){
        size_t sym = ifs.read(i)<<1UL | (i==0);
        rules_buffer.push_back(sym);
    }
    gram.r++;
    gram.g+=ifs.size();
    gram.c=ifs.size();
    ifs.close();

    assert(gram.g==rules_buffer.size());
    basic_store_vector_to_file(gram.rules_file, rules_buffer);
}

template<class par_type>
size_t par_phase_int(std::string& i_file, std::string& o_file, parsing_info& p_info,
                     size_t hbuff_size, size_t n_threads, lc_gram_buffer_t& gram, tmp_workspace& ws){
    size_t alph_size;
    if (n_threads > 1) {
        {
            using par_strat_t = mt_parse_strat_t<par_type, hash_functor, parse_functor>;
            par_strat_t p_strat(i_file, o_file, p_info, hbuff_size, n_threads);
            alph_size = par_round<par_strat_t>(p_strat, p_info, gram, ws);
        }
    } else {
        {
            using par_strat_t = st_parse_strat_t<par_type, hash_functor, parse_functor>;
            par_strat_t p_strat(i_file, o_file, p_info);
            alph_size = par_round<par_strat_t>(p_strat, p_info, gram, ws);
        }
    }
    return alph_size;
}

template<class sym_type>
size_t build_lc_grammar(std::string &i_file, std::string &gram_o_file, size_t n_threads, float hbuff_frac, tmp_workspace &ws) {

    std::cout<<"Reading the file"<<std::endl;
    str_collection str_coll = collection_stats<sym_type>(i_file);
    std::cout<<"Stats: "<<std::endl;
    std::cout<<"  Smallest symbol               : "<<str_coll.min_sym<<std::endl;
    std::cout<<"  Greatest symbol               : "<<str_coll.max_sym<<std::endl;
    std::cout<<"  Number of symbols in the file : "<<str_coll.n_syms<<std::endl;
    std::cout<<"  Number of strings             : "<<str_coll.n_strings<<std::endl;
    if constexpr (std::is_same<uint8_t, sym_type>::value){
        std::cout<<"  Max sym freq.                 : "<<str_coll.max_sym_freq<<std::endl;
    }

    auto hbuff_size = std::max<size_t>(64 * n_threads, size_t(ceil(float(str_coll.n_syms) * hbuff_frac)));

    std::cout << "Parsing the text:    " << std::endl;
    if(n_threads>1){
        std::cout<<"  Running with up to "<<n_threads<<" working threads "<<std::endl;
        std::cout<<"  Using "<<float(hbuff_size)/1000000<<" megabytes for the thread hash tables ("<< (float(hbuff_size)/1000000)/float(n_threads)<<" megabytes each)"<<std::endl;
    }

    std::string output_file = ws.get_file("tmp_output");
    std::string tmp_i_file = ws.get_file("tmp_input");

    parsing_info p_info;
    p_info.lms_phrases = str_coll.max_sym + 1;
    p_info.max_symbol = str_coll.max_sym;
    p_info.str_ptrs.swap(str_coll.str_ptrs);
    p_info.str_ptrs.push_back((long) str_coll.n_syms);
    p_info.str_ptrs.shrink_to_fit();
    p_info.longest_str = str_coll.longest_string;

    //TODO testing
    //lc_parser_t<i_file_stream<uint8_t>, true> lcp(i_file, p_info.str_ptrs, str_coll.min_sym, str_coll.max_sym);
    //lcp.partition_text(20);
    //
    //exit(0);

    std::string gram_file = ws.get_file("gram");
    lc_gram_buffer_t gram_buff(gram_file, str_coll.alphabet, p_info.str_ptrs, str_coll.sep_symbol);
    assert(gram_buff.rules_buffer.size()==(gram_buff.max_tsym+1));

    size_t iter = 1;
    size_t n_syms;
    using f_parser_t = lms_parsing<i_file_stream<sym_type>, true>;

    std::cout << "  Parsing round " << iter++ << std::endl;
    auto start = std::chrono::steady_clock::now();
    n_syms = par_phase_int<f_parser_t>(i_file, tmp_i_file, p_info, hbuff_size, n_threads, gram_buff, ws);
    auto end = std::chrono::steady_clock::now();

    report_time(start, end, 4);
    malloc_count_print_status();
    malloc_count_reset_peak();

    while (n_syms > 0) {
        start = std::chrono::steady_clock::now();
        std::cout << "  Parsing round " << iter++ << std::endl;
        size_t bps = sym_width(p_info.max_symbol);
        if(bps<=8){
            n_syms = par_phase_int<uint8t_parser_t>(tmp_i_file, output_file, p_info,
                                                    hbuff_size, n_threads, gram_buff, ws);
        }else if(bps<=16){
            n_syms = par_phase_int<uint16t_parser_t>(tmp_i_file, output_file, p_info,
                                                     hbuff_size, n_threads, gram_buff, ws);
        } else if(bps<=32){
            n_syms = par_phase_int<uint32t_parser_t>(tmp_i_file, output_file, p_info,
                                                     hbuff_size, n_threads, gram_buff, ws);
        } else{
            n_syms = par_phase_int<uint64t_parser_t>(tmp_i_file, output_file, p_info,
                                                     hbuff_size, n_threads, gram_buff, ws);
        }
        end = std::chrono::steady_clock::now();
        report_time(start, end, 4);

        remove(tmp_i_file.c_str());
        rename(output_file.c_str(), tmp_i_file.c_str());
        malloc_count_print_status();
        malloc_count_reset_peak();
    }
    gram_buff.rules_buffer.close();

    run_length_compress(gram_buff);

    size_t bps = sym_width(p_info.max_symbol);
    if(bps<=8){
        insert_comp_string<uint8_t>(gram_buff, tmp_i_file);
    }else if(bps<=16){
        insert_comp_string<uint16_t>(gram_buff, tmp_i_file);
    } else if(bps<=32){
        insert_comp_string<uint32_t>(gram_buff, tmp_i_file);
    } else{
        insert_comp_string<uint64_t>(gram_buff, tmp_i_file);
    }

    std::cout<<"Compacting and storing the grammar "<<std::endl;
    lc_gram_t gram(gram_buff);
    gram.print_stats();

    //TODO check
    check_plain_grammar(gram, i_file);
    //

    size_t written_bytes = store_to_file(gram_o_file, gram);
    std::cout<<"Grammar encoding "<<float(written_bytes)/1000000<<" MBs"<<std::endl;

    return iter - 2;
}

template<class parse_strategy_t>
size_t par_round(parse_strategy_t &p_strategy, parsing_info &p_info, lc_gram_buffer_t& gram, tmp_workspace &ws) {

#ifdef __linux__
    malloc_trim(0);
#endif
    std::cout << "    Computing the dictionary of phrases" << std::flush;
    auto start = std::chrono::steady_clock::now();
    auto res = p_strategy.get_phrases();
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 22);

    store_pl_vector(ws.get_file("str_ptr"), p_info.str_ptrs);
    std::vector<long>().swap(p_info.str_ptrs);

#ifdef __linux__
    malloc_trim(0);
#endif
    phrase_map_t &map = p_strategy.map;
    size_t psize;//<- for the iter stats
    assert(map.size() > 0);

    size_t dict_sym = res.first;
    {
        std::cout << "    Assigning metasymbols to the phrases" << std::flush;
        start = std::chrono::steady_clock::now();
        size_t j = gram.r;
        for (auto const &ptr: map) {
            size_t val = j++;
            map.insert_value_at(ptr, val);
        }
        create_lc_rules(gram, map, gram.r);
        end = std::chrono::steady_clock::now();
        report_time(start, end, 21);
    }

    std::cout<< "    Creating the parse of the text" << std::flush;
    start = std::chrono::steady_clock::now();
    load_pl_vector(ws.get_file("str_ptr"), p_info.str_ptrs);

    p_info.max_symbol = gram.r-1;

    size_t bps = sym_width(p_info.max_symbol);
    if(bps<=8){
        psize = p_strategy.template parse_text<uint8_t>();
    }else if(bps<=16){
        psize = p_strategy.template parse_text<uint16_t>();
    } else if(bps<=32){
        psize = p_strategy.template parse_text<uint32_t>();
    }else{
        psize = p_strategy.template parse_text<uint64_t>();
    }
    assert(psize>=map.size());//the parse can't be smaller than the number of phrases

    end = std::chrono::steady_clock::now();
    report_time(start, end, 27);

    p_info.lms_phrases = map.size();
    p_info.p_round++;

    std::cout << "    Stats:" << std::endl;
    std::cout << "      Parsing phrases:                  " << p_info.lms_phrases << std::endl;
    std::cout << "      Number of symbols in the phrases: " << dict_sym << std::endl;
    std::cout << "      Parse size:                       " << psize << std::endl;

    map.destroy_data();
    map.destroy_table();

#ifdef __linux__
    malloc_trim(0);
#endif
    return (p_info.str_ptrs.size()-1) == psize ? 0 : p_info.lms_phrases;
}

template unsigned long build_lc_grammar<uint8_t>(std::string &i_file, std::string& gram_o_file, size_t n_threads, float hbuff_frac, tmp_workspace &ws);
template unsigned long build_lc_grammar<uint16_t>(std::string &i_file, std::string& gram_o_file, size_t n_threads, float hbuff_frac, tmp_workspace &ws);
template unsigned long build_lc_grammar<uint32_t>(std::string &i_file, std::string& gram_o_file, size_t n_threads, float hbuff_frac, tmp_workspace &ws);
template unsigned long build_lc_grammar<uint64_t>(std::string &i_file, std::string& gram_o_file, size_t n_threads, float hbuff_frac, tmp_workspace &ws);
