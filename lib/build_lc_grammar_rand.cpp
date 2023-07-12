//
// Created by Diaz, Diego on 26.1.2023.
//

#include "malloc_count.h"
#include "build_lc_grammar_rand.hpp"
#include "parsing_strategies.h"
#include "grammar.h"

void create_lc_rules(lc_gram_buffer_t& gram, phrase_map_t& map){
    key_wrapper key_w{sym_width(gram.r), map.description_bits(), map.get_data()};
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


template<class sym_type, bool f_round>
size_t p_round_int(std::string& i_file, std::string& o_file, parsing_info& p_info, lc_gram_buffer_t& gram_buff, size_t n_tries, size_t n_threads,  tmp_workspace& ws){

    lc_parser_t<sym_type, f_round> lcp(i_file, p_info.str_ptrs, p_info.min_symbol, p_info.max_symbol, ws);
    std::pair<size_t, size_t> res= lcp.partition_text(n_tries);

    phrase_map_t &map = lcp.get_map();
    size_t j = gram_buff.r;
    for (auto const &ptr: map) {
        size_t val = j++;
        map.insert_value_at(ptr, val);
    }

    p_info.min_symbol = gram_buff.r;
    p_info.max_symbol = p_info.min_symbol+map.size()-1;

    create_lc_rules(gram_buff, map);
    size_t p_size = lcp.produce_next_string(o_file);

    std::cout << "    Stats:" << std::endl;
    std::cout << "      Parsing phrases:                  " << res.second << std::endl;
    std::cout << "      Number of symbols in the phrases: " << res.first << std::endl;
    std::cout << "      Parse size:                       " << p_size << std::endl;
    return p_size;
}

template<class sym_type>
size_t build_lc_grammar(std::string &i_file, std::string &gram_o_file, size_t n_tries, size_t n_threads, tmp_workspace &ws) {

    std::cout<<"Reading the file"<<std::endl;
    str_collection str_coll = collection_stats<sym_type>(i_file);
    std::cout<<"Stats: "<<std::endl;
    std::cout<<"  Effective alphabet            : "<<str_coll.alphabet.size()<<std::endl;
    std::cout<<"  Smallest symbol               : "<<str_coll.min_sym<<std::endl;
    std::cout<<"  Greatest symbol               : "<<str_coll.max_sym<<std::endl;
    std::cout<<"  Number of symbols in the file : "<<str_coll.n_syms<<std::endl;
    std::cout<<"  Number of strings             : "<<str_coll.n_strings<<std::endl;

    std::cout << "Parsing the text:    " << std::endl;
    std::string output_file = ws.get_file("tmp_output");
    std::string tmp_i_file = ws.get_file("tmp_input");

    parsing_info p_info;
    p_info.lms_phrases = str_coll.max_sym + 1;
    p_info.min_symbol = str_coll.min_sym;
    p_info.max_symbol = str_coll.max_sym;
    p_info.str_ptrs.swap(str_coll.str_ptrs);
    p_info.str_ptrs.push_back((long) str_coll.n_syms);
    p_info.str_ptrs.shrink_to_fit();
    p_info.longest_str = str_coll.longest_string;

    std::string gram_file = ws.get_file("gram");
    lc_gram_buffer_t gram_buff(gram_file, str_coll.alphabet, p_info.str_ptrs, str_coll.sep_symbol);
    assert(gram_buff.rules_buffer.size()==(gram_buff.max_tsym+1));

    size_t iter = 1;

    std::cout << "  Parsing round " << iter++ << std::endl;
    size_t n_syms = p_round_int<sym_type, true>(i_file, tmp_i_file, p_info, gram_buff, n_tries, n_threads, ws);

    while(n_syms!=p_info.str_ptrs.size()-1){
        std::cout << "  Parsing round " << iter++ << std::endl;
        size_t bps = sym_width(p_info.max_symbol);
        if(bps<=8){
            n_syms = p_round_int<uint8_t, false>(tmp_i_file, output_file, p_info, gram_buff, n_tries, n_threads, ws);
        }else if(bps<=16){
            n_syms = p_round_int<uint16_t, false>(tmp_i_file, output_file, p_info, gram_buff, n_tries, n_threads, ws);
        } else if(bps<=32){
            n_syms = p_round_int<uint32_t, false>(tmp_i_file, output_file, p_info, gram_buff, n_tries, n_threads, ws);
        } else{
            n_syms = p_round_int<uint64_t, false>(tmp_i_file, output_file, p_info, gram_buff, n_tries, n_threads, ws);
        }
        remove(tmp_i_file.c_str());
        rename(output_file.c_str(), tmp_i_file.c_str());
    }

    /*auto start = std::chrono::steady_clock::now();
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
    }*/

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

/*
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
*/
template unsigned long build_lc_grammar<uint8_t>(std::string &i_file, std::string& gram_o_file, size_t n_tries, size_t n_threads, tmp_workspace &ws);
template unsigned long build_lc_grammar<uint16_t>(std::string &i_file, std::string& gram_o_file, size_t n_tries, size_t n_threads, tmp_workspace &ws);
template unsigned long build_lc_grammar<uint32_t>(std::string &i_file, std::string& gram_o_file, size_t n_tries, size_t n_threads, tmp_workspace &ws);
template unsigned long build_lc_grammar<uint64_t>(std::string &i_file, std::string& gram_o_file, size_t n_tries, size_t n_threads, tmp_workspace &ws);
