//
// Created by Diaz, Diego on 26.1.2023.
//

#include "malloc_count.h"
#include "build_lc_grammar.hpp"
#include "grammar.h"

template<class sym_type>
void create_lc_rules(lc_gram_buffer_t& gram, std::vector<rand_order>& str_orders, std::vector<sym_type>& parsing_set){
    size_t sym, tot_syms=0, len;
    bool first;
    for(auto & str : str_orders){
        first = true;
        len = str.str_len;
        for(size_t j=0;j<len;j++){
            sym = parsing_set[str.str_ptr+j];
            assert(sym>0);
            sym = (sym<<1UL) | first;
            gram.rules_buffer.push_back(sym);
            first = false;
        }
        tot_syms+=len;
    }
    gram.r += str_orders.size();
    gram.g += tot_syms;
    gram.lvl_rules.push_back(str_orders.size());
    gram.lvl_size.push_back(tot_syms);
}

/*
template<class sym_type>
void assign_symbols_to_phrases_int(std::vector<rand_order>& str_orders, std::vector<sym_type>& parsing_set, size_t lvl,
                                   size_t start, size_t end, std::vector<hashing>& p_functions){

    if(p_functions.size()==lvl){
        hashing par_function;
        p_functions.push_back(par_function);
    }

    for(size_t i=start; i<=end; i++){
        str_orders[i].hash = p_functions[lvl].string_hash((char *)&parsing_set[str_orders[i].str_ptr], str_orders[i].str_len*sizeof(sym_type),64);
    }

    auto first = long(start);
    auto last = long(end)+1;

    std::sort(str_orders.begin()+first, str_orders.begin()+last, [&](auto a, auto b) -> bool{
        if(a.hash==b.hash && a.orig_order!=b.orig_order){
            std::cout<<" we had equal hash values "<<a.orig_order<<" "<<b.orig_order<<std::endl;
        }
        return a.hash<b.hash;
    });

    size_t prev_hash = str_orders[0].hash;
    size_t prev_pos=0;
    for(size_t i=start; i<=end; i++){
        if(prev_hash!=str_orders[i].hash){
            if((i-prev_pos)>1){
                assign_symbols_to_phrases_int<sym_type>(str_orders, parsing_set, lvl+1, prev_pos, i-1, p_functions);
            }
            prev_pos = i;
            prev_hash = str_orders[i].hash;
        }
    }
}*/

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

template<class sym_type, bool first_round>
void assign_symbols_to_phrases(phrase_map_t& map, hashing& p_func, lc_gram_buffer_t& gram_buff, parsing_info& p_info, tmp_workspace& ws){

    std::vector<rand_order> str_orders(map.size());
    std::vector<sym_type> parsing_set(p_info.par_symbols, 0);
    std::vector<size_t> hash_values_vector;

    if constexpr(first_round){
        hash_values_vector.resize(p_info.max_symbol+1);
        for(size_t i=0;i<=p_info.max_symbol;i++){
            hash_values_vector[i] = p_func.symbol_hash(i);
        }
    }else{
        std::string hash_values_file = ws.get_file("hash_values_round_"+std::to_string(p_info.p_round-1));
        load_pl_vector(hash_values_file, hash_values_vector);
        assert(!hash_values_vector.empty());
    }

    key_wrapper key_w{sym_width(p_info.max_symbol), map.description_bits(), map.get_data()};

    size_t k=0, u=0, len, sym;
    std::vector<size_t> tmp_vector;

    for (auto const &ptr: map) {
        len = key_w.size(ptr);
        str_orders[u] = {k, len, 0, u};
        tmp_vector.clear();
        for(size_t i=key_w.size(ptr);i-->0;){
            sym = key_w.read(ptr, i);
            parsing_set[k++] = sym;
            tmp_vector.push_back(hash_values_vector[sym-p_info.min_symbol]);
        }
        assert(tmp_vector.size()==len);
        str_orders[u].hash = p_func.string_hash((char *)tmp_vector.data(), tmp_vector.size()*sizeof(size_t),64);
        u++;
    }

    //sort the phrases according their hash values
    std::sort(str_orders.begin(), str_orders.end(), [&](auto a, auto b) -> bool{
        return a.hash<b.hash;
    });

    //resolve colliding phrases according the lexicographical rank of their random symbols
    size_t prev_hash = str_orders[0].hash;
    long prev_pos=0, tot_phrases=long(str_orders.size());
    for(long i=0; i<tot_phrases; i++){
        if(prev_hash!=str_orders[i].hash){
            if((i-prev_pos)>1){
                std::cout<<"Warning: we have "<<(i-prev_pos)<<" colliding phrases"<<std::endl;
                //sort the range [prev_pos..i-1]
                std::sort(str_orders.begin()+prev_pos, str_orders.begin()+i, [&](auto a, auto b) -> bool{
                    size_t len = std::min(a.str_len, b.str_len);
                    size_t j=0;
                    while(j<len && parsing_set[a.str_ptr+j]==parsing_set[b.str_ptr+j]) j++;
                    assert(j<len);
                    if constexpr (first_round){
                        return p_func.symbol_hash(parsing_set[a.str_ptr+j])<p_func.symbol_hash(parsing_set[b.str_ptr+j]);
                    }else{
                        return parsing_set[a.str_ptr+j]<parsing_set[b.str_ptr+j];
                    }
                });
            }
            prev_pos = i;
            prev_hash = str_orders[i].hash;
        }
    }

    hash_values_vector.resize(map.size());
    for(size_t i=0;i<map.size();i++){
        hash_values_vector[i] = str_orders[i].hash;
    }
    std::string hash_values_file = ws.get_file("hash_values_round_"+std::to_string(p_info.p_round));
    store_pl_vector(hash_values_file, hash_values_vector);

    std::vector<size_t> new_order(map.size());
    for(size_t i=0;i<map.size();i++){
        new_order[str_orders[i].orig_order] = i;
    }

    size_t j = 0;
    for (auto const &ptr: map) {
        size_t val = gram_buff.r+new_order[j++];
        map.insert_value_at(ptr, val);
    }
    p_info.min_symbol = gram_buff.r;
    p_info.max_symbol = p_info.min_symbol+map.size()-1;
    create_lc_rules<sym_type>(gram_buff, str_orders, parsing_set);
}

template<class sym_type, bool first_round>
size_t p_round_int(std::string& i_file, std::string& o_file, parsing_info& p_info, lc_gram_buffer_t& gram_buff,
                   size_t n_tries, size_t n_threads,  tmp_workspace& ws){

    lc_parser_t<sym_type, first_round> lcp(i_file, p_info, ws);
    std::tie(p_info.par_symbols, p_info.par_phrases) = lcp.partition_text(n_tries);

    assign_symbols_to_phrases<sym_type, first_round>(lcp.get_map(), lcp.get_par_function(), gram_buff, p_info, ws);
    size_t p_size = lcp.produce_next_string(o_file);
    gram_buff.par_functions.push_back(lcp.get_par_function());

    std::cout << "    Stats:" << std::endl;
    std::cout << "      Parsing phrases:                  " << p_info.par_phrases << std::endl;
    std::cout << "      Number of symbols in the phrases: " << p_info.par_symbols << std::endl;
    std::cout << "      Parse size:                       " << p_size << std::endl;

    return p_size;
}

template<class sym_type>
size_t build_lc_grammar(std::string &i_file, std::string& pf_file, std::string &gram_o_file, size_t n_tries, size_t n_threads, tmp_workspace &ws) {

    std::cout<<"Reading the file"<<std::endl;
    str_collection str_coll = collection_stats<sym_type>(i_file);
    std::cout<<"Stats: "<<std::endl;
    std::cout<<"  Effective alphabet            : "<<str_coll.alphabet.size()<<std::endl;
    std::cout<<"  Smallest symbol               : "<<str_coll.min_sym<<std::endl;
    std::cout<<"  Greatest symbol               : "<<str_coll.max_sym<<std::endl;
    std::cout<<"  Number of symbols in the file : "<<str_coll.n_syms<<std::endl;
    std::cout<<"  Number of strings             : "<<str_coll.n_strings<<std::endl;

    std::cout << "Creating the locally-consistent grammar:    " << std::endl;
    std::string output_file = ws.get_file("tmp_output");
    std::string tmp_i_file = ws.get_file("tmp_input");

    parsing_info p_info;
    p_info.par_phrases = str_coll.max_sym + 1;
    p_info.min_symbol = str_coll.min_sym;
    p_info.max_symbol = str_coll.max_sym;
    p_info.str_ptrs.swap(str_coll.str_ptrs);
    p_info.str_ptrs.push_back((long) str_coll.n_syms);
    p_info.str_ptrs.shrink_to_fit();
    p_info.longest_str = str_coll.longest_string;
    if(!pf_file.empty()){
        load_pl_vector(pf_file, p_info.p_functions);
        std::cout<<"Using the parsing functions in \""<<pf_file<<"\" for the first "<<p_info.p_functions.size()<<" parsing rounds "<<std::endl;
    }

    std::string gram_file = ws.get_file("gram");
    lc_gram_buffer_t gram_buff(gram_file, str_coll.alphabet, p_info.str_ptrs, p_info.longest_str, str_coll.sep_symbol);
    assert(gram_buff.rules_buffer.size()==(gram_buff.max_tsym+1));

    size_t iter = 1;

    std::cout << "  Parsing round " << iter++ << std::endl;
    size_t n_syms = p_round_int<sym_type, true>(i_file, tmp_i_file, p_info, gram_buff, n_tries, n_threads, ws);
    p_info.p_round++;

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
        p_info.p_round++;
        remove(tmp_i_file.c_str());
        rename(output_file.c_str(), tmp_i_file.c_str());
    }

    gram_buff.rules_buffer.close();


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

    lc_gram_t gram(gram_buff);
    store_to_file(gram_o_file, gram);

    return iter - 2;
}
template unsigned long build_lc_grammar<uint8_t>(std::string &i_file, std::string &pf_file, std::string& gram_o_file, size_t n_tries, size_t n_threads, tmp_workspace &ws);
template unsigned long build_lc_grammar<uint16_t>(std::string &i_file,std::string &pf_file, std::string& gram_o_file, size_t n_tries, size_t n_threads, tmp_workspace &ws);
template unsigned long build_lc_grammar<uint32_t>(std::string &i_file, std::string &pf_file, std::string& gram_o_file, size_t n_tries, size_t n_threads, tmp_workspace &ws);
template unsigned long build_lc_grammar<uint64_t>(std::string &i_file, std::string &pf_file, std::string& gram_o_file, size_t n_tries, size_t n_threads, tmp_workspace &ws);
