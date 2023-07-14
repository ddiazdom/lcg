//
// Created by Diaz, Diego on 14.7.2023.
//

#ifndef LCG_GRAMMAR_ALGORITHMS_H
#define LCG_GRAMMAR_ALGORITHMS_H

#include "grammar.h"
#include "build_lc_grammar.hpp"

size_t get_new_rl_rules(lc_gram_t& gram, phrase_map_t& ht) {

    size_t prev_sym, curr_sym, run_len, tmp_sym;
    string_t pair(2, sym_width(gram.rules.size()));
    size_t start_sym = gram.start_symbol();
    size_t new_id = start_sym;
    size_t new_size = gram.n_terminals();

    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++){
        auto range = gram.nt2phrase(sym);
        prev_sym = gram.rules[range.first];
        run_len = 1;

        for(size_t j=range.first+1;j<=range.second;j++){
            curr_sym = gram.rules[j];
            if(curr_sym!=prev_sym){
                if(run_len>1){
                    pair.write(0, prev_sym);
                    pair.write(1, run_len);
                    auto res = ht.insert(pair.data(), pair.n_bits(), 0);
                    if(res.second){
                        tmp_sym = new_id++;
                        ht.insert_value_at(res.first, tmp_sym);
                    }
                }
                new_size++;
                prev_sym = curr_sym;
                run_len=0;
            }
            run_len++;
        }

        if(run_len>1){
            pair.write(0, prev_sym);
            pair.write(1, run_len);
            auto res = ht.insert(pair.data(), pair.n_bits(), 0);
            if(res.second){
                tmp_sym = new_id++;
                ht.insert_value_at(res.first, tmp_sym);
            }
        }
        new_size++;
    }

    //deal with the strings
    for(size_t str=0;str<gram.n_strings();str++){
        auto range = gram.str2phrase(str);
        prev_sym = gram.rules[range.first];
        run_len = 1;

        for(size_t j=range.first+1;j<=range.second;j++){
            curr_sym = gram.rules[j];
            if(curr_sym!=prev_sym){
                if(run_len>1){
                    pair.write(0, prev_sym);
                    pair.write(1, run_len);
                    auto res = ht.insert(pair.data(), pair.n_bits(), 0);
                    if(res.second){
                        tmp_sym = new_id++;
                        ht.insert_value_at(res.first, tmp_sym);
                    }
                }
                new_size++;

                prev_sym = curr_sym;
                run_len=0;
            }
            run_len++;
        }

        if(run_len>1){
            pair.write(0, prev_sym);
            pair.write(1, run_len);
            auto res = ht.insert(pair.data(), pair.n_bits(), 0);
            if(res.second){
                tmp_sym = new_id;
                ht.insert_value_at(res.first, tmp_sym);
            }
        }
        new_size++;
    }

    return new_size+(ht.size()*2);
}

void run_length_compress(lc_gram_t& gram) {

    phrase_map_t ht;
    size_t new_size = get_new_rl_rules(gram, ht);

    int_array<size_t> new_rules(new_size, sym_width(gram.r+ht.size()));
    int_array<size_t> new_rl_ptrs(sym_width(gram.r+ht.size()), sym_width(new_size));
    size_t start_sym = gram.start_symbol(), prev_sym, run_len, curr_sym;
    string_t pair(2, sym_width(gram.rules.size()));

    for(size_t i=0;i<=gram.max_tsym;i++){
        new_rules.push_back(i);
    }

    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++) {

        new_rl_ptrs.push_back(new_rules.size());

        auto range = gram.nt2phrase(sym);
        prev_sym = gram.rules[range.first];
        run_len = 1;

        for(size_t j=range.first+1;j<=range.second;j++){

            curr_sym = gram.rules[j];
            if(curr_sym!=prev_sym){
                if(run_len>1){
                    pair.write(0, prev_sym);
                    pair.write(1, run_len);
                    curr_sym=0;
                    auto res = ht.find(pair.data(), pair.n_bits());
                    assert(res.second);
                    ht.get_value_from(res.first, curr_sym);
                }
                new_rules.push_back(curr_sym);
                prev_sym = curr_sym;
                run_len=0;
            }
            run_len++;
        }

        if(run_len>1){
            pair.write(0, prev_sym);
            pair.write(1, run_len);
            auto res = ht.insert(pair.data(), pair.n_bits(), 0);
            if(res.second){
                curr_sym = 0;
                ht.insert_value_at(res.first, curr_sym);
            }
        }
        new_rules.push_back(curr_sym);
    }

    //insert the new strings
    

    //deal with the strings
    new_rl_ptrs.push_back(new_rules.size());
    for(size_t str=0;str<gram.n_strings();str++){
        auto range = gram.str2phrase(str);
        gram.str_boundaries[str] = new_rules.size();

        prev_sym = gram.rules[range.first];
        run_len = 1;

        for(size_t j=range.first+1;j<=range.second;j++){
            curr_sym = gram.rules[j];
            if(curr_sym!=prev_sym){
                if(run_len>1){
                    pair.write(0, prev_sym);
                    pair.write(1, run_len);
                    curr_sym=0;
                    auto res = ht.find(pair.data(), pair.n_bits());
                    assert(res.second);
                    ht.get_value_from(res.first, curr_sym);
                }
                new_rules.push_back(curr_sym);
                prev_sym = curr_sym;
                run_len=0;
            }
            run_len++;
        }

        if(run_len>1){
            pair.write(0, prev_sym);
            pair.write(1, run_len);
            auto res = ht.insert(pair.data(), pair.n_bits(), 0);
            if(res.second){
                curr_sym = 0;
                ht.insert_value_at(res.first, curr_sym);
            }
        }
        new_rules.push_back(curr_sym);
    }
    new_rl_ptrs.push_back(new_rules.size());

    gram.run_len_nt.first = start_sym;
    gram.run_len_nt.second = ht.size();
    gram.r += ht.size();
    gram.g  = new_rules.size();
    gram.c = gram.g - gram.str_boundaries[0];

    std::cout<<"  Stats:"<<std::endl;
    std::cout<<"    Grammar size before:        "<<gram.rules.size()<<std::endl;
    std::cout<<"    Grammar size after:         "<<new_rules.size()<<std::endl;
    std::cout<<"    Number of new nonterminals: "<<ht.size()<<std::endl;
    std::cout<<"    Compression ratio:          "<<float(new_rules.size())/float(gram.rules.size())<<std::endl;
}

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
                decompression.push_back((char)gram.get_byte_ter(curr_sym));
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

std::pair<std::vector<uint8_t>, size_t> mark_disposable_symbols(const lc_gram_t& gram) {

    //compute which nonterminals are repeated and
    // which have a replacement of length 1
    std::vector<uint8_t> rep_nts(gram.r + 1, 0);
    for(size_t rule=gram.max_tsym+1;rule<gram.r;rule++){
        auto range = gram.nt2phrase(rule);
        if(gram.is_rl_sym((rule))){
            rep_nts[gram.rules[range.first]] = 2;
        }else{
            for(size_t i=range.first;i<=range.second;i++){
                rep_nts[gram.rules[i]]+=rep_nts[gram.rules[i]]<2;
            }
        }
    }

    std::vector<uint8_t> rem_nts(gram.r + 1, 0);
    //mark the rules to remove
    //1) nonterminals with freq 1
    //2) terminal symbols between [min_sym..max_sym] with
    // frequency zero: to compress the alphabet
    bool remove;
    size_t n_rem=0;
    size_t start_sym = gram.start_symbol();
    for(size_t sym=0;sym<start_sym;sym++){
        remove = rep_nts[sym]==0 || (rep_nts[sym]==1 && sym > gram.max_tsym && !gram.is_rl_sym(sym));
        rem_nts[sym] = remove;
        n_rem+=remove;
    }
    return {rem_nts, n_rem};
}

void simplify_grammar(lc_gram_t& gram) {

    assert(!gram.simplified);

    auto rem_syms = mark_disposable_symbols(gram);
    int_array<size_t> new_rules(gram.g-rem_syms.second, sym_width(gram.r-rem_syms.second));
    int_array<size_t> new_rl_ptrs(gram.r-rem_syms.second, sym_width((gram.g+1)-rem_syms.second));
    int_array<size_t> offsets(gram.r, sym_width(gram.r));

    size_t del_syms=0;
    for(size_t sym=0;sym<gram.r;sym++){
        offsets[sym] = del_syms;
        del_syms+=rem_syms.first[sym];
    }

    size_t n_ter=0, new_ter;
    for(size_t ter=0;ter<=gram.max_tsym;ter++){
        if(!rem_syms.first[ter]){
            new_ter = ter-offsets[ter];
            gram.terminals[new_ter] = ter;
            new_rules.push_back(new_ter);
            n_ter++;
        }
    }

    std::stack<size_t> stack;
    size_t start_sym = gram.start_symbol();
    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++){

        if(!rem_syms.first[sym]) {

            auto range = gram.nt2phrase(sym);
            new_rl_ptrs.push_back(new_rules.size());

            //TODO testing
            //std::cout<<sym<<" "<<sym-offsets[sym]<<" "<<sym-offsets[sym]-(sigma-1)<<" "<<(int)rem_syms.first[sym]<<" "<<new_rl_ptrs[sym-offsets[sym]-(sigma-1)]<<" "<<std::endl;
            //

            if(gram.is_rl_sym(sym)){
                new_rules.push_back(gram.rules[range.first]-offsets[gram.rules[range.first]]);
                new_rules.push_back(gram.rules[range.second]);
            }else{
                for(size_t j=range.first;j<=range.second;j++){
                    if(rem_syms.first[gram.rules[j]]) {
                        stack.push(gram.rules[j]);
                        while(!stack.empty()){
                            auto nt = stack.top();
                            stack.pop();
                            if(rem_syms.first[nt]){
                                assert(nt>gram.max_tsym);
                                assert(!gram.is_rl_sym(nt));
                                auto range2 = gram.nt2phrase(nt);
                                for(size_t k=range2.second+1;k-->range2.first;){
                                    stack.push(gram.rules[k]);
                                }
                            }else{
                                new_rules.push_back(nt - offsets[nt]);
                            }
                        }
                    }else{
                        new_rules.push_back(gram.rules[j] - offsets[gram.rules[j]]);
                    }
                }
            }
        }
    }

    //deal with the start symbol
    auto range = gram.nt2phrase(start_sym);
    size_t str=0;
    //std::cout<<start_sym<<" "<<start_sym-offsets[start_sym]<<" "<<new_rules.size()<<std::endl;
    new_rl_ptrs.push_back(new_rules.size());
    for(size_t j=range.first;j<=range.second;j++){
        gram.str_boundaries[str++] = new_rules.size();
        if(rem_syms.first[gram.rules[j]]) {
            stack.push(gram.rules[j]);
            while(!stack.empty()){
                auto nt = stack.top();
                stack.pop();
                if(rem_syms.first[nt]){
                    auto range2 = gram.nt2phrase(nt);
                    for(size_t k=range2.second+1;k-->range2.first;){
                        stack.push(gram.rules[k]);
                    }
                }else{
                    new_rules.push_back(nt - offsets[nt]);
                }
            }
        }else{
            new_rules.push_back(gram.rules[j] - offsets[gram.rules[j]]);
        }
    }
    new_rl_ptrs.push_back(new_rules.size());

    //TODO testing
    //size_t rm_nt=0;
    //for(size_t i=max_tsym+1;i<r;i++){
    //    if(rem_syms.first[i]) rm_nt++;
    //}
    //

    size_t del_nt = (rem_syms.second-(gram.max_tsym+1-n_ter));
    //size_t del_ter = rem_syms.second - del_nt;
    //std::cout<<del_ter<<std::endl;
    //std::cout<<new_rules.size()<<" "<<rules.size()<<" "<<rules.size()-new_rules.size()<<" "<<rem_syms.second<<std::endl;
    //std::cout<<new_rl_ptrs.size()<<" "<<rl_ptr.size()<<" "<<(rl_ptr.size()-new_rl_ptrs.size())<<" "<<rm_nt<<std::endl;

    assert(new_rules.size()==gram.rules.size()-rem_syms.second);
    assert(new_rl_ptrs.size()==(gram.rl_ptr.size()-del_nt));

    float rm_per = float(rem_syms.second)/float(gram.r)*100;
    float comp_rat = float(new_rules.size())/float(gram.rules.size());
    std::cout<<"  Stats:"<<std::endl;
    std::cout<<"    Grammar size before:  "<<gram.g<<std::endl;
    std::cout<<"    Grammar size after:   "<<new_rules.size()<<std::endl;
    std::cout<<"    Deleted nonterminals: "<<rem_syms.second<<" ("<<rm_per<<"%)"<<std::endl;
    std::cout<<"    Compression ratio:    "<<comp_rat<<std::endl;

    gram.c = new_rules.size()-gram.str_boundaries[0];
    gram.r -= rem_syms.second;
    gram.g = new_rules.size();

    gram.rules.swap(new_rules);
    gram.rl_ptr.swap(new_rl_ptrs);

    for(auto &sym : gram.lvl_rules){
        //std::cout<<sym<<" "<<sym-offsets[sym]<<" "<<offsets[sym]<<std::endl;
        sym -= offsets[sym];
    }

    gram.run_len_nt.first -= offsets[gram.run_len_nt.first];

    gram.max_tsym = n_ter-1;
    gram.simplified = true;

    //TODO test
    /*for(size_t sym=run_len_nt.first;sym<(run_len_nt.first+run_len_nt.second);sym++){
        auto range = nt2phrase(sym);
        std::cout<<sym<<" "<<range.first<<" "<<range.second<<std::endl;
        assert(is_rl_sym(sym));
        assert((range.second-range.first+1)==2);
    }*/
    //
}

/***
 *
 * @param i_file : input text file
 * @param n_threads : number of working threads
 * @param hbuff_size : buffer size for the hashing step
 */
template<class sym_type>
void gram_algo(std::string &i_file, std::string& o_file, tmp_workspace & tmp_ws, size_t n_tries, size_t n_threads){

    build_lc_grammar<sym_type>(i_file, o_file, n_tries, n_threads, tmp_ws);

    lc_gram_t gram;
    load_from_file(o_file, gram);

    std::cout<<"Run-length compressing the grammar"<<std::endl;
    run_length_compress(gram);

    //std::cout<<"Simplifying the grammar "<<std::endl;
    //simplify_grammar(gram);

    gram.print_stats();

    //optional check
    check_plain_grammar(gram, i_file);
    //

    store_to_file(o_file, gram);
    std::cout<<"The resulting grammar was stored in "<<o_file<<std::endl;
}
#endif //LCG_GRAMMAR_ALGORITHMS_H
