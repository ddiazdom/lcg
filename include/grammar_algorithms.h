//
// Created by Diaz, Diego on 14.7.2023.
//

#ifndef LCG_GRAMMAR_ALGORITHMS_H
#define LCG_GRAMMAR_ALGORITHMS_H

#include "se-build-strat/lc_parsing.h"
#include "lz-like-build-strat/lc_parsing.h"
#include "grammar.h"
#include "build_collage_system.h"

template<class gram_t>
void make_gram_fix_free(gram_t& gram){

    for(size_t i=0;i<gram.lvl_rules.size()-1; i++) {

        size_t n_lvl_rules = gram.lvl_rules[i+1]-gram.lvl_rules[i];//number of rules in the level

        if(n_lvl_rules>0) {
            std::vector<collage_data> level_rules;
            level_rules.resize(n_lvl_rules);
            size_t first_nt = gram.lvl_rules[i];
            size_t last_nt = gram.lvl_rules[i + 1] - 1;

            for(size_t nt=first_nt, j=0;nt<=last_nt;nt++,j++){
                level_rules[j].nt = nt;
                auto res = gram.nt2bitrange(nt);
                level_rules[j].position = res.first;
                level_rules[j].rule_size = res.second-res.first+1;
            }

            //sort the level rules in lexicographical order
            std::sort(level_rules.begin(), level_rules.end(), [&](const auto& a, const auto& b){
                size_t len = std::min(a.rule_size, b.rule_size);
                for(size_t i=0;i<len;i++){
                    if(gram.rules[a.position+i]!=gram.rules[b.position+i]){
                        return gram.rules[a.position+i] < gram.rules[b.position+i];
                    }
                }
                return a.rule_size<b.rule_size;
            });

            //get which rules are suffixes of others
            std::unordered_set<size_t> suff_nts;
            for(size_t k=1;k<level_rules.size();k++) {
                size_t u = 0;
                size_t min_size = std::min(level_rules[k - 1].rule_size, level_rules[k].rule_size);
                size_t pos_a = level_rules[k - 1].position;
                size_t pos_b = level_rules[k].position;

                while (u < min_size &&
                       gram.rules[pos_b + u] == gram.rules[pos_a + u]) {
                    u++;
                }

                if (min_size == u) {
                    assert(level_rules[k-1].rule_size == min_size);
                    suff_nts.insert(level_rules[k - 1].nt);
                }
            }
            std::cout<<"There are "<<suff_nts.size()<<" suffix rules "<<std::endl;

            //inspect the upper grammar levels to obtain the contexts of the suffix rules
            uint8_t width = gram.rules.width();
            std::unordered_set<size_t> ext_contexts;
            for(size_t l=i+1;l<gram.lvl_rules.size()-1;l++){
                first_nt = gram.lvl_rules[l];
                last_nt = gram.lvl_rules[l + 1] - 1;
                n_lvl_rules = last_nt - first_nt+1;
                if(n_lvl_rules>1) {
                    for (size_t nt = first_nt; nt <= last_nt; nt++) {
                        auto res = gram.nt2bitrange(nt);
                        for (size_t k = res.first; k <= res.second; k++) {
                            if (suff_nts.find(gram.rules[k]) != suff_nts.end()) {
                                //TODO hash the pair
                                if(k<res.second){
                                    size_t pair = (gram.rules[k]<<width) | gram.rules[k+1];
                                    ext_contexts.insert(pair);
                                }else{
                                    suff_nts.insert(nt);
                                }
                            }
                        }
                    }
                }
            }

            std::unordered_map<size_t, size_t> new_suff_nts;
            std::stack<size_t> left_branch, right_branch;
            //size_t new_nt = 0;
            for(auto const& context : ext_contexts){
                size_t left = context>>width;
                size_t right = context & ((1<<width)-1);

                size_t context_lvl = gram.parsing_level(left);

                //descend over the adjacent branches
                //size_t top_right = right;
                while(context_lvl>(i+1)){
                    auto res = gram.nt2bitrange(left);
                    left = gram.rules[res.second];
                    auto res2 = gram.nt2bitrange(right);
                    right = gram.rules[res2.first];

                    left_branch.push(left);
                    right_branch.push(right);
                    context_lvl--;
                }

                //create the new phrase gram.rules[res1.first..res1.second]{\cdot}right
                std::vector<uint32_t> new_phrase;
                auto res1 = gram.nt2bitrange(left);
                for(size_t u=res1.first;u<res1.second;u++){
                    new_phrase.push_back(gram.rules[u]);
                }
                auto res2 = gram.nt2bitrange(right);
                new_phrase.push_back(gram.rules[res2.first]);
                //TODO hash new_phrase and assign it a nonterminal value

                new_phrase.clear();
                //ascend back to do the renaming
                while(!right_branch.empty()){
                    left = left_branch.top();
                    left_branch.pop();

                    res1 = gram.nt2bitrange(left);
                    for(size_t u=res1.first; u<res1.second;u++){
                        new_phrase.push_back(gram.rules[u]);
                    }
                    //IMPORTANT: what happens if right is also modified? (I mean, its suffix?)
                    new_phrase.push_back(right);

                    //TODO hash the new phrase

                    right = right_branch.top();
                    right_branch.pop();
                    new_phrase.clear();
                }
                //TODO associate phrase,top_right with context in a map
            }
            //std::cout<<" there are "<<suff_nts.size()<<" distinct contexts for the suffix rules "<<std::endl;
        }
    }
}

template<class gram_t>
std::tuple<size_t, uint8_t, uint8_t> compute_exp_info(gram_t& gram, std::vector<uint64_t>& exp_len, std::vector<uint8_t>& str_width){

    size_t tmp_sym, len, exp_size;
    for(size_t i=0;i<=gram.max_tsym;i++) exp_len[i] = 1;

    size_t last_sym = gram.first_rl_sym(), rule_len, n_samples, samp_bits;

    //get the expansion length of each nonterminal in the grammar
    size_t max_r_samp_bits=0, r_samp_acc_bits=0;
    for(size_t sym=gram.max_tsym+1;sym<last_sym;sym++) {
        auto range = gram.nt2bitrange(sym);
        exp_size = 0;
        for(size_t i=range.first;i<=range.second;i+=gram.r_bits){
            tmp_sym = gram.rule_stream.read(i, i+gram.r_bits-1);
            len = 1;
            if(gram.is_rl_sym(tmp_sym)) {
                auto range2 = gram.nt2bitrange(tmp_sym);
                tmp_sym = gram.rule_stream.read(range2.first, range2.first + gram.r_bits-1);
                len = gram.rule_stream.read(range2.second, range2.second + gram.r_bits-1);
            }
            exp_size += exp_len[tmp_sym]*len;
        }
        exp_len[sym] = exp_size;

        rule_len = (range.second-range.first+gram.r_bits)/gram.r_bits;
        n_samples = (rule_len/gram.rl_samp_rate)+1;
        samp_bits = n_samples*sym_width(exp_size);
        r_samp_acc_bits += samp_bits;
        if(samp_bits>max_r_samp_bits) max_r_samp_bits = samp_bits;
    }

    last_sym = gram.last_rl_sym();
    for(size_t sym=gram.first_rl_sym();sym<=last_sym; sym++) {
        auto range = gram.nt2bitrange(sym);
        tmp_sym = gram.rule_stream.read(range.first, range.first + gram.r_bits-1);
        len = gram.rule_stream.read(range.second, range.second + gram.r_bits-1);
        exp_size = exp_len[tmp_sym]*len;
        exp_len[sym] = exp_size;
        samp_bits = sym_width(exp_size);
        r_samp_acc_bits+=samp_bits;
        if(samp_bits>max_r_samp_bits) max_r_samp_bits = samp_bits;
    }

    //compute the length of the strings
    size_t max_str_samp_bits=0, str_samp_acc_bits=0;
    for(size_t str=0;str<gram.n_strings();str++){
        auto range = gram.str2bitrange(str);
        exp_size = 0;
        for(size_t i=range.first;i<=range.second;i+=gram.r_bits){
            tmp_sym = gram.rule_stream.read(i, i+gram.r_bits-1);
            exp_size+=exp_len[tmp_sym];
        }
        rule_len = (range.second-range.first+gram.r_bits)/gram.r_bits;
        n_samples = rule_len/gram.str_samp_rate;
        str_width[str] = sym_width(exp_size);
        samp_bits = n_samples*str_width[str];
        str_samp_acc_bits += samp_bits;
        if(samp_bits>max_str_samp_bits) max_str_samp_bits = samp_bits;
    }

    size_t tot_samp_bits = r_samp_acc_bits + (sym_width(max_r_samp_bits))*(gram.r-gram.n_terminals()-1);
    tot_samp_bits += str_samp_acc_bits + (sym_width(max_str_samp_bits))*gram.n_strings();
    return {tot_samp_bits, sym_width(max_r_samp_bits), sym_width(max_str_samp_bits)};
}

template<class gram_t>
void add_random_access_support(gram_t& gram) {
    assert(!gram.has_cg_rules);

    std::vector<uint64_t> exp_tmp(gram.r, 0);
    std::vector<uint8_t> str_width(gram.s, 0);

    size_t samp_bits;
    uint8_t r_samp_bits, str_samp_bits;
    std::tie(samp_bits, r_samp_bits, str_samp_bits) = compute_exp_info(gram, exp_tmp, str_width);

    size_t samp_space = INT_CEIL(samp_bits, 8);
    std::cout<<"The sampling space is "<<report_space((off_t)samp_space)<<" "<<r_samp_bits<<" "<<str_samp_bits<<std::endl;

    bitstream<size_t> new_rule_stream;
    new_rule_stream.reserve_in_bits((gram.g*gram.r_bits)+samp_bits);
    int_array<size_t> new_rl_ptrs(gram.rl_ptr.size(), sym_width((gram.g*gram.r_bits)+samp_bits));

    //insert the terminals
    size_t bit_pos=0;
    for(size_t ter=0;ter<=gram.max_tsym;ter++){
        new_rule_stream.write(bit_pos, bit_pos+gram.r_bits-1, ter);
        bit_pos+=gram.r_bits;
    }

    //store the sampled elements of the rules
    size_t exp_acc, last_sym = gram.first_rl_sym(), tmp_sym, offset;
    size_t n_sampled, r_len, s_bits;
    for(size_t sym = gram.max_tsym+1; sym<last_sym; sym++){

        new_rl_ptrs.push_back(bit_pos);
        auto range = gram.nt2bitrange(sym);
        offset = bit_pos+(range.second-range.first+gram.r_bits);
        r_len = (range.second-range.first+gram.r_bits)/gram.r_bits;

        n_sampled = (r_len/gram.rl_samp_rate)+1;
        s_bits=sym_width(exp_tmp[sym]);

        exp_acc = 0;
        for(size_t i=range.first, j=1;i<=range.second;i+=gram.r_bits,j++){
            if(j%gram.rl_samp_rate==0){
                new_rule_stream.write(offset, offset+s_bits-1, exp_acc);
                offset+=s_bits;
            }
            tmp_sym = gram.rule_stream.read(i, i+gram.r_bits-1);
            new_rule_stream.write(bit_pos, bit_pos+gram.r_bits-1, tmp_sym);
            bit_pos+=gram.r_bits;
            exp_acc+=exp_tmp[tmp_sym];
        }

        assert(exp_acc==exp_tmp[sym]);
        new_rule_stream.write(offset, offset+s_bits-1, exp_acc);
        offset+=s_bits;
        new_rule_stream.write(offset, offset+r_samp_bits-1, s_bits*n_sampled);
        offset+=r_samp_bits;

        bit_pos=offset;
    }

    last_sym = gram.last_rl_sym();
    size_t tmp_len;
    for(size_t sym=gram.first_rl_sym();sym<=last_sym; sym++) {
        new_rl_ptrs.push_back(bit_pos);

        auto range = gram.nt2bitrange(sym);
        tmp_sym = gram.rule_stream.read(range.first, range.first+gram.r_bits-1);
        tmp_len = gram.rule_stream.read(range.second, range.second+gram.r_bits-1);

        new_rule_stream.write(bit_pos, bit_pos+gram.r_bits-1, tmp_sym);
        bit_pos+=gram.r_bits;
        new_rule_stream.write(bit_pos, bit_pos+gram.r_bits-1, tmp_len);
        bit_pos+=gram.r_bits;

        s_bits=sym_width(exp_tmp[sym]);
        new_rule_stream.write(bit_pos, bit_pos+s_bits-1, exp_tmp[sym]);
        bit_pos+=s_bits;
        new_rule_stream.write(bit_pos, bit_pos+r_samp_bits-1, s_bits);
        bit_pos+=r_samp_bits;
    }

    new_rl_ptrs.push_back(bit_pos);
    //compute the sample expansions for the strings
    for(size_t str=0;str<gram.n_strings();str++){

        auto range = gram.str2bitrange(str);
        gram.str_boundaries[str] = bit_pos;

        offset = bit_pos+(range.second-range.first+gram.r_bits);
        r_len = (range.second-range.first+gram.r_bits)/gram.r_bits;

        n_sampled = (r_len/gram.str_samp_rate);
        s_bits = str_width[str];

        exp_acc = 0;
        for(size_t i=range.first, j=1;i<=range.second;i+=gram.r_bits,j++){
            if(j%gram.str_samp_rate==0){
                new_rule_stream.write(offset, offset+s_bits-1, exp_acc);
                offset+=s_bits;
            }
            tmp_sym = gram.rule_stream.read(i, i+gram.r_bits-1);
            new_rule_stream.write(bit_pos, bit_pos+gram.r_bits-1, tmp_sym);
            bit_pos+=gram.r_bits;
            exp_acc+=exp_tmp[tmp_sym];
        }

        new_rule_stream.write(offset, offset+str_samp_bits-1, s_bits*n_sampled);
        offset+=str_samp_bits;

        bit_pos=offset;
    }

    new_rl_ptrs.push_back(bit_pos);

    gram.str_boundaries[gram.n_strings()] = bit_pos;
    gram.r_samp_bits = r_samp_bits;
    gram.str_samp_bits = str_samp_bits;
    assert(bit_pos==((gram.g*gram.r_bits)+samp_bits));

    new_rule_stream.swap(gram.rule_stream);
    new_rl_ptrs.swap(gram.rl_ptr);
}

//TODO implement this
/*
void merge_partial_grammars(lc_gram_t& gram_a, lc_gram_t& gram_b){

}*/

template<class gram_t>
size_t get_new_rl_rules(gram_t& gram, par_string_map<size_t>& ht) {

    size_t prev_sym, curr_sym, run_len, tmp_sym;
    size_t pair[2] = {0};
    size_t start_sym = gram.start_symbol();
    size_t new_id = start_sym;
    size_t new_size = gram.n_terminals();

    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++){
        auto range = gram.nt2bitrange(sym);
        //prev_sym = gram.rules[range.first];
        prev_sym = gram.rule_stream.read(range.first, range.first+gram.r_bits-1);
        run_len = 1;

        for(size_t j=range.first+gram.r_bits;j<=range.second;j+=gram.r_bits) {
            //curr_sym = gram.rules[j];
            curr_sym = gram.rule_stream.read(j, j+gram.r_bits-1);
            if(curr_sym!=prev_sym){
                if(run_len>1){
                    pair[0] = prev_sym;
                    pair[1] = run_len;
                    auto res = ht.insert((uint8_t *)&pair, sizeof(size_t)*2, 0);

                    if(res.second) {
                        tmp_sym = new_id++;
                        res.first = tmp_sym;
                    }
                }
                new_size++;
                prev_sym = curr_sym;
                run_len=0;
            }
            run_len++;
        }

        if(run_len>1){
            pair[0] = prev_sym;
            pair[1] = run_len;
            auto res = ht.insert((uint8_t *)&pair, sizeof(size_t)*2, 0);
            if(res.second){
                tmp_sym = new_id++;
                res.first = tmp_sym;
            }
        }
        new_size++;
    }

    //deal with the strings
    for(size_t str=0;str<gram.n_strings();str++){
        auto range = gram.str2bitrange(str);
        //prev_sym = gram.rules[range.first];
        prev_sym = gram.rule_stream.read(range.first, range.first+gram.r_bits-1);
        run_len = 1;

        for(size_t j=range.first+gram.r_bits;j<=range.second;j+=gram.r_bits){
            //curr_sym = gram.rules[j];
            curr_sym = gram.rule_stream.read(j, j+gram.r_bits-1);
            if(curr_sym!=prev_sym){

                if(run_len>1){
                    pair[0] = prev_sym;
                    pair[1] = run_len;
                    auto res = ht.insert((uint8_t *)&pair, sizeof(size_t)*2, 0);
                    if(res.second){
                        tmp_sym = new_id++;
                        res.first = tmp_sym;
                    }
                }
                new_size++;

                prev_sym = curr_sym;
                run_len=0;
            }
            run_len++;
        }

        if(run_len>1) {
            pair[0] = prev_sym;
            pair[1] = run_len;
            auto res = ht.insert((uint8_t *)&pair, sizeof(size_t)*2, 0);
            if(res.second){
                tmp_sym = new_id++;
                res.first = tmp_sym;
            }
        }
        new_size++;
    }
    return new_size+(ht.size()*2);
}

template<class gram_t>
void run_length_compress(gram_t& gram) {

    assert(gram.run_len_nt.second==0);
    par_string_map<size_t> ht(4, 0.6, 1);
    size_t new_size = get_new_rl_rules(gram, ht);
    uint8_t new_r_bits = sym_width(gram.r+ht.size());
    bitstream<size_t> new_rule_stream;
    new_rule_stream.reserve_in_bits(new_size*new_r_bits);

    int_array<size_t> new_rl_ptrs(gram.r+ht.size()+1, sym_width(new_size));
    size_t start_sym = gram.start_symbol(), prev_sym, run_len, curr_sym;
    size_t pair[2]={0};

    //insert terminals
    size_t new_bit_pos=0;
    for(size_t i=0;i<gram.n_terminals();i++){
        //new_rules.push_back(i);
        new_rule_stream.write(new_bit_pos, new_bit_pos+new_r_bits-1, i);
        new_bit_pos+=new_r_bits;
    }

    //insert regular rules
    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++) {

        new_rl_ptrs.push_back(new_bit_pos/new_r_bits);

        auto range = gram.nt2bitrange(sym);
        //prev_sym = gram.rules[range.first];
        prev_sym = gram.rule_stream.read(range.first, range.first+gram.r_bits-1);
        run_len = 1;

        for(size_t j=range.first+gram.r_bits;j<=range.second;j+=gram.r_bits){

            //curr_sym = gram.rules[j];
            curr_sym = gram.rule_stream.read(j, j+gram.r_bits-1);
            if(curr_sym!=prev_sym){
                if(run_len>1){
                    pair[0] = prev_sym;
                    pair[1] = run_len;
                    prev_sym = 0;
                    auto res = ht.find((uint8_t *)&pair, sizeof(size_t)*2, prev_sym);
                    assert(res);
                }

                //new_rule_stream.push_back(prev_sym);
                new_rule_stream.write(new_bit_pos, new_bit_pos+new_r_bits-1, prev_sym);
                new_bit_pos+=new_r_bits;
                prev_sym = curr_sym;
                run_len=0;
            }
            run_len++;
        }

        if(run_len>1){
            pair[0] = prev_sym;
            pair[1] = run_len;
            prev_sym = 0;
            auto res = ht.find((uint8_t *)&pair, sizeof(size_t)*2, prev_sym);
            assert(res);
            assert(prev_sym>=start_sym);
        }

        //new_rules.push_back(prev_sym);
        new_rule_stream.write(new_bit_pos, new_bit_pos+new_r_bits-1, prev_sym);
        new_bit_pos+=new_r_bits;
    }
    assert(new_rl_ptrs.size()==gram.n_nonterminals()-1);

    //insert the new rl rules
    for(auto const& phrase : ht){
        assert(phrase.first.second==sizeof(size_t)*2);
        //new_rl_ptrs.push_back(new_rules.size());
        new_rl_ptrs.push_back(new_bit_pos/new_r_bits);
        memcpy((char *)&pair, phrase.first.first, phrase.first.second);
        assert(pair[0]<start_sym);
        //new_rules.push_back(pair[0]);
        //new_rules.push_back(pair[1]);
        new_rule_stream.write(new_bit_pos, new_bit_pos+new_r_bits-1, pair[0]);
        new_bit_pos+=new_r_bits;
        new_rule_stream.write(new_bit_pos, new_bit_pos+new_r_bits-1, pair[1]);
        new_bit_pos+=new_r_bits;
    }

    //insert the compressed strings
    //new_rl_ptrs.push_back(new_rules.size());
    new_rl_ptrs.push_back(new_bit_pos/new_r_bits);

    for(size_t str=0;str<gram.n_strings();str++){

        auto range = gram.str2bitrange(str);
        //gram.str_boundaries[str] = new_rules.size();
        gram.str_boundaries[str] = new_bit_pos/new_r_bits;

        //prev_sym = gram.rules[range.first];
        prev_sym = gram.rule_stream.read(range.first, range.first+gram.r_bits-1);
        run_len = 1;

        for(size_t j=range.first+gram.r_bits;j<=range.second;j+=gram.r_bits){
            //curr_sym = gram.rules[j];
            curr_sym = gram.rule_stream.read(j, j+gram.r_bits-1);
            if(curr_sym!=prev_sym){
                if(run_len>1){
                    pair[0] = prev_sym;
                    pair[1] = run_len;
                    prev_sym=0;
                    auto res = ht.find((uint8_t *)&pair, sizeof(size_t)*2, prev_sym);
                    assert(res);
                }
                //new_rules.push_back(prev_sym);
                new_rule_stream.write(new_bit_pos, new_bit_pos+new_r_bits-1, prev_sym);
                new_bit_pos+=new_r_bits;
                prev_sym = curr_sym;
                run_len=0;
            }
            run_len++;
        }

        if(run_len>1){
            pair[0] = prev_sym;
            pair[1] = run_len;
            prev_sym=0;
            auto res = ht.find((uint8_t *)&pair, sizeof(size_t)*2, prev_sym);
            assert(res);
        }
        //new_rules.push_back(prev_sym);
        new_rule_stream.write(new_bit_pos, new_bit_pos+new_r_bits-1, prev_sym);
        new_bit_pos+=new_r_bits;
    }
    new_rl_ptrs.push_back(new_bit_pos/new_r_bits);
    gram.str_boundaries[gram.n_strings()]=new_bit_pos/new_r_bits;
    gram.r_bits = new_r_bits;

    assert((new_bit_pos/new_r_bits)==new_size);
    assert(gram.str2bitrange(gram.n_strings()-1).second==(new_bit_pos-new_r_bits));
    assert(new_rl_ptrs.size()==(gram.rl_ptr.size()+ht.size()));

    std::cout<<"  Stats:"<<std::endl;
    std::cout<<"    Grammar size before:        "<<gram.g<<std::endl;
    std::cout<<"    Grammar size after:         "<<new_size<<std::endl;
    std::cout<<"    Number of new nonterminals: "<<ht.size()<<std::endl;
    std::cout<<"    Compression ratio:          "<<float(new_size)/float(gram.g)<<std::endl;

    gram.run_len_nt.first = start_sym;
    gram.run_len_nt.second = ht.size();
    gram.r += ht.size();
    gram.g  = new_size;
    gram.c = gram.g - gram.str_boundaries[0];

    new_rule_stream.swap(gram.rule_stream);
    new_rl_ptrs.swap(gram.rl_ptr);
}

template<class gram_t>
void check_plain_grammar(gram_t& gram, std::string& uncomp_file) {

    std::cout<<"Checking the grammar produces the exact input string"<<std::endl;
    std::cout<<"  This step is optional and for debugging purposes"<<std::endl;
    std::cout<<"  Terminals:              "<<gram.n_terminals()<<std::endl;
    std::cout<<"  Number of nonterminals: "<<gram.n_nonterminals()<<std::endl;
    std::cout<<"  Compressed string:      "<<gram.comp_str_size()<<std::endl;

    i_file_stream<uint8_t> if_stream(uncomp_file, BUFFER_SIZE);
    std::stack<size_t> stack;
    size_t idx=0;
    std::string decompression;
    for(size_t str=0; str <gram.n_strings(); str++) {

        auto res = gram.str2bitrange(str);
        for(size_t i=res.first;i<=res.second;i+=gram.r_bits){
            stack.push(gram.bitpos2symbol(i));

            while(!stack.empty()) {
                auto curr_sym = stack.top() ;
                stack.pop();
                //std::cout<<curr_sym<<" "<<gram.r-gram.n_terminals()<<" "<<gram.rl_ptr.size()<<std::endl;
                if(gram.is_terminal(curr_sym)){
                    decompression.push_back((char)gram.get_byte_ter(curr_sym));
                }else{
                    auto res2 = gram.nt2bitrange(curr_sym);

                    if(gram.is_rl_sym(curr_sym)){
                        //std::cout<<curr_sym<<" "<<res2.first<<" "<<res2.second<<std::endl;
                        assert((res2.second-res2.first)==gram.r_bits);
                        size_t sym = gram.bitpos2symbol(res2.first);
                        size_t len = gram.bitpos2symbol(res2.second);
                        for(size_t j=0;j<len;j++) stack.emplace(sym);
                    }else{
                        for(off_t j=res2.second; j>=res2.first;j-=gram.r_bits){
                            stack.emplace(gram.bitpos2symbol(j));
                        }
                    }
                }
            }
        }

        for(char sym : decompression){
            if(sym!=(char)if_stream.read(idx)){
                std::cout<<"Error: decomp sym: "<<(int)sym<<", real sym: "<<if_stream.read(idx)<<", str: "<<str<<", position: "<<idx<<std::endl;
                assert(sym==(char)if_stream.read(idx));
            }
            idx++;
        }
        assert(if_stream.read(idx)==gram.sep_tsym);
        idx++;
        decompression.clear();
    }
    std::cout<<"Grammar is correct!!"<<std::endl;
}

template<class gram_t>
void estimate_alt_encodings(gram_t& gram){

    size_t max_sym=0;
    for(size_t i=0;i<gram.rules.size();i++){
        if(gram.rules[i]>max_sym){
            max_sym = gram.rules[i];
        }
    }
    std::vector<size_t> freqs(max_sym, 0);
    for(size_t i=0;i<gram.rules.size();i++){
        freqs[gram.rules[i]]++;
    }

    double h0 = 0;
    auto f = (double)gram.rules.size();
    for(unsigned long freq : freqs){
        if(freq!=0){
            auto f_c = (double)freq;
            h0 += (f_c/f)*log2(f/f_c);
        }
    }

    size_t lb_bytes= INT_CEIL(size_t(ceil(f*h0)), 8);
    std::cout<<"nH_0(G): "<<report_space((off_t)lb_bytes)<<std::endl;

    size_t curr_bytes = INT_CEIL(gram.rules.size()*gram.rules.width(), 8);
    std::cout<<"Current grammar space: "<<report_space((off_t)curr_bytes)<<std::endl;
    size_t tot_bits=0, n_bits;
    size_t new_sym=1;
    uint8_t code_len;
    uint8_t header_len = sym_width(sym_width(gram.r));
    for(auto & freq : freqs){
        code_len = sym_width(new_sym);
        n_bits = code_len + header_len;
        tot_bits+=n_bits*freq;
        new_sym++;
    }

    size_t tot_bytes = INT_CEIL(tot_bits, 8);
    std::cout<<"Delta codes without permutation: "<<report_space((off_t)tot_bytes)<<std::endl;

    tot_bytes=0;
    new_sym=1;
    for(auto & freq : freqs){
        size_t n_bytes = INT_CEIL(sym_width(new_sym), 8);
        n_bits = n_bytes + sym_width(new_sym);
        n_bytes = INT_CEIL(n_bits, 8);
        tot_bytes+=n_bytes*freq;
        new_sym++;
    }
    std::cout<<"Vbytes without permutation: "<<report_space((off_t)tot_bytes)<<std::endl;

    std::vector<std::pair<size_t, size_t>> sorted_freqs(gram.r);
    for(size_t i=0;i<freqs.size();i++){
        sorted_freqs[i] = {i, freqs[i]};
    }

    std::sort(sorted_freqs.begin(), sorted_freqs.end(), [](auto a, auto b){
        return a.second>b.second;
    });

    tot_bytes=0;
    new_sym=1;
    for(auto & sorted_freq : sorted_freqs){
        size_t n_bytes = INT_CEIL(sym_width(new_sym), 8);
        n_bits = n_bytes + sym_width(new_sym);
        n_bytes = INT_CEIL(n_bits, 8);
        tot_bytes+=n_bytes*sorted_freq.second;
        new_sym++;
    }
    std::cout<<"Vbytes with permutation: "<<report_space((off_t)tot_bytes)<<std::endl;
}

template<class gram_t>
std::pair<std::vector<uint8_t>, size_t> mark_disposable_symbols(const gram_t& gram) {

    //compute which nonterminals are repeated and
    // which have a replacement of length 1
    std::vector<uint8_t> rep_nts(gram.r + 1, 0);
    size_t nt;
    for(size_t rule=gram.max_tsym+1;rule<gram.r;rule++){
        auto range = gram.nt2bitrange(rule);
        if(gram.is_rl_sym((rule))){
            nt = gram.rule_stream.read(range.first, range.first+gram.r_bits-1);
            rep_nts[nt] = 2;
        }else{
            for(size_t i=range.first;i<=range.second;i+=gram.r_bits){
                nt = gram.rule_stream.read(i, i+gram.r_bits-1);
                rep_nts[nt]+=rep_nts[nt]<2;
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

template<class gram_t>
void balanced_simplification(gram_t& gram){

    assert(!gram.is_simplified);
    auto rem_syms = mark_disposable_symbols(gram);
    int_array<size_t> new_rules(gram.g-rem_syms.second, sym_width(gram.r-rem_syms.second));
    int_array<size_t> new_rl_ptrs(gram.r-rem_syms.second, sym_width((gram.g+1)-rem_syms.second));
    int_array<size_t> offsets(gram.r, sym_width(gram.r));

    size_t start_sym = gram.start_symbol();
    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++) {
        auto range = gram.nt2bitrange(sym);
        new_rl_ptrs.push_back(new_rules.size());
        if(!gram.is_rl_sym(sym)){
            for(size_t j=range.first;j<=range.second;j++){
                while(rem_syms.first[gram.rules[j]]){
                }
            }
        }
    }
}

template<class gram_t>
void simplify_grammar(gram_t& gram) {

    assert(!gram.is_simplified);

    auto rem_syms = mark_disposable_symbols(gram);
    //bitstream<size_t> new_rules(gram.g-rem_syms.second, sym_width(gram.r-rem_syms.second));

    uint8_t new_r_bits = sym_width(gram.r-rem_syms.second);
    bitstream<size_t> new_rules;
    new_rules.reserve_in_bits((gram.g-rem_syms.second)*new_r_bits);

    int_array<size_t> new_rl_ptrs(gram.r-rem_syms.second, sym_width(gram.g-rem_syms.second+1));
    int_array<size_t> offsets(gram.r, sym_width(gram.r));

    size_t del_syms=0;
    for(size_t sym=0;sym<gram.r;sym++){
        offsets[sym] = del_syms;
        del_syms+=rem_syms.first[sym];
    }

    size_t n_ter=0, new_ter, bit_pos=0;
    for(size_t ter=0;ter<=gram.max_tsym;ter++){
        if(!rem_syms.first[ter]){
            new_ter = ter-offsets[ter];
            gram.terminals[new_ter] = ter;

            new_rules.write(bit_pos, bit_pos+new_r_bits-1, new_ter);
            bit_pos+=new_r_bits;
            n_ter++;
        }
    }

    std::stack<size_t> stack;
    size_t start_sym = gram.start_symbol();
    auto tmp = gram.nt2bitrange(start_sym);
    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++) {

        if(!rem_syms.first[sym]) {

            auto range = gram.nt2bitrange(sym);
            new_rl_ptrs.push_back(bit_pos/new_r_bits);

            if(gram.is_rl_sym(sym)){
                //new_rules.push_back(gram.rules[range.first]-offsets[gram.rules[range.first]]);
                //new_rules.push_back(gram.rules[range.second]);

                size_t rl_sym = gram.rule_stream.read(range.first, range.first+gram.r_bits-1);
                size_t rl_len = gram.rule_stream.read(range.second, range.second+gram.r_bits-1);

                new_rules.write(bit_pos, bit_pos+new_r_bits-1, rl_sym-offsets[rl_sym]);
                bit_pos+=new_r_bits;
                new_rules.write(bit_pos, bit_pos+new_r_bits-1, rl_len);
                bit_pos+=new_r_bits;
            }else{

                for(size_t j=range.first;j<=range.second;j+=gram.r_bits){

                    size_t old_sym = gram.rule_stream.read(j, j+gram.r_bits-1);

                    if(rem_syms.first[old_sym]) {
                        stack.push(old_sym);
                        while(!stack.empty()){
                            auto nt = stack.top();
                            stack.pop();
                            if(rem_syms.first[nt]){
                                assert(nt>gram.max_tsym);
                                assert(!gram.is_rl_sym(nt));

                                auto range2 = gram.nt2bitrange(nt);
                                for(off_t k=range2.second;k>=range2.first;k-=gram.r_bits){
                                    stack.push(gram.rule_stream.read(k, k+gram.r_bits-1));
                                }
                            }else{
                                //new_rules.push_back(nt - offsets[nt]);
                                new_rules.write(bit_pos, bit_pos+new_r_bits-1, nt - offsets[nt]);
                                bit_pos+=new_r_bits;
                            }
                        }
                    }else{
                        //new_rules.push_back(gram.rules[j] - offsets[gram.rules[j]]);
                        new_rules.write(bit_pos, bit_pos+new_r_bits-1, old_sym - offsets[old_sym]);
                        bit_pos+=new_r_bits;
                    }
                }
            }
        }
    }

    //deal with the compressed strings
    auto range = gram.nt2bitrange(start_sym);
    size_t str=0;

    new_rl_ptrs.push_back(bit_pos/new_r_bits);
    for(size_t j=range.first;j<=range.second;j+=gram.r_bits){

        gram.str_boundaries[str++] = bit_pos/new_r_bits;

        size_t old_sym = gram.rule_stream.read(j, j+gram.r_bits-1);
        if(rem_syms.first[old_sym]) {
            stack.push(old_sym);

            while(!stack.empty()){
                auto nt = stack.top();
                stack.pop();
                if(rem_syms.first[nt]){
                    auto range2 = gram.nt2bitrange(nt);

                    for(off_t k=range2.second;k>=range2.first;k-=gram.r_bits){
                        stack.push(gram.rule_stream.read(k, k+gram.r_bits-1));
                    }
                }else{
                    //new_rules.push_back(nt - offsets[nt]);
                    new_rules.write(bit_pos, bit_pos+new_r_bits-1, nt - offsets[nt]);
                    bit_pos+=new_r_bits;
                }
            }
        }else{
            //new_rules.push_back(gram.rules[j] - offsets[gram.rules[j]]);
            new_rules.write(bit_pos, bit_pos+new_r_bits-1, old_sym - offsets[old_sym]);
            bit_pos+=new_r_bits;
        }
    }

    new_rl_ptrs.push_back(bit_pos/new_r_bits);
    gram.str_boundaries[gram.n_strings()]=bit_pos/new_r_bits;

    size_t del_nt = (rem_syms.second-(gram.max_tsym+1-n_ter));
    gram.r_bits = new_r_bits;

    assert((bit_pos/new_r_bits)==(gram.g-rem_syms.second));
    assert(gram.str2bitrange(gram.n_strings()-1).second==(bit_pos-new_r_bits));
    assert(new_rl_ptrs.size()==(gram.rl_ptr.size()-del_nt));

    size_t new_size = bit_pos/new_r_bits;
    float rm_per = float(rem_syms.second)/float(gram.r)*100;
    float comp_rat = float(new_size)/float(gram.g);

    std::cout<<"  Stats:"<<std::endl;
    std::cout<<"    Grammar size before:  "<<gram.g<<std::endl;
    std::cout<<"    Grammar size after:   "<<new_size<<std::endl;
    std::cout<<"    Deleted nonterminals: "<<rem_syms.second<<" ("<<rm_per<<"%)"<<std::endl;
    std::cout<<"    Compression ratio:    "<<comp_rat<<std::endl;

    gram.g = new_size;
    gram.r -= rem_syms.second;
    gram.c = gram.g - gram.str_boundaries[0];

    gram.rule_stream.swap(new_rules);
    gram.rl_ptr.swap(new_rl_ptrs);

    for(auto &sym : gram.lvl_rules){
        sym -= offsets[sym];
    }

    gram.run_len_nt.first -= offsets[gram.run_len_nt.first];

    gram.max_tsym = n_ter-1;
    gram.is_simplified = true;
}

void print_metadata(std::string& gram_file){
    lc_gram_t gram;
    std::ifstream ifs(gram_file, std::ios::binary);
    gram.load_metadata(ifs);
    std::cout<<"Metadata of "<<std::filesystem::path(gram_file).filename()<<std::endl;
    gram.stats(2);
}

void get_par_seed(std::string& gram_file){
    std::cout<<"Extracting the parsing functions from "<<gram_file<<std::endl;
    lc_gram_t gram;
    std::ifstream ifs(gram_file, std::ios::binary);
    gram.load_metadata(ifs);
    std::cout<<"This grammar was built using the seed "<<gram.par_seed<<std::endl;
}

/***
 *
 * @param i_file : input text file
 * @param n_threads : number of working threads
 */
template<class sym_type, class gram_type>
void build_gram(std::string &i_file, std::string& o_file, tmp_workspace & tmp_ws, size_t n_threads,
                size_t n_chunks, off_t chunk_size, size_t par_seed, bool se_build, bool skip_simp) {

    // the grammar encoding with random access support works differently, so this is a hack
    using tmp_gram_type = lc_gram_t<gram_type::has_cg_rules, gram_type::has_rl_rules, false>;

    std::cout<<"Building a locally-consistent grammar"<<std::endl;
    auto start = std::chrono::steady_clock::now();
    if(se_build){
        lc_parsing_algo<sym_type, tmp_gram_type>(i_file, o_file, tmp_ws, n_threads, n_chunks, chunk_size, par_seed);
    }else{
        lz_like_strat::lc_parsing_algo<sym_type, tmp_gram_type>(i_file, o_file, tmp_ws, n_threads, n_chunks, chunk_size, par_seed);
    }

    tmp_gram_type gram;
    load_from_file(o_file, gram);
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 2);

    /*if(gram_type::has_cg_rules){
        std::cout<<"Transforming the grammar into a college system"<<std::endl;
        start = std::chrono::steady_clock::now();
        compute_exp_len(gram);
        make_gram_fix_free(gram);
        //make_collage_system(gram);
        end = std::chrono::steady_clock::now();
        report_time(start, end, 2);
    }*/

    if(gram_type::has_rl_rules){
        std::cout<<"Run-length compressing the grammar"<<std::endl;
        start = std::chrono::steady_clock::now();
        run_length_compress(gram);
        end = std::chrono::steady_clock::now();
        report_time(start, end, 2);
    }

    if(!skip_simp){
        std::cout<<"Simplifying the grammar"<<std::endl;
        start = std::chrono::steady_clock::now();
        simplify_grammar(gram);
        end = std::chrono::steady_clock::now();
        report_time(start, end, 2);
    }

    if(gram_type::has_rand_access){
        std::cout<<"Adding random access support"<<std::endl;
        start = std::chrono::steady_clock::now();
        add_random_access_support(gram);
        end = std::chrono::steady_clock::now();
        report_time(start, end, 2);
    }

    // we need to do this final swap because the grammar
    // encoding with random access support works differently
    gram_type final_gram;
    final_gram.swap(gram);

    //optional check
    check_plain_grammar(final_gram, i_file);
    //

    std::cout<<"Stats for the final grammar:"<<std::endl;
    final_gram.breakdown(2);
    size_t written_bytes = store_to_file(o_file, final_gram);
    std::cout<<"The resulting grammar uses "+ report_space((off_t)written_bytes)+" and was stored in "<<o_file<<std::endl;
}


template<class gram_type>
void se_rand_access(std::string& gram_file, size_t str_id, size_t start, size_t end){
    gram_type gram;
    std::ifstream ifs(gram_file, std::ios::binary);
    gram.load_metadata(ifs);
    gram.load_pointers(ifs);
}
#endif //LCG_GRAMMAR_ALGORITHMS_H
