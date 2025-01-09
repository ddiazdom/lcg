//
// Created by Diaz, Diego on 14.7.2023.
//

#ifndef LCG_GRAMMAR_ALGORITHMS_H
#define LCG_GRAMMAR_ALGORITHMS_H

#include "compact_gram.h"
#include "build_gram/lc_parsing.h"
#include "cds/file_streams.hpp"
#include "cds/ts_string_map.h"

/*
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
*/

template<class gram_t>
std::tuple<size_t, uint8_t, uint8_t> compute_exp_info(gram_t& gram, std::vector<uint64_t>& exp_len, std::vector<uint8_t>& str_width){

    size_t tmp_sym, len, exp_size;
    for(size_t i=0;i<=gram.max_tsym;i++) exp_len[i] = 1;

    size_t last_sym = gram.n_lc_syms(), rule_len, n_samples, samp_bits;

    //get the expansion length of each nonterminal in the grammar
    size_t max_r_samp_bits=0, r_samp_acc_bits=0;
    for(size_t sym=gram.max_tsym+1;sym<last_sym;sym++) {
        auto range = gram.nt2bitrange(sym);
        exp_size = 0;
        for(off_t i=range.first;i<=range.second;i+=gram.r_bits){
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

    if(gram.has_run_length_rules()){
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
    }

    //compute the length of the strings
    size_t max_str_samp_bits=0, str_samp_acc_bits=0;
    for(size_t str=0;str<gram.n_strings();str++){
        auto range = gram.str2bitrange(str);
        exp_size = 0;
        for(off_t i=range.first;i<=range.second;i+=gram.r_bits){
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

    max_r_samp_bits = sym_width(max_r_samp_bits);
    // When all the string are compressed to one symbol,
    // then the strings have zero expansion samples and hence max_str_samp_bits is 0.
    // This situation happens when the grammar is not simplified, and produces inconsistency
    max_str_samp_bits = std::max<uint8_t>(1, sym_width(max_str_samp_bits));

    size_t tot_samp_bits = r_samp_acc_bits + max_r_samp_bits*(gram.r-gram.n_terminals()-1);
    tot_samp_bits += str_samp_acc_bits + max_str_samp_bits*gram.n_strings();

    return {tot_samp_bits, max_r_samp_bits, max_str_samp_bits};
}

template<log_lvl lvl_msg=INFO, class gram_t>
void add_random_access_support(gram_t& gram) {

    logger<lvl_msg>::info("Adding random access support");

    auto start = std::chrono::steady_clock::now();
    assert(!gram.has_cg_rules);

    std::vector<uint64_t> exp_tmp(gram.r, 0);
    std::vector<uint8_t> str_width(gram.s, 0);

    size_t samp_bits;
    uint8_t r_samp_bits, str_samp_bits;
    std::tie(samp_bits, r_samp_bits, str_samp_bits) = compute_exp_info(gram, exp_tmp, str_width);

    //size_t samp_space = INT_CEIL(samp_bits, 8);
    //std::cout<<"The sampling space is "<<report_space((off_t)samp_space)<<" "<<r_samp_bits<<" "<<str_samp_bits<<std::endl;

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
    size_t exp_acc, last_sym = gram.n_lc_syms(), tmp_sym, offset;
    size_t n_sampled, r_len, s_bits;
    for(size_t sym = gram.max_tsym+1; sym<last_sym; sym++){

        new_rl_ptrs.push_back(bit_pos);
        auto range = gram.nt2bitrange(sym);
        offset = bit_pos+(range.second-range.first+gram.r_bits);
        r_len = (range.second-range.first+gram.r_bits)/gram.r_bits;

        n_sampled = (r_len/gram.rl_samp_rate)+1;
        s_bits=sym_width(exp_tmp[sym]);

        exp_acc = 0;

        for(off_t i=range.first, j=1;i<=range.second;i+=gram.r_bits,j++){
            tmp_sym = gram.rule_stream.read(i, i+gram.r_bits-1);
            new_rule_stream.write(bit_pos, bit_pos+gram.r_bits-1, tmp_sym);
            bit_pos+=gram.r_bits;
            exp_acc+=exp_tmp[tmp_sym];
            if(j%gram.rl_samp_rate==0){
                new_rule_stream.write(offset, offset+s_bits-1, exp_acc);
                offset+=s_bits;
            }
        }

        assert(exp_acc==exp_tmp[sym]);
        new_rule_stream.write(offset, offset+s_bits-1, exp_acc);
        offset+=s_bits;
        new_rule_stream.write(offset, offset+r_samp_bits-1, s_bits*n_sampled);
        offset+=r_samp_bits;

        bit_pos=offset;
    }

    if(gram.has_run_length_rules()){
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

        for(off_t i=range.first, j=1;i<=range.second;i+=gram.r_bits,j++){
            tmp_sym = gram.rule_stream.read(i, i+gram.r_bits-1);
            new_rule_stream.write(bit_pos, bit_pos+gram.r_bits-1, tmp_sym);
            bit_pos+=gram.r_bits;
            exp_acc+=exp_tmp[tmp_sym];

            if(j%gram.str_samp_rate==0){
                new_rule_stream.write(offset, offset+s_bits-1, exp_acc);
                offset+=s_bits;
            }
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
    auto end = std::chrono::steady_clock::now();
    logger<lvl_msg>::info(report_time(start, end, 2));
}

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

        for(off_t j=range.first+gram.r_bits;j<=range.second;j+=gram.r_bits) {
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

        for(off_t j=range.first+gram.r_bits;j<=range.second;j+=gram.r_bits){
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

template<log_lvl msg_lvl=INFO, class gram_t>
void run_length_compress(gram_t& gram) {

    logger<msg_lvl>::info("Run-length compressing the grammar");

    auto start = std::chrono::steady_clock::now();
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

        for(off_t j=range.first+gram.r_bits;j<=range.second;j+=gram.r_bits){

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

        for(off_t j=range.first+gram.r_bits;j<=range.second;j+=gram.r_bits){
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
    assert(gram.str2bitrange(gram.n_strings()-1).second==off_t(new_bit_pos-new_r_bits));
    assert(new_rl_ptrs.size()==(gram.rl_ptr.size()+ht.size()));

    [[maybe_unused]] size_t old_size = gram.g;

    gram.run_len_nt.first = start_sym;
    gram.run_len_nt.second = ht.size();
    gram.r += ht.size();
    gram.g  = new_size;
    gram.c = gram.g - gram.str_boundaries[0];

    new_rule_stream.swap(gram.rule_stream);
    new_rl_ptrs.swap(gram.rl_ptr);
    auto end = std::chrono::steady_clock::now();

    std::string msg ="  Stats:";
    msg+="\n    Grammar size before:        "+std::to_string(old_size);
    msg+="\n    Grammar size after:         "+std::to_string(new_size);
    msg+="\n    Number of new nonterminals: "+std::to_string(ht.size());
    msg+="\n    Compression ratio:          "+std::to_string(float(new_size)/float(old_size));

    logger<msg_lvl>::debug(msg);
    logger<msg_lvl>::info(report_time(start, end, 0));
}

template<class gram_t>
void check_plain_grammar(gram_t& gram, std::vector<std::string>& uncomp_files, std::vector<text_format>& i_fmts) {

    std::cout<<"Checking the grammar produces the exact input string"<<std::endl;
    std::cout<<"  This step is optional and for debugging purposes"<<std::endl;
    std::cout<<"  Terminals:              "<<gram.n_terminals()<<std::endl;
    std::cout<<"  Number of nonterminals: "<<gram.n_nonterminals()<<std::endl;
    std::cout<<"  Tot. strings:           "<<gram.n_strings()<<std::endl;

    size_t str=0, file_pos, cmp_pos=0;
    std::stack<size_t> stack;
    std::string decompression;

    for(size_t i=0;i<uncomp_files.size();i++){

        assert(i_fmts[i]==PLAIN);
        i_file_stream<uint8_t> file_buffer(uncomp_files[i], BUFFER_SIZE);
        file_pos=0;

        while(file_pos<file_buffer.size()){

            auto res = gram.str2bitrange(cmp_pos++);

            for(off_t i=res.first;i<=res.second;i+=gram.r_bits){
                stack.push(gram.bitpos2symbol(i));
                while(!stack.empty()) {
                    auto curr_sym = stack.top() ;
                    stack.pop();
                    if(gram.is_terminal(curr_sym)){
                        decompression.push_back((char)gram.get_byte_ter(curr_sym));
                    }else{
                        auto res2 = gram.nt2bitrange(curr_sym);

                        if(gram.is_rl_sym(curr_sym)){
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

            assert(file_pos+decompression.size()<=file_buffer.size());

            size_t k=0;
            for(char sym : decompression){
                if(gram_t::has_rand_access){
                    size_t sym2 = gram.im_sym_rand_access(str, k);
                    if(sym2!=size_t(sym)){
                        std::cout<<"Error: decomp sym: "<<(int)sym<<", accessed sym: "<<sym2<<", real sym: "<<file_buffer.read(file_pos)<<", str: "<<str<<", file"<<uncomp_files[i]<<", file position"<<file_pos<<", string position: "<<k<<std::endl;
                        gram.print_parse_tree(str, true);
                        assert(size_t(sym)==sym2);
                    }
                }

                if(sym!=(char)file_buffer.read(file_pos)){
                    std::cerr<<"Error: decomp sym: "<<(int)sym<<", real sym: "<<file_buffer.read(file_pos)<<", str: "<<str<<", file: "<<uncomp_files[i]<<", file position: "<<file_pos<<" string position: "<<k<<std::endl;
                    gram.print_parse_tree(str, true);
                    exit(1);
                }
                file_pos++;
                k++;
            }

            if(file_buffer.read(file_pos)==gram.sep_tsym){
                file_pos++;
                std::cout<<"\rCorrect strings: "<<(++str)<<"/"<<gram.n_strings()<<std::flush;
            }
            decompression.clear();
        }
    }
    assert(cmp_pos==gram.n_strings());
    std::cout<<"\nGrammar is correct!!"<<std::endl;
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
            for(off_t i=range.first;i<=range.second;i+=gram.r_bits){
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

template<log_lvl lvl_msg=INFO, class gram_t>
void simplify_grammar(gram_t& gram) {

    logger<lvl_msg>::info("Simplifying the grammar");

    assert(!gram.is_simplified);
    [[maybe_unused]] auto start = std::chrono::steady_clock::now();

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
    //auto tmp = gram.nt2bitrange(start_sym);
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

                for(off_t j=range.first;j<=range.second;j+=gram.r_bits){

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
    for(off_t j=range.first;j<=range.second;j+=gram.r_bits){

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
    assert(gram.str2bitrange(gram.n_strings()-1).second==off_t(bit_pos-new_r_bits));
    assert(new_rl_ptrs.size()==(gram.rl_ptr.size()-del_nt));

    size_t new_size = bit_pos/new_r_bits;
    [[maybe_unused]] size_t old_size = gram.g;

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
    [[maybe_unused]] auto end = std::chrono::steady_clock::now();

    float rm_per = float(rem_syms.second)/float(gram.r)*100;
    float comp_rat = float(new_size)/float(old_size);

    std::string msg ="  Stats:";
    msg+="\n    Grammar size before:  "+std::to_string(old_size);
    msg+="\n    Grammar size after:   "+std::to_string(new_size);
    msg+="\n    Deleted nonterminals: "+std::to_string(rem_syms.second)+" ("+std::to_string(rm_per)+"%)";
    msg+="\n    Compression ratio:    "+std::to_string(comp_rat);

    logger<lvl_msg>::debug(msg);
    logger<lvl_msg>::info(report_time(start, end, 0));
}

std::tuple<bool, bool, bool> read_grammar_flags(std::string& file){
    std::ifstream ifs(file, std::ios::binary);
    bool tmp_has_rl_rules, tmp_has_cg_rules, tmp_has_rand_rules;
    load_elm(ifs, tmp_has_rl_rules);
    load_elm(ifs, tmp_has_cg_rules);
    load_elm(ifs, tmp_has_rand_rules);
    return {tmp_has_rl_rules, tmp_has_cg_rules, tmp_has_rand_rules};
}

void print_metadata(std::string& gram_file){

    bool has_rl_rules, has_cg_rules, has_rand_access;
    std::tie(has_rl_rules, has_cg_rules, has_rand_access) = read_grammar_flags(gram_file);
    std::ifstream ifs(gram_file, std::ios::binary);

    if(has_cg_rules){
        if(has_rl_rules){
            lc_gram_t<true, true, true> gram;
            gram.load_metadata(ifs);
            std::cout<<"Metadata of "<<std::filesystem::path(gram_file).filename()<<std::endl;
            std::cout<<gram.stats(2)<<std::endl;
        }else{
            lc_gram_t<true, false, true> gram;
            gram.load_metadata(ifs);
            std::cout<<"Metadata of "<<std::filesystem::path(gram_file).filename()<<std::endl;
            std::cout<<gram.stats(2)<<std::endl;
        }
    }else{
        if(has_rl_rules){
            lc_gram_t<false, true, true> gram;
            gram.load_metadata(ifs);
            std::cout<<"Metadata of "<<std::filesystem::path(gram_file).filename()<<std::endl;
            std::cout<<gram.stats(2)<<std::endl;
        }else{
            lc_gram_t<false, false, true> gram;
            gram.load_metadata(ifs);
            std::cout<<"Metadata of "<<std::filesystem::path(gram_file).filename()<<std::endl;
            std::cout<<gram.stats(2)<<std::endl;
        }
    }
}

void get_par_seed(std::string& gram_file){
    std::cout<<"Extracting the parsing functions from "<<gram_file<<std::endl;
    lc_gram_t gram;
    std::ifstream ifs(gram_file, std::ios::binary);
    gram.load_metadata(ifs);
    std::cout<<"This grammar was built using the seed "<<gram.par_seed<<std::endl;
}

template<log_lvl msg_lvl=INFO, class gram_type>
void complete_and_pack_grammar(plain_gram& p_gram, gram_type& new_gram){

    //TODO Reorder the elements in the stream so random access is more disk-friendly
    //This order is better: compressed string, and then the rules in decreasing level
    //Remember to pad some bits to byte-align each level, this will simplify random access
    logger<msg_lvl>::info("Compacting the grammar");

    auto start = std::chrono::steady_clock::now();
    new_gram.n = p_gram.txt_size();
    new_gram.r = p_gram.gram_alphabet();
    new_gram.g = p_gram.gram_size();
    new_gram.s = p_gram.tot_strings();
    new_gram.c = p_gram.tot_strings();

    new_gram.max_tsym = size_t(std::numeric_limits<uint8_t>::max());
    new_gram.sep_tsym = p_gram.sep_sym();

    size_t gram_lvl = p_gram.tot_rounds();
    new_gram.lvl_rules.reserve(gram_lvl);
    size_t acc=new_gram.max_tsym+1, tmp;

    for(size_t i=0;i<gram_lvl;i++){
        tmp = p_gram.phrases_round(i);//the rounds are zero-based in p_gram
        new_gram.lvl_rules.push_back(acc);
        acc+=tmp;
    }
    assert((acc+1)==new_gram.r);//the +1 is because of the grammar's start symbol
    new_gram.lvl_rules.push_back(acc);

    size_t bit_pos = 0;
    off_t r_bits = sym_width(new_gram.r);
    new_gram.r_bits = r_bits;

    // insert the terminals into the stream.
    // I don't remember why, but there was a reason for this.
    new_gram.terminals.resize(new_gram.max_tsym+1);
    new_gram.rule_stream.reserve_in_bits(r_bits*(new_gram.max_tsym+1));
    for(size_t ter=0;ter<=new_gram.max_tsym;ter++){
        new_gram.rule_stream.write(bit_pos, bit_pos+r_bits-1, ter);
        bit_pos+=r_bits;
        new_gram.terminals[ter] = ter;
    }

    size_t sym, rule=0, rule_start_ptr=bit_pos;

    //insert the first level
    auto it = p_gram.ter_dict.begin();
    auto it_end = p_gram.ter_dict.end();
    new_gram.rule_stream.reserve_in_bits(new_gram.rule_stream.capacity_in_bits() +
                                         (p_gram.ter_dict.tot_symbols()*r_bits));
    new_gram.rl_ptr.set_width(sym_width(new_gram.g));
    new_gram.rl_ptr.resize(p_gram.ter_dict.size());
    while(it!=it_end){
        auto phr = *it;
        for(size_t i=0;i<phr.len;i++){
            new_gram.rule_stream.write(bit_pos, bit_pos+r_bits-1, phr.phrase[i]);
            bit_pos+=r_bits;
        }
        new_gram.rl_ptr.write(rule++, rule_start_ptr/r_bits);
        rule_start_ptr = bit_pos;
        ++it;
    }
    size_t min_sym = new_gram.max_tsym;
    size_t max_sym = min_sym+p_gram.ter_dict.size();
    p_gram.ter_dict.destroy();

    uint8_t lvl=0;
    for(auto &nt_dict : p_gram.nt_dicts) {

        auto it_nt = nt_dict.begin();
        auto it_nt_end = nt_dict.end();

        new_gram.rule_stream.reserve_in_bits(new_gram.rule_stream.capacity_in_bits() +
                                             (nt_dict.tot_symbols()*r_bits));
        new_gram.rl_ptr.resize(new_gram.rl_ptr.size()+nt_dict.size());

        while(it_nt!=it_nt_end){
            auto phr = *it_nt;
            for(size_t i=0;i<phr.len;i++){
                sym = min_sym+phr.phrase[i];
                assert(sym==(phr.phrase[i]+new_gram.lvl_rules[lvl]-1));
                new_gram.rule_stream.write(bit_pos, bit_pos+r_bits-1, sym);
                bit_pos+=r_bits;
            }
            new_gram.rl_ptr.write(rule++, rule_start_ptr/r_bits);
            rule_start_ptr = bit_pos;
            ++it_nt;
        }
        min_sym=max_sym;
        max_sym+=nt_dict.size();
        nt_dict.destroy();
        lvl++;
    }

    //allocate space for the sequence of compressed strings
    new_gram.rule_stream.reserve_in_bits(new_gram.rule_stream.capacity_in_bits() + (p_gram.comp_string.size()*r_bits));

    //the +2 is because we need to record the start and end+1
    // where the sequence of compressed strings lies in the string
    new_gram.rl_ptr.resize(new_gram.rl_ptr.size()+2);

    //insert the sequence of compressed strings. It is a bit
    // tricky because in the current grammar's construction
    // algorithm, symbols for the compressed strings can have
    // different levels.
    size_t pos =0;
    for(auto const& str_subset: p_gram.str_orders){
        for(size_t j=0;j<str_subset.n_strings;j++){
            sym = new_gram.lvl_rules[str_subset.lvl-1]+p_gram.comp_string[pos]-1;
            new_gram.rule_stream.write(bit_pos, bit_pos+r_bits-1, sym);
            bit_pos+=r_bits;
            pos++;
        }
    }
    new_gram.rl_ptr.write(rule++, rule_start_ptr/r_bits);
    destroy(p_gram.comp_string);
    destroy(p_gram.str_orders);
    //

    assert((bit_pos/r_bits)==new_gram.g);
    assert(rule==(new_gram.r-(new_gram.max_tsym+1)));
    new_gram.rl_ptr.write(rule, bit_pos/r_bits);

    //store a pointer to the set of strings. I might be useless to store a pointer for each string given
    // their compressed size is one, but we do it this way because later we can simplify the grammar. In that case,
    // the compressed representation of each string is a substring of variable length
    size_t offset = new_gram.g-new_gram.c;
    new_gram.str_boundaries.resize(new_gram.s+1);
    for(size_t str=0;str<=new_gram.s;str++){
        new_gram.str_boundaries[str] = offset+str;
    }
    assert(new_gram.str_boundaries[0]==offset);
    assert(new_gram.str_boundaries.back()==(bit_pos/r_bits));
    auto end = std::chrono::steady_clock::now();
    logger<msg_lvl>::info(report_time(start, end, 0));
}

template<class gram_type, log_lvl lvl_msg=INFO>
void build_gram(std::vector<std::string> i_files, std::vector<text_format>& i_fmts, std::string& o_file,
                size_t n_threads, off_t chunk_size, float i_frac, bool skip_simp, bool par_gram, bool check_gram,
                bool semi_ext, std::string& parent_tmp_folder) {

    // create a temporary folder inside parent_tmp_folder
    // the folder will be created only if semi_ext is true, otherwise
    // it will be just ignored
    tmp_workspace wrk_space(parent_tmp_folder, true, "lcg", semi_ext);

    if(semi_ext){
        logger<lvl_msg>::info("\nRunning in semi external mode. Temporary folder is \""+wrk_space.tmp_folder+"\"");
    }

    logger<lvl_msg>::info("\nCompressing "+std::to_string(i_files.size())+" files");
    // We set the number of grammar level to 40.
    // This limit is enough to process strings of up to 4TB in length.
    // In any case, the combined length of the strings in the collection can be arbitrarily large.
    plain_gram p_gram(40, '\n');

    build_lc_gram<lvl_msg>(i_files, i_fmts, p_gram, n_threads, chunk_size, i_frac, semi_ext, wrk_space.tmp_folder);

    if(par_gram){
        //TODO uncomment the line below
        //store_to_file(o_file, sink_gram);
        logger<lvl_msg>::info("The resulting grammar was stored in "+o_file);
        return;
    }

    // the grammar encoding with random access support works differently, so this hack
    using tmp_gram_type = lc_gram_t<gram_type::has_cg_rules, gram_type::has_rl_rules, false>;
    tmp_gram_type compact_gram;

    complete_and_pack_grammar<lvl_msg>(p_gram, compact_gram);

    //gram.print_parse_tree(0, true);
    /*if(gram_type::has_cg_rules){
        compute_exp_len(gram);
        make_gram_fix_free(gram);
        //make_collage_system(gram);
    }*/

    if(gram_type::has_rl_rules){
        run_length_compress<lvl_msg>(compact_gram);
    }

    if(!skip_simp){
        simplify_grammar<lvl_msg>(compact_gram);
    }

    if(gram_type::has_rand_access){
        add_random_access_support<lvl_msg>(compact_gram);
    }

    // we need to do this final swap because the grammar
    // encoding with random access support works differently
    gram_type final_gram;
    final_gram.swap(compact_gram);

    //optional check
    if(check_gram) {
        check_plain_grammar(final_gram, i_files, i_fmts);
    }
    //

    logger<lvl_msg>::debug("Stats for the final grammar:\n"+final_gram.breakdown(2));

    size_t written_bytes = store_to_file(o_file, final_gram);
    logger<lvl_msg>::info("The resulting grammar uses "+report_space((off_t)written_bytes)+" and was stored in "+o_file);
}
#endif //LCG_GRAMMAR_ALGORITHMS_H
