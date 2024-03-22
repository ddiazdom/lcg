//
// Created by Diaz, Diego on 14.7.2023.
//

#ifndef LCG_GRAMMAR_ALGORITHMS_H
#define LCG_GRAMMAR_ALGORITHMS_H

#include "semi-external-strategy/lc_parsing.h"
#include "lz-like-strategy/lc_parsing.h"
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
                auto res = gram.nt2phrase(nt);
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
                        auto res = gram.nt2phrase(nt);
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
            size_t new_nt = 0;
            for(auto const& context : ext_contexts){
                size_t left = context>>width;
                size_t right = context & ((1<<width)-1);

                size_t context_lvl = gram.parsing_level(left);

                //descend over the adjacent branches
                size_t top_right = right;
                while(context_lvl>(i+1)){
                    auto res = gram.nt2phrase(left);
                    left = gram.rules[res.second];
                    auto res2 = gram.nt2phrase(right);
                    right = gram.rules[res2.first];

                    left_branch.push(left);
                    right_branch.push(right);
                    context_lvl--;
                }

                //create the new phrase gram.rules[res1.first..res1.second]{\cdot}right
                std::vector<uint32_t> new_phrase;
                auto res1 = gram.nt2phrase(left);
                for(size_t u=res1.first;u<res1.second;u++){
                    new_phrase.push_back(gram.rules[u]);
                }
                auto res2 = gram.nt2phrase(right);
                new_phrase.push_back(gram.rules[res2.first]);
                //TODO hash new_phrase and assign it a nonterminal value

                new_phrase.clear();
                //ascend back to do the renaming
                while(!right_branch.empty()){
                    left = left_branch.top();
                    left_branch.pop();

                    res1 = gram.nt2phrase(left);
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
void compute_exp_len(gram_t& gram){

    int_array<size_t> exp(gram.r, sym_width(gram.longest_str));
    size_t tmp_sym, len, exp_size;
    size_t last_sym;

    for(size_t i=0;i<=gram.max_tsym;i++){
        exp[i] = 1;
    }

    if constexpr (gram_t::has_rl_rules){
        last_sym = gram.first_rl_sym();
    }else{
        last_sym = gram.r-1;
    }

    //get the expansion length of each nonterminal in the grammar
    for(size_t sym=gram.max_tsym+1;sym<last_sym;sym++) {
        auto range = gram.nt2phrase(sym);
        exp_size = 0;
        for(size_t i=range.first;i<=range.second;i++){
            tmp_sym = gram.rules.read(i);
            len = 1;
            if(gram.is_rl_sym(tmp_sym)) {
                auto range2 = gram.nt2phrase(tmp_sym);
                tmp_sym = gram.rules.read(range2.first);
                len = gram.rules.read(range2.second);
            }

            if(gram.is_terminal(tmp_sym)){
                exp_size += len;
            } else{
                exp_size += exp.read(tmp_sym)*len;
            }
        }
        exp.write(sym, exp_size);
    }

    if constexpr (gram_t::has_rl_rules){
        last_sym = gram.last_rl_sym();
        for(size_t sym=gram.first_rl_sym();sym<=last_sym; sym++) {
            auto range = gram.nt2phrase(sym);
            tmp_sym = gram.rules.read(range.first);
            len = gram.rules.read(range.second);
            if(gram.is_terminal(tmp_sym)){
                exp_size = len;
            }else{
                exp_size = exp.read(tmp_sym)*len;
            }
            exp.write(sym, exp_size);
        }
    }
    gram.rule_exp.swap(exp);
}

template<class gram_t>
void add_random_access_support(gram_t& gram){

    int_array<size_t> exp(gram.n_nonterminals(), sym_width(gram.longest_str));
    size_t tmp_sym, len, exp_size;
    size_t last_sym = gram.first_rl_sym();

    //get the expansion length of each nonterminal in the grammar
    for(size_t sym=gram.max_tsym+1;sym<last_sym;sym++) {

        auto range = gram.nt2phrase(sym);
        exp_size = 0;
        for(size_t i=range.first;i<=range.second;i++){
            tmp_sym = gram.rules.read(i);
            len = 1;
            if(gram.is_rl_sym(tmp_sym)) {
                auto range2 = gram.nt2phrase(tmp_sym);
                tmp_sym = gram.rules.read(range2.first);
                len = gram.rules.read(range2.second);
            }

            if(gram.is_terminal(tmp_sym)){
                exp_size += len;
            } else{
                exp_size += exp.read(tmp_sym-gram.max_tsym)*len;
            }
        }
        exp.write(sym-gram.max_tsym, exp_size);
    }

    last_sym = gram.last_rl_sym();
    for(size_t sym=gram.first_rl_sym();sym<=last_sym; sym++) {
        auto range = gram.nt2phrase(sym);
        tmp_sym = gram.rules.read(range.first);
        len = gram.rules.read(range.second);
        if(gram.is_terminal(tmp_sym)){
            exp_size = len;
        }else{
            exp_size = exp.read(tmp_sym-gram.max_tsym)*len;
        }
        exp.write(sym-gram.max_tsym, exp_size);
    }

    size_t n_samples = INT_CEIL(gram.g-gram.max_tsym, gram.samp_rate);
    int_array<size_t> samp_exp(n_samples, sym_width(gram.longest_str));
    samp_exp.resize(n_samples);

    //store the sampled elements
    size_t first_sym = gram.max_tsym+1, exp_acc;
    last_sym = gram.last_rl_sym();
    for(size_t sym = first_sym; sym<=last_sym; sym++){
        auto range = gram.nt2phrase(sym);
        exp_acc = 0;
        for(size_t i=range.first, j=1;i<range.second;i++,j++){
            if(j%gram.samp_rate==0){
                samp_exp.push_back(exp_acc);
            }
            tmp_sym = gram.rules.read(i);
            if(gram.is_terminal(tmp_sym)){
                exp_acc+=1;
            }else{
                exp_acc+=exp.read(tmp_sym-gram.max_tsym);
            }
        }
    }

    for(size_t str=0;str<gram.n_strings();str++){
        auto range = gram.str2phrase(str);
        exp_acc = 0;
        for(size_t i=range.first, j=1;i<range.second;i++,j++){
            if(j%gram.samp_rate==0){
                samp_exp.push_back(exp_acc);
            }
            tmp_sym = gram.rules.read(i);
            if(gram.is_terminal(tmp_sym)){
                exp_acc+=1;
            }else{
                exp_acc+=exp.read(tmp_sym-gram.max_tsym);
            }
        }
    }

    exp.resize(exp.size());
    samp_exp.resize(samp_exp.size());

    gram.rule_exp.swap(exp);
    gram.sampled_exp.swap(samp_exp);
    gram.has_rand_access = true;
}

//TODO implement this
/*
void merge_grammars(lc_gram_t& gram_a, lc_gram_t& gram_b){

}*/

template<class gram_t>
size_t get_new_rl_rules(gram_t& gram, par_string_map<size_t>& ht) {

    size_t prev_sym, curr_sym, run_len, tmp_sym;
    size_t pair[2] = {0};
    size_t start_sym = gram.start_symbol();
    size_t new_id = start_sym;
    size_t new_size = gram.n_terminals();

    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++){
        auto range = gram.nt2phrase(sym);
        prev_sym = gram.rules[range.first];
        run_len = 1;

        for(size_t j=range.first+1;j<=range.second;j++) {
            curr_sym = gram.rules[j];
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
        auto range = gram.str2phrase(str);
        prev_sym = gram.rules[range.first];
        run_len = 1;

        for(size_t j=range.first+1;j<=range.second;j++){
            curr_sym = gram.rules[j];
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

    int_array<size_t> new_rules(new_size, sym_width(gram.r+ht.size()));
    int_array<size_t> new_rl_ptrs(sym_width(gram.r+ht.size()), sym_width(new_size));
    size_t start_sym = gram.start_symbol(), prev_sym, run_len, curr_sym;
    size_t pair[2]={0};

    //insert terminals
    for(size_t i=0;i<gram.n_terminals();i++){
        new_rules.push_back(i);
    }

    //insert regular rules
    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++) {

        new_rl_ptrs.push_back(new_rules.size());

        auto range = gram.nt2phrase(sym);
        prev_sym = gram.rules[range.first];
        run_len = 1;

        for(size_t j=range.first+1;j<=range.second;j++){

            curr_sym = gram.rules[j];
            if(curr_sym!=prev_sym){
                if(run_len>1){
                    pair[0] = prev_sym;
                    pair[1] = run_len;
                    prev_sym = 0;
                    auto res = ht.find((uint8_t *)&pair, sizeof(size_t)*2, prev_sym);
                    assert(res);
                }

                new_rules.push_back(prev_sym);
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
        new_rules.push_back(prev_sym);
    }
    assert(new_rl_ptrs.size()==gram.n_nonterminals()-1);

    //insert the new rl rules
    for(auto const& phrase : ht){
        assert(phrase.first.second==sizeof(size_t)*2);
        new_rl_ptrs.push_back(new_rules.size());
        memcpy((char *)&pair, phrase.first.first, phrase.first.second);
        assert(pair[0]<start_sym);
        new_rules.push_back(pair[0]);
        new_rules.push_back(pair[1]);
    }

    //insert the compressed strings
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
                    pair[0] = prev_sym;
                    pair[1] = run_len;
                    prev_sym=0;
                    auto res = ht.find((uint8_t *)&pair, sizeof(size_t)*2, prev_sym);
                    assert(res);
                }
                new_rules.push_back(prev_sym);
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
        new_rules.push_back(prev_sym);
    }
    new_rl_ptrs.push_back(new_rules.size());
    gram.str_boundaries[gram.n_strings()]=new_rules.size();

    assert(gram.str2phrase(gram.n_strings()-1).second==new_rules.size()-1);
    assert(new_rules.size()==new_size);
    assert(new_rl_ptrs.size()==(gram.rl_ptr.size()+ht.size()));

    std::cout<<"  Stats:"<<std::endl;
    std::cout<<"    Grammar size before:        "<<gram.rules.size()<<std::endl;
    std::cout<<"    Grammar size after:         "<<new_rules.size()<<std::endl;
    std::cout<<"    Number of new nonterminals: "<<ht.size()<<std::endl;
    std::cout<<"    Compression ratio:          "<<float(new_rules.size())/float(gram.rules.size())<<std::endl;

    gram.run_len_nt.first = start_sym;
    gram.run_len_nt.second = ht.size();
    gram.r += ht.size();
    gram.g  = new_rules.size();
    gram.c = gram.g - gram.str_boundaries[0];

    new_rules.swap(gram.rules);
    new_rl_ptrs.swap(gram.rl_ptr);

    if(gram.has_rand_access){
        add_random_access_support(gram);
    }
}

template<class gram_t>
void check_plain_grammar(gram_t& gram, std::string& uncomp_file) {

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
                std::string rand_file = "error_pf_functions.pf";
                store_pl_vector(rand_file, gram.par_functions);
                std::cout<<"Error: decomp sym: "<<(int)sym<<", real sym: "<<if_stream.read(idx)<<", str: "<<str<<", position: "<<idx<<std::endl;
                std::cout<<"The parsing functions were stored in the file "<<rand_file<<" for debugging"<<std::endl;
                assert(sym==(char)if_stream.read(idx));
            }
            idx++;
        }
        if(gram.str_boundaries[str+1]==(i+1)){
            idx++;
            str++;
        }
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

    /*auto res1 = gram.nt2phrase(gram.first_rl_sym());
    auto res2 = gram.nt2phrase(gram.last_rl_sym());
    for(size_t i=0;i<res1.first;i++){
        freqs[gram.rules[i]]++;
    }

    for(size_t i=res1.first;i<=res2.second;i+=2){
        freqs[gram.rules[i]]++;
    }

    for(size_t i=res2.second+1;i<gram.rules.size();i++){
        freqs[gram.rules[i]]++;
    }*/

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

template<class gram_t>
void balanced_simplification(gram_t& gram){

    assert(!gram.is_simplified);
    auto rem_syms = mark_disposable_symbols(gram);
    int_array<size_t> new_rules(gram.g-rem_syms.second, sym_width(gram.r-rem_syms.second));
    int_array<size_t> new_rl_ptrs(gram.r-rem_syms.second, sym_width((gram.g+1)-rem_syms.second));
    int_array<size_t> offsets(gram.r, sym_width(gram.r));

    size_t start_sym = gram.start_symbol();
    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++) {
        auto range = gram.nt2phrase(sym);
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
    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++) {

        if(!rem_syms.first[sym]) {

            auto range = gram.nt2phrase(sym);
            new_rl_ptrs.push_back(new_rules.size());

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

    //deal with the compressed strings
    auto range = gram.nt2phrase(start_sym);
    size_t str=0;
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
    gram.str_boundaries[gram.n_strings()]=new_rules.size();

    size_t del_nt = (rem_syms.second-(gram.max_tsym+1-n_ter));
    assert(new_rules.size()==gram.rules.size()-rem_syms.second);
    assert(new_rl_ptrs.size()==(gram.rl_ptr.size()-del_nt));
    assert(gram.str2phrase(gram.n_strings()-1).second==new_rules.size()-1);

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
        sym -= offsets[sym];
    }

    gram.run_len_nt.first -= offsets[gram.run_len_nt.first];

    gram.max_tsym = n_ter-1;
    gram.is_simplified = true;

    if(gram.has_rand_access){
        add_random_access_support(gram);
    }
}

void print_metadata(std::string& gram_file){
    lc_gram_t gram;
    std::ifstream ifs(gram_file, std::ios::binary);
    gram.load_metadata(ifs);
    gram.stats(2);
}

void get_par_functions(std::string& gram_file, std::string& output_file){
    std::cout<<"Extracting the parsing functions from "<<gram_file<<std::endl;
    lc_gram_t gram;
    std::ifstream ifs(gram_file, std::ios::binary);
    gram.load_metadata(ifs);
    store_pl_vector(output_file, gram.par_functions);
    std::cout<<gram.par_functions.size()<<" functions were extracted from the grammar"<<std::endl;
    std::cout<<"The functions (PF format) are in "<<gram_file<<std::endl;
}

/***
 *
 * @param i_file : input text file
 * @param n_threads : number of working threads
 */
template<class sym_type, class gram_type>
void gram_algo(std::string &i_file, std::string& pf_file, std::string& o_file, tmp_workspace & tmp_ws, size_t n_threads, size_t n_chunks, off_t chunk_size, bool se_p_rounds){

    std::cout<<"Building a locally-consistent grammar"<<std::endl;
    auto start = std::chrono::steady_clock::now();
    if(se_p_rounds){
        lc_parsing_algo<sym_type, gram_type>(i_file, pf_file, o_file, tmp_ws, n_threads, n_chunks, chunk_size);
    }else{
        lz_like_strat::lc_parsing_algo<sym_type>(i_file, pf_file, o_file, tmp_ws, n_threads, n_chunks, chunk_size);
    }

    gram_type gram;
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

    //if(gram_type::has_rl_rules){
    //    std::cout<<"Run-length compressing the grammar"<<std::endl;
    //    start = std::chrono::steady_clock::now();
    //    run_length_compress(gram);
    //    end = std::chrono::steady_clock::now();
    //    report_time(start, end, 2);
    //}

    //std::cout<<"Simplifying the grammar"<<std::endl;
    //start = std::chrono::steady_clock::now();
    //simplify_grammar(gram);
    //end = std::chrono::steady_clock::now();
    //report_time(start, end, 2);

    //estimate_alt_encodings(gram);
    //std::cout<<"Adding random access support"<<std::endl;
    //start = std::chrono::steady_clock::now();
    //add_random_access_support(gram);
    //end = std::chrono::steady_clock::now();
    //report_time(start, end, 2);

    //optional check
    //check_plain_grammar(gram, i_file);
    //

    std::cout<<"Stats for the final grammar:"<<std::endl;
    gram.breakdown(2);
    size_t written_bytes = store_to_file(o_file, gram);
    std::cout<<"The resulting grammar uses "+ report_space((off_t)written_bytes)+" and was stored in "<<o_file<<std::endl;
}
#endif //LCG_GRAMMAR_ALGORITHMS_H
