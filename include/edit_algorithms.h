//
// Created by Diaz, Diego on 11.5.2024.
//

#ifndef LCG_EDITION_ALGORITHMS_H
#define LCG_EDITION_ALGORITHMS_H

#include "grammar.h"
#include "grammar_algorithms.h"
#include "unordered_map"

struct rule_type{
    off_t nt{};
    uint8_t lvl{};
    bool exp_branch=false;
    std::vector<size_t> rhs;

    rule_type()=default;
    rule_type(rule_type&& other) noexcept {
        std::swap(nt, other.nt);
        std::swap(lvl, other.lvl);
        std::swap(exp_branch, other.exp_branch);
        rhs.swap(other.rhs);
    }

    rule_type(off_t nt_, uint8_t lvl_, bool exp_branch_, std::vector<size_t>& rhs_)  {
        nt = nt_;
        lvl = lvl_;
        exp_branch = exp_branch_;
        rhs = rhs_;
    }

    rule_type(off_t nt_, uint8_t lvl_, bool exp_branch_, std::vector<size_t>&& rhs_)  {
        nt = nt_;
        lvl = lvl_;
        exp_branch = exp_branch_;
        rhs.swap(rhs_);
    }
};

struct new_rule_type{
    off_t nt{};
    uint64_t fp{};
    uint8_t lvl{};
    std::vector<size_t> rhs;

    new_rule_type()=default;

    new_rule_type(new_rule_type&& other) noexcept {
        std::swap(nt, other.nt);
        std::swap(fp, other.fp);
        std::swap(lvl, other.lvl);
        rhs.swap(other.rhs);
    }

    new_rule_type(off_t nt_, uint64_t& fp_, uint8_t lvl_, std::vector<size_t>& rhs_)  {
        nt = nt_;
        fp = fp_;
        lvl = lvl_;
        rhs = rhs_;
    }

    new_rule_type(off_t nt_, uint64_t fp_, uint8_t lvl_, std::vector<size_t>&& rhs_)  {
        nt = nt_;
        fp = fp_;
        lvl = lvl_;
        rhs.swap(rhs_);
    }

    new_rule_type& operator=(new_rule_type&& other) noexcept {
        std::swap(nt, other.nt);
        std::swap(fp, other.fp);
        std::swap(lvl, other.lvl);
        rhs.swap(other.rhs);
        return *this;
    }
};

template<class gram_type>
void retrieve_left_offset(std::stack<std::tuple<off_t, off_t, bool>>& off_stack, gram_type& gram, uint8_t r_bits, std::vector<size_t>& offset_seq){
    offset_seq.clear();
    auto left_side = off_stack.top();
    off_stack.pop();
    if(std::get<2>(left_side)){//the sequence is in run-length format
        size_t len = (std::get<1>(left_side)-std::get<0>(left_side)+r_bits)/r_bits;
        size_t tmp_sym = gram.bitpos2symbol(std::get<0>(left_side));
        for(size_t i=0;i<len;i++){
            offset_seq.push_back(tmp_sym);
        }
        if(!off_stack.empty()){
            left_side = off_stack.top();
            off_stack.pop();
            assert(!std::get<2>(left_side));
            gram.get_rhs_suffix(std::get<0>(left_side), std::get<1>(left_side), true, offset_seq);
        }
    } else {
        gram.get_rhs_suffix(std::get<0>(left_side), std::get<1>(left_side), true, offset_seq);
    }
    std::reverse(offset_seq.begin(), offset_seq.end());
}

template<class gram_type>
void retrieve_right_offset(std::stack<std::tuple<off_t, off_t, bool>>& off_stack, gram_type& gram, uint8_t r_bits, std::vector<size_t>& offset_seq){
    offset_seq.clear();
    auto right_side = off_stack.top();
    off_stack.pop();
    if(std::get<2>(right_side)){
        size_t len = std::get<1>(right_side)/r_bits;
        size_t tmp_sym = gram.rule_stream.read(std::get<0>(right_side), std::get<0>(right_side)+r_bits-1);
        for(size_t i=0;i<len;i++){
            offset_seq.push_back(tmp_sym);
        }
        if(!off_stack.empty()){
            right_side = off_stack.top();
            off_stack.pop();
            assert(!std::get<2>(right_side));
            gram.get_rhs_prefix(std::get<0>(right_side), std::get<1>(right_side), true, offset_seq);
        }
    } else {
        gram.get_rhs_prefix(std::get<0>(right_side), std::get<1>(right_side), true, offset_seq);
    }
}

off_t find_parsing_break (std::vector<size_t> &left_side, std::vector<size_t> &right_side,
                          std::vector<uint64_t> &all_fps,
                          std::unordered_map<uint64_t, std::pair<std::vector<size_t>, uint64_t>> &new_gram_rules) {

    if (left_side.empty () || right_side.empty ()) return -1;

    size_t i = 1;
    size_t lm_sym = right_side[ 0 ];
    while (i < right_side.size () && lm_sym == right_side[ i ]) i++;

    uint64_t lm_fp;
    uint64_t lm_fp_next;

    if (lm_sym > all_fps.size ()) {
        lm_fp = new_gram_rules[ lm_sym ].second;
    } else {
        lm_fp = all_fps[ lm_sym ];
    }
    lm_fp_next = lm_fp - 1;
    if (i < right_side.size ()) {
        if (right_side[ i ] >= all_fps.size ()) {
            lm_fp_next = new_gram_rules[ right_side[ i ]].second;
        } else {
            lm_fp_next = all_fps[ right_side[ i ]];
        }
    }

    uint64_t rm_fp;
    if (left_side.back () >= all_fps.size ()) {
        rm_fp = new_gram_rules[ left_side.back () ].second;
    } else {
        rm_fp = all_fps[ left_side.back () ];
    }
    bool is_local_minimum = rm_fp > lm_fp && lm_fp < lm_fp_next;
    return (-1) + is_local_minimum;
}

template<class gram_t>
void insert_left_offset(std::stack<rule_type>& offset_stack, exp_data& rhs_data, gram_t& gram, uint8_t lvl, bool& lm_branch){
    std::vector<size_t> rhs;

    if(rhs_data.is_rl){
        size_t sym = gram.bitpos2symbol(rhs_data.bp_rhs_s);
        off_t len = (rhs_data.bp_exp_s-rhs_data.bp_rhs_s)/gram.r_bits;
        for(off_t i=0;i<len;i++){
            rhs.push_back(sym);
        }
    }else{
        for(off_t i=rhs_data.bp_rhs_s;i<=(rhs_data.bp_exp_s-gram.r_bits);i+=gram.r_bits){
            rhs.push_back(gram.bitpos2symbol(i));
        }
    }

    lm_branch = lm_branch && rhs_data.bp_rhs_s==rhs_data.bp_exp_s;

    if(!rhs.empty() || lm_branch){
        offset_stack.emplace(-1, lvl, !rhs.empty(), rhs);
    }
}

template<class gram_t>
void insert_right_offset(std::stack<rule_type>& offset_stack, exp_data& rhs_data, gram_t& gram, uint8_t lvl, bool& rm_branch){
    std::vector<size_t> rhs;
    if(rhs_data.is_rl){
        size_t sym = gram.bitpos2symbol(rhs_data.bp_rhs_s);
        off_t len = ((rhs_data.bp_rhs_e-rhs_data.bp_exp_e)-rhs_data.bp_rhs_s+gram.r_bits)/gram.r_bits;
        for(off_t i=0;i<len;i++){
            rhs.push_back(sym);
        }
    }else{
        for(off_t i=rhs_data.bp_exp_e+gram.r_bits;i<=rhs_data.bp_rhs_e;i+=gram.r_bits){
            rhs.push_back(gram.bitpos2symbol(i));
        }
    }

    rm_branch = rm_branch && rhs_data.bp_rhs_e==rhs_data.bp_exp_e;

    if(!rhs.empty() || rm_branch){
        offset_stack.emplace(-1, lvl, !rhs.empty(), rhs);
    }
}

template<class gram_t>
off_t parse_seq(size_t* text, off_t txt_size, gram_t& gram,
                std::unordered_map<uint64_t, new_rule_type>& new_gram_rules,
                size_t &next_av_nt, std::vector<uint64_t>& all_fps, uint64_t pf_seed, uint8_t p_level) {

    size_t mt_sym;
    size_t prev_sym, curr_sym, next_sym, end_sym = std::numeric_limits<size_t>::max();
    off_t txt_pos = 0, phrase_len, lb, rb;
    lz_like_map<size_t> map(text);

    uint64_t curr_hash, prev_hash, next_hash;

    bool inserted;
    off_t sym_bytes = sizeof(size_t);

    lb = txt_pos;
    prev_sym = text[txt_pos++];
    prev_hash = prev_sym<gram.r ? all_fps[prev_sym] : new_gram_rules[prev_sym].fp;

    curr_sym = text[txt_pos++];
    while(curr_sym==prev_sym) curr_sym = text[txt_pos++];
    rb = txt_pos-1;

    if(curr_sym==end_sym){
        phrase_len = rb-lb;
        mt_sym = map.insert(lb, phrase_len, inserted);
        if(!inserted){//we can not replace the first phrase occurrence as we use it as source for the dictionary
            assert(text[lb]!=end_sym);
            text[lb] = mt_sym+next_av_nt;//store the metasymbol in the first phrase position
            memset(&text[lb+1], (int)end_sym, sym_bytes*(phrase_len-1));//pad the rest of the phrase with dummy symbols
        }
    }else{
        curr_hash = curr_sym<gram.r ? all_fps[curr_sym] : new_gram_rules[curr_sym].fp;

        next_sym = text[txt_pos++];
        while(next_sym==curr_sym) next_sym = text[txt_pos++];

        while(next_sym!=end_sym){

            next_hash = next_sym<gram.r ? all_fps[next_sym] : new_gram_rules[next_sym].fp;
            if(prev_hash>curr_hash && curr_hash<next_hash){//local minimum
                phrase_len = rb-lb;
                mt_sym = map.insert(lb, phrase_len, inserted);
                if(!inserted){//we can not replace the first phrase occurrence as we use it as source for the dictionary
                    assert(text[lb]!=end_sym);
                    text[lb] = mt_sym+next_av_nt;
                    memset(&text[lb+1], (int)end_sym, sym_bytes*(phrase_len-1));
                }
                lb = rb;
            }

            rb = txt_pos-1;

            prev_hash = curr_hash;

            curr_hash = next_hash;
            curr_sym = next_sym;
            next_sym = text[txt_pos++];
            while(next_sym==curr_sym) next_sym = text[txt_pos++];
        }

        phrase_len = txt_pos-1-lb;
        assert(text[lb+phrase_len]==end_sym);
        mt_sym = map.insert(lb, phrase_len, inserted);
        if(!inserted){
            //we can not replace the first phrase occurrence as we use it as source for the dictionary
            assert(text[lb]!=end_sym);
            text[lb] = mt_sym+next_av_nt;//store the metasymbol in the first phrase position
            memset(&text[lb+1], (int)end_sym, sym_bytes*(phrase_len-1));//pad the rest of the phrase with dummy symbols
        }
        map.shrink_to_fit();
        map.destroy_table();
    }

    uint64_t fp;
    uint32_t source, end;
    std::vector<uint64_t> fp_sequence;
    std::vector<size_t> new_phrase;

    // create the parse in place
    map.insert_dummy_entry({uint32_t(txt_size-1), 0});
    size_t tot_phrases = map.phrase_set.size()-1;//do not count the dummy
    mt_sym = 0, lb = 0;
    off_t i=0, parse_size=0;

    while(mt_sym<tot_phrases) {
        assert(i==lb);
        source = map.phrase_set[mt_sym].source;
        end = source + map.phrase_set[mt_sym].len;
        for(size_t j=source;j<end;j++){
            fp_sequence.push_back(text[j]<gram.r ? all_fps[text[j]]: new_gram_rules[text[j]].fp);
            new_phrase.push_back(text[j]);
        }
        fp = XXH64(fp_sequence.data(), fp_sequence.size()*sizeof(uint64_t), pf_seed);
        new_gram_rules[mt_sym+next_av_nt] = {0, fp, p_level, new_phrase};
        fp_sequence.clear();
        new_phrase.clear();

        text[parse_size++] = mt_sym+next_av_nt;
        i+= map.phrase_set[mt_sym].len;//move out of the phrase boundary

        mt_sym++;
        lb = map.phrase_set[mt_sym].source;//position for the next phrase
        while(i<lb){//process the text area between consecutive phrases
            text[parse_size++] = text[i++];
            while(text[i]==end_sym && i<lb) i++;
        }
    }
    next_av_nt+=tot_phrases;
    return parse_size;
}

template<class gram_type>
void rem_txt_from_gram_int(gram_type& gram, std::vector<std::tuple<size_t, off_t, off_t>>& rem_coordinates){

    size_t sym;
    off_t d_seq_s, d_seq_e;
    uint8_t r_bits = gram.r_bits;
    std::vector<uint64_t> p_seeds = gram.get_parsing_seeds();
    size_t next_av_nt = gram.r, end_sym = std::numeric_limits<size_t>::max();

    std::vector<uint64_t> all_fps = gram.get_all_fps();

    std::stack<rule_type> left_off_stack;
    std::stack<rule_type> right_off_stack;
    std::unordered_map<uint64_t, new_rule_type> new_gram_rules;

    gram.breakdown(2);

    for(auto & coord : rem_coordinates){

        //we will delete exp(sym)[d_seq_s..d_seq_e] from the grammar,
        // where sym is a string id in the range [0..n_strings-1]
        sym = std::get<0>(coord);
        d_seq_s = std::get<1>(coord);
        d_seq_e = std::get<2>(coord);

        gram.print_parse_tree(sym, true);

        size_t lm_sym, rm_sym;
        uint8_t left_lvl, right_lvl;
        bool str_lm_branch=true, str_rm_branch=true;

        exp_data rhs_data{};
        gram.template exp_search_range<STR_EXP>(sym, d_seq_s, d_seq_e, rhs_data);

        left_lvl = gram.lc_par_tree_height();
        insert_left_offset(left_off_stack, rhs_data, gram, left_lvl, str_lm_branch);
        lm_sym = gram.bitpos2symbol(rhs_data.bp_exp_s);

        right_lvl = gram.lc_par_tree_height();
        insert_right_offset(right_off_stack, rhs_data, gram, right_lvl, str_rm_branch);
        rm_sym = gram.bitpos2symbol(rhs_data.bp_exp_e);

        off_t k = (rhs_data.bp_exp_e-rhs_data.bp_exp_s+r_bits) / r_bits;

        //left and rightmost branches of exp(sym)[d_seq_e..d_seq_e]'s parse tree descend together in
        // the parse tree of exp(sym) while k=1
        while(k == 1) {
            left_lvl = gram.parsing_level(lm_sym);
            gram.template exp_search_range<RULE_EXP>(lm_sym, d_seq_s, d_seq_e, rhs_data);
            insert_left_offset(left_off_stack, rhs_data, gram, left_lvl, str_lm_branch);
            insert_right_offset(right_off_stack, rhs_data, gram, left_lvl, str_rm_branch);

            lm_sym = gram.bitpos2symbol(rhs_data.bp_exp_s);
            k = (rhs_data.bp_exp_e-rhs_data.bp_exp_s+r_bits) / r_bits;
        }

        //with k>1, exp(sym)[d_seq_s] and exp(sym)[d_seq_e] are in different symbols of the right-hand side of sym's rule
        //rm_sym and lm_sym, respectively.
        rm_sym = rhs_data.is_rl? lm_sym: gram.bitpos2symbol(rhs_data.bp_exp_e);

        //descend over the leftmost branch of exp(sym)[d_seq_s..d_seq_e]'s parse tree
        while (lm_sym > gram.max_tsym) {
            left_lvl = gram.parsing_level(lm_sym);
            gram.template exp_search<RULE_EXP>(d_seq_s, lm_sym, rhs_data);
            insert_left_offset(left_off_stack, rhs_data, gram, left_lvl, str_lm_branch);
        }

        //descend over the rightmost branch of exp(sym)[d_seq_s..d_seq_e]'s parse tree
        while (rm_sym > gram.max_tsym) {
            right_lvl = gram.parsing_level(lm_sym);
            gram.template exp_search<RULE_EXP>(d_seq_e, rm_sym, rhs_data);
            insert_right_offset(right_off_stack, rhs_data, gram, right_lvl, str_rm_branch);
        }

        off_t g_level=0, last_level = gram.lc_par_tree_height()-1;
        std::vector<size_t> rhs;
        std::vector<size_t> new_seq;

        while(g_level<last_level){

            rule_type * left_offset = &left_off_stack.top();
            rule_type * right_offset = &right_off_stack.top();

            if(left_offset->exp_branch){
                rm_sym = left_offset->rhs.back();
                left_offset->rhs.pop_back();
                while(gram.parsing_level(rm_sym)>g_level){
                    left_offset->exp_branch = false;
                    gram.get_rhs(rm_sym, true, rhs);
                    size_t tmp_sym = rhs.back();
                    rhs.pop_back();
                    left_off_stack.emplace(rm_sym, 0, false, rhs);
                    rm_sym = tmp_sym;
                    left_offset = &left_off_stack.top();
                }
                left_offset->rhs.push_back(rm_sym);
            }

            if(right_offset->exp_branch){
                lm_sym = right_offset->rhs[0];
                right_offset->rhs.erase(right_offset->rhs.begin());
                while(gram.parsing_level(lm_sym)>g_level){
                    right_offset->exp_branch = false;
                    gram.get_rhs(lm_sym, true, rhs);
                    size_t tmp_sym = rhs[0];
                    rhs.erase(rhs.begin());
                    right_off_stack.emplace(lm_sym, 0, false, rhs);
                    lm_sym = tmp_sym;
                    right_offset = &right_off_stack.top();
                }
                right_offset->rhs.insert(right_offset->rhs.begin(), lm_sym);
            }

            //create the alternative sequence
            for(unsigned long s : new_seq){
                left_offset->rhs.push_back(s);
            }
            for(unsigned long s : right_offset->rhs){
                left_offset->rhs.push_back(s);
            }
            left_offset->rhs.push_back(end_sym);
            off_t p_size = parse_seq(left_offset->rhs.data(), left_offset->rhs.size(), gram, new_gram_rules, next_av_nt, all_fps, p_seeds[g_level+1], g_level);
            std::swap(left_offset->rhs, new_seq);
            new_seq.resize(p_size);

            left_off_stack.pop();
            right_off_stack.pop();

            /*off_t parsing_break = find_parsing_break(left_offset->rhs, right_offset->rhs, all_fps, new_gram_rules);
            if(parsing_break>0){//there is a break

                rm_sym = left_offset->nt;
                if(left_offset->nt==-1){
                    std::vector<uint64_t> fp_seq;
                    for(unsigned long & s : left_offset->rhs){
                        if(s>gram.r){
                            fp_seq.push_back(new_gram_rules[s].second);
                        }else{
                            fp_seq.push_back(all_fps[s]);
                        }
                    }
                    uint64_t new_fp = XXH64(fp_seq.data(), sizeof(uint64_t)*fp_seq.size(), p_seeds[g_level+1]);
                    rm_sym = next_av_nt++;
                    new_gram_rules[rm_sym] = {left_offset->rhs, new_fp};
                }

                lm_sym = right_offset->nt;
                if(right_offset->nt==-1){
                    std::vector<uint64_t> fp_seq;
                    for(unsigned long & s : right_offset->rhs){
                        if(s>gram.r){
                            fp_seq.push_back(new_gram_rules[s].second);
                        }else{
                            fp_seq.push_back(all_fps[s]);
                        }
                    }
                    uint64_t new_fp = XXH64(fp_seq.data(), sizeof(uint64_t)*fp_seq.size(), p_seeds[g_level+1]);
                    lm_sym = next_av_nt++;
                    new_gram_rules[lm_sym] = {right_offset->rhs, new_fp};
                }
                right_off_stack.pop();
                right_off_stack.top().rhs.push_back(lm_sym);
                right_off_stack.top().exp_branch = false;

            } else {//the phrases are merged

                std::vector<uint64_t> fp_seq;
                std::vector<size_t> new_seq;
                for(unsigned long s : left_offset->rhs){
                    fp_seq.push_back(all_fps[s]);
                    new_seq.push_back(s);
                }
                for(unsigned long & s : right_offset->rhs){
                    fp_seq.push_back(all_fps[s]);
                    new_seq.push_back(s);
                }
                uint64_t new_fp = XXH64(fp_seq.data(), sizeof(uint64_t)*fp_seq.size(), p_seeds[g_level+1]);
                rm_sym = next_av_nt++;
                new_gram_rules[rm_sym] = {new_seq, new_fp};

                right_off_stack.pop();
            }
            left_off_stack.pop();
            left_off_stack.top().rhs.push_back(rm_sym);
            left_off_stack.top().exp_branch = false;*/

            g_level++;
        }

        std::vector<new_rule_type*> sorted_new_rules(new_gram_rules.size());
        size_t i=0;
        for(auto & rule : new_gram_rules){
            sorted_new_rules[i++] = &rule.second;
        }

        std::sort(sorted_new_rules.begin(), sorted_new_rules.end(), [](auto const& a, auto const& b){
            if(a->lvl!=b->lvl){
                return a->lvl<b->lvl;
            }
            return a->fp < b->fp;
        });

        for(auto const &rule : sorted_new_rules){
            std::cout<<"lvl: "<<int(rule->lvl)<<", fp: "<<rule->fp<<", seq: ";
            for(auto const& s : rule->rhs){
                std::cout<<s<<" ";
            }
            std::cout<<""<<std::endl;
        }

    }
}

void rem_txt_from_gram(std::string& input_gram, std::vector<std::tuple<size_t, off_t, off_t>>& rem_coordinates){

    bool has_rl_rules, has_cg_rules, has_rand_access;
    std::tie(has_rl_rules, has_cg_rules, has_rand_access) = read_grammar_flags(input_gram);
    assert(has_rand_access);
    if(has_cg_rules){
        if(has_rl_rules){
            lc_gram_t<true, true, true> gram;
            load_from_file(input_gram, gram);
            rem_txt_from_gram_int(gram, rem_coordinates);
        }else{
            lc_gram_t<true, false, true> gram;
            load_from_file(input_gram, gram);
            rem_txt_from_gram_int(gram, rem_coordinates);
        }
    }else{
        if(has_rl_rules){
            lc_gram_t<false, true, true> gram;
            load_from_file(input_gram, gram);
            rem_txt_from_gram_int(gram, rem_coordinates);
        }else{
            lc_gram_t<false, false, true> gram;
            load_from_file(input_gram, gram);
            rem_txt_from_gram_int(gram, rem_coordinates);
        }
    }
}
#endif //LCG_EDITION_ALGORITHMS_H
