//
// Created by Diaz, Diego on 11.5.2024.
//

#ifndef LCG_EDITION_ALGORITHMS_H
#define LCG_EDITION_ALGORITHMS_H

#include "grammar.h"
#include "grammar_algorithms.h"
#include "unordered_map"

#define LEFT_CONTEXT 0
#define RIGHT_CONTEXT 1
#define BOTH_CONTEXTS 2

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
    uint64_t exp_len{};
    bool exist{};
    std::vector<size_t> rhs;

    new_rule_type()=default;

    new_rule_type(new_rule_type&& other) noexcept {
        std::swap(nt, other.nt);
        std::swap(fp, other.fp);
        std::swap(lvl, other.lvl);
        std::swap(exp_len, other.exp_len);
        std::swap(exist, other.exist);
        rhs.swap(other.rhs);
    }

    new_rule_type(off_t nt_, uint64_t& fp_, uint8_t lvl_, uint64_t exp_len_, bool exist_, std::vector<size_t>& rhs_)  {
        nt = nt_;
        fp = fp_;
        lvl = lvl_;
        exp_len = exp_len_;
        exist = exist_;
        rhs = rhs_;
    }

    new_rule_type(off_t nt_, uint64_t fp_, uint8_t lvl_, uint64_t exp_len_, bool exist_, std::vector<size_t>&& rhs_)  {
        nt = nt_;
        fp = fp_;
        lvl = lvl_;
        exp_len = exp_len_;
        exist = exist_;
        rhs.swap(rhs_);
    }

    new_rule_type& operator=(new_rule_type&& other) noexcept {
        std::swap(nt, other.nt);
        std::swap(fp, other.fp);
        std::swap(lvl, other.lvl);
        std::swap(exp_len, other.exp_len);
        std::swap(exist, other.exist);
        rhs.swap(other.rhs);
        return *this;
    }
};

template<class gram_t, class vector_type>
void insert_left_offset(std::stack<rule_type>& offset_stack, exp_data& rhs_data, gram_t& gram, uint8_t lvl, bool& lm_branch, vector_type& nt_freqs){
    std::vector<size_t> rhs;

    if(rhs_data.is_rl){
        size_t sym = gram.bitpos2symbol(rhs_data.bp_rhs_s);
        off_t len = (rhs_data.bp_exp_s-rhs_data.bp_rhs_s)/gram.r_bits;
        nt_freqs[sym]-=len;
        for(off_t i=0;i<len;i++){
            rhs.push_back(sym);
        }
    }else{
        for(off_t i=rhs_data.bp_rhs_s;i<=(rhs_data.bp_exp_s-gram.r_bits);i+=gram.r_bits){
            size_t sym = gram.bitpos2symbol(i);
            rhs.push_back(sym);
            nt_freqs[sym]--;
            //std::cout<<sym<<std::endl;
        }
    }

    lm_branch = lm_branch && rhs_data.bp_rhs_s==rhs_data.bp_exp_s;

    if(!rhs.empty() || lm_branch){
        offset_stack.emplace(-1, lvl, !rhs.empty(), rhs);
    }
}

template<class gram_t, class vector_type>
void insert_right_offset(std::stack<rule_type>& offset_stack, exp_data& rhs_data, gram_t& gram, uint8_t lvl, bool& rm_branch, vector_type& nt_freqs){
    std::vector<size_t> rhs;
    if(rhs_data.is_rl){
        size_t sym = gram.bitpos2symbol(rhs_data.bp_rhs_s);
        off_t len = ((rhs_data.bp_rhs_e-rhs_data.bp_exp_e)-rhs_data.bp_rhs_s+gram.r_bits)/gram.r_bits;
        nt_freqs[sym]-=len;
        for(off_t i=0;i<len;i++){
            rhs.push_back(sym);
        }
    }else{
        for(off_t i=rhs_data.bp_exp_e+gram.r_bits;i<=rhs_data.bp_rhs_e;i+=gram.r_bits){
            size_t sym = gram.bitpos2symbol(i);
            rhs.push_back(sym);
            nt_freqs[sym]--;
            //std::cout<<sym<<std::endl;
        }
    }

    rm_branch = rm_branch && rhs_data.bp_rhs_e==rhs_data.bp_exp_e;

    if(!rhs.empty() || rm_branch){
        offset_stack.emplace(-1, lvl, !rhs.empty(), rhs);
    }
}

template<class gram_t>
off_t parse_seq(size_t* text, off_t txt_size, gram_t& gram,
                std::vector<std::vector<new_rule_type>>& new_gram_rules,
                std::vector<uint64_t>& all_fps, uint64_t pf_seed, uint8_t p_level) {

    size_t mt_sym, next_av_nt=gram.r;
    size_t prev_sym, curr_sym, next_sym, end_sym = std::numeric_limits<size_t>::max();
    off_t txt_pos = 0, phrase_len, lb, rb;
    lz_like_map<size_t> map(text);

    uint64_t curr_hash, prev_hash, next_hash;

    bool inserted;
    off_t sym_bytes = sizeof(size_t);

    lb = txt_pos;
    prev_sym = text[txt_pos++];
    prev_hash = prev_sym<gram.r ? all_fps[prev_sym] : new_gram_rules[p_level-1][prev_sym-gram.r].fp;

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
        curr_hash = curr_sym<gram.r ? all_fps[curr_sym] : new_gram_rules[p_level-1][curr_sym-gram.r].fp;

        next_sym = text[txt_pos++];
        while(next_sym==curr_sym) next_sym = text[txt_pos++];

        while(next_sym!=end_sym){

            next_hash = next_sym<gram.r ? all_fps[next_sym] : new_gram_rules[p_level-1][next_sym-gram.r].fp;
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
    uint64_t exp_len=0;

    while(mt_sym<tot_phrases) {

        assert(i==lb);
        source = map.phrase_set[mt_sym].source;
        end = source + map.phrase_set[mt_sym].len;

        exp_len=0;
        for(size_t j=source;j<end;j++){
            fp_sequence.push_back(text[j]<gram.r ? all_fps[text[j]]: new_gram_rules[p_level-1][text[j]-gram.r].fp);
            new_phrase.push_back(text[j]);
            exp_len += text[j]<gram.r?  gram.exp_len(text[j]) : new_gram_rules[p_level-1][text[j]-gram.r].exp_len;
        }
        fp = XXH64(fp_sequence.data(), fp_sequence.size()*sizeof(uint64_t), pf_seed);
        new_gram_rules[p_level].emplace_back(0, fp, p_level, exp_len, false, new_phrase);
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

    return parse_size;
}

template<class gram_type>
void recompute_exp_samples(std::vector<size_t>& rhs, gram_type& gram, uint8_t lvl,
                           std::vector<std::vector<new_rule_type>>& edited_rules,
                           size_t sampling_rate){

    size_t acc =0, j=1, pos=0;
    for(size_t i=0;i<rhs.size();i++){
        if(rhs[i]<gram.r){
            acc += gram.exp_len(rhs[i]);
        }else{
            acc += edited_rules[lvl-1][rhs[i]-gram.r].exp_len;
        }

        if((j % sampling_rate)==0){
            rhs[pos++] = acc;
        }
        j++;
    }
    rhs.resize(pos);
    rhs.push_back(acc);
}


template<class gram_type, class vector_type>
void find_grammar_places(std::vector<std::vector<new_rule_type>>& edited_rules,
                         std::vector<std::pair<off_t, std::vector<size_t>>>& edited_strings,
                         gram_type& gram, vector_type& nt_freqs){
    size_t lvl=0;
    size_t eff_new_rules = 0;
    size_t del_rules = 0;

    for(;lvl<(gram.lvl_rules.size()-1);lvl++) {

        std::vector<new_rule_type *> perm(edited_rules[lvl].size());
        size_t i = 0;
        for (auto &new_rule: edited_rules[lvl]) {
            perm[i++] = &new_rule;
        }
        std::sort(perm.begin(), perm.end(), [](auto const &a, auto const &b) {
            //TODO handle collisions
            return a->fp < b->fp;
        });

        off_t middle;
        off_t left = gram.lvl_rules[lvl];
        off_t right = gram.lvl_rules[lvl + 1] - 1;
        off_t end_nt = gram.lvl_rules[lvl + 1];

        uint64_t fp_first, fp_second;

        for (auto const &p: perm) {

            //binary search of the rule's fingerprint
            while (left <= right) {
                middle = left + (right - left) / 2;
                fp_first = gram.get_fp(middle);
                if(fp_first > p->fp) {
                    right = middle - 1;
                    continue;
                }
                if((middle+1)==end_nt){
                    /*for(auto const& s : p.second->rhs){
                        std::cout<<s<<" ";
                    }
                    std::cout<<" | "<<middle<<" exist: "<<(p.second->fp==fp_first? "yes" : "no")<<" -> "<<gram.get_fp(middle)<<" "<<p.second->fp<<" exp:"<<p.second->exp_len<<std::endl;*/
                    break;
                }
                fp_second = gram.get_fp(middle+1);
                if(p->fp<fp_second) {
                    /*for(auto const& s : p.second->rhs){
                        std::cout<<s<<" ";
                    }
                    std::cout<<" | "<<middle<<" exist: "<<(p.second->fp==fp_first? "yes" : "no")<<" -> "<<gram.get_fp(middle)<<" "<<p.second->fp<<" "<<gram.get_fp(middle+1)<<" exp:"<<p.second->exp_len<<std::endl;*/
                    break;
                }
                left = middle + 1;
            }

            std::cout<<"new seq: ";
            for(auto &sym : p->rhs){
                if(sym>=gram.r && edited_rules[lvl-1][sym-gram.r].exist){
                    sym = edited_rules[lvl-1][sym-gram.r].nt;
                }
                if(sym<gram.r) nt_freqs[sym]++;
                std::cout<<sym<<" ";
            }
            std::cout<<" exist? "<<(p->fp==fp_first)<<std::endl;

            p->nt = middle;
            eff_new_rules +=p->fp!=fp_first;
            if(p->fp==fp_first){
                p->exist = true;
                //TODO handle collisions
            }
            left = middle+1;
            right = gram.lvl_rules[lvl+1]-1;
        }
    }

    for(auto & pair: edited_strings){
        std::cout<<"the resulting sequence for string "<<pair.first<<" is ";
        for(auto &s: pair.second){
            if(s>=gram.r && edited_rules[lvl-1][s-gram.r].exist){
                s = edited_rules[lvl-1][s-gram.r].nt;
            }
            std::cout<<s<<" ";
        }
        std::cout<<""<<std::endl;
    }

    /*std::stack<size_t> stack;
    size_t copies;
    for(size_t nt=gram.max_tsym+1;nt<nt_freqs.size();nt++){
        if(nt_freqs[nt]==0){
            stack.push(nt);
            while(!stack.empty()){
                size_t sym = stack.top();
                stack.pop();
                auto range = gram.nt2bitrange(sym);
                for(off_t bp=range.second;bp>=range.second;bp-=gram.r_bits){

                    sym = gram.bitpos2symbol(bp);
                    copies=1;
                    if(gram.is_rl_sym(sym)){
                        range = gram.nt2bitrange(sym);
                        sym = gram.bitpos2symbol(range.first);
                        copies = gram.bitpos2symbol(range.second);
                    }

                    if(nt_freqs[sym]>0){//for the case when we have a path that was removed in the tree
                        nt_freqs[sym]-=copies;
                        if(nt_freqs[sym]==0){
                            stack.push(sym);
                        }
                    }
                }
            }
        }
    }*/
}

template<class gram_type>
void insert_edited_rules(std::vector<std::vector<new_rule_type>>& edited_rules,
                         std::vector<std::pair<off_t, std::vector<size_t>>>& edited_strings,
                         gram_type& gram){


    off_t nt_offset=0;
    int_array<size_t> nt_offsets(gram.r, sym_width(edited_rules.size()));
    o_file_stream<size_t> edited_gram("pepe.txt", BUFFER_SIZE, std::ios::binary | std::ios::out);

    for(size_t i=0;i<=gram.max_tsym;i++){
        nt_offsets[i] = 0;
    }

    size_t lvl=0;
    for(;lvl<(gram.lvl_rules.size()-1);lvl++) {

        std::vector<new_rule_type *> perm(edited_rules[lvl].size());
        size_t i=0;
        for(auto & new_rule : edited_rules[lvl]){
            perm[i++] =&new_rule;
        }
        std::sort(perm.begin(), perm.end(), [](auto const& a, auto const& b){
            return a->fp < b->fp;
        });

        off_t middle;
        off_t left = gram.lvl_rules[lvl];
        off_t right = gram.lvl_rules[lvl+1]-1;

        off_t first_nt = gram.lvl_rules[lvl];
        off_t end_nt = gram.lvl_rules[lvl+1];

        uint64_t fp_first, fp_second;

        for(auto const &p : perm){

            //binary search of the rule's fingerprint
            while (left <= right) {
                middle = left + (right - left) / 2;
                fp_first = gram.get_fp(middle);
                if(fp_first > p->fp) {
                    right = middle - 1;
                    continue;
                }
                if((middle+1)==end_nt){
                    for(auto const& s : p->rhs){
                        std::cout<<s<<" ";
                    }
                    std::cout<<" | "<<middle<<"+"<<(p->fp==fp_first? nt_offset : nt_offset+1)<<" -> "<<gram.get_fp(middle)<<" "<<p->fp<<" exp:"<<p->exp_len<<std::endl;
                    break;
                }
                fp_second = gram.get_fp(middle+1);
                if(p->fp<fp_second) {
                    for(auto const& s : p->rhs){
                        std::cout<<s<<" ";
                    }
                    std::cout<<" | "<<middle<<"+"<<(p->fp==fp_first? nt_offset : nt_offset+1)<<" -> "<<gram.get_fp(middle)<<" "<<p->fp<<" "<<gram.get_fp(middle+1)<<" exp:"<<p->exp_len<<std::endl;
                    break;
                }
                left = middle + 1;
            }

            for(off_t nt=first_nt;nt<=middle;nt++){
                nt_offsets[nt] = nt_offset;
                auto res = gram.nt2bitrange(nt);
                for(off_t j=res.first;j<=res.second;j+=gram.r_bits){
                    size_t sym = gram.bitpos2symbol(j);
                    sym+=nt_offsets[sym];
                    edited_gram.push_back((sym<<1UL) | (j==res.second));
                }
                auto e_data = gram.template nt2expdata<RULE_EXP>(nt);
                for(off_t j = std::get<0>(e_data);j<=std::get<1>(e_data);j+=std::get<2>(e_data)){
                    size_t exp = gram.rule_stream.read(j, j+std::get<2>(e_data)-1);
                    edited_gram.push_back((exp<<1UL) | (j==std::get<1>(e_data)));
                }
            }

            if(p->fp==fp_first){
                //TODO handle collisions
                p->nt = middle+nt_offset;
            } else{
                //handle the new rule
                size_t last = p->rhs.size()-1, j=0;
                for(auto sym : p->rhs){
                    if(sym<gram.r){
                        sym +=nt_offsets[sym];
                    }else{
                        sym = edited_rules[lvl-1][sym-gram.r].nt;
                    }
                    //std::cout<<sym<<", ";
                    edited_gram.push_back(sym<<1UL | (j==last));
                    j++;
                }
                //std::cout<<""<<std::endl;

                recompute_exp_samples(p->rhs, gram, lvl, edited_rules, gram.rl_samp_rate);
                for(size_t u=0;u<p->rhs.size();u++){
                    edited_gram.push_back((p->rhs[u]<<1UL) | (u==(p->rhs.size()-1)));
                }
                nt_offset++;
                p->nt = middle+nt_offset;
            }
            first_nt = middle+1;
            left = middle+1;
            right = gram.lvl_rules[lvl+1]-1;
        }
    }

    //TODO in case the edition create more levels in the grammar
    for(;lvl<edited_rules.size();lvl++){
    }

    //insert the strings in the grammar (it is simpler than the rules)
    size_t e_str_pos=0;
    size_t next_edited_str=edited_strings[e_str_pos].first;
    for(size_t str=0;str<gram.n_strings();str++){
        if(str==next_edited_str){

            size_t last = edited_strings[e_str_pos].second.size()-1, j=0;
            for(auto sym : edited_strings[e_str_pos].second){
                if(sym<gram.r){
                    sym +=nt_offsets[sym];
                }else{
                    sym = edited_rules[lvl-1][sym-gram.r].nt;
                }
                edited_gram.push_back(sym<<1UL | (j==last));
                j++;
            }
            recompute_exp_samples(edited_strings[e_str_pos].second, gram, lvl, edited_rules, gram.rl_samp_rate);
            last = edited_strings[e_str_pos].second.size()-1;
            size_t u=0;
            for(auto e : edited_strings[e_str_pos].second){
                edited_gram.push_back((e<<1UL) | (u==last));
                u++;
            }
            next_edited_str=edited_strings[++e_str_pos].first;
        }else{
            auto res = gram.str2bitrange(str);
            for(off_t j=res.first;j<=res.second;j+=gram.r_bits){
                size_t sym = gram.bitpos2symbol(j);
                sym+=nt_offsets[sym];
                edited_gram.push_back((sym<<1UL) | (j==res.second));
            }
            auto e_data = gram.template nt2expdata<STR_EXP>(str);
            for(off_t j = std::get<0>(e_data);j<=std::get<1>(e_data);j+=std::get<2>(e_data)){
                size_t exp = gram.rule_stream.read(j, j+std::get<2>(e_data)-1);
                edited_gram.push_back((exp<<1UL) | (j==std::get<1>(e_data)));
            }
        }
    }
}

/*template<class gram_t, class vector_t>
void subtract_context(gram_t& gram, vector_t& nt_freqs, exp_data& exp_data, uint8_t context){
    switch (context) {
        case BOTH_CONTEXTS:
            for(off_t bp = exp_data.bp_rhs_s; bp<exp_data.bp_exp_s; bp+=gram.r_bits){
                size_t symbol = gram.bitpos2symbol(bp);
                nt_freqs[symbol]--;
            }
            for(off_t bp = exp_data.bp_exp_e+gram.r_bits; bp<exp_data.bp_rhs_e; bp+=gram.r_bits){
                size_t symbol = gram.bitpos2symbol(bp);
                nt_freqs[symbol]--;
            }
            break;
        case LEFT_CONTEXT:
            for(off_t bp = exp_data.bp_rhs_s; bp<exp_data.bp_exp_s; bp+=gram.r_bits){
                size_t symbol = gram.bitpos2symbol(bp);
                nt_freqs[symbol]--;
            }
            break;
        case RIGHT_CONTEXT:
            for(off_t bp = exp_data.bp_exp_e+gram.r_bits; bp<exp_data.bp_rhs_e; bp+=gram.r_bits){
                size_t symbol = gram.bitpos2symbol(bp);
                nt_freqs[symbol]--;
            }
            break;
        default:
            std::cout<<"There is bug"<<std::endl;
            exit(1);
    }
}*/

template<class gram_t, class vector_t>
void subtract_int_par_trees(gram_t& gram, vector_t& nt_freqs, std::vector<std::pair<off_t, off_t>>& internal_par_trees){

    std::stack<size_t> stack;
    size_t copies;
    for(auto const& i_par_tree: internal_par_trees){

        for(off_t bp=i_par_tree.first;bp<=i_par_tree.second;bp+=gram.r_bits){
            size_t sym = gram.bitpos2symbol(bp);
            stack.push(sym);
        }

        while(!stack.empty()){
            size_t sym = stack.top();
            stack.pop();
            if(sym>gram.max_tsym){
                auto range = gram.nt2bitrange(sym);
                for(off_t bp=range.second;bp>=range.second;bp-=gram.r_bits){
                    sym = gram.bitpos2symbol(bp);
                    copies=1;
                    if(gram.is_rl_sym(sym)){
                        range = gram.nt2bitrange(sym);
                        sym = gram.bitpos2symbol(range.first);
                        copies = gram.bitpos2symbol(range.second);
                    }
                    if(nt_freqs[sym]>0){//for the case when we have a path that was removed in the tree
                        nt_freqs[sym]-=copies;
                        if(nt_freqs[sym]==0){
                            stack.push(sym);
                        }
                    }
                }
            }
        }
    }

}


template<class gram_type>
void rem_txt_from_gram_int(gram_type& gram, std::vector<str_coord_type>& coordinates){

    size_t sym;
    off_t d_seq_s, d_seq_e;
    uint8_t r_bits = gram.r_bits;
    std::vector<uint64_t> p_seeds = gram.get_parsing_seeds();
    size_t end_sym = std::numeric_limits<size_t>::max();

    std::vector<uint64_t> all_fps = gram.get_all_fps();
    std::vector<size_t> nt_freqs(gram.r-1, 0);
    gram.get_nt_freqs(nt_freqs);

    std::stack<rule_type> left_off_stack;
    std::stack<rule_type> right_off_stack;

    std::vector<std::vector<new_rule_type>> edited_rules(gram.lc_par_tree_height());
    std::vector<std::pair<off_t, std::vector<size_t>>> edited_strings;
    std::vector<std::pair<off_t, off_t>> int_par_trees;

    for(auto & coord : coordinates){

        //we will delete exp(sym)[d_seq_s..d_seq_e] from the grammar,
        // where sym is a string id in the range [0..n_strings-1]
        sym = coord.str;
        d_seq_s = coord.start;
        d_seq_e = coord.end;

        gram.print_parse_tree(sym, true);

        size_t lm_sym, rm_sym;
        bool str_lm_branch=true, str_rm_branch=true;
        exp_data rhs_data{};

        gram.template exp_search_range<STR_EXP>(sym, d_seq_s, d_seq_e, rhs_data);
        insert_left_offset(left_off_stack, rhs_data, gram, gram.lc_par_tree_height(), str_lm_branch, nt_freqs);
        insert_right_offset(right_off_stack, rhs_data, gram, gram.lc_par_tree_height(), str_rm_branch, nt_freqs);
        lm_sym = gram.bitpos2symbol(rhs_data.bp_exp_s);
        nt_freqs[lm_sym]--;
        off_t k = (rhs_data.bp_exp_e-rhs_data.bp_exp_s+r_bits) / r_bits;

        //left and rightmost branches of exp(sym)[d_seq_e..d_seq_e]'s parse tree descend together in
        // the parse tree of exp(sym) while k=1
        while(k == 1 && lm_sym>gram.max_tsym) {
            uint8_t lvl = gram.parsing_level(lm_sym);
            gram.template exp_search_range<RULE_EXP>(lm_sym, d_seq_s, d_seq_e, rhs_data);
            insert_left_offset(left_off_stack, rhs_data, gram, lvl, str_lm_branch, nt_freqs);
            insert_right_offset(right_off_stack, rhs_data, gram, lvl, str_rm_branch, nt_freqs);
            lm_sym = gram.bitpos2symbol(rhs_data.bp_exp_s);
            nt_freqs[lm_sym]--;
            k = (rhs_data.bp_exp_e-rhs_data.bp_exp_s+r_bits) / r_bits;
        }

        rm_sym = lm_sym;
        //with k>1, exp(sym)[d_seq_s] and exp(sym)[d_seq_e] are in different symbols of the right-hand side of sym's rule
        //rm_sym and lm_sym, respectively.
        if(k>1){//k=1 is a corner case: when we are removing only one symbol
            rm_sym = rhs_data.is_rl? lm_sym: gram.bitpos2symbol(rhs_data.bp_exp_e);
            nt_freqs[rm_sym]--;
            int_par_trees.emplace_back(rhs_data.bp_exp_s+gram.r_bits, rhs_data.bp_exp_e-gram.r_bits);
        }

        //descend over the rightmost branch of exp(sym)[d_seq_s..d_seq_e]'s parse tree
        while (rm_sym > gram.max_tsym) {
            uint8_t lvl = gram.parsing_level(lm_sym);
            gram.template exp_search<RULE_EXP>(d_seq_e, rm_sym, rhs_data);
            insert_right_offset(right_off_stack, rhs_data, gram, lvl, str_rm_branch, nt_freqs);
            int_par_trees.emplace_back(rhs_data.bp_rhs_s, rhs_data.bp_exp_s-gram.r_bits);
            nt_freqs[rm_sym]--;
        }

        //descend over the leftmost branch of exp(sym)[d_seq_s..d_seq_e]'s parse tree
        while (lm_sym > gram.max_tsym) {
            uint8_t lvl = gram.parsing_level(lm_sym);
            gram.template exp_search<RULE_EXP>(d_seq_s, lm_sym, rhs_data);
            insert_left_offset(left_off_stack, rhs_data, gram, lvl, str_lm_branch, nt_freqs);
            int_par_trees.emplace_back(rhs_data.bp_exp_e+gram.r_bits, rhs_data.bp_rhs_e);
            nt_freqs[lm_sym]--;
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

                    for(auto const& s : rhs){
                        //std::cout<<s<<std::endl;
                        nt_freqs[s]--;
                    }

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
                    for(auto const& s : rhs){
                        //std::cout<<s<<std::endl;
                        nt_freqs[s]--;
                    }

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
            off_t p_size = parse_seq(left_offset->rhs.data(), left_offset->rhs.size(), gram, edited_rules, all_fps, p_seeds[g_level+1], g_level);
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
        edited_strings.emplace_back(coord.str, new_seq);
        find_grammar_places(edited_rules, edited_strings, gram, nt_freqs);
        subtract_int_par_trees(gram, nt_freqs, int_par_trees);
    }
    //insert_edited_rules(edited_rules, edited_strings, gram);
}

void rem_txt_from_gram(std::string& input_gram, std::vector<str_coord_type>& rem_coordinates){

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
