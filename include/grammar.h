//
// Created by Diaz, Diego on 3.7.2023.
//

#ifndef LCG_GRAMMAR_H
#define LCG_GRAMMAR_H
#include <stack>
#include <queue>

#include "cds/cdt_common.hpp"
#include "cds/int_array.h"
#include "cds/file_streams.hpp"
#include "hashing.h"
#include "cds/utils.h"
#include "lz-like-strategy/partial_gram.h"

struct lc_gram_buffer_t{

    size_t                             n{}; //n: number of symbols
    size_t                             r{}; //r: number of grammar symbols (nter + ter)
    size_t                             prev_r{}; //r: number of grammar symbols (nter + ter)

    size_t                             c{}; //c: length of the right-hand of the start symbol
    size_t                             g{}; //g: sum of the rules' right-hand sides
    size_t                             max_tsym{}; //highest terminal symbol
    size_t                             long_str{};
    uint8_t                            sep_tsym{}; //separator symbol in the collection. This is 0 if the text is a single string
    std::pair<size_t, size_t>          rl_rules={0,0}; //mark the range of run-length compressed rules

    std::string                        rules_file; // rules are concatenated in this array
    std::vector<uint8_t>               terminals;
    std::vector<hashing>               par_functions;

    o_file_stream<size_t>              rules_buffer;
    std::vector<off_t>&                str_boundaries;//
    std::vector<size_t>                lvl_rules;     // number of rules generated in every parsing round
    std::vector<size_t>                lvl_size;      // number of rules generated in every parsing round


    lc_gram_buffer_t(std::string& rules_f,
                     std::vector<uint8_t>& alphabet,
                     std::vector<off_t>& str_bd_,
                     size_t long_str_,
                     uint8_t sep_symbol_): sep_tsym(sep_symbol_),
                                           rules_file(rules_f),
                                           rules_buffer(rules_file, BUFFER_SIZE, std::ios::out),
                                           str_boundaries(str_bd_){

        terminals = alphabet;
        n = str_boundaries.back();
        max_tsym = terminals.back();
        prev_r = 0;
        r = max_tsym + 1;
        g = r;
        long_str = long_str_;

        for(size_t i=0;i<=max_tsym; i++){
            rules_buffer.push_back((i<<1UL) | 1UL);
        }
    }

    ~lc_gram_buffer_t(){
        rules_buffer.close();
    }

    template<class vector_type>
    void create_lc_rules(std::vector<rand_order>& str_orders, vector_type& parsing_set){
        size_t sym, tot_syms=0, len;
        bool first;
        for(auto & str : str_orders){
            first = true;
            len = str.str_len;
            for(size_t j=0;j<len;j++){
                sym = prev_r+parsing_set[str.str_ptr+j];
                //assert(sym>0);
                sym = (sym<<1UL) | first;
                rules_buffer.push_back(sym);
                first = false;
            }
            tot_syms+=len;
        }
        prev_r = r;
        r += str_orders.size();
        g += tot_syms;
        lvl_rules.push_back(str_orders.size());
        lvl_size.push_back(tot_syms);
    }

    template<class sym_type>
    void insert_comp_string(std::string& input_file){

        i_file_stream<sym_type> ifs(input_file, BUFFER_SIZE);
        for(size_t i=0;i<ifs.size();i++){
            str_boundaries[i] = i;
            size_t sym = prev_r+ifs.read(i);
            sym = sym<<1UL | (i==0);
            rules_buffer.push_back(sym);
        }
        str_boundaries.back() = ifs.size();
        r++;
        g+=ifs.size();
        c=ifs.size();
        ifs.close();
        assert(g==rules_buffer.size());
    }
};

template<bool is_cg=false, bool is_rl=false>
struct lc_gram_t {

    size_t                             n{}; //n: number of symbols in the original text
    size_t                             r{}; //r: number of grammar symbols (nter + ter)
    size_t                             c{}; //c: length of the right-hand of the start symbol
    size_t                             g{}; //g: sum of the rules' right-hand sides
    size_t                             s{}; //s: number of strings
    size_t                             longest_str{}; //length of the longest string encoded in the grammar
    size_t                             max_tsym{}; //highest terminal symbol
    size_t                             min_nter{}; //smallest nonterminal symbol
    size_t                             cg_mark{};
    size_t                             n_cg_rules{};

    uint8_t                            sep_tsym{}; //separator symbol in the collection. This is 0 if the text is a single string
    bool                               is_simplified=false;
    const static bool                  has_rl_rules=is_rl;
    const static bool                  has_cg_rules=is_cg;
    bool                               has_rand_access=false;

    std::vector<uint8_t>               terminals; //set of terminals
    std::vector<hashing>               par_functions; //list of hash functions from which the grammar was constructed
    std::vector<size_t>                str_boundaries; // start position of every string in the compressed string
    std::vector<size_t>                lvl_rules; //number of rules generated in every round of locally-consistent parsing

    int_array<size_t>                  rules; //concatenated set of grammar rules
    int_array<size_t>                  rl_ptr; //pointer in "rules" to the leftmost symbol of each rule

    int_array<size_t>                  rule_exp;// length of the nt expansions
    int_array<size_t>                  sampled_exp;// sampled nt expansions in "rules"
    size_t                             samp_rate=4;//sampling rate to sample nt expansions in "rules"

    std::pair<size_t, size_t>          run_len_nt{0,0};//first run-length rule and total number of run-length rules

    lc_gram_t()= default;

    explicit lc_gram_t(std::string& p_gram_file, std::vector<uint64_t>& fp_seeds_){

        partial_gram<uint8_t> p_gram;
        std::ifstream ifs(p_gram_file, std::ios::binary);

        p_gram.load_metadata(ifs);

        n = p_gram.txt_size();
        r = p_gram.tot_gram_symbols();
        g = p_gram.gram_size();
        s = p_gram.tot_strings();
        c = p_gram.tot_strings();

        max_tsym = size_t(p_gram.max_terminal_symbol());
        sep_tsym = p_gram.separator_symbol();

        //add the set of seeds we use for the hash function
        //TODO add the set of seeds we used for the hash functions
        //fp_seeds = p_gram.par_functions;

        lvl_rules.reserve(p_gram.metadata.size());
        size_t acc=max_tsym+1, tmp;
        for(size_t i=1;i<p_gram.metadata.size();i++){
            tmp = p_gram.metadata[i].n_rules;
            lvl_rules.push_back(acc);
            acc+=tmp;
        }
        lvl_rules.push_back(acc);

        rules.set_width(sym_width(r));
        rules.resize(g);

        rl_ptr.set_width(sym_width(g));
        rl_ptr.resize(r-max_tsym);

        terminals.resize(max_tsym+1);
        for(size_t j=0;j<=max_tsym;j++){
            rules.write(j, j);
            terminals[j] = j;
        }

        bool last;
        size_t rule=0, acc_rules=0, pos, width, n_bits, j=max_tsym+1, rule_start_ptr=j, sym;
        bitstream<size_t> rules_buffer;

        for(size_t i=0;i<p_gram.lvl;i++){

            p_gram.load_next_rule_set(ifs, i, rules_buffer);
            pos = 0;
            width = p_gram.metadata[i+1].sym_width;
            n_bits = p_gram.metadata[i+1].n_bits();

            while(pos<n_bits){
                sym = rules_buffer.read(pos, pos+width-1);
                last = sym & 1UL;
                sym = acc_rules+(sym>>1);
                rules.write(j, sym);
                pos+=width;
                j++;

                if(last){
                    rl_ptr.write(rule, rule_start_ptr);
                    rule_start_ptr = j;
                    rule++;
                }
            }
            acc_rules+=p_gram.metadata[i+1].n_rules;
            assert(acc_rules==rule);
            assert(pos==n_bits);
        }
        assert(rules.size()==g);
        assert(rule==(r-(max_tsym+1)));
        rl_ptr.write(rule, g);

        size_t offset = g-c;
        str_boundaries.resize(s+1);
        for(size_t str=0;str<=s;str++){
            str_boundaries[str] = offset+str;
        }
        assert(str_boundaries[0]==offset);
        assert(str_boundaries.back()==rules.size());
    }

    explicit lc_gram_t(lc_gram_buffer_t& gram_buff) {

        i_file_stream<size_t> rules_buffer(gram_buff.rules_file, BUFFER_SIZE);

        n = gram_buff.n;
        r = gram_buff.r;
        g = gram_buff.g;
        c = gram_buff.c;
        s = gram_buff.str_boundaries.size()-1;
        longest_str = gram_buff.long_str;

        max_tsym = gram_buff.max_tsym;
        sep_tsym = gram_buff.sep_tsym;

        terminals = gram_buff.terminals;
        par_functions = gram_buff.par_functions;
        lvl_rules = gram_buff.lvl_rules;

        run_len_nt = gram_buff.rl_rules;

        rules.set_width(sym_width(r));
        rules.resize(g);

        rl_ptr.set_width(sym_width(g));
        rl_ptr.resize(r-max_tsym);

        for(size_t i=0;i<=max_tsym;i++){
            rules.write(i, rules_buffer.read(i)>>1UL);
        }

        size_t rule=0;
        for(size_t i=max_tsym+1;i<rules.size();i++){
            size_t sym = rules_buffer.read(i);
            bool first = sym & 1UL;

            rules.write(i, sym>>1UL);
            if(first){
                rl_ptr.write(rule, i);
                rule++;
            }
        }
        assert(rules.size()==g);
        assert(rule==(r-(max_tsym+1)));

        rl_ptr.write(rule, g);

        size_t offset = g-c;
        str_boundaries.resize(gram_buff.str_boundaries.size());
        size_t str=0;
        for(off_t & str_boundary : gram_buff.str_boundaries){
            str_boundaries[str++] = str_boundary + offset;
        }
        assert(str_boundaries[0]==offset);
        assert(str_boundaries.back()==rules.size());

        size_t acc=max_tsym+1, tmp;
        for(unsigned long & lvl_rule : lvl_rules){
            tmp = lvl_rule;
            lvl_rule = acc;
            acc+=tmp;
        }
        lvl_rules.push_back(acc);
        rules_buffer.close();
    }

    size_t serialize(std::ofstream &ofs){
        size_t written_bytes=0;
        written_bytes +=serialize_elm(ofs, n);
        written_bytes +=serialize_elm(ofs, r);
        written_bytes +=serialize_elm(ofs, g);
        written_bytes +=serialize_elm(ofs, c);
        written_bytes +=serialize_elm(ofs, s);
        written_bytes +=serialize_elm(ofs, samp_rate);
        written_bytes +=serialize_elm(ofs, longest_str);
        written_bytes +=serialize_elm(ofs, max_tsym);
        written_bytes +=serialize_elm(ofs, min_nter);
        written_bytes +=serialize_elm(ofs, sep_tsym);
        written_bytes +=serialize_elm(ofs, is_simplified);
        written_bytes +=serialize_elm(ofs, has_rl_rules);
        written_bytes +=serialize_elm(ofs, has_cg_rules);
        written_bytes +=serialize_elm(ofs, cg_mark);
        written_bytes +=serialize_elm(ofs, n_cg_rules);
        written_bytes +=serialize_elm(ofs, has_rand_access);
        written_bytes +=serialize_elm(ofs, run_len_nt.first);
        written_bytes +=serialize_elm(ofs, run_len_nt.second);
        written_bytes +=serialize_plain_vector(ofs, lvl_rules);
        written_bytes +=serialize_plain_vector(ofs, par_functions);

        written_bytes +=serialize_plain_vector(ofs, terminals);
        written_bytes +=serialize_plain_vector(ofs, str_boundaries);
        written_bytes +=rules.serialize(ofs);
        written_bytes +=rl_ptr.serialize(ofs);

        //written_bytes +=rule_exp.serialize(ofs);
        //written_bytes +=sampled_exp.serialize(ofs);
        return written_bytes;
    }

    void load_metadata(std::ifstream &ifs){
        load_elm(ifs, n);
        load_elm(ifs, r);
        load_elm(ifs, g);
        load_elm(ifs, c);
        load_elm(ifs, s);
        load_elm(ifs, samp_rate);
        load_elm(ifs, longest_str);
        load_elm(ifs, max_tsym);
        load_elm(ifs, min_nter);
        load_elm(ifs, sep_tsym);
        load_elm(ifs, is_simplified);
        bool tmp;
        load_elm(ifs, tmp);
        assert(has_rl_rules==is_rl);
        load_elm(ifs, tmp);
        assert(has_cg_rules==is_cg);
        load_elm(ifs, cg_mark);
        load_elm(ifs, n_cg_rules);
        load_elm(ifs, has_rand_access);
        load_elm(ifs, run_len_nt.first);
        load_elm(ifs, run_len_nt.second);
        load_plain_vector(ifs, lvl_rules);
        load_plain_vector(ifs, par_functions);
    }

    void load(std::ifstream &ifs){
        load_metadata(ifs);
        load_plain_vector(ifs, terminals);
        load_plain_vector(ifs, str_boundaries);
        rules.load(ifs);
        rl_ptr.load(ifs);
        //rule_exp.load(ifs);
    }

    [[nodiscard]] inline std::pair<size_t, size_t> nt2phrase(size_t sym) const {
        assert(sym>max_tsym);
        size_t pos = sym - max_tsym-1;
        return {rl_ptr[pos], rl_ptr[pos+1]-1};
    }

    [[nodiscard]] inline bool is_terminal(size_t sym) const {
        return sym<=max_tsym;
    }

    [[nodiscard]] inline bool is_rl_sym(size_t symbol) const{
        return run_len_nt.first<=symbol && symbol<(run_len_nt.first+run_len_nt.second);
    }

    [[nodiscard]] inline bool is_cg_nt(size_t sym) const{
        auto res = nt2phrase(sym);
        return rules[res.first]>=cg_mark;
    }

    [[nodiscard]] inline size_t first_rl_sym() const{
        return run_len_nt.first;
    }

    [[nodiscard]] inline size_t last_rl_sym() const{
        return run_len_nt.first+run_len_nt.second-1;
    }

    [[nodiscard]] inline size_t n_terminals() const {
        return (max_tsym+1);
    }

    [[nodiscard]] inline size_t n_nonterminals() const {
        return r-n_terminals();
    }

    [[nodiscard]] inline size_t comp_str_size() const {
        return c;
    }

    inline size_t get_byte_ter(size_t sym){
        assert(is_terminal(sym));
        if(is_simplified){
            return terminals[sym];
        }else{
            return sym;
        }
    }

    [[nodiscard]] inline off_t parsing_level(size_t symbol) const {
        if(symbol <= max_tsym) return 0;
        for(off_t i=0;i<(off_t)lvl_rules.size();i++){
            if(lvl_rules[i]<=symbol && symbol<lvl_rules[i+1]){
                return i+1;
            }
        }
        return -1;
    }

    [[nodiscard]] inline size_t pos2symbol(size_t idx) const{
        assert(idx<rules.size());
        return rules.read(idx);
    }

    [[nodiscard]] inline size_t start_symbol() const {
        return r-1;
    }

    [[nodiscard]] inline size_t n_strings() const {
        return str_boundaries.size()-1;
    }

    [[nodiscard]] inline std::pair<size_t, size_t> str2phrase(size_t str) const {
        return {str_boundaries[str], str_boundaries[str+1]-1};
    }

    void stats(size_t pad) const {
        size_t n_ter = n_terminals();
        size_t n_nter = n_nonterminals();

        auto pt_bytes = INT_CEIL((r-n_ter)*sym_width(g), 8);//space of the pointers for the nonterminals
        auto g_bytes = INT_CEIL(g*sym_width(r), 8); //space of the expansions
        auto pt_str_bytes = (s+1)*sizeof(size_t);  //space of the pointers to the strings

        float comp_ratio = float(n)/float(pt_bytes+g_bytes+pt_str_bytes);

        std::string pad_string(pad,' ');

        std::cout<<pad_string<<"Number of compressed symbols:   "<<n<<std::endl;
        std::cout<<pad_string<<"Number of compressed strings:   "<<s<<" ("<<report_space((off_t)pt_str_bytes)<<" in pointers)"<<std::endl;
        std::cout<<pad_string<<"Separator symbol:               "<<(int)sep_tsym<<std::endl;
        std::cout<<pad_string<<"Longest string:                 "<<longest_str<<std::endl;
        std::cout<<pad_string<<"Number of terminals:            "<<n_ter<<std::endl;
        std::cout<<pad_string<<"Number of non-terminals:        "<<n_nter<<" ("<<report_space((off_t)pt_bytes)<<" in pointers)"<<std::endl;
        std::cout<<pad_string<<"Grammar size:                   "<<g<<" ("<<report_space((off_t)g_bytes)<<")"<<std::endl;
        std::cout<<pad_string<<"Length of the comp. collection: "<<c<<std::endl;
        std::cout<<pad_string<<"Approx. compression ratio:      "<<comp_ratio<<std::endl;
        std::cout<<pad_string<<"Simplified:                     "<<(is_simplified ? "yes" : "no")<<std::endl;
        std::cout<<pad_string<<"Run-len rules:                  "<<(has_rl_rules ? "yes" : "no")<<std::endl;
        std::cout<<pad_string<<"Collage rules:                  "<<(has_cg_rules ? "yes" : "no")<<std::endl;
        std::cout<<pad_string<<"Random access support:          "<<(has_rand_access? "yes" : "no");
        if(has_rand_access){
            auto ras_bytes = INT_CEIL((rule_exp.size()*rule_exp.width()+ sampled_exp.size()+sampled_exp.width()), 8);
            std::cout<<" ("<<report_space((off_t)ras_bytes)<<" in data)"<<std::endl;
            std::cout<<pad_string<<"  Samp. rate for non.ter exps: 1/"<<samp_rate<<std::endl;
        }else{
            std::cout<<""<<std::endl;
        }
    }

    void breakdown(size_t pad){
        assert(g==rules.size());
        assert(r-(max_tsym+1)==(rl_ptr.size()-1));
        assert(s==str_boundaries.size()-1);
        stats(pad);

        std::string pad_string(pad,' ');

        size_t tot_sym, n_rules;
        std::cout<<pad_string<<"Grammar rules per level"<<std::endl;
        for(size_t i=0;i<lvl_rules.size()-1; i++){
            n_rules = lvl_rules[i+1]-lvl_rules[i];//number of rules in the level
            if(n_rules>0){
                auto res1 = nt2phrase(lvl_rules[i]);
                auto res2 = nt2phrase(lvl_rules[i + 1] - 1);
                tot_sym = res2.second-res1.first+1;
                std::cout<<pad_string<<"  Level " << (i + 1) << ": number of rules: " << n_rules << ",  number of symbols: " << tot_sym << std::endl;
            }
        }

        if(has_cg_rules){
            std::cout <<pad_string<<"Number of collage rules: " << n_cg_rules <<std::endl;
        }

        if(has_rl_rules){
            std::cout <<pad_string<<"Number of run-length rules: " << run_len_nt.second <<std::endl;
        }

        std::cout <<pad_string<<"Length of the compressed sequence (start symbol's rule): " << c<<std::endl;
    }

    void set_samp_rate(size_t new_samp_rate){
        samp_rate = new_samp_rate;
        if(has_rand_access){
        }
    }

    void get_prefix(size_t sym, size_t len, std::vector<size_t>& dc_string){

    }

    void get_suffix(size_t sym, size_t len, std::vector<size_t>& dc_string){

    }

    [[nodiscard]] uint8_t access_pos(size_t sym, size_t idx) const {

        std::stack<size_t> stack;
        stack.push(sym);
        size_t offset, u;
        assert(idx<rule_exp[sym]);

        while(!stack.empty()){

            sym = stack.top();
            stack.pop();

            auto range = nt2phrase(sym);
            u=0, offset=0;
            while(offset+rule_exp[rules[range.first+u]]<=idx){
                offset += rule_exp[rules[range.first+u]];
                u++;
            }
            idx -=offset;

            sym = rules[range.first+u];

            if(sym>max_tsym){
                stack.push(sym);
            }
        }
        return (uint8_t) sym;
    }

    //returns longest common prefix between exp(sym_a)[idx..] and exp(sym_b)[idx..]
    [[nodiscard]] long longest_com_pref(size_t sym_a, size_t sym_b, size_t idx) const {

        if(rule_exp[sym_a]<=idx || rule_exp[sym_b]<=idx){
            return -1;
        }

        std::stack<size_t> stack_a;
        stack_a.push(sym_a);
        size_t idx_a=idx, u, offset;
        while(idx_a>0){
            sym_a = stack_a.top();
            stack_a.pop();

            auto range = nt2phrase(sym_a);
            u=0, offset=0;
            while(offset+rule_exp[rules[range.first+u]]<=idx_a){
                offset += rule_exp[rules[range.first+u]];
                u++;
            }
            idx_a -=offset;

            size_t first = range.first+u;
            for(size_t j = range.second;j>=first;j--){
                stack_a.push(rules[j]);
            }
        }

        std::stack<size_t> stack_b;
        stack_b.push(sym_b);
        size_t idx_b=idx;
        while(idx_b>0){
            sym_b = stack_b.top();
            stack_b.pop();

            auto range = nt2phrase(sym_b);
            u=0, offset=0;
            while(offset+rule_exp[rules[range.first+u]]<=idx_b){
                offset += rule_exp[rules[range.first+u]];
                u++;
            }
            idx_b -=offset;

            size_t first = range.first+u;
            for(size_t j = range.second;j>=first;j--){
                stack_b.push(rules[j]);
            }
        }

        size_t p_level_a, p_level_b;
        long n_eq=0;

        while(!stack_a.empty() && !stack_b.empty()){

            sym_a = stack_a.top();
            sym_b = stack_b.top();
            p_level_a = parsing_level(stack_a.top());
            p_level_b = parsing_level(stack_b.top());

            if(sym_a==sym_b){
                n_eq += sym_a<=max_tsym;
                stack_a.pop();
                stack_b.pop();
            }else {

                if(sym_a<=max_tsym && sym_b<=max_tsym){
                    break;
                }

                if(p_level_b<=p_level_a){
                    stack_a.pop();
                    auto range= nt2phrase(sym_a);
                    for(size_t j=range.second; j>=range.first;j--){
                        stack_a.push(rules[j]);
                    }
                }

                if(p_level_a<=p_level_b){
                    stack_b.pop();
                    auto range= nt2phrase(sym_b);
                    for(size_t j=range.second; j>=range.first;j--){
                        stack_b.push(rules[j]);
                    }
                }
            }
        }
        return n_eq;
    }

    //returns true if exp(sym_a)[idx..] <_{lex} exp(sym_b)[idx..], or false otherwise
    [[nodiscard]] bool substring_lex_comp(size_t sym_a, size_t sym_b, size_t idx) const {

        if(rule_exp[sym_a]<=idx || rule_exp[sym_b]<=idx){
            return rule_exp[sym_a]<rule_exp[sym_b];
        }

        std::stack<size_t> stack_a;
        stack_a.push(sym_a);
        size_t idx_a=idx, u, offset;
        while(idx_a>0){
            sym_a = stack_a.top();
            stack_a.pop();

            auto range = nt2phrase(sym_a);
            u=0, offset=0;
            while(offset+rule_exp[rules[range.first+u]]<=idx_a){
                offset += rule_exp[rules[range.first+u]];
                u++;
            }
            idx_a -=offset;

            size_t first = range.first+u;
            for(size_t j = range.second;j>=first;j--){
                stack_a.push(rules[j]);
            }
        }

        std::stack<size_t> stack_b;
        stack_b.push(sym_b);
        size_t idx_b=idx;
        while(idx_b>0){
            sym_b = stack_b.top();
            stack_b.pop();

            auto range = nt2phrase(sym_b);
            u=0, offset=0;
            while(offset+rule_exp[rules[range.first+u]]<=idx_b){
                offset += rule_exp[rules[range.first+u]];
                u++;
            }
            idx_b -=offset;

            size_t first = range.first+u;
            for(size_t j = range.second;j>=first;j--){
                stack_b.push(rules[j]);
            }
        }

        //std::cout<<"We will compare "<<stack_a.top()<<" "<<stack_b.top()<<std::endl;
        size_t p_level_a, p_level_b;

        while(!stack_a.empty() && !stack_b.empty()){

            sym_a = stack_a.top();
            sym_b = stack_b.top();
            p_level_a = parsing_level(stack_a.top());
            p_level_b = parsing_level(stack_b.top());

            if(sym_a==sym_b){
                stack_a.pop();
                stack_b.pop();
            }else {

                if(sym_a<=max_tsym && sym_b<=max_tsym){
                    //std::cout<<"They differ in "<<sym_a<<" "<<sym_b<<std::endl;
                    return sym_a<sym_b;
                }

                if(p_level_b<=p_level_a){
                    stack_a.pop();
                    auto range= nt2phrase(sym_a);
                    for(size_t j=range.second; j>=range.first;j--){
                        stack_a.push(rules[j]);
                    }
                }

                if(p_level_a<=p_level_b){
                    stack_b.pop();
                    auto range= nt2phrase(sym_b);
                    for(size_t j=range.second; j>=range.first;j--){
                        stack_b.push(rules[j]);
                    }
                }
            }
        }
        return rule_exp[sym_a]<rule_exp[sym_b];
    }

    void print_parse_tree(size_t sym) const {

        std::queue<std::tuple<size_t, bool, bool>> queue;
        std::tuple<size_t, bool, bool> dc_sym;
        queue.emplace(sym, true, true);

        while(!queue.empty()) {

            dc_sym = queue.front();
            queue.pop();

            sym = std::get<0>(dc_sym);
            bool is_last_sib = std::get<1>(dc_sym);
            bool is_left_branch = std::get<2>(dc_sym);
            std::cout<<sym<<" ";

            if(is_last_sib){
                std::cout<<"| ";
            }

            if(is_left_branch){
                std::cout<<" "<<std::endl;
            }

            if(!is_terminal(sym)){
                auto range = nt2phrase(sym);
                for(size_t j=range.first; j<=range.second;j++){
                    queue.emplace(pos2symbol(j), j==range.second, (j==range.second && is_left_branch));
                }
            }
        }
    }

    size_t in_memory_decompression(size_t sym, std::string& dc_string) const {

        std::stack<size_t> stack;
        stack.push(sym);
        size_t exp_len = 0;

        while(!stack.empty()) {

            sym = stack.top();
            stack.pop();

            if(is_terminal(sym)){
                dc_string.push_back((char)sym);
                exp_len++;
            }else{
                auto range = nt2phrase(sym);

                if constexpr (is_rl) {
                    if(is_rl_sym(sym)) {
                        assert(range.second - range.first + 1 == 2);
                        size_t len = pos2symbol(range.second);
                        for (size_t j = 0; j < len; j++) {
                            stack.emplace(pos2symbol(range.first));
                        }
                        continue;
                    }
                }

                /*if constexpr (is_cg) {
                    if(rules[range.first]>=cg_mark){
                        assert(range.second-range.first+1==3);
                        if(rules[range.first] == cg_mark){
                            get_prefix(rules[range.first+1], rules[range.first+2], dc_string);
                        }else{
                            get_suffix(rules[range.first+1], rules[range.first+2], dc_string);
                        }
                        continue;
                    }
                }*/

                for(size_t j=range.second+1; j-->range.first;){
                    stack.emplace(pos2symbol(j));
                }
            }
        }
        return exp_len;
    }

    [[nodiscard]] size_t exp_size(size_t sym) const {
        if(is_terminal(sym)) return 1;
        auto range = nt2phrase(sym);

        size_t pos = range.second/samp_rate;
        size_t exp = sampled_exp.read(pos);
        size_t last_sampled = pos*samp_rate;

        for(size_t i=last_sampled;i<=range.second;i++){
            exp += exp_size(rules.read(i));
        }
        return exp;
    }
};

#endif //LCG_GRAMMAR_H
