//
// Created by Diaz, Diego on 3.7.2023.
//

#ifndef LCG_GRAMMAR_H
#define LCG_GRAMMAR_H
#include "cds/cdt_common.hpp"
#include "cds/int_array.h"
#include "cds/file_streams.hpp"
#include "hashing.h"
#include "cds/utils.h"

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
    std::vector<size_t>                lvl_rules; // number of rules generated in every parsing round
    std::vector<size_t>                lvl_size; // number of rules generated in every parsing round



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

struct lc_gram_t {

    size_t                             n{}; //n: number of symbols in the original text
    size_t                             r{}; //r: number of grammar symbols (nter + ter)
    size_t                             c{}; //c: length of the right-hand of the start symbol
    size_t                             g{}; //g: sum of the rules' right-hand sides
    size_t                             s{}; //s: number of strings
    size_t                             longest_str{}; //length of the longest string encoded in the grammar
    size_t                             max_tsym{}; //highest terminal symbol
    uint8_t                            sep_tsym{}; //separator symbol in the collection. This is 0 if the text is a single string
    bool                               is_simplified=false;
    bool                               has_rl_rules=false;
    bool                               has_rand_access=false;

    std::vector<uint8_t>               terminals; //set of terminals
    std::vector<hashing>               par_functions; //list of hash functions from which the grammar was constructed
    std::vector<size_t>                str_boundaries; // start position of every string in the compressed string
    std::vector<size_t>                lvl_rules; //number of rules generated in every round of locally-consistent parsing

    int_array<size_t>                  rules; //concatenated set of grammar rules
    int_array<size_t>                  rl_ptr; //pointer in "rules" to the leftmost symbol of each rule
    int_array<size_t>                  rule_exp;// length of the nt expansions
    int_array<size_t>                  sampled_exp;// sampled nt expansions in "rules"
    size_t                             samp_rate=4; //sampling rate to sample nt expansions in "rules"

    std::pair<size_t, size_t>          run_len_nt; //first run-length rule and total number of run-length rules

    lc_gram_t()= default;

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
        written_bytes +=serialize_elm(ofs, sep_tsym);
        written_bytes +=serialize_elm(ofs, is_simplified);
        written_bytes +=serialize_elm(ofs, has_rl_rules);
        written_bytes +=serialize_elm(ofs, has_rand_access);
        written_bytes +=serialize_elm(ofs, run_len_nt.first);
        written_bytes +=serialize_elm(ofs, run_len_nt.second);
        written_bytes +=serialize_plain_vector(ofs, lvl_rules);
        written_bytes +=serialize_plain_vector(ofs, par_functions);

        written_bytes +=serialize_plain_vector(ofs, terminals);
        written_bytes +=serialize_plain_vector(ofs, str_boundaries);
        written_bytes +=rules.serialize(ofs);
        written_bytes +=rl_ptr.serialize(ofs);
        written_bytes +=rule_exp.serialize(ofs);
        written_bytes +=sampled_exp.serialize(ofs);
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
        load_elm(ifs, sep_tsym);
        load_elm(ifs, is_simplified);
        load_elm(ifs, has_rl_rules);
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
        rule_exp.load(ifs);
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

        std::cout<<pad_string<<"Number of compressed symbols: "<<n<<std::endl;
        std::cout<<pad_string<<"Number of compressed strings: "<<s<<" ("<<report_space((off_t)pt_str_bytes)<<" in pointers)"<<std::endl;
        std::cout<<pad_string<<"Separator symbol:             "<<(int)sep_tsym<<std::endl;
        std::cout<<pad_string<<"Longest string:               "<<longest_str<<std::endl;
        std::cout<<pad_string<<"Number of terminals:          "<<n_ter<<std::endl;
        std::cout<<pad_string<<"Number of non-terminals:      "<<n_nter<<" ("<<report_space((off_t)pt_bytes)<<" in pointers)"<<std::endl;
        std::cout<<pad_string<<"Grammar size:                 "<<g<<" ("<<report_space((off_t)g_bytes)<<")"<<std::endl;
        std::cout<<pad_string<<"Approx. compression ratio:    "<<comp_ratio<<std::endl;
        std::cout<<pad_string<<"Simplified:                   "<<(is_simplified ? "yes" : "no")<<std::endl;
        std::cout<<pad_string<<"Run-len rules:                "<<(has_rl_rules ? "yes" : "no")<<std::endl;
        std::cout<<pad_string<<"Random access support:        "<<(has_rand_access? "yes" : "no");
        if(has_rand_access){
            auto ras_bytes = INT_CEIL((rule_exp.size()*rule_exp.width()+ sampled_exp.size()+sampled_exp.width()), 8);
            std::cout<<" ("<<report_space((off_t)ras_bytes)<<" in data)"<<std::endl;
            std::cout<<pad_string<<"  sampling rate for the expansions: 1/"<<samp_rate<<std::endl;
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

    //TODO implement this
    /*
    void access(size_t str_id, size_t start, size_t end, std::string& output){
    }*/
};

#endif //LCG_GRAMMAR_H
