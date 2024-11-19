//
// Created by Diaz, Diego on 3.7.2023.
//

#ifndef LCG_GRAMMAR_H
#define LCG_GRAMMAR_H
#include <stack>
#include <queue>
#include <random>

#include "cds/cdt_common.hpp"
#include "cds/int_array.h"
#include "cds/utils.h"
#include "xxhash.h"

#define STR_EXP true
#define RULE_EXP false

/** exp_data is a struct storing expansion information. Let B -> A_1 A_2 ... A_x be a grammar rule, and
 * let exp(B)[exp_s..exp_e] be a substring of exp(B) covered by the substring A_i A_{i+1} .. A_{i+k} of the
 * right-hand side (rhs) of B's rule. Here, "covered" means that if we remove A_i or A_{i+k}, the remaining substring
 * no longer contains exp(X)[exp_s..exp_e]. The meaning of each field in the struct is described below. Additionally,
 * when B -> A_1 A_2 ... A_x is a sequence of x copies of some symbol A, we encode the rule as B -> (A, x) in the grammar.
 * In that case, the fields indicating bit positions are relative to the uncompressed rule B -> A_1=A A_2=A ... A_x=A.
 * The flag is_rl indicates when B's rule is run-length encoded.
 **/
struct exp_data{
    off_t bp_rhs_s; //leftmost bit of A_{1} within the rhs of X's rule
    off_t bp_rhs_e; //leftmost bit of A_{x} within the rhs of X's rule
    off_t bp_exp_s; //leftmost bit of A_{i} within the rhs of X's rule
    off_t bp_exp_e; //leftmost bit of A_{i+k} within the rhs of X's rule
    bool is_rl;
};

struct str_coord_type{
    size_t str;
    off_t start;
    off_t end;
    str_coord_type(size_t str_, off_t start_, off_t end_): str(str_), start(start_), end(end_){};
};

template<bool is_cg=false, bool is_rl=false, bool has_ra=false>
struct lc_gram_t {

    size_t n{}; //n: number of symbols in the original text
    size_t r{}; //r: number of grammar symbols (nter + ter)
    size_t c{}; //c: length of the right-hand of the start symbol
    size_t g{}; //g: sum of the rules' right-hand sides
    size_t s{}; //s: number of strings
    size_t max_tsym{}; //highest terminal symbol
    size_t min_nter{}; //smallest nonterminal symbol
    size_t cg_mark{};
    size_t n_cg_rules{};
    uint64_t par_seed{}; //list of hash functions from which the grammar was constructed
    uint8_t r_bits = 0;
    uint8_t r_samp_bits = 0;
    uint8_t str_samp_bits = 0;

    uint8_t sep_tsym{}; //separator symbol in the collection. This is 0 if the text is a single string
    bool is_simplified = false;
    const static bool has_rl_rules = is_rl;
    const static bool has_cg_rules = is_cg;
    const static bool has_rand_access = has_ra;

    std::vector<uint8_t> terminals; //set of terminals
    std::vector<size_t> str_boundaries; // start position of every string in the compressed string
    std::vector<size_t> lvl_rules; //number of rules generated in every round of locally-consistent parsing

    bitstream<size_t> rule_stream;
    int_array<size_t> rl_ptr; //pointer in "rules" to the leftmost symbol of each rule
    size_t rl_samp_rate = 4;//sampling rate to sample nt expansions in "rules"
    size_t str_samp_rate = 8;//sampling rate to sample nt expansions in "rules"
    std::pair<size_t, size_t> run_len_nt{0, 0};//first run-length rule and total number of run-length rules

    lc_gram_t() = default;

    size_t serialize(std::ofstream &ofs) {
        size_t written_bytes = 0;
        written_bytes += serialize_elm(ofs, has_rl_rules);
        written_bytes += serialize_elm(ofs, has_cg_rules);
        written_bytes += serialize_elm(ofs, has_rand_access);
        written_bytes += serialize_elm(ofs, n);
        written_bytes += serialize_elm(ofs, r);
        written_bytes += serialize_elm(ofs, g);
        written_bytes += serialize_elm(ofs, c);
        written_bytes += serialize_elm(ofs, s);
        written_bytes += serialize_elm(ofs, rl_samp_rate);
        written_bytes += serialize_elm(ofs, str_samp_rate);
        written_bytes += serialize_elm(ofs, r_bits);
        written_bytes += serialize_elm(ofs, r_samp_bits);
        written_bytes += serialize_elm(ofs, str_samp_bits);
        written_bytes += serialize_elm(ofs, max_tsym);
        written_bytes += serialize_elm(ofs, min_nter);
        written_bytes += serialize_elm(ofs, sep_tsym);
        written_bytes += serialize_elm(ofs, is_simplified);
        written_bytes += serialize_elm(ofs, cg_mark);
        written_bytes += serialize_elm(ofs, n_cg_rules);
        written_bytes += serialize_elm(ofs, run_len_nt.first);
        written_bytes += serialize_elm(ofs, run_len_nt.second);
        written_bytes += serialize_elm(ofs, par_seed);
        written_bytes += serialize_plain_vector(ofs, lvl_rules);

        written_bytes += serialize_plain_vector(ofs, terminals);
        written_bytes += serialize_plain_vector(ofs, str_boundaries);
        written_bytes += rl_ptr.serialize(ofs);
        written_bytes += rule_stream.serialize(ofs);

        return written_bytes;
    }

    void load_metadata(std::ifstream &ifs) {
        bool tmp;
        load_elm(ifs, tmp);
        assert(tmp == has_rl_rules);
        load_elm(ifs, tmp);
        assert(tmp == has_cg_rules);
        load_elm(ifs, tmp);
        assert(tmp == has_rand_access);

        load_elm(ifs, n);
        load_elm(ifs, r);
        load_elm(ifs, g);
        load_elm(ifs, c);
        load_elm(ifs, s);
        load_elm(ifs, rl_samp_rate);
        load_elm(ifs, str_samp_rate);
        load_elm(ifs, r_bits);
        load_elm(ifs, r_samp_bits);
        load_elm(ifs, str_samp_bits);
        load_elm(ifs, max_tsym);
        load_elm(ifs, min_nter);
        load_elm(ifs, sep_tsym);
        load_elm(ifs, is_simplified);
        load_elm(ifs, cg_mark);
        load_elm(ifs, n_cg_rules);
        load_elm(ifs, run_len_nt.first);
        load_elm(ifs, run_len_nt.second);
        load_elm(ifs, par_seed);
        load_plain_vector(ifs, lvl_rules);
    }

    template<class gram_type>
    void swap(gram_type &other) {
        std::swap(n, other.n);
        std::swap(r, other.r);
        std::swap(g, other.g);
        std::swap(c, other.c);
        std::swap(s, other.s);
        std::swap(rl_samp_rate, other.rl_samp_rate);
        std::swap(str_samp_rate, other.str_samp_rate);
        std::swap(r_bits, other.r_bits);
        std::swap(r_samp_bits, other.r_samp_bits);
        std::swap(str_samp_bits, other.str_samp_bits);
        std::swap(max_tsym, other.max_tsym);
        std::swap(min_nter, other.min_nter);
        std::swap(sep_tsym, other.sep_tsym);
        std::swap(is_simplified, other.is_simplified);
        std::swap(cg_mark, other.cg_mark);
        std::swap(n_cg_rules, other.n_cg_rules);
        std::swap(run_len_nt, other.run_len_nt);
        std::swap(par_seed, other.par_seed);
        std::swap(lvl_rules, other.lvl_rules);
        terminals.swap(other.terminals);
        str_boundaries.swap(other.str_boundaries);
        rl_ptr.swap(other.rl_ptr);
        rule_stream.swap(other.rule_stream);
    }

    void load_pointers(std::ifstream &ifs) {
        load_plain_vector(ifs, terminals);
        load_plain_vector(ifs, str_boundaries);
        rl_ptr.load(ifs);
    }

    void load(std::ifstream &ifs) {
        load_metadata(ifs);
        load_pointers(ifs);
        rule_stream.load(ifs);
    }

    [[nodiscard]] inline std::pair<off_t, off_t> nt2bitrange(size_t sym) const {
        assert(sym > max_tsym);
        size_t nt = sym - max_tsym - 1;
        if constexpr (has_ra) {
            size_t start = rl_ptr.read(nt);
            size_t end = rl_ptr.read(nt + 1);
            size_t ra_bits = rule_stream.read(end - r_samp_bits, end - 1);
            end -= r_samp_bits + ra_bits;
            return {start, end - r_bits};
        } else {
            return {rl_ptr.read(nt) * r_bits, (rl_ptr.read(nt + 1) - 1) * r_bits};
        }
    }

    template<bool is_str>
    [[nodiscard]] inline std::tuple<off_t, off_t, off_t> nt2expdata(size_t sym) const {

        assert(has_ra);
        off_t rhs_start, exp_start, ra_bits, rhs_len, n_samples, bps;

        if constexpr (is_str) {
            rhs_start = str_boundaries[sym];
            exp_start = str_boundaries[sym + 1];
            ra_bits = rule_stream.read(exp_start - str_samp_bits, exp_start - 1);
            exp_start -= str_samp_bits + ra_bits;

            rhs_len = (exp_start - rhs_start) / r_bits;
            n_samples = rhs_len / str_samp_rate;
            bps = ra_bits / n_samples;
        } else {
            sym -= max_tsym + 1;

            rhs_start = rl_ptr.read(sym);
            exp_start = rl_ptr.read(sym + 1);
            ra_bits = rule_stream.read(exp_start - r_samp_bits, exp_start - 1);
            exp_start -= r_samp_bits + ra_bits;

            rhs_len = (exp_start - rhs_start) / r_bits;
            n_samples = rhs_len / rl_samp_rate;
            bps = ra_bits / (n_samples + 1);// +1 because the rules also include the value of |exp(sym)| as a sample
        }

        return {exp_start, exp_start+ra_bits, bps};
    }

    [[nodiscard]] inline bool is_terminal(size_t sym) const {
        return sym <= max_tsym;
    }

    [[nodiscard]] inline size_t n_lc_syms() const {
        return lvl_rules.back();
    }

    [[nodiscard]] inline bool has_run_length_rules() const {
        return has_rl_rules;
    }

    [[nodiscard]] inline bool is_rl_sym(size_t symbol) const {
        return run_len_nt.first <= symbol && symbol < (run_len_nt.first + run_len_nt.second);
    }

    [[nodiscard]] inline size_t first_rl_sym() const {
        return run_len_nt.first;
    }

    [[nodiscard]] inline size_t last_rl_sym() const {
        return run_len_nt.first + run_len_nt.second - 1;
    }

    [[nodiscard]] inline size_t n_terminals() const {
        return (max_tsym + 1);
    }

    [[nodiscard]] inline size_t n_nonterminals() const {
        return r - n_terminals();
    }

    [[nodiscard]] inline size_t comp_str_size() const {
        return c;
    }

    [[nodiscard]] inline off_t lc_par_tree_height() const {
        return lvl_rules.size();
    }

    inline size_t get_byte_ter(size_t sym) {
        assert(is_terminal(sym));
        if (is_simplified) {
            return terminals[sym];
        } else {
            return sym;
        }
    }

    [[nodiscard]] inline off_t parsing_level(size_t symbol) const {
        if (symbol <= max_tsym) return 0;
        for (off_t i = 0; i < (off_t) lvl_rules.size(); i++) {
            if (lvl_rules[i] <= symbol && symbol < lvl_rules[i + 1]) {
                return i + 1;
            }
        }
        return -1;
    }

    [[nodiscard]] inline size_t bitpos2symbol(size_t bit_pos) const {
        return rule_stream.read(bit_pos, bit_pos + r_bits - 1);
    }

    [[nodiscard]] inline size_t start_symbol() const {
        return r - 1;
    }

    [[nodiscard]] inline size_t n_strings() const {
        return str_boundaries.size() - 1;
    }

    [[nodiscard]] inline std::pair<off_t, off_t> str2bitrange(size_t str) const {
        if constexpr (has_ra) {
            size_t start = str_boundaries[str];
            size_t end = str_boundaries[str + 1];
            size_t ra_bits = rule_stream.read(end - str_samp_bits, end - 1);
            end -= str_samp_bits + ra_bits;
            return {start, end - r_bits};
        } else {
            return {str_boundaries[str] * r_bits, (str_boundaries[str + 1] - 1) * r_bits};
        }
    }

    void stats(size_t pad) const {

        size_t n_ter = n_terminals();
        size_t n_nter = n_nonterminals();

        auto pt_bytes = INT_CEIL((r * sym_width(g)), 8);//space of the pointers for the nonterminals
        auto g_bytes = INT_CEIL((g * r_bits), 8); //space of the expansions
        auto pt_str_bytes = (s + 1) * sizeof(size_t);  //space of the pointers to the strings

        float comp_ratio = float(n) / float(pt_bytes + g_bytes + pt_str_bytes);

        std::string pad_string(pad, ' ');
        std::cout << pad_string << "Number of compressed symbols:   " << n << " (" << report_space((off_t) n) << ")"
                  << std::endl;
        std::cout << pad_string << "Number of compressed strings:   " << s << " (" << report_space((off_t) pt_str_bytes)
                  << " in pointers)" << std::endl;
        std::cout << pad_string << "Separator symbol:               " << (int) sep_tsym << std::endl;
        std::cout << pad_string << "Number of terminals:            " << n_ter << std::endl;
        std::cout << pad_string << "Number of nonterminals:         " << n_nter << " ("
                  << report_space((off_t) pt_bytes) << " in pointers)" << std::endl;
        std::cout << pad_string << "Grammar size:                   " << g << " (" << report_space((off_t) g_bytes)
                  << ")" << std::endl;
        std::cout << pad_string << "Length of the comp. collection: " << c << std::endl;
        std::cout << pad_string << "Approx. compression ratio:      " << comp_ratio << std::endl;
        std::cout << pad_string << "Simplified:                     " << (is_simplified ? "yes" : "no") << std::endl;
        std::cout << pad_string << "Run-length rules:               " << (has_rl_rules ? "yes" : "no") << std::endl;
        std::cout << pad_string << "Collage system rules:           " << (has_cg_rules ? "yes" : "no") << std::endl;
        std::cout << pad_string << "Random access support:          " << (has_rand_access ? "yes" : "no") << std::endl;
        if (has_rand_access) {
            auto ras_bytes = rule_stream.bit_capacity() - (g * r_bits);//the cost of the expansion samples
            ras_bytes += rl_ptr.n_bits() - (r * sym_width(g));//the cost of expanding the rules' pointers
            ras_bytes = INT_CEIL(ras_bytes, 8);
            std::cout << pad_string << "  Samp. rate for non.ter exps:  1/" << rl_samp_rate << std::endl;
            std::cout << pad_string << "  Samp. rate for string exps:   1/" << str_samp_rate << std::endl;

            if(rule_stream.bit_capacity()!=0){//for the cases when we are reading from a file, and we do not have the grammar loaded
                std::cout << pad_string << "  Space overhead:               " << report_space((off_t) ras_bytes)<< std::endl;
            }
        }
    }

    void breakdown(size_t pad) {
        //assert(g==rules.size());
        //assert(r-(max_tsym+1)==(rl_ptr.size()-1));
        assert(s == str_boundaries.size() - 1);
        stats(pad);

        std::string pad_string(pad, ' ');

        size_t n_rules;
        std::cout << pad_string << "Grammar rules per level" << std::endl;
        for (size_t i = 0; i < lvl_rules.size() - 1; i++) {
            n_rules = lvl_rules[i + 1] - lvl_rules[i];//number of rules in the level
            if (n_rules > 0) {
                std::cout << pad_string << "  Level " << (i + 1) << ": number of rules: " << n_rules << std::endl;
            }
        }

        if (has_cg_rules) {
            std::cout << pad_string << "Number of collage rules: " << n_cg_rules << std::endl;
        }

        if (has_rl_rules) {
            std::cout << pad_string << "Number of run-length rules: " << run_len_nt.second << std::endl;
        }
        std::cout << pad_string << "Length of the compressed sequence (start symbol's rule): " << c << std::endl;
    }

    /*
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
     */

    void print_parse_tree(size_t sym, bool is_str=false) const {
        std::queue<std::pair<size_t, bool>> queue;
        std::pair<size_t, bool> dc_sym;
        bool next_break_known;
        if(is_str){
            auto range = str2bitrange(sym);
            for (off_t j = range.first; j <= range.second; j += r_bits) {
                bool flag = false;
                if (!next_break_known && bitpos2symbol(j) > max_tsym) {
                    next_break_known = true;
                    flag = true;
                }
                std::cout << bitpos2symbol(j) << " " << std::flush;
                queue.emplace(bitpos2symbol(j), flag);
            }
            std::cout << "| " << std::flush;
        }else{
            queue.emplace(sym, true);
            next_break_known=true;
            std::cout << sym << " |" << std::flush;
        }


        while (!queue.empty()) {

            dc_sym = queue.front();
            queue.pop();

            sym = std::get<0>(dc_sym);

            if (sym > max_tsym) {
                auto range = nt2bitrange(sym);

                if (is_rl_sym(sym)) {
                    assert(range.second - range.first == r_bits);
                    size_t rl_sym = bitpos2symbol(range.first);
                    size_t len = bitpos2symbol(range.second);

                    if (std::get<1>(dc_sym)) {
                        next_break_known = false;
                        std::cout << " " << std::endl;
                    }
                    for (size_t j = 0; j < len; j++) {
                        bool flag = false;
                        if (!next_break_known && rl_sym > max_tsym) {
                            next_break_known = true;
                            flag = true;
                        }
                        std::cout << rl_sym << " " << std::flush;
                        queue.emplace(rl_sym, flag);
                    }
                    std::cout << "| ";
                } else {
                    if (std::get<1>(dc_sym)) {
                        next_break_known = false;
                        std::cout << " " << std::endl;
                    }
                    for (off_t j = range.first; j <= range.second; j += r_bits) {
                        bool flag = false;
                        if (!next_break_known && bitpos2symbol(j) > max_tsym) {
                            next_break_known = true;
                            flag = true;
                        }
                        std::cout << bitpos2symbol(j) << " " << std::flush;
                        queue.emplace(bitpos2symbol(j), flag);
                    }
                    std::cout << "| " << std::flush;
                }
            }
        }
        std::cout<<""<<std::endl;
    }

    void im_nt_decomp(size_t sym, std::string &dc_string) const {

        std::stack<size_t> stack;
        stack.push(sym);

        while (!stack.empty()) {
            sym = stack.top();
            stack.pop();
            if (sym <= max_tsym) {
                dc_string.push_back(terminals[sym]);
            } else {
                auto range = nt2bitrange(sym);
                if constexpr (is_rl) {
                    if (is_rl_sym(sym)) {
                        assert(range.second - range.first == r_bits);
                        size_t len = bitpos2symbol(range.second);
                        for (size_t j = 0; j < len; j++) {
                            stack.emplace(bitpos2symbol(range.first));
                        }
                        continue;
                    }
                }
                for (off_t j = range.second; j >= range.first; j -= r_bits) {
                    stack.emplace(bitpos2symbol(j));
                }
            }
        }
    }

    void im_str_decomp(size_t str, std::string &dc_string) const {
        assert(str < n_strings());
        size_t sym;

        auto range = str2bitrange(str);
        std::stack<size_t> stack;
        for (off_t j = range.second; j >= range.first; j -= r_bits) {
            stack.emplace(bitpos2symbol(j));
        }

        while (!stack.empty()) {
            sym = stack.top();
            stack.pop();
            if (sym <= max_tsym) {
                dc_string.push_back(terminals[sym]);
            } else {
                range = nt2bitrange(sym);
                if constexpr (is_rl) {
                    if (is_rl_sym(sym)) {
                        assert(range.second - range.first == r_bits);
                        size_t len = bitpos2symbol(range.second);
                        for (size_t j = 0; j < len; j++) {
                            stack.emplace(bitpos2symbol(range.first));
                        }
                        continue;
                    }
                }
                for (off_t j = range.second; j >= range.first; j -= r_bits) {
                    stack.emplace(bitpos2symbol(j));
                }
            }
        }
    }

    void im_str_rand_access(size_t sym, off_t exp_start, off_t exp_end, std::string &dc_string) const {

        assert(exp_start < exp_end);
        dc_string.clear();

        exp_data rhs_data{};
        exp_search_range<STR_EXP>(sym, exp_start, exp_end, rhs_data);
        size_t lm_sym = rule_stream.read(rhs_data.bp_exp_s, rhs_data.bp_exp_s + r_bits - 1);
        off_t k = (rhs_data.bp_exp_e - rhs_data.bp_exp_s + r_bits) / r_bits;
        while (k == 1) {
            exp_search_range<RULE_EXP>(lm_sym, exp_start, exp_end, rhs_data);
            lm_sym = rule_stream.read(rhs_data.bp_exp_s, rhs_data.bp_exp_s + r_bits - 1);
            k = (rhs_data.bp_exp_e - rhs_data.bp_exp_s + r_bits) / r_bits;
        }

        size_t rm_sym;
        std::stack<uint64_t> exp_stack;
        if(rhs_data.is_rl) {
            for (off_t u = rhs_data.bp_exp_e - r_bits; u > rhs_data.bp_exp_s; u -= r_bits) {
                exp_stack.push(lm_sym);
            }
            rm_sym = lm_sym;
        } else {
            for (off_t u = rhs_data.bp_exp_e - r_bits; u > rhs_data.bp_exp_s; u -= r_bits) {
                exp_stack.push(rule_stream.read(u, u + r_bits - 1));
            }
            rm_sym = rule_stream.read(rhs_data.bp_exp_e, rhs_data.bp_exp_e + r_bits - 1);
        }

        while (lm_sym > max_tsym) {
            exp_search<RULE_EXP>(exp_start, lm_sym, rhs_data);
            if (rhs_data.is_rl) {
                for (off_t u = rhs_data.bp_rhs_e; u > rhs_data.bp_exp_s; u -= r_bits) {
                    exp_stack.push(lm_sym);
                }
            } else {
                for (off_t u = rhs_data.bp_rhs_e; u > rhs_data.bp_exp_s; u -= r_bits) {
                    exp_stack.push(rule_stream.read(u, u + r_bits - 1));
                }
            }
        }
        dc_string.push_back(terminals[lm_sym]);

        while (!exp_stack.empty()) {
            sym = exp_stack.top();
            exp_stack.pop();

            if (sym <= max_tsym) {
                dc_string.push_back(terminals[sym]);
            } else {
                auto range = nt2bitrange(sym);
                if (is_rl_sym(sym)) {
                    assert(range.second - range.first == r_bits);
                    size_t len = bitpos2symbol(range.second);
                    for (size_t j = 0; j < len; j++) {
                        exp_stack.emplace(bitpos2symbol(range.first));
                    }
                } else {
                    for (off_t j = range.second; j >= range.first; j -= r_bits) {
                        exp_stack.emplace(bitpos2symbol(j));
                    }
                }
            }
        }

        size_t stack_sym;
        while (rm_sym > max_tsym) {
            exp_search<RULE_EXP>(exp_end, rm_sym, rhs_data);
            if(rhs_data.is_rl) {
                for (off_t u = rhs_data.bp_exp_s - r_bits; u >= rhs_data.bp_rhs_s; u -= r_bits) {
                    exp_stack.push(rm_sym);
                }
            } else {
                for (off_t u = rhs_data.bp_exp_s - r_bits; u >= rhs_data.bp_rhs_s; u -= r_bits) {
                    exp_stack.push(rule_stream.read(u, u + r_bits - 1));
                }
            }

            while (!exp_stack.empty()) {
                stack_sym = exp_stack.top();
                exp_stack.pop();
                if (stack_sym <= max_tsym) {
                    dc_string.push_back(terminals[stack_sym]);
                } else {
                    auto range = nt2bitrange(stack_sym);
                    if (is_rl_sym(stack_sym)) {
                        assert(range.second - range.first == r_bits);
                        size_t len = bitpos2symbol(range.second);
                        for (size_t j = 0; j < len; j++) {
                            exp_stack.emplace(bitpos2symbol(range.first));
                        }
                        continue;
                    }
                    for (off_t j = range.second; j >= range.first; j -= r_bits) {
                        exp_stack.emplace(bitpos2symbol(j));
                    }
                }
            }
        }
        dc_string.push_back(terminals[rm_sym]);
    }

    [[nodiscard]] size_t im_sym_rand_access(size_t sym, off_t pos) const {
        assert(has_rand_access);
        exp_data rhs{};
        exp_search<STR_EXP>(pos, sym, rhs);
        while (sym > max_tsym) {
            exp_search<RULE_EXP>(pos, sym, rhs);
        }
        return terminals[sym];
    }

    /**
     * sym is a grammar symbol (nonterminal or string id)
     * Given the rule sym -> A_1 A_2 ... A_x, this function returns the minimum substring A_i A_{i+1} .. A_{i+k} in
     * the right-hand side of that rule covering exp(sym)[exp_s..exp_e] (exp_e is inclusive).
     * The function produces the following output:
     *      exp_s: the relative position of exp(sym)[exp_s] within exp(A_i)
     *      exp_e: the relative position of exp(sym)[exp_e] within exp(A_{i+k})
     *      rhs_data: a struct of the type "exp_data" storing relevant information about the expansion
    **/
    template<bool is_str>
    inline void exp_search_range(size_t sym, off_t &exp_s, off_t &exp_e, exp_data& rhs_data) const {

        assert(exp_s <= exp_e);
        off_t rhs_start, exp_start, ra_bits, last_rhs_idx, n_samples, acc_exp, rhs_idx_l, rhs_idx_r, bps;

        if constexpr (is_str) {
            rhs_start = str_boundaries[sym];
            exp_start = str_boundaries[sym + 1];
            ra_bits = rule_stream.read(exp_start - str_samp_bits, exp_start - 1);
            exp_start -= str_samp_bits + ra_bits;
            last_rhs_idx = (exp_start - rhs_start) / r_bits;
            n_samples = last_rhs_idx / str_samp_rate;
            last_rhs_idx--;
            bps = ra_bits / n_samples;
        } else {

            bool is_rl_rule = is_rl_sym(sym);
            sym -= max_tsym + 1;

            rhs_start = rl_ptr.read(sym);
            exp_start = rl_ptr.read(sym + 1);
            ra_bits = rule_stream.read(exp_start - r_samp_bits, exp_start - 1);
            exp_start -= r_samp_bits + ra_bits;

            if (is_rl_rule) {//the right-hand side A_1 A_2 ... A_x = A^{l} is a sequence of l copies of symbol A
                off_t rl_len = rule_stream.read(rhs_start + r_bits, rhs_start + (2 * r_bits) - 1);//read the value of l
                acc_exp = rule_stream.read(exp_start, exp_start + ra_bits - 1);//read the value of |exp(A)|*l
                acc_exp /= rl_len;// read |exp(A)|
                //compute how many copies of "A" we need to cover exp(sym)[exp_s..exp_e]
                off_t n_copies = (exp_e / acc_exp) - (exp_s / acc_exp) + 1;
                exp_s %= acc_exp;
                exp_e %= acc_exp;

                rhs_data.bp_rhs_s = rhs_start;
                rhs_data.bp_exp_s = rhs_start;
                rhs_data.bp_exp_e = rhs_start + (n_copies - 1) * r_bits;
                rhs_data.bp_rhs_e = rhs_start + (rl_len-1) * r_bits;
                rhs_data.is_rl = true;
                return;
            }

            last_rhs_idx = (exp_start - rhs_start) / r_bits;
            n_samples = last_rhs_idx / rl_samp_rate;
            last_rhs_idx--;
            bps = ra_bits / (n_samples + 1);// +1 because the rules also include the value of |exp(sym)| as a sample
        }

        //binary search of exp_s in the expansion samples
        off_t left = 0, exp_m = 0, exp_n, bit_pos, middle = -1;
        auto right = (off_t) n_samples - 1;
        while (left <= right) {
            middle = left + (right - left) / 2;
            bit_pos = exp_start + middle * bps;
            exp_m = rule_stream.read(bit_pos, bit_pos + bps - 1);
            if (exp_m > exp_s) {
                right = middle - 1;
                middle = -1;
                exp_m = 0;
                continue;
            }

            //exp_m<=pos
            bit_pos += bps;
            exp_n = rule_stream.read(bit_pos, bit_pos + bps - 1);
            if ((middle + 1) == n_samples || exp_s < exp_n) {
                //exp_m<=pos && pos<exp_n
                break;
            }
            //exp_n=<pos
            left = middle + 1;
        }

        acc_exp = exp_m;
        if constexpr (is_str) {
            rhs_idx_l = (middle + 1) * str_samp_rate;
        } else {
            rhs_idx_l = (middle + 1) * rl_samp_rate;
        }

        size_t rhs_bit_pos = rhs_start + rhs_idx_l * r_bits;
        sym = rule_stream.read(rhs_bit_pos, rhs_bit_pos + r_bits - 1);
        off_t e_len = exp_len(sym);
        acc_exp += e_len;

        while (acc_exp <= exp_s && rhs_idx_l < last_rhs_idx) {
            rhs_idx_l++;
            rhs_bit_pos += r_bits;
            sym = rule_stream.read(rhs_bit_pos, rhs_bit_pos + r_bits - 1);
            e_len = exp_len(sym);
            acc_exp += e_len;
        }
        assert(acc_exp > exp_s);
        exp_s -= acc_exp - e_len;

        //binary search of exp_e in the expansion samples
        if constexpr (is_str) {
            left = rhs_idx_l / str_samp_rate;
        } else {
            left = rhs_idx_l / rl_samp_rate;
        }
        right = (off_t) n_samples - 1;
        while (left <= right) {
            middle = left + (right - left) / 2;
            bit_pos = exp_start + middle * bps;
            exp_m = rule_stream.read(bit_pos, bit_pos + bps - 1);
            if (exp_m > exp_e) {
                right = middle - 1;
                middle = -1;
                exp_m = 0;
                continue;
            }
            //exp_m<=pos
            bit_pos += bps;
            exp_n = rule_stream.read(bit_pos, bit_pos + bps - 1);
            if ((middle + 1) == n_samples || exp_e < exp_n) {
                //exp_m<=pos && pos<exp_n
                break;
            }
            //exp_n=<pos
            left = middle + 1;
        }
        acc_exp = exp_m;
        if constexpr (is_str) {
            rhs_idx_r = (middle + 1) * str_samp_rate;
        } else {
            rhs_idx_r = (middle + 1) * rl_samp_rate;
        }
        rhs_bit_pos = rhs_start + rhs_idx_r * r_bits;
        sym = rule_stream.read(rhs_bit_pos, rhs_bit_pos + r_bits - 1);
        e_len = exp_len(sym);
        acc_exp += e_len;

        while (acc_exp <= exp_e && rhs_idx_r < last_rhs_idx) {
            rhs_idx_r++;
            rhs_bit_pos += r_bits;
            sym = rule_stream.read(rhs_bit_pos, rhs_bit_pos + r_bits - 1);
            e_len = exp_len(sym);
            acc_exp += e_len;
        }
        assert(acc_exp > exp_e);
        exp_e -= acc_exp - e_len;

        rhs_data.bp_rhs_s = rhs_start;
        rhs_data.bp_exp_s = rhs_start + rhs_idx_l * r_bits;
        rhs_data.bp_exp_e = rhs_start + rhs_idx_r * r_bits;
        rhs_data.bp_rhs_e = exp_start - r_bits;
        rhs_data.is_rl = false;
    }


    /**
     * Input:
     *      sym: a grammar symbol (nonterminal or string id)
     *      exp_pos: the position for the symbol exp(sym)[exp_pos]
     * Given the rule sym -> A_1 A_2 ... A_x, this function returns symbol A_i in
     * the right-hand side of that rule containing exp(sym)[exp_pos]
     * Output:
     *      exp_pos: the relative position of exp(sym)[exp_pos] within exp(A_i)
     *      sym: the value of A_i
     *      bit_pos_s: leftmost bit of A_1 within the grammar encoding
     *      bit_pos_exp: leftmost bit of A_i within the grammar encoding
     *      bit_pos_e: leftmost bit of A_x within the grammar encoding
     * Caveat: when A_1 A_2 ... A_x is a run A^{l} of l copies of some symbol A, then this rule is encoded as
     * A l. Therefore, bit_pos_s=bit_pos_exp=bit_pos_e is the leftmost bit of A
     * within (A l) in the grammar encoding
    **/
    template<bool is_str>
    inline void exp_search(off_t &exp_pos, size_t &sym, exp_data& rhs_data) const {

        off_t rhs_start, exp_start, ra_bits, last_rhs_idx, n_samples, acc_exp, rhs_idx, bps;

        if constexpr (is_str) {
            rhs_start = str_boundaries[sym];
            exp_start = str_boundaries[sym + 1];
            ra_bits = rule_stream.read(exp_start - str_samp_bits, exp_start - 1);
            exp_start -= str_samp_bits + ra_bits;
            last_rhs_idx = (exp_start - rhs_start) / r_bits;
            n_samples = last_rhs_idx / str_samp_rate;
            last_rhs_idx--;
            bps = ra_bits / n_samples;
        } else {
            bool is_rl_rule = is_rl_sym(sym);
            sym -= max_tsym + 1;

            rhs_start = rl_ptr.read(sym);
            exp_start = rl_ptr.read(sym + 1);
            ra_bits = rule_stream.read(exp_start - r_samp_bits, exp_start - 1);
            exp_start -= r_samp_bits + ra_bits;

            if (is_rl_rule) {
                sym = rule_stream.read(rhs_start, rhs_start + r_bits - 1);
                off_t rl_len = rule_stream.read(rhs_start + r_bits, rhs_start + (2 * r_bits) - 1);
                acc_exp = rule_stream.read(exp_start, exp_start + ra_bits - 1);

                acc_exp /= rl_len;

                rhs_data.bp_rhs_s = rhs_start;
                rhs_data.bp_exp_s = rhs_start + (exp_pos / acc_exp) * r_bits;
                rhs_data.bp_exp_e = rhs_data.bp_exp_s;
                rhs_data.bp_rhs_e = rhs_start + (rl_len - 1) * r_bits;
                rhs_data.is_rl = true;
                exp_pos %= acc_exp;
                return;
            }

            last_rhs_idx = (exp_start - rhs_start) / r_bits;
            n_samples = last_rhs_idx / rl_samp_rate;
            last_rhs_idx--;
            bps = ra_bits / (n_samples + 1);// the rules also include the expansion size as a sample
        }

        off_t left = 0, exp_m = 0, exp_n, bit_pos, middle = -1;
        auto right = (off_t) n_samples - 1;

        while (left <= right) {
            middle = left + (right - left) / 2;
            bit_pos = exp_start + middle * bps;
            exp_m = rule_stream.read(bit_pos, bit_pos + bps - 1);
            if (exp_m > exp_pos) {
                right = middle - 1;
                middle = -1;
                exp_m = 0;
                continue;
            }

            //exp_m<=pos
            bit_pos += bps;
            exp_n = rule_stream.read(bit_pos, bit_pos + bps - 1);
            if ((middle + 1) == n_samples || exp_pos < exp_n) {
                //exp_m<=pos && pos<exp_n
                break;
            }
            //exp_n=<pos
            left = middle + 1;
        }

        acc_exp = exp_m;
        if constexpr (is_str) {
            rhs_idx = (middle + 1) * str_samp_rate;
        } else {
            rhs_idx = (middle + 1) * rl_samp_rate;
        }

        off_t rhs_bit_pos = rhs_start + rhs_idx * r_bits;
        sym = rule_stream.read(rhs_bit_pos, rhs_bit_pos + r_bits - 1);
        off_t e_len = exp_len(sym);
        acc_exp += e_len;

        while (acc_exp <= exp_pos && rhs_idx < last_rhs_idx) {
            rhs_idx++;
            rhs_bit_pos += r_bits;
            sym = rule_stream.read(rhs_bit_pos, rhs_bit_pos + r_bits - 1);
            e_len = exp_len(sym);
            acc_exp += e_len;
        }
        exp_pos -= acc_exp - e_len;

        rhs_data.bp_rhs_s = rhs_start;
        rhs_data.bp_exp_s = rhs_bit_pos;
        rhs_data.bp_exp_e = rhs_bit_pos;
        rhs_data.bp_rhs_e = exp_start - r_bits;
        rhs_data.is_rl = false;
    }

    //returns |exp(sym)|
    [[nodiscard]] inline size_t exp_len(size_t sym) const {
        assert(has_rand_access);

        if (sym <= max_tsym) return 1;

        sym -= max_tsym + 1;
        size_t rl_start = rl_ptr.read(sym);
        size_t end = rl_ptr.read(sym + 1);
        size_t ra_bits = rule_stream.read(end - r_samp_bits, end - 1);
        end -= r_samp_bits;
        size_t rl_end = end - ra_bits;

        size_t rule_len = (rl_end - rl_start) / r_bits;
        size_t n_samples = (rule_len / rl_samp_rate) + 1;
        size_t bps = ra_bits / n_samples;

        return rule_stream.read(end - bps, end - 1);
    }

    uint64_t fps2fps(std::vector<uint64_t> &buffer, size_t pos, size_t e_pos) const {

        //corner case: the buffer has only one symbol
        if((e_pos-pos)==1){
            buffer[0] = XXH3_64bits(buffer.data(), sizeof(uint64_t));
            return 1;
        }

        uint64_t lb, rb, prev_sym, curr_sym, next_sym, ps_pos = pos, first_pos=pos, phrase_len;

        lb = pos;
        prev_sym = buffer[pos++];
        curr_sym = buffer[pos++];
        while (curr_sym == prev_sym && pos < e_pos) curr_sym = buffer[pos++];
        rb = pos - 1;

        //corner case: the buffer is a sequence of one symbol
        if (pos == e_pos) {
            phrase_len = rb - lb + 1;
            buffer[ps_pos++] = XXH3_64bits(buffer.data() + lb, phrase_len * sizeof(uint64_t));
            return ps_pos-first_pos;
        }

        while (pos < e_pos) {
            next_sym = buffer[pos++];
            while (next_sym == curr_sym && pos < e_pos) next_sym = buffer[pos++];

            if (prev_sym > curr_sym && curr_sym < next_sym) {
                phrase_len = rb - lb;
                buffer[ps_pos++] = XXH3_64bits(buffer.data() + lb, phrase_len * sizeof(uint64_t));
                lb = rb;
            }

            rb = pos - 1;
            prev_sym = curr_sym;
            curr_sym = next_sym;
        }

        phrase_len = pos - lb;
        buffer[ps_pos++] = XXH3_64bits(buffer.data() + lb, phrase_len * sizeof(uint64_t));

        return ps_pos-first_pos;
    }

    void get_rhs(size_t sym, bool dc_rl, std::vector<size_t>& rhs){
        rhs.clear();
        auto range = nt2bitrange(sym);
        size_t tmp_sym, len;
        for(off_t u=range.first;u<=range.second;u+=r_bits){
            tmp_sym = bitpos2symbol(u);
            if(dc_rl && is_rl_sym(tmp_sym)){
                auto range2 = nt2bitrange(tmp_sym);
                tmp_sym = bitpos2symbol(range2.first);
                len = bitpos2symbol(range2.second);
                for(size_t j=0;j<len;j++){
                    rhs.push_back(tmp_sym);
                }
            }else{
                rhs.push_back(tmp_sym);
            }
        }
    }

    void get_rhs_prefix(off_t start, off_t end, bool dc_rl, std::vector<size_t>& output){
        size_t tmp_sym, len;
        for(off_t i=start;i<=end;i+=r_bits){
            tmp_sym = bitpos2symbol(i);
            if(dc_rl && is_rl_sym(tmp_sym)){
                auto range = nt2bitrange(tmp_sym);
                tmp_sym = bitpos2symbol(range.first);
                len = bitpos2symbol(range.second);
                for(size_t j=0;j<len;j++){
                    output.push_back(tmp_sym);
                }
            }else{
                output.push_back(tmp_sym);
            }
        }
    }

    void get_rhs_suffix(off_t start, off_t end, bool dc_rl, std::vector<size_t>& output){

        if(start>end) return;
        size_t tmp_sym, len;
        for(off_t i=end;i>=start;i-=r_bits){
            tmp_sym = bitpos2symbol(i);
            if(dc_rl && is_rl_sym(tmp_sym)){
                auto range = nt2bitrange(tmp_sym);
                tmp_sym = bitpos2symbol(range.first);
                len = bitpos2symbol(range.second);
                for(size_t j=0;j<len;j++){
                    output.push_back(tmp_sym);
                }
            }else{
                output.push_back(tmp_sym);
            }
        }
    }

    uint64_t get_fp_int(size_t sym, std::vector<uint64_t> &p_seeds) const {

        assert(sym < r);
        assert(!is_rl_sym(sym));

        if (sym <= max_tsym) {
            return XXH3_64bits(&terminals[sym], sizeof(uint8_t));
        }

        uint8_t g_level = parsing_level(sym);
        auto range = nt2bitrange(sym);
        size_t rule_len = (range.second - range.first + r_bits) / r_bits;
        std::vector<uint64_t> buff_fp_seq;
        std::vector<uint8_t> buff_fp_lvl;
        buff_fp_seq.reserve(rule_len);
        buff_fp_lvl.reserve(rule_len);
        uint8_t min_g_level=255;

        for (off_t j = range.first; j <= range.second; j += r_bits) {
            sym = rule_stream.read(j, j + r_bits - 1);
            if (is_rl_sym(sym)) {
                auto range2 = nt2bitrange(sym);
                sym = rule_stream.read(range2.first, range2.first + r_bits - 1);
                rule_len = rule_stream.read(range2.second, range2.second + r_bits - 1);
                uint64_t fp = get_fp_int(sym, p_seeds);
                uint8_t p_lvl = parsing_level(sym);
                if(p_lvl<min_g_level) min_g_level = p_lvl;
                for (size_t i = 0; i < rule_len; i++) {
                    buff_fp_seq.emplace_back(fp);
                    buff_fp_lvl.emplace_back(p_lvl);
                }
            } else {
                uint8_t p_lvl = parsing_level(sym);
                if(p_lvl<min_g_level) min_g_level = p_lvl;
                buff_fp_seq.emplace_back(get_fp_int(sym, p_seeds));
                buff_fp_lvl.emplace_back(p_lvl);
            }
        }

        size_t start, n_fps, comp_pos;
        bool in_min_p_level;
        while(min_g_level<g_level) {
            comp_pos = 0, start=0;
            for (size_t i = 1; i < buff_fp_seq.size(); i++) {
                if(buff_fp_lvl[start] != buff_fp_lvl[i]) {
                    in_min_p_level =  buff_fp_lvl[start] == min_g_level;
                    n_fps = i - start;
                    if(in_min_p_level) {
                        n_fps = fps2fps(buff_fp_seq, start, i);
                    }
                    for (size_t u = start; u < start + n_fps; u++) {
                        buff_fp_seq[comp_pos] = buff_fp_seq[u];
                        buff_fp_lvl[comp_pos++] = buff_fp_lvl[u] + in_min_p_level;
                    }
                    start = i;
                }
            }

            in_min_p_level =  buff_fp_lvl[start] == min_g_level;
            n_fps = buff_fp_seq.size() - start;
            if (in_min_p_level) {
                n_fps = fps2fps(buff_fp_seq, start, buff_fp_seq.size());
            }

            for (size_t u = start; u < start + n_fps; u++) {
                buff_fp_seq[comp_pos] = buff_fp_seq[u];
                buff_fp_lvl[comp_pos++] = buff_fp_lvl[u] + in_min_p_level;
            }
            buff_fp_seq.resize(comp_pos);
            buff_fp_lvl.resize(comp_pos);
            min_g_level++;
        }
        assert(buff_fp_seq.size()==1);

        return buff_fp_seq[0];
    }

    [[nodiscard]] uint64_t get_fp(size_t sym) const {
        std::vector<uint64_t> p_seeds = get_parsing_seeds();
        return get_fp_int(sym, p_seeds);
    }

    [[nodiscard]] std::vector<uint64_t> get_parsing_seeds() const {
        std::vector<uint64_t> p_seeds;
        p_seeds.resize(lvl_rules.size()+1);
        std::random_device rd;
        std::mt19937 gen(par_seed);
        std::uniform_int_distribution<uint64_t> distrib(1, std::numeric_limits<uint64_t>::max());
        for(uint64_t & seed : p_seeds){
            seed = distrib(gen);
        }
        return p_seeds;
    }

    std::vector<uint64_t> get_all_fps(){
        std::vector<uint64_t> p_seeds = get_parsing_seeds();
        std::vector<uint64_t> fps(r, 0);

        for(size_t i=0;i<=max_tsym;i++){
            fps[i] = XXH3_64bits(&i, sizeof(uint8_t));
        }

        std::vector<uint64_t> fp_seq;
        for (off_t i = 0; i < (off_t) lvl_rules.size()-1; i++) {
            uint64_t prev=0;
            for(size_t nt=lvl_rules[i];nt<lvl_rules[i+1];nt++){
                auto range = nt2bitrange(nt);
                for(off_t j=range.first;j<=range.second;j+=r_bits){
                    size_t sym = bitpos2symbol(j);
                    assert(sym<nt);
                    fp_seq.push_back(fps[sym]);
                }
                fps[nt] = XXH3_64bits(fp_seq.data(), sizeof(uint64_t)*fp_seq.size());
                //std::cout<<nt<<" "<<max_tsym<<" "<<fps[nt]<<" "<<prev<<std::endl;
                assert(fps[nt]>=prev);
                prev=fps[nt];
                fp_seq.clear();
            }
        }
        return fps;
    }

    template<class vector_type>
    void get_nt_freqs(vector_type& nt_freq){
        size_t last_sym = n_lc_syms();
        for(size_t nt = max_tsym+1; nt<last_sym;nt++){
            auto range = nt2bitrange(nt);
            for(off_t bp=range.first;bp<=range.second;bp+=r_bits){
                size_t sym = bitpos2symbol(bp);
                nt_freq[sym]++;
            }
        }

        if constexpr (is_rl){
            last_sym = last_rl_sym();
            for(size_t nt=first_rl_sym(); nt<=last_sym;nt++){
            }
        }

        for(size_t str=0;str<n_strings();str++){
            auto range = str2bitrange(str);
            for(off_t bp=range.first;bp<=range.second;bp+=r_bits){
                size_t sym = bitpos2symbol(bp);
                nt_freq[sym]++;
            }
        }
    }
};
#endif //LCG_GRAMMAR_H
