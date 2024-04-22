//
// Created by Diaz, Diego on 19.4.2024.
//

#ifndef SE_STRAT_GRAM_BUFFER_H
#define SE_STRAT_GRAM_BUFFER_H

#include "cds/file_streams.hpp"
#include "../old/hashing.h"

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
    std::vector<uint64_t>              par_functions;

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
                //std::cout<<parsing_set[str.str_ptr+j]<<" ";
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
            str_boundaries[i] = (off_t)i;
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

    template<class gram_type>
    void make_gram(gram_type& new_gram, uint64_t par_seed_) {

        rules_buffer.close();
        i_file_stream<size_t> r_buffer(rules_file, BUFFER_SIZE);

        new_gram.n = n;
        new_gram.r = r;
        new_gram.g = g;
        new_gram.c = c;
        new_gram.s = str_boundaries.size()-1;
        //longest_str = gram_buff.long_str;

        new_gram.max_tsym = max_tsym;
        new_gram.sep_tsym = sep_tsym;

        new_gram.terminals = terminals;
        new_gram.par_seed = par_seed_;
        new_gram.lvl_rules = lvl_rules;

        new_gram.run_len_nt = rl_rules;

        uint8_t r_bits = sym_width(new_gram.r);
        new_gram.r_bits = r_bits;
        //new_gram.rules.set_width(sym_width(new_gram.r));
        //new_gram.rules.resize(new_gram.g);
        new_gram.rule_stream.reserve_in_bits(new_gram.r_bits*new_gram.g);


        new_gram.rl_ptr.set_width(sym_width(new_gram.g*new_gram.r_bits));
        new_gram.rl_ptr.resize(new_gram.r-new_gram.max_tsym);

        size_t bit_pos = 0;
        for(size_t ter=0;ter<=new_gram.max_tsym;ter++){
            //new_gram.rules.write(i, r_buffer.read(i)>>1UL);
            new_gram.rule_stream.write(bit_pos, bit_pos+r_bits, r_buffer.read(ter)>>1UL);
            bit_pos+=r_bits;
        }

        size_t rule=0;
        for(size_t i=new_gram.max_tsym+1;i<r_buffer.size();i++){
            size_t sym = r_buffer.read(i);
            bool first = sym & 1UL;

            //new_gram.rules.write(i, sym>>1UL);
            new_gram.rule_stream.write(bit_pos, bit_pos+r_bits-1, sym>>1UL);
            if(first){
                new_gram.rl_ptr.write(rule, bit_pos);
                rule++;
            }
            bit_pos+=r_bits;
        }
        assert((bit_pos/new_gram.r_bits)==new_gram.g);
        assert(rule==(new_gram.r-(new_gram.max_tsym+1)));

        new_gram.rl_ptr.write(rule, bit_pos);

        size_t offset = new_gram.r_bits*(new_gram.g-new_gram.c);
        new_gram.str_boundaries.resize(str_boundaries.size());
        size_t str=0;
        for(off_t & str_boundary : str_boundaries){
            new_gram.str_boundaries[str++] = (str_boundary*r_bits) + offset;
        }
        assert(new_gram.str_boundaries[0]==offset);
        assert(new_gram.str_boundaries.back()==bit_pos);


        size_t acc=new_gram.max_tsym+1, tmp;
        for(unsigned long & lvl_rule : new_gram.lvl_rules){
            tmp = lvl_rule;
            lvl_rule = acc;
            acc+=tmp;
        }
        new_gram.lvl_rules.push_back(acc);
        r_buffer.close();
    }
};


#endif //SE_STRAT_GRAM_BUFFER_H
