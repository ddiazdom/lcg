//
// Created by Diaz, Diego on 16.1.2024.
//

#ifndef LCG_BUILD_COLLAGE_SYSTEM_H
#define LCG_BUILD_COLLAGE_SYSTEM_H
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include "cds/int_array.h"
#include "sorting_strings.h"

struct collage_data{
    size_t nt{};
    std::pair<long, long> pref_range;
    size_t position;
    size_t rule_size{};
    uint8_t f_sym;
    uint8_t mm_sym;

    long ref_nt =-1;
    bool side; //suffix or prefix
    size_t lcp{};
    std::string dc_string;
};

size_t get_cs_left_rules(std::vector<collage_data>& dc_rules) {

    std::sort(dc_rules.begin(), dc_rules.end(), [&](auto const& a, auto const& b) -> bool{
        size_t len = std::min(a.dc_string.size(), b.dc_string.size());
        for(size_t k=0;k<len;k++){
            if(a.dc_string[k]!=b.dc_string[k]){
                return a.dc_string[k]<b.dc_string[k];
            }
        }
        return a.dc_string.size()<b.dc_string.size();
    });

    size_t nt_compressed=0;
    dc_rules[0].lcp = 0;
    for(size_t k=1;k<dc_rules.size();k++){
        size_t u = 0;
        size_t tmp = std::min(dc_rules[k-1].dc_string.size(), dc_rules[k].dc_string.size());

        while(u<tmp &&
              dc_rules[k].dc_string[u]==dc_rules[k-1].dc_string[u]){
            u++;
        }
        dc_rules[k].lcp = u;
    }

    for(size_t k=0;k<dc_rules.size()-1;k++) {
        size_t len = dc_rules[k].dc_string.size();
        if(dc_rules[k+1].lcp<len) continue;

        size_t prev_len = len;
        size_t ptr_rl = k+1;
        while(ptr_rl<dc_rules.size() && dc_rules[ptr_rl].lcp>=prev_len){
            prev_len = dc_rules[ptr_rl].lcp;
            ptr_rl++;
        }

        if(dc_rules[k].ref_nt<0 && dc_rules[k].rule_size>3){
            dc_rules[k].ref_nt = (long)dc_rules[ptr_rl-1].nt;
            dc_rules[k].side = false;
            //auto res =std::mismatch(dc_rules[k].dc_string.begin(), dc_rules[k].dc_string.end(), dc_rules[ptr_rl-1].dc_string.begin());
            //assert(res.first==dc_rules[k].dc_string.end());
            nt_compressed++;
        }
    }
    return nt_compressed;
}

size_t get_cs_right_rules(std::vector<collage_data>& dc_rules) {

    std::sort(dc_rules.begin(), dc_rules.end(), [&](auto const& a, auto const& b) -> bool{
        size_t len = std::min(a.dc_string.size(), b.dc_string.size());
        size_t pos_a = a.dc_string.size()-1;
        size_t pos_b = b.dc_string.size()-1;

        for(size_t k=0;k<len;k++){
            if(a.dc_string[pos_a-k]!=b.dc_string[pos_b-k]){
                return a.dc_string[pos_a-k]<b.dc_string[pos_b-k];
            }
        }
        return a.dc_string.size()<b.dc_string.size();
    });

    size_t prev_len = dc_rules[0].dc_string.size()-1;
    dc_rules[0].lcp = 0;

    for(size_t k=1;k<dc_rules.size();k++){
        size_t u = 0;
        size_t len = dc_rules[k].dc_string.size()-1;

        size_t tmp = std::min(len, prev_len);
        while(u<=tmp &&
              dc_rules[k].dc_string[len-u]==dc_rules[k-1].dc_string[prev_len-u]){
            u++;
        }
        dc_rules[k].lcp = u;
        prev_len = len;
    }

    size_t nt_compressed=0;
    for(size_t k=0;k<dc_rules.size()-1;k++) {
        size_t len = dc_rules[k].dc_string.size();
        if(dc_rules[k+1].lcp<len) continue;

        prev_len = len;
        size_t ptr_rl = k+1;
        while(ptr_rl<dc_rules.size() && dc_rules[ptr_rl].lcp>=prev_len){
            prev_len = dc_rules[ptr_rl].lcp;
            ptr_rl++;
        }

        if(dc_rules[k].ref_nt<0 && dc_rules[k].rule_size>3){
            dc_rules[k].ref_nt = (long)dc_rules[ptr_rl-1].nt;
            dc_rules[k].side = true;
            //std::string& str = dc_rules[dc_rules[k].ref_nt].dc_string;
            //std::string& suffix = dc_rules[k].dc_string;
            //assert(str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0);
            nt_compressed++;
        }
    }
    return nt_compressed;
}

template<class gram_t>
std::pair<size_t, size_t> compute_cg_rules(gram_t& gram, std::unordered_map<size_t, std::tuple<size_t, size_t, bool>>& ht){

    size_t max_len = 0, new_g_size=0;
    for(size_t i=0;i<gram.lvl_rules.size()-1; i++){
        size_t lvl_rules = gram.lvl_rules[i+1]-gram.lvl_rules[i];//number of rules in the level

        std::cout<<"Finding colleague rules in level "<<i+1<<std::endl;
        if(lvl_rules>0) {

            std::vector<collage_data> dc_rules;
            dc_rules.resize(lvl_rules);

            size_t first_nt = gram.lvl_rules[i];
            size_t last_nt = gram.lvl_rules[i + 1] - 1;
            std::cout<<"  Decompressing the strings"<<std::endl;
            for(size_t nt=first_nt, j=0;nt<=last_nt;nt++,j++){
                dc_rules[j].nt = nt;
                auto res = gram.nt2phrase(nt);
                size_t rule_size = res.second-res.first+1;
                gram.in_memory_decompression(nt, dc_rules[j].dc_string);
                dc_rules[j].rule_size = rule_size;
            }
            std::cout<<"  Computing suffix and prefix rules"<<std::endl;
            get_cs_left_rules(dc_rules);
            get_cs_right_rules(dc_rules);

            size_t comp_rules=0;
            for(auto const& rule : dc_rules){
                if(rule.ref_nt>=0){
                    if(rule.dc_string.size()>max_len) max_len = rule.dc_string.size();
                    ht.insert({rule.nt, {rule.ref_nt, rule.dc_string.size(), rule.side}});
                    new_g_size+=3;
                    comp_rules++;
                }else{
                    new_g_size+=rule.rule_size;
                }
            }
            std::cout<<"  There are "<<comp_rules<<" compressible rules in this level"<<std::endl;
            if(comp_rules==0) break;
        }
    }
    return {new_g_size, max_len};
}

template<class gram_t>
size_t get_cs_left_rules_new(std::vector<collage_data>& level_rules,  const gram_t& gram,
                             std::vector<uint32_t>& g_lcp, std::vector<uint32_t>& sym_perm,
                             size_t& lvl_offset, size_t first_lvl_nt, size_t lvl, std::unordered_set<size_t>& suff_nts) {

    std::sort(level_rules.begin(), level_rules.end(), [&](auto const& a, auto const& b){
        size_t len = std::min(a.rule_size, b.rule_size);
        for(size_t i=0;i<len;i++){
            if(gram.rules[a.position+i]!=gram.rules[b.position+i]){
                size_t sym_a = sym_perm[gram.rules[a.position+i]-lvl_offset];
                size_t sym_b = sym_perm[gram.rules[b.position+i]-lvl_offset];
                //bool is_suffix = (sym_a+1)<=sym_b && sym_b<=suff_arr[sym_a] || (sym_b+1) <= suff_arr[sym_b];
                //return gram.substring_lex_comp(a.nt, b.nt, exp);

                return sym_a < sym_b;
            }
        }
        return a.rule_size<b.rule_size;
    });

    std::vector<uint32_t> new_g_lcp(level_rules.size(), 0);
    std::vector<uint32_t> new_sym_perm(level_rules.size(), 0);

    new_sym_perm[level_rules[0].nt-first_lvl_nt] = 0;

    for(size_t k=1;k<level_rules.size();k++){
        size_t u = 0;
        size_t min_size = std::min(level_rules[k-1].rule_size, level_rules[k].rule_size);
        size_t pos_a = level_rules[k-1].position;
        size_t pos_b = level_rules[k].position;

        while(u<min_size &&
              gram.rules[pos_b+u]==gram.rules[pos_a+u]){
            new_g_lcp[k] += gram.rule_exp[gram.rules[pos_b+u]];
            u++;
        }

        if(min_size==u){
            suff_nts.insert(level_rules[k-1].nt);
        }

        //if(u==min_size){
        //    std::string tmp;
        //    gram.in_memory_decompression(level_rules[k-1].nt, tmp);
        //    std::cout<<"This is prefix: "<<k-1<<"\t"<<tmp<<std::endl;
        //}

        //if(lvl>0 && u<level_rules[k-1].rule_size && u<level_rules[k].rule_size){
            //new_g_lcp[k]+=gram.longest_com_pref(level_rules[k].nt, level_rules[k-1].nt, new_g_lcp[k]);
            //size_t sym_a = sym_perm[gram.rules[pos_a+u]-lvl_offset];
            //size_t sym_b = sym_perm[gram.rules[pos_b+u]-lvl_offset];
            /*if(sym_a>sym_b) std::swap(sym_a, sym_b);
            //TODO a dummy way to simulate a rmq
            size_t pref_match = std::numeric_limits<size_t>::max();
            for(size_t j=sym_a+1;j<=sym_b;j++){
                if(g_lcp[j]<pref_match){
                    pref_match = g_lcp[j];
                }
            }
            //
            new_g_lcp[k] +=pref_match;
            if(u<(level_rules[k].rule_size-1) && pref_match==gram.rule_exp[gram.rules[pos_a+u]]){
                //gram.print_parse_tree(level_rules[k-1].nt);
                //gram.print_parse_tree(level_rules[k].nt);
                //std::cout<<"accessing "<<new_g_lcp[k]<<" from "<<level_rules[k].nt<<" "<<gram.rule_exp[level_rules[k].nt]<<std::endl;
                //std::cout<<"We got: "<<gram.access_pos(level_rules[k].nt, new_g_lcp[k])<<" in level "<<lvl<<std::endl;
                //if(lvl>=4){
                //    std::cout<<gram.substring_lex_comp(level_rules[k].nt, level_rules[k-1].nt, new_g_lcp[k])<<std::endl;
                //}
                //std::cout<<"Ups : "<<pref_match<<"\n"<<std::endl;
                //TODO do something else
            }*/
        //}
        new_sym_perm[level_rules[k].nt-first_lvl_nt] = k;
    }

    /*size_t k=0;
    for(auto & rule : level_rules){
        std::cout<<new_g_lcp[k]<<"\t";
        //size_t end = rule.position+rule.rule_size-1;
        for(size_t u=0;u<rule.rule_size;u++){
            std::cout<<gram.rules[rule.position+u]<<" ";
        }
        std::cout<<" "<<std::endl;
        k++;
        std::string tmp;
        gram.in_memory_decompression(rule.nt, tmp);
        std::cout<<tmp<<std::endl;
    }*/

    std::swap(new_g_lcp, g_lcp);
    std::swap(new_sym_perm, sym_perm);

    size_t nt_compressed=0;
    /*for(size_t k=0;k<level_rules.size()-1;k++) {
        size_t len = level_rules[k].dc_string.size();
        if(level_rules[k+1].lcp<len) continue;

        size_t prev_len = len;
        size_t ptr_rl = k+1;
        while(ptr_rl<level_rules.size() && level_rules[ptr_rl].lcp>=prev_len){
            prev_len = level_rules[ptr_rl].lcp;
            ptr_rl++;
        }

        if(level_rules[k].ref_nt<0 && level_rules[k].rule_size>3){
            level_rules[k].ref_nt = (long)level_rules[ptr_rl-1].nt;
            level_rules[k].side = false;
            nt_compressed++;
        }
    }*/
    return nt_compressed;
}

template<class gram_t>
std::pair<size_t, size_t> compute_cg_rules_new(gram_t& gram, std::unordered_map<size_t, std::tuple<size_t, size_t, bool>>& ht){

    size_t max_len = 0, new_g_size=0, lvl_offset=0;
    std::vector<uint32_t> g_lcp(gram.max_tsym+1, 0);
    std::vector<uint32_t> sym_perm(gram.max_tsym+1, 0);

    for(size_t i=0;i<=gram.max_tsym;i++){
        sym_perm[i] = i;
    }

    std::unordered_set<size_t> suff_rules;

    for(size_t i=0;i<gram.lvl_rules.size()-1; i++){

        size_t n_lvl_rules = gram.lvl_rules[i+1]-gram.lvl_rules[i];//number of rules in the level

        std::cout<<"Finding colleague rules in level "<<i+1<<std::endl;
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

            get_cs_left_rules_new(level_rules, gram, g_lcp, sym_perm, lvl_offset, first_nt, i, suff_rules);
            lvl_offset = first_nt;

            /*std::sort(level_rules.begin(), level_rules.end(), [&](auto const& a, auto const& b){
                size_t len = std::min(a.rule_size, b.rule_size);
                size_t end_a = a.position+a.rule_size-1;
                size_t end_b = b.position+b.rule_size-1;
                for(size_t i=0;i<len;i++){
                    if(gram.rules[end_a-i]!=gram.rules[end_b-i]){
                        return gram.rules[end_a-i]<gram.rules[end_b-i];
                    }
                }
                return a.rule_size<b.rule_size;
            });*/


            /*std::cout<<"  Computing suffix and prefix rules"<<std::endl;
            get_cs_left_rules(level_rules);
            get_cs_right_rules(level_rules);

            size_t comp_rules=0;
            for(auto const& rule : level_rules){
                if(rule.ref_nt>=0){
                    if(rule.dc_string.size()>max_len) max_len = rule.dc_string.size();
                    ht.insert({rule.nt, {rule.ref_nt, rule.dc_string.size(), rule.side}});
                    new_g_size+=3;
                    comp_rules++;
                }else{
                    new_g_size+=rule.rule_size;
                }
            }
            std::cout<<"  There are "<<comp_rules<<" compressible rules in this level"<<std::endl;
            if(comp_rules==0) break;*/
        }
    }
    std::cout<<"We have "<<suff_rules.size()<<" rules "<<std::endl;
    return {new_g_size, max_len};
}

template<class gram_t>
void make_collage_system(gram_t& gram){

    std::unordered_map<size_t, std::tuple<size_t, size_t, bool>> ht;
    auto cs_res = compute_cg_rules_new(gram, ht);

    /*int_array<size_t> new_rules(cs_res.first, sym_width(std::max(gram.r, cs_res.second)));
    int_array<size_t> new_rl_ptrs(sym_width(gram.r), sym_width(cs_res.first));

    //insert terminals
    for(size_t i=0;i<gram.n_terminals();i++){
        new_rules.push_back(i);
    }

    //insert regular rules
    size_t start_sym = gram.start_symbol();
    for(size_t sym=gram.max_tsym+1;sym<start_sym;sym++) {
        new_rl_ptrs.push_back(new_rules.size());

        auto res = ht.find(sym);
        if(res!=ht.end()){

            //TODO test
            /*auto r1 = gram.nt2phrase(sym);
            auto r2 = gram.nt2phrase(std::get<0>(res->second));

            size_t matches=0;
            size_t len = r1.second-r1.first+1;
            if(!std::get<2>(res->second)){
                for(size_t k=0;k<len;k++){
                    matches += gram.rules[r1.first+k] == gram.rules[r2.first+k];
                }
            }else{
                for(size_t k=0;k<len;k++){
                    matches += gram.rules[r1.second-k] == gram.rules[r2.second-k];
                }
            }

            if(matches<(len-1)){
                std::cout<<" comp rule :"<<sym<<" -> ";
                for(size_t j=r1.first;j<=r1.second;j++){
                    std::cout<<gram.rules[j]<<" ";
                }
                std::cout<<" ref rule :"<<std::get<0>(res->second)<<" -> ";
                for(size_t j=r2.first;j<=r2.second;j++){
                    std::cout<<gram.rules[j]<<" ";
                }
                std::cout<<" "<<len<<" "<<matches<<std::endl;
                std::string a;
                std::string b;
                gram.in_memory_decompression(sym, a);
                gram.in_memory_decompression(std::get<0>(res->second), b);
                for(size_t i=0;i<a.size();i++){
                    std::cout<<a[i]<<" ";
                }
                std::cout<<""<<std::endl;
                gram.print_parse_tree(sym);
                for(size_t i=0;i<b.size();i++){
                    std::cout<<b[i]<<" ";
                }
                gram.print_parse_tree(std::get<0>(res->second));
                std::cout<<"\n\n"<<std::endl;

            } * /
            //
            new_rules.push_back(std::get<2>(res->second));
            new_rules.push_back(std::get<0>(res->second));
            new_rules.push_back(std::get<1>(res->second));
        }else{
            auto range = gram.nt2phrase(sym);
            for(size_t j=range.first;j<=range.second;j++){
                new_rules.push_back(gram.rules[j]);
            }
        }
    }

    //insert the compressed strings
    new_rl_ptrs.push_back(new_rules.size());
    for(size_t str=0;str<gram.n_strings();str++){
        auto range = gram.str2phrase(str);
        gram.str_boundaries[str] = new_rules.size();
        for(size_t j=range.first;j<=range.second;j++){
            new_rules.push_back(gram.rules[j]);
        }
    }
    new_rl_ptrs.push_back(new_rules.size());
    gram.str_boundaries[gram.n_strings()]=new_rules.size();

    std::cout<<"  Stats:"<<std::endl;
    std::cout<<"    Grammar size before: "<<gram.rules.size()<<std::endl;
    std::cout<<"    Grammar size after:  "<<new_rules.size()<<std::endl;
    std::cout<<"    Number of cs rules:  "<<ht.size()<<" ("<<float(ht.size())/float(gram.r)*100<<"%)"<<std::endl;
    std::cout<<"    Compression ratio:   "<<float(new_rules.size())/float(gram.rules.size())<<std::endl;

    gram.g  = new_rules.size();
    new_rules.swap(gram.rules);
    new_rl_ptrs.swap(gram.rl_ptr);

    //TODO test
    start_sym = gram.start_symbol();
    std::vector<uint8_t> freq(gram.r+1, 0);
    for(size_t sym=gram.max_tsym+1;sym<=start_sym;sym++) {
        auto range = gram.nt2phrase(sym);
        if(gram.rules[range.first]==0){
            freq[gram.rules[range.first+1]++];
        }else{
            for(size_t j=range.first;j<=range.second;j++){
                freq[gram.rules[j]]++;
            }
        }
    }

    size_t n_rems = 0;
    for(size_t i=gram.max_tsym+1;i<gram.r;i++){
        n_rems+=freq[i]==0;
    }
    std::cout<<"We can remove "<<n_rems<<std::endl;*/
}
#endif //LCG_BUILD_COLLAGE_SYSTEM_H
