//
// Created by Diaz, Diego on 6.5.2024.
//
#include "include/grammar.h"
#include "include/grammar_algorithms.h"

#include <fstream>

std::vector<std::string> file2strcoll(std::string& input_file){
    std::ifstream infile(input_file);
    std::vector<std::string> str_coll;
    std::string line;
    while (std::getline(infile, line)) {
        line.pop_back();
        str_coll.push_back(line);
    }
    return str_coll;
}

template<class gram_t>
void check_random_access_symbol(std::vector<std::string>& str_coll, gram_t& gram){
    for(size_t i=0;i<str_coll.size();i++){
        for(size_t j=0;j<str_coll[i].size();j++){
            uint8_t sym = gram.im_sym_rand_access(i, j);
            if(sym!=str_coll[i][j]){
                std::cout<<"Mismatch: str:"<<i<<" pos:"<<j<<" real:"<<str_coll[i][j]<<" retrieved:"<<sym<<std::endl;
            }
            assert(str_coll[i][j]==sym);
        }
    }
    std::cout<<"All symbols were retrieved from the grammar correctly"<<std::endl;
}

template<class gram_t>
void check_random_access_string(std::vector<std::string>& str_coll, gram_t& gram){
    std::string dc_string;
    for(size_t i=0;i<str_coll.size();i++){
        for(size_t j=1;j<str_coll[i].size();j++){
            for(size_t u=0;u<j;u++){
                gram.im_str_rand_access(i, u, j, dc_string);
                for(size_t a=u, b=0;a<=j;a++,b++){
                    assert(str_coll[i][a]==dc_string[b]);
                }
                dc_string.clear();
            }
        }
    }
    std::cout<<"All possible substrings were retrieved from the grammar correctly"<<std::endl;
}

void check_delete_seq(std::vector<std::string>& dc_coll){
}

void check_add_seq(std::vector<std::string>& dc_coll){
}

void check_modify_seq(std::vector<std::string>& dc_coll){
}

void check_get_fp(std::vector<std::string>& dc_coll){
}

int main(int argc, char** argv) {

    std::string str_coll_file = "/Users/ddiaz/CLionProjects/lcg-2023/data/sample_file.txt";
    std::string gram_file = "/Users/ddiaz/CLionProjects/lcg-2023/cmake-build-debug/sample_file.lcg";

    std::vector<std::string> str_coll = file2strcoll(str_coll_file);
    bool has_rl_rules, has_cg_rules, has_rand_access;
    std::tie(has_rl_rules, has_cg_rules, has_rand_access) = read_grammar_flags(gram_file);

    if(has_rand_access){
        if(has_cg_rules){
            if(has_rl_rules){
                lc_gram_t<true, true, true> gram;
                load_from_file(gram_file, gram);
                check_random_access_symbol(str_coll, gram);
                check_random_access_string(str_coll, gram);
            }else{
                lc_gram_t<true, false, true> gram;
                load_from_file(gram_file, gram);
                check_random_access_symbol(str_coll, gram);
                check_random_access_string(str_coll, gram);
            }
        }else{
            if(has_rl_rules){
                lc_gram_t<false, true, true> gram;
                load_from_file(gram_file, gram);
                check_random_access_symbol(str_coll, gram);
                check_random_access_string(str_coll, gram);
            }else{
                lc_gram_t<false, false, true> gram;
                load_from_file(gram_file, gram);
                check_random_access_symbol(str_coll, gram);
                check_random_access_string(str_coll, gram);
            }
        }
    }else{
        if(has_cg_rules){
            if(has_rl_rules){
                lc_gram_t<true, true, false> gram;
                load_from_file(gram_file, gram);
            }else{
                lc_gram_t<true, false, false> gram;
                load_from_file(gram_file, gram);
            }
        }else{
            if(has_rl_rules){
                lc_gram_t<false, true, false> gram;
                load_from_file(gram_file, gram);
            }else{
                lc_gram_t<false, false, false> gram;
                load_from_file(gram_file, gram);
            }
        }
    }

    return 0;
}
