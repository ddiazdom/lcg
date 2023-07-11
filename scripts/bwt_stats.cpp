//
// Created by Diaz, Diego on 5.2.2022.
//
#include <iostream>
#include <vector>
#include "bwt_io.h"
#include <cmath>

int main(int argc, char** argv){

    if(argc!=2){
        std::cout<<"usage: ./stats file.rlbwt\n"
                   "file.rlbwt is the bcr bwt file\n"<<std::endl;
        exit(0);
    }

    errno = 0;
    std::string input_file = std::string(argv[1]);
    size_t min_run=std::numeric_limits<size_t>::max();
    size_t max_run=0;

    bwt_buff_reader bwt_reader(input_file);
    size_t sym=0, len=0, acc=0;
    std::vector<uint32_t> bwt_lens(bwt_reader.size(), 0);
    size_t one_byte=0;
    size_t two_bytes=0;
    size_t three_bytes=0;
    size_t alphabet[256]={0};
    size_t freq[256]={0};
    uint8_t min_sym = 255, max_sym=0;

    std::cout<<"Reading the input BWT"<<std::endl;
    for(size_t i=0;i<bwt_reader.size();i++){
        bwt_reader.read_run(i, sym, len);
        if(len<min_run) min_run = len;
        if(len>max_run) max_run = len;
        alphabet[sym]++;
        freq[sym]+=len;

        acc+=len;
        bwt_lens[i] = len;
        if(len<=255){
            one_byte++;
        }else if(len<=65535){
            two_bytes++;
        }else{
            three_bytes++;
        }
        if(sym<min_sym) min_sym = sym;
        if(sym>max_sym) max_sym = sym;
    }

    size_t sigma=0;
    for(unsigned char i : alphabet){
        if(i!=0){
            sigma++;
        }
    }

    std::sort(bwt_lens.begin(), bwt_lens.end());

    auto l = (double)bwt_reader.size();
    double prop=0.1;

    std::vector<std::tuple<uint8_t, size_t, size_t>> rl_stats;
    for(size_t i=0;i<256;i++){
        if(alphabet[i]!=0){
            rl_stats.emplace_back(i, alphabet[i], freq[i]);
        }
    }
    std::sort(rl_stats.begin(), rl_stats.end(), [&](auto a, auto b){ return std::get<1>(a)<std::get<1>(b);});

    std::cout<<"Number of runs: "<<bwt_reader.size()<<std::endl;
    std::cout<<"Text alphabet: "<<sigma<<std::endl;
    std::cout<<"Run stats:"<<std::endl;

    size_t k=1;
    //size_t n_bits = sizeof(unsigned int)*8 - __builtin_clz(acc);
    size_t space_acc=0, space;
    for(auto const& tuple : rl_stats){
        space=32*std::get<1>(tuple)*2;
        std::cout<<k++<<") Symbol:"<<(int)std::get<0>(tuple)<<"\t\tnumber of runs in the BWT:"<<std::get<1>(tuple)<<"\t\ttext frequency:"<<std::get<2>(tuple)<<" | "<<space<<" "<<space_acc<<std::endl;
        space_acc+=space;
    }
    std::cout<<"Text size: "<<acc<<std::endl;
    std::cout<<"n/r: "<<double(acc)/double(bwt_reader.size())<<std::endl;
    std::cout<<(double(one_byte)/l)*100<<"% of the run lenghts fit 1 byte"<<std::endl;
    std::cout<<(double(two_bytes)/l)*100<<"% of the run lenghts fit 2 bytes"<<std::endl;
    std::cout<<(double(three_bytes)/l)*100<<"% of the run lenghts fit 3 or more bytes"<<std::endl;
    std::cout<<"Deciles: "<<std::endl;
    for(size_t i=0;i<9;i++){
        size_t q = ceil(l*prop);
        std::cout<<"  q"<<(i+1)<<": "<<bwt_lens[q]<<std::endl;
        prop+=0.1;
    }
    std::cout<<"Longest run: "<<max_run<<std::endl;
    std::cout<<"Shortest run: "<<min_run<<std::endl;

}

