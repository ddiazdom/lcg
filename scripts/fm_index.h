//
// Created by Diaz, Diego on 3.2.2022.
//

#ifndef TEST_RL_BCR_BWT_FM_INDEX_H
#define TEST_RL_BCR_BWT_FM_INDEX_H

#include "bwt_io.h"
#include "sdsl/wavelet_trees.hpp"

struct fm_index{

    sdsl::wt_huff<> bwt;
    sdsl::int_vector<64> C;
    sdsl::int_vector<8> sym_map;
    uint8_t dummy=0;

    explicit fm_index(std::string& input_file){

        std::string plain_bwt_file = "tmp_file";
        std::ofstream ofs(plain_bwt_file, std::ios::out | std::ios::binary);
        uint8_t buffer[1024]={0};
        bwt_buff_reader bwt_reader(input_file);
        size_t sym, freq, k=0;
        size_t sym_freqs[256]={0};
        for(size_t i=0;i<bwt_reader.size();i++){
            bwt_reader.read_run(i, sym, freq);
            for(size_t j=0;j<freq;j++){
                buffer[k++] = sym;
                if(k==1024){
                    ofs.write((char *)buffer, 1024);
                    k=0;
                }
            }
            sym_freqs[sym]+=freq;
        }
        if(k!=0){
            ofs.write((char *)buffer, (std::streamsize)k);
        }
        bwt_reader.close();
        ofs.close();
        sdsl::construct(bwt, plain_bwt_file, 1);
        remove(plain_bwt_file.c_str());

        k=0;
        std::vector<size_t> C_tmp;
        sym_map.resize(256);
        sdsl::util::set_to_value(sym_map, 0);
        for(size_t i=0;i<256;i++){
            if(sym_freqs[i]!=0){
                if(k==0) dummy = i;
                sym_map[i] = k++;
                C_tmp.push_back(sym_freqs[i]);
            }
        }

        size_t acc=0, tmp;
        for(size_t i=0;i<k;i++){
            tmp = C_tmp[i];
            C_tmp[i] = acc;
            acc+=tmp;
        }
        C_tmp[k] = acc;

        C.resize(C_tmp.size());
        for(size_t i=0;i<C.size();i++){
            C[i] = C_tmp[i];
        }
    }

    [[nodiscard]] size_t size() const {
        return bwt.size();
    }

    [[nodiscard]] size_t alphabet() const {
        return C.size();
    }

    [[nodiscard]] std::pair<uint8_t, size_t> lf(size_t idx) const {
        auto res = bwt.inverse_select(idx);
        size_t next = C[sym_map[res.second]] + res.first;
        return {res.second, next};
    }

    [[nodiscard]] inline uint8_t get_dummy() const {
        return dummy;
    }

    [[nodiscard]] inline size_t tot_strings() const {
        return C[sym_map[dummy]+1]-C[sym_map[dummy]];
    }
};
#endif //TEST_RL_BCR_BWT_FM_INDEX_H
