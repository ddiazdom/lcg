//
// Created by Diaz, Diego on 28.3.2022.
//

#ifndef CDT_BWT_DT_H
#define CDT_BWT_DT_H

#include <sdsl/wavelet_trees.hpp>
#include "utils.h"

template<class input_bwt_type>
class bwt_dt{
    typedef sdsl::wt_huff<> wt_type;
    wt_type                 m_bwt;
    std::vector<size_t>     m_C;
    std::vector<uint8_t>    m_inv_sym_map;
    std::vector<uint8_t>    m_sym_map;

    void move(bwt_dt&& other){
        m_bwt.swap(other.m_bwt);
        m_C.swap(other.m_C);
        m_inv_sym_map.swap(other.m_inv_sym_map);
        m_sym_map.swap(other.m_sym_map);
    }

    void copy(const bwt_dt& other){
        m_bwt = other.m_bwt;
        m_C = other.m_C;
        m_inv_sym_map = other.m_inv_sym_map;
        m_sym_map = other.m_sym_map;
    }

public:

    const wt_type & bwt =                     m_bwt;
    const std::vector<size_t>& C =            m_C;
    const std::vector<uint8_t>& inv_sym_map = m_inv_sym_map;
    const std::vector<uint8_t>& sym_map =     m_sym_map;
    typedef wt_type::value_type value_type;
    typedef wt_type::size_type  size_type;

    bwt_dt() = default;

    bwt_dt(std::string& plain_bwt_file, tmp_workspace& ws): m_inv_sym_map(255, 0) {

        input_bwt_type bwt_io(plain_bwt_file);
        size_t sym, run_len, k=0;

        std::string tmp_file = ws.get_file("uc_bwt");
        uint8_t buffer[8192]={0};
        size_t sym_freqs[255]={0};
        std::ofstream ofs(tmp_file, std::ios::binary);

        for(size_t i=0;i<bwt_io.size();i++){
            bwt_io.read_run(i, sym, run_len);
            for(size_t j=0;j<run_len;j++){
                buffer[k++] = (uint8_t) sym;
                if(k == 8192) {
                    ofs.write((char *)buffer, 8192);
                    k = 0;
                }
            }
            sym_freqs[sym]+=run_len;
        }
        ofs.write((char *)buffer, (std::streamsize)k);
        ofs.close();

        sdsl::construct(m_bwt, tmp_file, 1);
        ws.remove_file("uc_bwt");

        size_t max_sym, c_sym=0;
        for(size_t i=0;i<255;i++){
            if(sym_freqs[i]>0){
                m_C.push_back(sym_freqs[i]);
                m_sym_map.push_back(i);
                m_inv_sym_map[i] = c_sym++;
                max_sym=i;
            }
        }
        m_sym_map.shrink_to_fit();
        m_inv_sym_map.resize(max_sym+1);
        m_inv_sym_map.shrink_to_fit();

        size_t acc=0, tmp;
        for(unsigned long & freq : m_C){
            tmp = freq;
            freq = acc;
            acc +=tmp;
        }
        m_C.push_back(acc);
        m_C.shrink_to_fit();
    }

    bwt_dt(const bwt_dt& other) {
        copy(other);
    }

    bwt_dt(bwt_dt&& other)  noexcept {
        move(std::forward<bwt_dt>(other));
    }

    inline value_type operator[](size_t idx) const {
        return m_bwt[idx];
    }

    [[nodiscard]] inline size_t size() const {
        return m_bwt.size();
    }

    [[nodiscard]] inline size_t rank(size_t idx, value_type sym) const {
        return m_bwt.rank(idx, sym);
    }

    [[nodiscard]] inline size_t select(size_t idx, value_type sym) const {
        return m_bwt.select(idx, sym);
    }

    [[nodiscard]] inline std::pair<value_type, size_t> lf(size_t idx) const {
        auto res = m_bwt.inverse_select(idx);
        return {res.second, m_C[m_inv_sym_map[res.second]] + res.first};
    }

    [[nodiscard]] inline std::pair<value_type, size_t> inv_lf(size_t idx) const {
        value_type i=0;
        while(i<C.size() && idx>=C[i+1]) i++;
        size_t rank = idx-C[i]+1;
        return{sym_map[i], m_bwt.select(rank, sym_map[i])};
    }

    [[nodiscard]] inline std::vector<uint8_t> label(size_t idx, size_t len) const {
        assert(len>0);
        std::vector<uint8_t> label;
        auto res = inv_lf(idx);
        label.push_back(res.first);
        idx = res.second;
        while(label.size()<len && label.back()!=sym_map[0]){
            res = inv_lf(idx);
            label.push_back(res.first);
            idx = res.second;
        }
        return label;
    }

    [[nodiscard]] size_t alphabet() const {
        return m_C.size()-1;
    }

    inline void interval_symbols(size_t i, size_t j, size_type& k, std::vector<uint8_t>& symbols,
                                 std::vector<size_type>& l_ranks, std::vector<size_type>& r_ranks) const {
        m_bwt.interval_symbols(i, j, k , symbols, l_ranks, r_ranks);
    }

    inline std::pair<size_t, size_t> backward_search(std::vector<uint8_t>& pattern) const {
        std::pair<size_t, size_t> range;
        range.first=0, range.second=m_bwt.size()-1;
        value_type sym;
        for(size_t i=pattern.size();i-->0;){
            sym = m_inv_sym_map[pattern[i]];
            range.first = m_C[sym] + m_bwt.rank(range.first, pattern[i]);
            range.second = m_C[sym] + (m_bwt.rank(range.second+1, pattern[i])-1);
            if(range.second<range.first) return range;
        }
        return range;
    }


    bwt_dt& operator=(const bwt_dt& other){
        if(this!=&other){
            copy(other);
        }
        return *this;
    }

    bwt_dt& operator=(bwt_dt&& other) noexcept {
        if(this!=&other){
            move(std::forward<bwt_dt>(other));
        }
        return *this;
    }

    [[nodiscard]] inline uint8_t min_byte_symbol() const {
        return m_sym_map[0];
    }

    [[nodiscard]] inline uint8_t max_byte_symbol() const {
        return m_sym_map.back();
    }

    size_t serialize(std::ostream &out) const{
        size_t written_bytes = m_bwt.serialize(out, nullptr, "bwt_dt");
        written_bytes += serialize_plain_vector(out, m_C);
        written_bytes += serialize_plain_vector(out, m_inv_sym_map);
        written_bytes += serialize_plain_vector(out, m_sym_map);
        return written_bytes;
    }

    void load(std::istream &in){
        m_bwt.load(in);
        load_plain_vector(in, m_C);
        load_plain_vector(in, m_inv_sym_map);
        load_plain_vector(in, m_sym_map);
    }

    void swap(bwt_dt& other){
        move(std::move(other));
    }
};

#endif //CDT_BWT_DT_H
