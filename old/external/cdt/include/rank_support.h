//
// Created by Diaz, Diego on 22.11.2022.
//

#ifndef CDT_RANK_SUPPORT_H
#define CDT_RANK_SUPPORT_H
#include <vector>
#include <iostream>
#include "int_array.h"

constexpr size_t W = std::numeric_limits<size_t>::digits;

//TODO k>W is supported just if k is a power of two (we can support Rank for any K but I have to modify rank operation)
template<size_t k=W, size_t s=k*4>
class rank_support {
private:
    static_assert(k >= W && k%W==0 && k<s && s%k==0);
    constexpr static uint8_t shift = __builtin_ctz(W);
    const bit_array *m_bv = nullptr;
    size_t words_per_block = k/W;
    std::vector<size_t> R;
    int_array<size_t> r;

private:
    void populate_rank_samples(){
        size_t  R_acc=0, r_acc=0, pop_count, bits=0;

        R.reserve(m_bv->size()/s + 1);
        r.set_width(sym_width(s-k));
        r.reserve(m_bv->size()/k + 1);

        R.push_back(0);
        r.push_back(0);

        for(size_t i=0;i<m_bv->bits.stream_size;i++){

            pop_count = __builtin_popcountll(m_bv->bits.stream[i]);
            R_acc+=pop_count;
            r_acc+=pop_count;
            bits+= bit_array::stream_t::word_bits;

            if(bits%s==0){
                R.push_back(R_acc);
                r_acc = 0;
            }

            if(bits%k==0){
                r.push_back(r_acc);
            }
        }
    };

    void copy(const rank_support &other){
        m_bv = other.m_bv;
        R = other.R;
        r = other.r;
    };

public:
    explicit rank_support(const bit_array* bv): m_bv(bv),
                                                r(0, sym_width(s-k)){
        populate_rank_samples();
    };

    rank_support()= default;

    void set_bit_vector(const bit_array *bv){
        m_bv = bv;
        populate_rank_samples();
    };

    size_t operator()(size_t pos) const{
        assert(pos<=m_bv->size());
        if(pos==0) return 0;

        pos--;
        size_t index = pos/k;
        size_t n_words = (pos%k)>>shift;
        size_t start = index*words_per_block;
        size_t end = start + n_words;

        size_t acc=0;
        for(size_t i=start;i<end;i++){
            acc += __builtin_popcountll(m_bv->bits.stream[i]);
        }
        return R[pos/s] + r[index] + acc + m_bv->bits.pop_count((index*k) + (n_words<<shift), pos);
    };

    rank_support& operator=(const rank_support& other){
        if(this!=&other){
            copy(other);
        }
        return *this;
    };
};
#endif //CDT_RANK_SUPPORT_H
