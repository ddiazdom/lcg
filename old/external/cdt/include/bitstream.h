//
// Created by diego on 02-09-20.
//

#ifndef LPG_COMPRESSOR_BITSTREAM_H
#define LPG_COMPRESSOR_BITSTREAM_H

#include<iostream>
#include <limits>
#include "macros.h"
#include "cdt_common.hpp"

template<class word_t,
         uint8_t max_dist=std::numeric_limits<word_t>::digits>
struct bitstream{

    constexpr static uint8_t word_bits = std::numeric_limits<word_t>::digits;
    constexpr static uint8_t word_shift = __builtin_ctz(word_bits);
    const static size_t masks[65];
    using size_type = size_t;

    word_t *stream=nullptr;
    size_t stream_size=0;

    bitstream(): stream(nullptr), stream_size(0){};

    [[nodiscard]] inline size_t n_bits() const {
        return stream_size<<word_shift;
    }

    inline bitstream<word_t>& swap(bitstream<word_t>& other) {
        std::swap(stream, other.stream);
        std::swap(stream_size, other.stream_size);
        return *this;
    }

    inline bitstream<word_t>& operator=(bitstream<word_t> const& other){
        if(&other!=this){
            if(stream_size!=other.stream_size){
                stream = reinterpret_cast<word_t *>(realloc(stream, other.stream_size*sizeof(word_t)));
                stream_size = other.stream_size;
            }
            memcpy(stream, other.stream, stream_size*sizeof(word_t));
        }
    }

    inline void write(size_t i, size_t j, size_t value){
        size_t cell_i = i >> word_shift;
        size_t i_pos = (i & (word_bits-1UL));

        if constexpr (max_dist==1){
            stream[cell_i] &= ~(1UL << i_pos);
            stream[cell_i] |= value << i_pos;
        }else{
            size_t cell_j = j >> word_shift;
            if(cell_i==cell_j){
                stream[cell_i] &= ~(masks[(j - i + 1UL)] << i_pos);
                stream[cell_i] |= value << i_pos;
            }else{
                size_t right = word_bits - i_pos;
                size_t left = 1+(j & (word_bits - 1UL));
                stream[cell_i] = (stream[cell_i] & ~(masks[right] << i_pos)) | (value << i_pos);
                stream[cell_j] = (stream[cell_j] & ~masks[left]) | (value >> right);
            }
        }
    }

    inline void write_chunk(const void* source, size_t i, size_t j){
        size_t tot_bits = j-i+1;
        size_t n_words = INT_CEIL(tot_bits, word_bits);
        size_t left = i & (word_bits - 1UL);
        size_t right = word_bits - left;
        size_t cell_i = i >> word_shift;

        auto tmp_src = reinterpret_cast<const word_t *>(source);

        //TODO use SIMD instructions for this segment
        for(size_t k=0; k < n_words - 1; k++){
            stream[cell_i] = (stream[cell_i] & ~(masks[right] << left)) | (tmp_src[k] << left);
            cell_i++;
            stream[cell_i] = (stream[cell_i] & ~masks[left]) | (tmp_src[k] >> right);
        }
        //

        size_t read_bits = ((n_words - 1) << word_shift);

        write(i + read_bits, j, (tmp_src[n_words-1] & masks[tot_bits-read_bits]));
    }

    [[nodiscard]] inline size_t read(size_t i, size_t j) const{
        if constexpr (max_dist==1){
            return (stream[i>>word_shift] >> (i & (word_bits - 1UL))) & 1UL;
        }else{
            size_t cell_i = i >> word_shift;
            size_t i_pos = (i & (word_bits - 1UL));
            size_t cell_j = j >> word_shift;
            if(cell_i == cell_j){
                return (stream[cell_i] >> i_pos) & masks[(j - i + 1UL)];
            }else{
                size_t right = word_bits-i_pos;
                size_t left = 1+(j & (word_bits - 1UL));
                return ((stream[cell_j] & masks[left]) << right) | ((stream[cell_i] >> i_pos) & masks[right]);
            }
        }
    }

    [[nodiscard]] inline size_t pop_count(size_t i, size_t j) const{
        size_t cell_i = i >> word_shift;
        size_t i_pos = (i & (word_bits - 1UL));
        size_t cell_j = j >> word_shift;
        size_t val;
        if(cell_i == cell_j){
            val = (stream[cell_i] >> i_pos) & masks[(j - i + 1UL)];
        }else{
            size_t right = word_bits-i_pos;
            size_t left = 1+(j & (word_bits - 1UL));
            val = ((stream[cell_j] & masks[left]) << right) | ((stream[cell_i] >> i_pos) & masks[right]);
        }
        return __builtin_popcountll(val);
    }

    inline void read_chunk(void* dst, size_t i, size_t j) const{
        size_t tot_bits = j-i+1;
        size_t n_words = INT_CEIL(tot_bits, word_bits);
        size_t left = i & (word_bits - 1UL);
        size_t right = word_bits - left;
        size_t cell_i = i >> word_shift;

        auto tmp_dst = reinterpret_cast<word_t *>(dst);

        //TODO use SIMD instructions for this segment
        for(size_t k=0; k < n_words - 1; k++){
            tmp_dst[k] = (stream[cell_i] >> left) & masks[right];
            tmp_dst[k] |= (stream[++cell_i] & masks[left]) << right;
        }
        //

        size_t read_bits = ((n_words - 1) << word_shift);
        tmp_dst[n_words-1] &= ~masks[tot_bits-read_bits];
        tmp_dst[n_words-1] |= read(i + read_bits, j);
    }

    //compare a segment of the stream with an external source of bits
    inline bool compare_chunk(const void* input, size_t i, size_t bits) const {

        size_t n_words = INT_CEIL(bits, word_bits);
        size_t left = i & (word_bits - 1UL);
        size_t right = word_bits - left;
        size_t cell_i = i >> word_shift;

        auto tmp_in = reinterpret_cast<const word_t *>(input);

        //TODO use SIMD instructions for this segment
        size_t tmp_data;
        for(size_t k=0; k < n_words - 1; k++){
            tmp_data = (stream[cell_i] >> left) & masks[right];
            tmp_data |= (stream[++cell_i] & masks[left]) << right;
            if(tmp_data != tmp_in[k]) return false;
        }
        //
        size_t read_bits = ((n_words - 1) << word_shift);
        return (tmp_in[n_words - 1] & masks[(bits-read_bits)]) == read(i + read_bits, i+bits-1);
    }

    //compare a segment of the stream with an external source of bits
    /*inline bool compare_segment(const uint8_t* input, size_t i, size_t bits) const {

        size_t n_words = INT_CEIL(bits, word_bits);
        size_t left = i & (word_bits - 1UL);
        size_t right = word_bits - left;

        for(size_t k=0, cell_i=(i>>word_shift); k < n_words - 1; k+=8, cell_i++){
            const size_t tmp_data = ((stream[cell_i] >> left) & masks[right]) | ((stream[cell_i+1] & masks[left]) << right);
            if(memcmp(input+k, &tmp_data, sizeof(size_t))!=0) return false;
        }

        //size_t cell_i = i >> word_shift;
        size_t read_bits = ((n_words - 1) << word_shift);
        //return (tmp_in[n_words - 1] & masks[(bits-read_bits)]) == read(i + read_bits, i+bits-1);
    }*/

    //compare the segment ]a-bits..a] with the segment ]b-bits..b+bits]
    //return the bit_pos (0-based) of the rightmost different bit (return len if the segments are equal)
    inline size_t inv_com_segments(size_t a, size_t b, size_t& bits) const {
        size_t n_words = INT_CEIL(bits, word_bits);
        size_t rem_bits, data_a, data_b, read_bits=0;

        if(n_words>1){

            size_t rs_a = a & (word_bits - 1UL);
            size_t ls_a = word_bits - rs_a - 1;
            size_t cell_a = a >> word_shift;

            size_t rs_b = b & (word_bits - 1UL);
            size_t ls_b = word_bits - rs_b - 1;
            size_t cell_b = b >> word_shift;

            for(size_t k=n_words; k-->1;){
                //the mask helps to deal with corner cases
                // (if the shift is 64-bits long, for instance)
                data_a = stream[cell_a] << ls_a;
                data_a |= (stream[--cell_a] & (masks[ls_a]<<(rs_a+1UL))) >> (rs_a+1UL);

                data_b = stream[cell_b] << ls_b;
                data_b |= (stream[--cell_b] & (masks[ls_b]<<(rs_b+1UL))) >> (rs_b+1UL);

                if(data_a != data_b){
                    return read_bits + __builtin_clzl(data_a^data_b);
                }
                read_bits+=word_bits;
            }
        }

        rem_bits = bits-read_bits;
        data_a = read(a-bits+1, a-read_bits);
        data_b = read(b-bits+1, b-read_bits);

        data_a^=data_b;
        if(data_a==0){
            return bits;
        }else{
            return read_bits + (rem_bits- ((8*sizeof(unsigned long) - __builtin_clzl(data_a))));
        }
    }

    size_type serialize(std::ostream &out) const{
        size_t written_bytes = serialize_elm(out, stream_size);
        out.write((char *)stream, sizeof(word_t)*stream_size);
        return written_bytes + (sizeof(word_t)*stream_size);
    }

    void load(std::istream &in){
        load_elm(in, stream_size);
        if(stream==nullptr){
            stream = (word_t *) malloc(sizeof(word_t)*stream_size);
        }else{
            stream = (word_t *) realloc(stream, sizeof(word_t)*stream_size);
        }

        in.read((char *)stream, sizeof(word_t)*stream_size);
    }
};

template<class word_t, uint8_t max_dist>
const size_t bitstream<word_t, max_dist>::masks[65]={0x0,
                               0x1,0x3, 0x7,0xF,
                               0x1F,0x3F, 0x7F,0xFF,
                               0x1FF,0x3FF, 0x7FF,0xFFF,
                               0x1FFF,0x3FFF, 0x7FFF,0xFFFF,
                               0x1FFFF,0x3FFFF, 0x7FFFF,0xFFFFF,
                               0x1FFFFF,0x3FFFFF, 0x7FFFFF,0xFFFFFF,
                               0x1FFFFFF,0x3FFFFFF, 0x7FFFFFF,0xFFFFFFF,
                               0x1FFFFFFF,0x3FFFFFFF, 0x7FFFFFFF,0xFFFFFFFF,
                               0x1FFFFFFFF,0x3FFFFFFFF, 0x7FFFFFFFF,0xFFFFFFFFF,
                               0x1FFFFFFFFF,0x3FFFFFFFFF, 0x7FFFFFFFFF,0xFFFFFFFFFF,
                               0x1FFFFFFFFFF,0x3FFFFFFFFFF, 0x7FFFFFFFFFF,0xFFFFFFFFFFF,
                               0x1FFFFFFFFFFF,0x3FFFFFFFFFFF, 0x7FFFFFFFFFFF,0xFFFFFFFFFFFF,
                               0x1FFFFFFFFFFFF,0x3FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,0xFFFFFFFFFFFFF,
                               0x1FFFFFFFFFFFFF,0x3FFFFFFFFFFFFF, 0x7FFFFFFFFFFFFF,0xFFFFFFFFFFFFFF,
                               0x1FFFFFFFFFFFFFF,0x3FFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFF,
                               0x1FFFFFFFFFFFFFFF,0x3FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF};

#endif //LPG_COMPRESSOR_BITSTREAM_H
