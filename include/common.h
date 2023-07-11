//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef GRLBWT_COMMON_H
#define GRLBWT_COMMON_H

//#include <sdsl/bit_vectors.hpp>
#include "hash_table.hpp"
#include "int_array.h"
#include "rank_support.h"
#include "file_streams.hpp"

#define BUFFER_SIZE 8388608 //8MB of buffer

//typedef sdsl::bit_vector                        bv_t;
//typedef sdsl::bit_vector::rank_1_type           bv_rs_t;
typedef int_array<size_t>                       vector_t;
typedef int_array<size_t>                       string_t;
typedef hash_table<size_t, 44>                  phrase_map_t;
typedef buffered_hash_table<size_t, 44>         buffered_map_t;
typedef typename phrase_map_t::buff_t           ht_buff_t;

// the phrases are stored in a bit-compressed hash table:
// this wrapper reinterprets the bits back as phrases
struct key_wrapper{
    size_t width;
    size_t d_bits;//bits used to describe the string
    const bitstream<ht_buff_t>& stream;

    //offset is the bit where the key description starts
    [[nodiscard]] inline size_t read(size_t offset, size_t idx) const {
        return stream.read(offset + d_bits + (idx * width),
                           offset + d_bits + ((idx + 1) * width - 1));
    }

    //offset is the bit where the key description starts
    [[nodiscard]] inline size_t size(size_t offset) const{
        return stream.read(offset, offset + d_bits - 1) / width;
    }

    //compare two phrases are different positions of the bit stream
    [[nodiscard]] inline bool compare(size_t a, size_t b) const{

        size_t a_bits = stream.read(a, a + d_bits - 1);
        size_t b_bits = stream.read(b, b + d_bits - 1);
        size_t min_bits = std::min<size_t>(a_bits, b_bits);

        size_t a_pos = a+d_bits+a_bits-1;
        size_t b_pos = b+d_bits+b_bits-1;
        size_t rm_diff = stream.inv_com_segments(a_pos, b_pos, min_bits);

        if(rm_diff < min_bits){
            a_pos = a+d_bits+(((a_bits - rm_diff-1) / width) * width);
            b_pos = b+d_bits+(((b_bits - rm_diff-1) / width) * width);
            size_t sym_a = stream.read(a_pos, a_pos+width-1);
            size_t sym_b = stream.read(b_pos, b_pos+width-1);
            return sym_a<sym_b;
        }else{
            return a_bits>b_bits;
        }
    }
};
#endif //GRLBWT_COMMON_H
