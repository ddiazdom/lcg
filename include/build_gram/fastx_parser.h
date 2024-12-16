//
// Created by Diaz, Diego on 9.12.2024.
//

#ifndef LCG_FASTX_PARSER_H
#define LCG_FASTX_PARSER_H

#include <iostream>

static inline size_t fasta_parsing_scalar(uint8_t *stream, size_t size) {

    size_t i=0, pos=0, seq_len=0, n_empty=0;
    int in_header = 0, prev_in_header;
    bool h2s_tran, emp_entry;

    while(i<size){
        uint8_t c = stream[i++];
        stream[pos] = c;

        prev_in_header = in_header;
        //in_header=1 after the operations below means with are in a header area
        in_header += c=='>';
        in_header -= (in_header>0 && c=='\n');

        h2s_tran = (prev_in_header>0 && in_header==0);//transition from header to sequence
        emp_entry = h2s_tran && seq_len==0;
        n_empty+= emp_entry;
        seq_len*= !h2s_tran;

        pos+= (c!='\n' && in_header==0) || (h2s_tran && !emp_entry);
        seq_len += (c!='\n' && in_header==0);
    }
    stream[pos++]='\n';
    std::cout<<"There are "<<n_empty-1<<" entries "<<std::endl;
    return pos;
}

#ifdef __ARM_NEON__
#include <arm_neon.h>
static const uint8_t __attribute__((aligned(16))) mask_shuffle[256*8] = {
        0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 1,0,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,
        2,0,0,0,0,0,0,0, 0,2,0,0,0,0,0,0, 1,2,0,0,0,0,0,0, 0,1,2,0,0,0,0,0,
        3,0,0,0,0,0,0,0, 0,3,0,0,0,0,0,0, 1,3,0,0,0,0,0,0, 0,1,3,0,0,0,0,0,
        2,3,0,0,0,0,0,0, 0,2,3,0,0,0,0,0, 1,2,3,0,0,0,0,0, 0,1,2,3,0,0,0,0,
        4,0,0,0,0,0,0,0, 0,4,0,0,0,0,0,0, 1,4,0,0,0,0,0,0, 0,1,4,0,0,0,0,0,
        2,4,0,0,0,0,0,0, 0,2,4,0,0,0,0,0, 1,2,4,0,0,0,0,0, 0,1,2,4,0,0,0,0,
        3,4,0,0,0,0,0,0, 0,3,4,0,0,0,0,0, 1,3,4,0,0,0,0,0, 0,1,3,4,0,0,0,0,
        2,3,4,0,0,0,0,0, 0,2,3,4,0,0,0,0, 1,2,3,4,0,0,0,0, 0,1,2,3,4,0,0,0,
        5,0,0,0,0,0,0,0, 0,5,0,0,0,0,0,0, 1,5,0,0,0,0,0,0, 0,1,5,0,0,0,0,0,
        2,5,0,0,0,0,0,0, 0,2,5,0,0,0,0,0, 1,2,5,0,0,0,0,0, 0,1,2,5,0,0,0,0,
        3,5,0,0,0,0,0,0, 0,3,5,0,0,0,0,0, 1,3,5,0,0,0,0,0, 0,1,3,5,0,0,0,0,
        2,3,5,0,0,0,0,0, 0,2,3,5,0,0,0,0, 1,2,3,5,0,0,0,0, 0,1,2,3,5,0,0,0,
        4,5,0,0,0,0,0,0, 0,4,5,0,0,0,0,0, 1,4,5,0,0,0,0,0, 0,1,4,5,0,0,0,0,
        2,4,5,0,0,0,0,0, 0,2,4,5,0,0,0,0, 1,2,4,5,0,0,0,0, 0,1,2,4,5,0,0,0,
        3,4,5,0,0,0,0,0, 0,3,4,5,0,0,0,0, 1,3,4,5,0,0,0,0, 0,1,3,4,5,0,0,0,
        2,3,4,5,0,0,0,0, 0,2,3,4,5,0,0,0, 1,2,3,4,5,0,0,0, 0,1,2,3,4,5,0,0,
        6,0,0,0,0,0,0,0, 0,6,0,0,0,0,0,0, 1,6,0,0,0,0,0,0, 0,1,6,0,0,0,0,0,
        2,6,0,0,0,0,0,0, 0,2,6,0,0,0,0,0, 1,2,6,0,0,0,0,0, 0,1,2,6,0,0,0,0,
        3,6,0,0,0,0,0,0, 0,3,6,0,0,0,0,0, 1,3,6,0,0,0,0,0, 0,1,3,6,0,0,0,0,
        2,3,6,0,0,0,0,0, 0,2,3,6,0,0,0,0, 1,2,3,6,0,0,0,0, 0,1,2,3,6,0,0,0,
        4,6,0,0,0,0,0,0, 0,4,6,0,0,0,0,0, 1,4,6,0,0,0,0,0, 0,1,4,6,0,0,0,0,
        2,4,6,0,0,0,0,0, 0,2,4,6,0,0,0,0, 1,2,4,6,0,0,0,0, 0,1,2,4,6,0,0,0,
        3,4,6,0,0,0,0,0, 0,3,4,6,0,0,0,0, 1,3,4,6,0,0,0,0, 0,1,3,4,6,0,0,0,
        2,3,4,6,0,0,0,0, 0,2,3,4,6,0,0,0, 1,2,3,4,6,0,0,0, 0,1,2,3,4,6,0,0,
        5,6,0,0,0,0,0,0, 0,5,6,0,0,0,0,0, 1,5,6,0,0,0,0,0, 0,1,5,6,0,0,0,0,
        2,5,6,0,0,0,0,0, 0,2,5,6,0,0,0,0, 1,2,5,6,0,0,0,0, 0,1,2,5,6,0,0,0,
        3,5,6,0,0,0,0,0, 0,3,5,6,0,0,0,0, 1,3,5,6,0,0,0,0, 0,1,3,5,6,0,0,0,
        2,3,5,6,0,0,0,0, 0,2,3,5,6,0,0,0, 1,2,3,5,6,0,0,0, 0,1,2,3,5,6,0,0,
        4,5,6,0,0,0,0,0, 0,4,5,6,0,0,0,0, 1,4,5,6,0,0,0,0, 0,1,4,5,6,0,0,0,
        2,4,5,6,0,0,0,0, 0,2,4,5,6,0,0,0, 1,2,4,5,6,0,0,0, 0,1,2,4,5,6,0,0,
        3,4,5,6,0,0,0,0, 0,3,4,5,6,0,0,0, 1,3,4,5,6,0,0,0, 0,1,3,4,5,6,0,0,
        2,3,4,5,6,0,0,0, 0,2,3,4,5,6,0,0, 1,2,3,4,5,6,0,0, 0,1,2,3,4,5,6,0,
        7,0,0,0,0,0,0,0, 0,7,0,0,0,0,0,0, 1,7,0,0,0,0,0,0, 0,1,7,0,0,0,0,0,
        2,7,0,0,0,0,0,0, 0,2,7,0,0,0,0,0, 1,2,7,0,0,0,0,0, 0,1,2,7,0,0,0,0,
        3,7,0,0,0,0,0,0, 0,3,7,0,0,0,0,0, 1,3,7,0,0,0,0,0, 0,1,3,7,0,0,0,0,
        2,3,7,0,0,0,0,0, 0,2,3,7,0,0,0,0, 1,2,3,7,0,0,0,0, 0,1,2,3,7,0,0,0,
        4,7,0,0,0,0,0,0, 0,4,7,0,0,0,0,0, 1,4,7,0,0,0,0,0, 0,1,4,7,0,0,0,0,
        2,4,7,0,0,0,0,0, 0,2,4,7,0,0,0,0, 1,2,4,7,0,0,0,0, 0,1,2,4,7,0,0,0,
        3,4,7,0,0,0,0,0, 0,3,4,7,0,0,0,0, 1,3,4,7,0,0,0,0, 0,1,3,4,7,0,0,0,
        2,3,4,7,0,0,0,0, 0,2,3,4,7,0,0,0, 1,2,3,4,7,0,0,0, 0,1,2,3,4,7,0,0,
        5,7,0,0,0,0,0,0, 0,5,7,0,0,0,0,0, 1,5,7,0,0,0,0,0, 0,1,5,7,0,0,0,0,
        2,5,7,0,0,0,0,0, 0,2,5,7,0,0,0,0, 1,2,5,7,0,0,0,0, 0,1,2,5,7,0,0,0,
        3,5,7,0,0,0,0,0, 0,3,5,7,0,0,0,0, 1,3,5,7,0,0,0,0, 0,1,3,5,7,0,0,0,
        2,3,5,7,0,0,0,0, 0,2,3,5,7,0,0,0, 1,2,3,5,7,0,0,0, 0,1,2,3,5,7,0,0,
        4,5,7,0,0,0,0,0, 0,4,5,7,0,0,0,0, 1,4,5,7,0,0,0,0, 0,1,4,5,7,0,0,0,
        2,4,5,7,0,0,0,0, 0,2,4,5,7,0,0,0, 1,2,4,5,7,0,0,0, 0,1,2,4,5,7,0,0,
        3,4,5,7,0,0,0,0, 0,3,4,5,7,0,0,0, 1,3,4,5,7,0,0,0, 0,1,3,4,5,7,0,0,
        2,3,4,5,7,0,0,0, 0,2,3,4,5,7,0,0, 1,2,3,4,5,7,0,0, 0,1,2,3,4,5,7,0,
        6,7,0,0,0,0,0,0, 0,6,7,0,0,0,0,0, 1,6,7,0,0,0,0,0, 0,1,6,7,0,0,0,0,
        2,6,7,0,0,0,0,0, 0,2,6,7,0,0,0,0, 1,2,6,7,0,0,0,0, 0,1,2,6,7,0,0,0,
        3,6,7,0,0,0,0,0, 0,3,6,7,0,0,0,0, 1,3,6,7,0,0,0,0, 0,1,3,6,7,0,0,0,
        2,3,6,7,0,0,0,0, 0,2,3,6,7,0,0,0, 1,2,3,6,7,0,0,0, 0,1,2,3,6,7,0,0,
        4,6,7,0,0,0,0,0, 0,4,6,7,0,0,0,0, 1,4,6,7,0,0,0,0, 0,1,4,6,7,0,0,0,
        2,4,6,7,0,0,0,0, 0,2,4,6,7,0,0,0, 1,2,4,6,7,0,0,0, 0,1,2,4,6,7,0,0,
        3,4,6,7,0,0,0,0, 0,3,4,6,7,0,0,0, 1,3,4,6,7,0,0,0, 0,1,3,4,6,7,0,0,
        2,3,4,6,7,0,0,0, 0,2,3,4,6,7,0,0, 1,2,3,4,6,7,0,0, 0,1,2,3,4,6,7,0,
        5,6,7,0,0,0,0,0, 0,5,6,7,0,0,0,0, 1,5,6,7,0,0,0,0, 0,1,5,6,7,0,0,0,
        2,5,6,7,0,0,0,0, 0,2,5,6,7,0,0,0, 1,2,5,6,7,0,0,0, 0,1,2,5,6,7,0,0,
        3,5,6,7,0,0,0,0, 0,3,5,6,7,0,0,0, 1,3,5,6,7,0,0,0, 0,1,3,5,6,7,0,0,
        2,3,5,6,7,0,0,0, 0,2,3,5,6,7,0,0, 1,2,3,5,6,7,0,0, 0,1,2,3,5,6,7,0,
        4,5,6,7,0,0,0,0, 0,4,5,6,7,0,0,0, 1,4,5,6,7,0,0,0, 0,1,4,5,6,7,0,0,
        2,4,5,6,7,0,0,0, 0,2,4,5,6,7,0,0, 1,2,4,5,6,7,0,0, 0,1,2,4,5,6,7,0,
        3,4,5,6,7,0,0,0, 0,3,4,5,6,7,0,0, 1,3,4,5,6,7,0,0, 0,1,3,4,5,6,7,0,
        2,3,4,5,6,7,0,0, 0,2,3,4,5,6,7,0, 1,2,3,4,5,6,7,0, 0,1,2,3,4,5,6,7,
};

// credit: Martins Mozeiko
static inline size_t fasta_parsing_neon(uint8_t *stream, size_t size) {

    assert(stream[0]=='>');

    size_t i = 0, pos = 0, seq_len=0, n_empty=0, prev_pos;
    uint8x16_t bitmask = {1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128 };
    uint8x16_t new_line = vdupq_n_u8('\n');
    uint8x16_t new_entry = vdupq_n_u8('>');
    int in_header = 0, prev_in_header;
    bool h2s_tran, emp_entry;

    while(i + 16 <= size) {

        uint8x16_t vec = vld1q_u8(stream + i);

        uint8x16_t has_new_entry = vceqq_u8(vec, new_entry);
        if(vmaxvq_u8(has_new_entry)){
            size_t end = std::min(i+15, size);
            while(end<(size-1) && stream[end]!='\n') end++;

            while(i<=end){
                uint8_t c = stream[i++];
                stream[pos] = c;

                prev_in_header = in_header;
                //in_header=1 after the operations below means with are in a header area
                in_header += c=='>';
                in_header -= (in_header>0 && c=='\n');

                h2s_tran = (prev_in_header>0 && in_header==0);//transition from header to sequence
                emp_entry = h2s_tran && seq_len==0;
                n_empty+= emp_entry;
                seq_len*= !h2s_tran;

                pos+= (c!='\n' && in_header==0) || (h2s_tran && !emp_entry);
                seq_len += (c!='\n' && in_header==0);
            }
            continue;
        }

        //remove new lines
        uint8x16_t cmp = vmvnq_u8(vceqq_u8(vec, new_line));
        uint64x2_t mask = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(vandq_u8(cmp, bitmask))));

        uint8_t mlow = vgetq_lane_u8(vreinterpretq_u8_u64(mask), 0);
        uint8_t mhigh = vgetq_lane_u8(vreinterpretq_u8_u64(mask), 8);

        uint8x8_t low = vtbl1_u8(vget_low_u8(vec), vld1_u8(mask_shuffle + mlow*8));
        uint8x8_t high = vtbl1_u8(vget_high_u8(vec), vld1_u8(mask_shuffle + mhigh*8));

        prev_pos = pos;
        vst1_u8(stream + pos, low);
        pos += __builtin_popcount(mlow);

        vst1_u8(stream + pos, high);
        pos += __builtin_popcount(mhigh);
        seq_len +=pos-prev_pos;
        i+=16;
    }

    while(i<size){
        uint8_t c = stream[i++];
        stream[pos] = c;

        prev_in_header = in_header;
        //in_header=1 after the operations below means with are in a header area
        in_header += c=='>';
        in_header -= (in_header>0 && c=='\n');

        h2s_tran = (prev_in_header>0 && in_header==0);//transition from header to sequence
        emp_entry = h2s_tran && seq_len==0;
        n_empty+= emp_entry;
        seq_len*= !h2s_tran;

        pos+= (c!='\n' && in_header==0) || (h2s_tran && !emp_entry);
        seq_len += (c!='\n' && in_header==0);
    }
    stream[pos++]='\n';
    std::cout<<"There are "<<n_empty-1<<" entries "<<std::endl;
    return pos;
}
#endif

#ifdef __SSE4_2__

#include <x86intrin.h>
#include "simd_tables.h"

static inline size_t fasta_parsing_sse42(uint8_t *stream, size_t size) {

    assert(stream[0]=='>');

    size_t i = 0, pos = 0, seq_len=0, n_empty=0, valid_syms;
    __m128i new_entry = _mm_set1_epi8('>');
    __m128i new_line = _mm_set1_epi8('\n');

    int in_header = 0, prev_in_header;
    bool h2s_tran, emp_entry;

    while(i + 16 <= size) {

        //load 16 symbols starting from string[i]
        __m128i x = _mm_loadu_si128((const __m128i *)(stream + i));

        //create a mask with the position equal to a '>' symbol
        int new_entry_mask = _mm_movemask_epi8(_mm_cmpeq_epi8(x, new_entry));

        //handle the header area
        if(new_entry_mask){
            size_t end = std::min(i+15, size);
            while(end<(size-1) && stream[end]!='\n') end++;
            while(i<=end){
                uint8_t c = stream[i++];
                stream[pos] = c;

                prev_in_header = in_header;
                //in_header=1 after the operations below means with are in a header area
                in_header += c=='>';
                in_header -= (in_header>0 && c=='\n');

                h2s_tran = (prev_in_header>0 && in_header==0);//transition from header to sequence
                emp_entry = h2s_tran && seq_len==0;
                n_empty+= emp_entry;
                seq_len*= !h2s_tran;

                pos+= (c!='\n' && in_header==0) || (h2s_tran && !emp_entry);
                seq_len += (c!='\n' && in_header==0);
            }
            continue;
        }

        //remove new lines
        int newlinemask = _mm_movemask_epi8(_mm_cmpeq_epi8(x, new_line));
        x = _mm_shuffle_epi8(x, _mm_loadu_si128((const __m128i *)despace_mask16 + (newlinemask & 0x7fff)));
        _mm_storeu_si128((__m128i *)(stream + pos), x);
        valid_syms = 16 - _mm_popcnt_u32(newlinemask);
        pos += valid_syms;
        seq_len+ = valid_syms;
        i+=16;
    }

    while(i<size){
        uint8_t c = stream[i++];
        stream[pos] = c;

        prev_in_header = in_header;
        //in_header=1 after the operations below means with are in a header area
        in_header += c=='>';
        in_header -= (in_header>0 && c=='\n');

        h2s_tran = (prev_in_header>0 && in_header==0);//transition from header to sequence
        emp_entry = h2s_tran && seq_len==0;
        n_empty+= emp_entry;
        seq_len*= !h2s_tran;

        pos+= (c!='\n' && in_header==0) || (h2s_tran && !emp_entry);
        seq_len += (c!='\n' && in_header==0);
    }
    stream[pos++]='\n';
    std::cout<<"There are "<<n_empty-1<<" entries "<<std::endl;
    return pos;
}
#endif

#ifdef __AVX2__
#include <x86intrin.h>
#include "simd_tables.h"
#include <immintrin.h>

/* GCC <10.1 doesn't include the split store/load intrinsics so define them here. */
#if defined(__GNUC__) &&  __GNUC__ < 10
static inline void __attribute__((__always_inline__))
_mm256_storeu2_m128i(__m128i* const hi, __m128i* const lo, const __m256i a) {
    _mm_storeu_si128(lo, _mm256_castsi256_si128(a));
    _mm_storeu_si128(hi, _mm256_extracti128_si256(a, 1));
}

static inline __m256i
_mm256_loadu2_m128i(__m128i const* __addr_hi, __m128i const* __addr_lo) {
    __m256i __v256 = _mm256_castsi128_si256(_mm_loadu_si128(__addr_lo));
    return _mm256_insertf128_si256(__v256, _mm_loadu_si128(__addr_hi), 1);
}
#endif /* defined(__GNUC__) */

static inline size_t fasta_parsing_avx2(uint8_t *stream, size_t size) {

    assert(stream[0]=='>');

    size_t i = 0, pos = 0, seq_len=0, n_empty=0, n_valid;
    __m256i newline = _mm256_set1_epi8('\n');
    __m256i new_entry = _mm256_set1_epi8('>');
    int in_header = 0, prev_in_header;

    while(i + 32 <= size) {

        //load 32 symbols starting from stream[i]
        __m256i x = _mm256_loadu_si256((const __m256i *)(stream + i));

        //create a mask with the position equal to a '>' symbol
        int new_entry_mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(x, new_entry));

        //handle the header area
        if(new_entry_mask){
            size_t end = std::min(i+31, size);
            while(end<(size-1) && stream[end]!='\n') end++;
            while(i<=end){
                uint8_t c = stream[i++];
                stream[pos] = c;

                prev_in_header = in_header;
                //in_header=1 after the operations below means with are in a header area
                in_header += c=='>';
                in_header -= (in_header>0 && c=='\n');

                h2s_tran = (prev_in_header>0 && in_header==0);//transition from header to sequence
                emp_entry = h2s_tran && seq_len==0;
                n_empty+= emp_entry;
                seq_len*= !h2s_tran;

                pos+= (c!='\n' && in_header==0) || (h2s_tran && !emp_entry);
                seq_len += (c!='\n' && in_header==0);
            }
            continue;
        }

        //remove new lines
        uint32_t newlinemask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(x, newline));
        if(!newlinemask) { // no newline
            _mm256_storeu_si256((__m256i *)(stream + pos), x);
            n_valid=32;
        } else {
            unsigned int maskhigh = (newlinemask) >> 16;
            unsigned int masklow = (newlinemask)&0xFFFF;

            __m256i mask = _mm256_loadu2_m128i((const __m128i *)despace_mask16 + (maskhigh & 0x7fff),
                                               (const __m128i *)despace_mask16 + (masklow & 0x7fff));

            x = _mm256_shuffle_epi8(x, mask);
            int offset1 = 16 - _mm_popcnt_u32(masklow);
            int offset2 = 16 - _mm_popcnt_u32(maskhigh);
            _mm256_storeu2_m128i((__m128i *)(stream + pos + offset1),
                                 (__m128i *)(stream + pos), x);
            n_valid=offset1+offset2;
        }

        pos += n_valid;
        seq_len += n_valid;
        i+=32;
    }

    while(i<size){
        uint8_t c = stream[i++];
        stream[pos] = c;

        prev_in_header = in_header;
        //in_header=1 after the operations below means with are in a header area
        in_header += c=='>';
        in_header -= (in_header>0 && c=='\n');

        h2s_tran = (prev_in_header>0 && in_header==0);//transition from header to sequence
        emp_entry = h2s_tran && seq_len==0;
        n_empty+= emp_entry;
        seq_len*= !h2s_tran;

        pos+= (c!='\n' && in_header==0) || (h2s_tran && !emp_entry);
        seq_len += (c!='\n' && in_header==0);
    }
    stream[pos++]='\n';
    std::cout<<"There are "<<n_empty-1<<" entries "<<std::endl;
    return pos;
}
#endif

#if defined(__AVX2__)
#define PARSE_FASTA(param1, param2) fasta_parsing_avx2(param1, param2)
#elif defined(__SSE4_2__)
#define PARSE_FASTA(param1, param2) fasta_parsing_sse42(param1, param2)
#elif defined(__ARM_NEON__) || defined(__ARM_NEON)
#define PARSE_FASTA(param1, param2) fasta_parsing_neon(param1, param2)
#else
#define PARSE_FASTA(param1, param2) fasta_parsing_scalar(param1, param2)
#endif

#endif
