//
// Created by Diaz, Diego on 9.12.2024.
//

#ifndef LCG_FASTX_PARSER_H
#define LCG_FASTX_PARSER_H

#include <iostream>


static inline size_t fasta_parsing_scalar(uint8_t *stream, size_t size) {

    size_t i=0, pos=0;
    int in_header = 0;
    int prev_in_header;

    while(i<size){
        uint8_t c = stream[i++];
        stream[pos] = c;

        prev_in_header = in_header;
        //in_header=1 after the operations below means with are in a header area
        in_header += c=='>';
        in_header -= (in_header>0 && c=='\n');

        pos+= (c!='\n' && in_header==0) || //new lines in the middle of the sequences
              (prev_in_header>0 && in_header==0 && pos!=0);//transition from header to sequence
    }
    stream[pos++]='\n';
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

    size_t i = 0, pos = 0;
    uint8x16_t bitmask = {1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128 };
    uint8x16_t new_line = vdupq_n_u8('\n');
    uint8x16_t new_entry = vdupq_n_u8('>');
    int in_header = 0;
    int prev_header_area;

    while(i + 16 <= size) {

        uint8x16_t vec = vld1q_u8(stream + i);

        uint8x16_t has_new_entry = vceqq_u8(vec, new_entry);
        if(vmaxvq_u8(has_new_entry)){
            size_t end = std::min(i+15, size);
            while(end<(size-1) && stream[end]!='\n') end++;
            while(i<=end){
                uint8_t c = stream[i++];
                stream[pos] = c;

                prev_header_area = in_header;

                //in_header=1 after the operations below means with are in a header area
                in_header += c=='>';
                in_header -= (in_header>0 && c=='\n');

                pos+= (c!='\n' && in_header==0) || //new lines in the middle of the sequences
                      (prev_header_area>0 && in_header==0 && pos!=0);//transition from header to sequence
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

        vst1_u8(stream + pos, low);
        pos += __builtin_popcount(mlow);

        vst1_u8(stream + pos, high);
        pos += __builtin_popcount(mhigh);
        i+=16;
    }

    while(i<size){
        uint8_t c = stream[i++];
        stream[pos] = c;
        prev_header_area = in_header;

        //in_header=1 after the operations below means with are in a header area
        in_header += c=='>';
        in_header -= (in_header>0 && c=='\n');

        pos+= (c!='\n' && in_header==0) || //new lines in the middle of the sequences
              (prev_header_area>0 && in_header==0 && pos!=0);//transition from header to sequence
    }
    stream[pos++]='\n';
    return pos;
}
#endif

#ifdef __SSE4_2__

#include <x86intrin.h>
#include "simd_tables.h"

static inline size_t fasta_parsing_sse42(uint8_t *stream, size_t size) {

    assert(stream[0]=='>');
    size_t pos = 0, i=0;
    __m128i new_entry = _mm_set1_epi8('>');
    __m128i new_line = _mm_set1_epi8('\n');

    int in_header = 0;
    int prev_header_area;

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

                prev_header_area = in_header;

                //in_header=1 after the operations below means with are in a header area
                in_header += c=='>';
                in_header -= (in_header>0 && c=='\n');

                pos+= (c!='\n' && in_header==0) || //new lines in the middle of the sequences
                      (prev_header_area>0 && in_header==0 && pos!=0);//transition from header to sequence
            }
            continue;
        }

        //remove new lines
        int newlinemask = _mm_movemask_epi8(_mm_cmpeq_epi8(x, new_line));
        x = _mm_shuffle_epi8(x, _mm_loadu_si128((const __m128i *)despace_mask16 + (newlinemask & 0x7fff)));
        _mm_storeu_si128((__m128i *)(stream + pos), x);
        pos += 16 - _mm_popcnt_u32(newlinemask);
        i+=16;
    }

    while(i<size){
        uint8_t c = stream[i++];
        stream[pos] = c;
        prev_header_area = in_header;

        //in_header=1 after the operations below means with are in a header area
        in_header += c=='>';
        in_header -= (in_header>0 && c=='\n');

        pos+= (c!='\n' && in_header==0) || //new lines in the middle of the sequences
              (prev_header_area>0 && in_header==0 && pos!=0);//transition from header to sequence
    }
    stream[pos++]='\n';
    return pos;
}
#endif

#ifdef __AVX2__
#include <x86intrin.h>
#include "simd_tables.h"
#include <immintrin.h>

static inline size_t fasta_parsing_avx2(uint8_t *stream, size_t size) {

    assert(stream[0]=='>');
    size_t pos = 0, i=0;

    __m256i newline = _mm256_set1_epi8('\n');
    __m256i new_entry = _mm256_set1_epi8('>');
    int in_header = 0;
    int prev_header_area;

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

                prev_header_area = in_header;

                //in_header=1 after the operations below means with are in a header area
                in_header += c=='>';
                in_header -= (in_header>0 && c=='\n');

                pos+= (c!='\n' && in_header==0) || //new lines in the middle of the sequences
                      (prev_header_area>0 && in_header==0 && pos!=0);//transition from header to sequence
            }
            continue;
        }

        //remove new lines
        uint32_t newlinemask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(x, newline));
        if(!newlinemask) { // no newline
            _mm256_storeu_si256((__m256i *)(stream + pos), x);
            pos += 32;
        } else {
            unsigned int maskhigh = (newlinemask) >> 16;
            unsigned int masklow = (newlinemask)&0xFFFF;

            //TODO the load below does not work with g++ 9.4 (probably with older versions either)
            __m256i mask = _mm256_loadu2_m128i((const __m128i *)despace_mask16 + (maskhigh & 0x7fff),
                                               (const __m128i *)despace_mask16 + (masklow & 0x7fff));
            x = _mm256_shuffle_epi8(x, mask);
            int offset1 = 16 - _mm_popcnt_u32(masklow);
            int offset2 = 16 - _mm_popcnt_u32(maskhigh);
            _mm256_storeu2_m128i((__m128i *)(stream + pos + offset1),
                                 (__m128i *)(stream + pos), x);
            pos += offset1 + offset2;
        }
        i+=32;
    }

    while(i<size){
        uint8_t c = stream[i++];
        stream[pos] = c;
        prev_header_area = in_header;

        //in_header=1 after the operations below means with are in a header area
        in_header += c=='>';
        in_header -= (in_header>0 && c=='\n');

        pos+= (c!='\n' && in_header==0) || //new lines in the middle of the sequences
              (prev_header_area>0 && in_header==0 && pos!=0);//transition from header to sequence
    }
    stream[pos++]='\n';
    return pos;
}
#endif


// Detect ARM NEON
#if defined(USE_AVX2)
#define PARSE_FASTA(param1, param2) fasta_parsing_avx2(param1, param2)
#elif defined(USE_SSE4_2)
#define PARSE_FASTA(param1, param2) fasta_parsing_sse42(param1, param2)
#elif defined(USE_NEON)
#define PARSE_FASTA(param1, param2) fasta_parsing_neon(param1, param2)
#else
#define PARSE_FASTA(param1, param2) fasta_parsing_scalar(param1, param2)
#endif

#endif
