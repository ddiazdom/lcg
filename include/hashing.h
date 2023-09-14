//
// Created by Diaz, Diego on 19.4.2023.
//

#ifndef SIMPLE_EXAMPLE_HASHING_H
#define SIMPLE_EXAMPLE_HASHING_H
#include <random>
#include <cinttypes>

struct rand_order{
    size_t str_ptr;
    size_t str_len;
    size_t hash;
    size_t orig_order;
};

constexpr unsigned __int128 ultra_long_mersenne_number() {
    __uint128_t tmp = 1;
    tmp<<=89;
    tmp--;
    return tmp;
}

struct hashing {

    uint32_t str_chunk[64]={0};
    uint64_t a1[65]={0};
    uint64_t a2[65]={0};

    uint8_t mersenne_primes_uint8[4] = {3, 7, 31, 127};
    uint16_t mersenne_primes_uint16[1] = {8191};
    uint32_t mersenne_primes_uint32[3] = {131071, 524287, 2147483647};
    uint64_t mersenne_primes_uint64[1] = {2305843009213693951UL};

    unsigned __int128 mersenne_prime_uint128;
    unsigned __int128 mp_ul_a;
    unsigned __int128 mp_ul_b;
    unsigned __int128 mp_ul_c;

    uint64_t powers[65] = {0ULL, 1ULL, 3ULL, 7ULL, 15ULL, 31ULL, 63ULL, 127ULL, 255ULL, 511ULL, 1023ULL, 2047ULL,
                           4095ULL, 8191ULL, 16383ULL, 32767ULL, 65535ULL, 131071ULL, 262143ULL, 524287ULL, 1048575ULL,
                           2097151ULL, 4194303ULL, 8388607ULL, 16777215ULL, 33554431ULL, 67108863ULL, 134217727ULL,
                           268435455ULL, 536870911ULL, 1073741823ULL, 2147483647ULL, 4294967295ULL, 8589934591ULL,
                           17179869183ULL, 34359738367ULL, 68719476735ULL, 137438953471ULL, 274877906943ULL, 549755813887ULL,
                           1099511627775ULL, 2199023255551ULL, 4398046511103ULL, 8796093022207ULL, 17592186044415ULL,
                           35184372088831ULL, 70368744177663ULL, 140737488355327ULL, 281474976710655ULL, 562949953421311ULL,
                           1125899906842623ULL, 2251799813685247ULL, 4503599627370495ULL, 9007199254740991ULL, 18014398509481983ULL,
                           36028797018963967ULL, 72057594037927935ULL, 144115188075855871ULL, 288230376151711743ULL,
                           576460752303423487ULL, 1152921504606846975ULL, 2305843009213693951ULL, 4611686018427387903ULL,
                           9223372036854775807ULL, 18446744073709551615ULL};

    void generate_random_numbers(){
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<uint64_t> distrib(1, std::numeric_limits<uint64_t>::max());

        for(size_t i=0;i<65;i++){
            a1[i] = distrib(gen);
            while(!(a1[i] & 1UL)) a1[i] = distrib(gen);
            a2[i] = distrib(gen);
            while(!(a2[i] & 1UL)) a2[i] = distrib(gen);
        }

        std::mt19937 gen2(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<unsigned __int128> distrib2(1, mersenne_prime_uint128-1);
        mp_ul_a=distrib2(gen2);
        mp_ul_b=distrib2(gen2);
        mp_ul_c=distrib2(gen2);
    }

    hashing(): mersenne_prime_uint128(ultra_long_mersenne_number()),
               mp_ul_a{0},
               mp_ul_b{0},
               mp_ul_c{0}{
        generate_random_numbers();
    }

    static inline uint32_t pair_multiply_shift(const uint32_t *x, size_t len, uint32_t l, const uint64_t* a) {
        uint64_t acc = 0;
        for(size_t i=0;i<len;i+=2){
            acc+=(a[i]+x[i+1])*(a[i+1]+x[i]);
        }
        return (acc+a[len]) >> (64-l);
    }

    inline uint64_t short_string_short_hash(const char *string, size_t len, size_t l){
        assert(l<=32);
        size_t n_words = (len+3)>>2; // fast version of ceil(len, 4);
        str_chunk[n_words-1] = 0;
        memcpy(str_chunk, string, len);
        return pair_multiply_shift(str_chunk, n_words, l, a1);
    }

    inline uint64_t short_string_long_hash(const char *string, size_t len, size_t l){
        assert(l>32);
        size_t n_words = (len+3)>>2; // fast version of ceil(len, 4);
        str_chunk[n_words-1] = 0;
        memcpy(str_chunk, string, len);
        uint64_t hash1 = pair_multiply_shift(str_chunk, n_words, 32, a1);
        uint64_t hash2 = pair_multiply_shift(str_chunk, n_words, l-32, a2);
        return hash2<<32UL | hash1;
    }

    [[nodiscard]] inline uint64_t symbol_hash(uint64_t x, uint64_t l=64) const {
        return (a1[0]*x) >> (64-l);
    }

    //the resulting hash will be l-bits long
    uint64_t string_hash(const char * string, size_t len, size_t l){

        if(len<=256){//short string
            if(l<=32){
                return short_string_short_hash(string, len, l);
            }else{
                return short_string_long_hash(string, len, l);
            }
        } else{//long string
            size_t n_chunks = (len+255) >> 8;//fast version of ceil(len, 256);
            uint64_t xi;
            unsigned __int128 hash = 0;
            size_t chunk_len;

            for(size_t i = 0, j=0; i<n_chunks; i++, j+=256){
                chunk_len = std::min<size_t>(256, len);
                xi = short_string_long_hash(&string[j], chunk_len, 64);
                // (c*hash +xi) % p = 2^{89} -1
                hash = (mp_ul_c*hash + xi);
                hash = (hash & mersenne_prime_uint128) + (hash >> 89); // hash % p = x % 2^{89} + floor(hash/2^{89}) = (hash & (2^{89}-1)) + (x>>89)
                if(hash >= mersenne_prime_uint128) hash-=mersenne_prime_uint128;
                len-=chunk_len;
            }

            //composition with a multiple-mod-prime universal hash function
            hash = (mp_ul_a * hash) + mp_ul_b;
            hash = (hash & mersenne_prime_uint128) + (hash >> 89);
            if(hash >= mersenne_prime_uint128) hash-=mersenne_prime_uint128;

            hash = hash & powers[l];
            return (uint64_t) hash;
        }
    }
};
#endif //SIMPLE_EXAMPLE_HASHING_H
