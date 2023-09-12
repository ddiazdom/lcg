//
// Created by Diaz, Diego on 20.7.2023.
//

#ifndef TEXT_PARSER_INT_DECODERS_H
#define TEXT_PARSER_INT_DECODERS_H
#include <cstring>

template<class sym_t>
struct plain_decoder{

    typedef sym_t sym_type;
    static const size_t bytes = sizeof(sym_t);
    static const size_t shift = __builtin_ctz(bytes);

    inline static off_t forward(__attribute__((unused)) uint8_t *& ptr, off_t n_pos) {
        return bytes*n_pos;
    }

    inline static off_t backward(__attribute__((unused)) uint8_t * ptr, __attribute__((unused)) const uint8_t *& l_boundary, off_t n_pos) {
        return bytes*n_pos;
    }

    //rightmost byte of the rightmost valid value before ptr:
    // For instance, 0 0 0 1* 0 0* 0 1 (1 is the last byte of each symbol).
    // The function returns 3 with pos = 5 (assuming uint32_t symbols of sizeof(uint32_t)= 4 bytes)
    inline static off_t mov_to_prev_valid_byte(uint8_t *&ptr, off_t& pos) {
        off_t rem = (pos+1) % bytes;
        ptr-=rem;
        return rem;
    }

    inline static off_t write_forward(uint8_t * ptr, sym_type& sym) {
        memcpy(ptr, &sym, bytes);
        return bytes;
    }

    inline static off_t read_forward(const uint8_t * ptr, sym_type& sym) {
        memcpy(&sym, ptr, bytes);
        return bytes;
    }

    //this assumes ptr points to the rightmost byte of a symbol in the stream
    //notice this function changes the memory address
    inline static off_t read_backwards(uint8_t *&ptr, __attribute__((unused))  const uint8_t *& boundary, sym_type& sym) {
        ptr -= bytes-1;
        memcpy(&sym, ptr, bytes);
        return bytes-1;
    }

    //get the rightmost symbol that differs from the symbol in ptr
    //it scans the bytes preceding ptr to look for the first mismatching symbol
    //the function assumes ptr and boundary are addresses within the same memory area
    inline static off_t mov_to_prev_diff_sym(uint8_t *& ptr, const uint8_t *boundary, sym_type& sym) {

        //off_t n_comp = (ptr-boundary-1)/bytes;
        off_t n_comp = (ptr-boundary-1)>>shift;
        sym = 0;
        auto *tmp = ptr;
        ptr-=bytes;
        while(n_comp>0 && memcmp(tmp, ptr, bytes)==0){
            ptr-=bytes;
            n_comp--;
        }

        if(n_comp>0){
            assert(ptr!=boundary);
            memcpy(&sym, ptr, bytes);
        }

        return tmp-ptr;
    }

    //get the leftmost symbol that differs from the symbol in ptr
    //it will stop checking once the scan reaches boundary address
    //the function assumes ptr and boundary are addresses within the same memory area
    inline static off_t mov_to_next_diff_sym(uint8_t *& ptr, const uint8_t* boundary, sym_type& sym) {
        off_t n_comp = ((boundary-1)-ptr)>>shift;
        sym = 0;
        auto *tmp = ptr;
        ptr+=bytes;
        while(n_comp>0 && memcmp(tmp, ptr, bytes)==0){
            ptr+=bytes;
            n_comp--;
        }

        if(n_comp>0){
            assert(ptr!=boundary);
            memcpy(&sym, ptr, bytes);
        }
        return ptr-tmp;
    }
};
#endif //TEXT_PARSER_INT_DECODERS_H
