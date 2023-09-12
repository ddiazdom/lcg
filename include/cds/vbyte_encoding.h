//
// Created by Diaz, Diego on 28.4.2023.
//

#ifndef SIMPLE_EXAMPLE_VBYTE_ENCODING_H
#define SIMPLE_EXAMPLE_VBYTE_ENCODING_H
#include<iostream>

/*
//maximum n-bytes integer we can encode in vbyte format
size_t max_vbyte_symbol(size_t n_bytes){
    assert(n_bytes>0 && n_bytes<=8);
    size_t n_bits = (n_bytes*8)-n_bytes;
    return (1UL<<n_bits)-1UL;
}

void vbyte_encode(size_t x, size_t& comp_x, size_t& len){
    size_t c=0, b=0;
    while(x>=128){
        c += (x & 127)<<b;
        b+=8;
        x>>=7;
    }
    comp_x = c+((x+128)<<b);
    len = (b+8);
}

template<class stream_t>
off_t vbyte_encode_stream(const stream_t* stream, size_t len, uint8_t * comp_stream){
    size_t tot_bytes = 0, code_len;
    for(size_t i=0;i<len;i++){
        size_t c=0, b=0, x = stream[i], comp_x=0;
        while(x>=128){
            c += (x & 127)<<b;
            b+=8;
            x>>=7;
        }
        comp_x = c+((x+128)<<b);
        code_len = (b+8)>>3;

        memcpy((char *)(comp_stream+tot_bytes), (void *)&comp_x,  code_len);
        tot_bytes += code_len;
    }
    return (off_t)tot_bytes;
}

uint8_t vbyte_len(size_t value){
    if(value==0) return 1;
    uint8_t n_bits = sym_width(value);
    n_bits+=((n_bits+7)/8);
    return (n_bits+7)/8;
}

off_t next_valid_pos(const uint8_t* stream, off_t stream_len, off_t byte_pos) {
    while(byte_pos>0 && byte_pos < stream_len && stream[byte_pos-1]<128) byte_pos++;
    return byte_pos<stream_len? byte_pos : -1;
}

inline static off_t mov_to_prev_valid_byte(uint8_t *&ptr, off_t& pos) {
    while(pos>=0 && *ptr<128) pos--;
    return pos>0 ? pos : -1;
}

template<typename int_type>
off_t write_forward(uint8_t *& stream, int_type x){
    uint8_t len;
    size_t comp_x;
    vbyte_encode(x, comp_x, len);
    memcpy(stream, &comp_x, len);
    return len;
}

template<typename int_type>
off_t read_forward(const uint8_t * stream, int_type& x){

    off_t i=0;
    uint8_t max_code_len = 8;
    while(max_code_len>0){
        if(stream[i]>=128) break;
        max_code_len--;
        i++;
    }
    off_t len = i+1;

    switch (len) {
        case 1:
            x = stream[i]-128;
            break;
        case 2:
            x = stream[i--]-128;
            x = (x<<7) + stream[i];
            break;
        case 3:
            x = stream[i--]-128;
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i];
            break;
        case 4:
            x = stream[i--]-128;
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i];
            break;
        case 5:
            x = stream[i--]-128;
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i];
            break;
        case 6:
            x = stream[i--]-128;
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i];
            break;
        case 7:
            x = stream[i--]-128;
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i];
            break;
        case 8 :
            x = stream[i--]-128;
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i--];
            x = (x<<7) + stream[i];
            break;
        default:
            exit(1);
    }
    return len;
}

template<typename int_type>
off_t vbyte_decode_inverse(const uint8_t * stream, off_t& pos, int_type& x){
    assert(pos>=0 && stream[pos]>=128);
    off_t len=1;
    x = stream[pos--]-128;
    while(pos>=0 && stream[pos]<128){
        x = (x<<7) + stream[pos--];
        len++;
    }
    return len;
}

template<class stream_t>
off_t vbyte_decode_stream(const uint8_t* comp_stream, off_t comp_stream_len, stream_t* decomp_stream, off_t& decomp_len){
    off_t tail = comp_stream_len-1;
    off_t consumed_bytes=0;
    while(tail>=0 && comp_stream[tail]<128) tail--;
    decomp_len = 0;
    assert(tail>=0 && comp_stream[tail]>=128);
    while(consumed_bytes<=tail){
        consumed_bytes +=vbyte_decode(&comp_stream[consumed_bytes], decomp_stream[decomp_len++]);
    }
    return consumed_bytes;
}

template<typename int_type>
struct vbyte_forward_iterator {

    off_t byte_pos;
    off_t stream_len;
    off_t len{};
    uint8_t const* stream = nullptr;
    int_type curr_val = std::numeric_limits<int_type>::max();

    vbyte_forward_iterator(const uint8_t *_stream, off_t _stream_len, off_t _byte_pos) : byte_pos(_byte_pos),
                                                                                         stream_len(_stream_len),
                                                                                         stream(_stream) {

        assert(_byte_pos>=0);
        while(byte_pos>0 && byte_pos < stream_len && stream[byte_pos-1]<128){
            byte_pos++;
        }

        if(byte_pos<stream_len){
            len = read_forward<int_type>(&stream[byte_pos], curr_val);
        }
    }

    void operator++(){
        byte_pos+=len;
        if(byte_pos<stream_len){
            len = read_forward<int_type>(&stream[byte_pos], curr_val);
        }else{
            byte_pos=stream_len;
        }
    }

    inline bool operator==(vbyte_forward_iterator& other) const{
        return byte_pos==other.byte_pos;
    }

    inline bool operator!=(vbyte_forward_iterator& other) const{
        return byte_pos!=other.byte_pos;
    }

    inline bool operator<(vbyte_forward_iterator& other) const{
        return byte_pos<other.byte_pos;
    }

    inline bool operator<=(vbyte_forward_iterator& other) const{
        return byte_pos<=other.byte_pos;
    }

    [[nodiscard]] inline off_t get_byte_pos() const{
        return byte_pos;
    }

    inline int_type operator*() const{
        return curr_val;
    }

    [[nodiscard]] inline off_t code_bytes() const{
        return len;
    }
};

template<typename int_type>
struct vbyte_reverse_iterator{

    off_t byte_pos;
    off_t stream_len;
    off_t len{};
    uint8_t const* stream = nullptr;
    int_type curr_val = std::numeric_limits<int_type>::max();

    vbyte_reverse_iterator(const uint8_t *_stream, off_t _stream_len, off_t _byte_pos) : byte_pos(_byte_pos),
                                                                                         stream_len(_stream_len),
                                                                                         stream(_stream) {

        assert(byte_pos<stream_len);
        while(byte_pos>=0 && stream[byte_pos]<128) byte_pos--;
        if(byte_pos>=0) {
            len = vbyte_decode_inverse(stream, byte_pos, curr_val);
        }
    }

    void operator--() {
        if(byte_pos>0){
            len = vbyte_decode_inverse(stream, byte_pos, curr_val);
        } else {
            byte_pos=-1;
        }
    }

    inline bool operator==(vbyte_reverse_iterator& other) const{
        return byte_pos==other.byte_pos;
    }

    inline bool operator!=(vbyte_reverse_iterator& other) const{
        return byte_pos!=other.byte_pos;
    }

    inline bool operator<(vbyte_reverse_iterator& other) const{
        return byte_pos<other.byte_pos;
    }

    inline bool operator<=(vbyte_reverse_iterator& other) const{
        return byte_pos<=other.byte_pos;
    }

    [[nodiscard]] inline off_t get_byte_pos() const{
        return byte_pos+1;
    }

    inline int_type operator*() const{
        return curr_val;
    }

    [[nodiscard]] inline off_t code_bytes() const{
        return len;
    }
};
*/

uint8_t vbyte_len(size_t value){
    if(value==0) return 1;
    uint8_t n_bits = sym_width(value);
    n_bits+=((n_bits+7)/8);
    return (n_bits+7)/8;
}

template<class sym_t>
struct vbyte_decoder{

    typedef sym_t sym_type;

    inline static off_t forward(const uint8_t * ptr, off_t n_pos) {
        off_t n_bytes=0;
        for(off_t i=0;i<n_pos;i++){
            while(ptr[n_bytes]<128) n_bytes++;
            n_bytes++;
        }
        return n_bytes;
    }

    inline static off_t backward(uint8_t * ptr, const uint8_t *& l_boundary, off_t n_pos) {

        if(ptr==l_boundary) return 0;

        uint8_t *tmp = ptr-1;
        //std::cout<<(int)*(ptr-3)<<" "<<(int)*(ptr-2)<<" "<<(int)*(ptr-1)<<" "<<(int)*ptr<<std::endl;
        while(tmp!=l_boundary && n_pos>0){
            tmp--;
            while(tmp!=l_boundary && *tmp<128){
                --tmp;
            }
            n_pos--;
        }
        tmp++;
        assert(tmp!=l_boundary);
        return ptr-tmp;
    }

    //rightmost byte of the rightmost valid value before ptr:
    // For instance, 0 0 0 1* 0 0* 0 1 (1 is the last byte of each symbol).
    // The function returns 3 with pos = 5 (assuming uint32_t symbols of sizeof(uint32_t)= 4 bytes)
    inline static off_t mov_to_prev_valid_byte(uint8_t *&ptr, __attribute__((unused)) off_t& pos) {
        uint8_t * tmp = ptr;
        while(*ptr<128) --ptr;
        return tmp-ptr;
    }

    inline static size_t write_forward(uint8_t * stream, sym_type& x){
        size_t len;
        size_t comp_x;
        size_t c=0, b=0;
        while(x>=128){
            c += (x & 127)<<b;
            b+=8;
            x>>=7;
        }
        comp_x = c+((x+128)<<b);
        len = (b+8)>>3UL;//in bytes
        memcpy(stream, &comp_x, len);
        return len;
    }

    static off_t read_forward(const uint8_t * ptr, sym_type& x){
        off_t i=0;
        uint8_t max_code_len = 8;
        while(max_code_len>0){
            if(ptr[i]>=128) break;
            max_code_len--;
            i++;
        }
        off_t len = i+1;

        switch (len) {
            case 1:
                x = ptr[i]-128;
                break;
            case 2:
                x = ptr[i--]-128;
                x = (x<<7) + ptr[i];
                break;
            case 3:
                x = ptr[i--]-128;
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i];
                break;
            case 4:
                x = ptr[i--]-128;
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i];
                break;
            case 5:
                x = ptr[i--]-128;
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i];
                break;
            case 6:
                x = ptr[i--]-128;
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i];
                break;
            case 7:
                x = ptr[i--]-128;
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i];
                break;
            case 8 :
                x = ptr[i--]-128;
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i--];
                x = (x<<7) + ptr[i];
                break;
            default:
                std::cerr<<"\nIncorrect vbyte length "<<len<<std::endl;
                exit(1);
        }
        return len;
    }

    //this assumes ptr points to the rightmost byte of a symbol in the stream
    //notice this function changes the memory address
    inline static off_t read_backwards(uint8_t *&ptr, const uint8_t *& l_boundary, sym_type& sym) {
        assert(ptr!=l_boundary);
        assert(*ptr>=128);

        uint8_t * tmp = ptr;
        sym = *ptr-128;
        ptr--;
        while(ptr!=l_boundary && *ptr<128){
            sym = (sym<<7) + *ptr;
            --ptr;
        }
        ptr++;
        return tmp-ptr;
    }

    //get the rightmost symbol that differs from the symbol in ptr
    //it scans prefix number of bytes preceding ptr to look for the mismatching symbol
    inline static off_t mov_to_prev_diff_sym(uint8_t *& ptr, const uint8_t *l_boundary, sym_type& sym) {
        //get the length in bytes of the current symbol
        uint8_t * tmp = ptr;
        while(tmp!=l_boundary && *tmp<128) tmp++;
        off_t bytes = tmp-ptr+1;
        sym = 0;
        off_t dist = ptr-l_boundary-1;//number of bytes we can visit

        off_t n_comp = dist/bytes;
        tmp = ptr;
        ptr-=bytes;
        while(n_comp>0 && memcmp(tmp, ptr, bytes)==0){
            ptr -= bytes;
            n_comp--;
        }

        ptr+=bytes-1;
        if(ptr!=l_boundary){
            //this extra loop is to cover a weird corner case.
            // E.g., 92 131, 131, 157.
            // In this case, the comparison of rightmost occurrence of 131 will finish in 92,
            // but the backward read has to start from the leftmost occurrence of 131
            while(*ptr<128) ptr++;
            assert(tmp-ptr>0);
            read_backwards(ptr, l_boundary, sym);
        }
        return tmp-ptr;
    }

    //get the leftmost symbol that differs from the symbol in ptr
    //it will stop checking once the scan reaches boundary address
    //the function assumes ptr and boundary are addresses within the same memory area
    inline static off_t mov_to_next_diff_sym(uint8_t *& ptr, const uint8_t* r_boundary, sym_type& sym) {

        //get the length in bytes of the current symbol
        uint8_t * tmp = ptr;
        while(*tmp<128) tmp++;
        off_t bytes = tmp-ptr+1;

        sym = 0;
        tmp = ptr;
        off_t dist = r_boundary-ptr-1;//number of bytes we can visit
        off_t n_comp = dist/bytes;
        ptr+=bytes;
        while(n_comp>0 && memcmp(tmp, ptr, bytes)==0){
            ptr+=bytes;
            n_comp--;
        }

        if(ptr!=r_boundary){
            read_forward(ptr, sym);
        }
        return ptr-tmp;
    }
};
#endif //SIMPLE_EXAMPLE_VBYTE_ENCODING_H
