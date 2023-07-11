//
// Created by diego on 27-08-20.
//

#ifndef LPG_COMPRESSOR_INT_ARRAY_H
#define LPG_COMPRESSOR_INT_ARRAY_H

#include <iostream>
#include <limits>
#include <cassert>
#include <utility>
#include "memory_handler.hpp"
#include "bitstream.h"

template<class word_t,
         uint8_t w_bits=std::numeric_limits<word_t>::digits>
struct int_array{

    static_assert(w_bits<=std::numeric_limits<word_t>::digits);
    typedef size_t size_type;
    typedef word_t value_type;
    //typedef typename std::conditional<w_bits==1, bool, typename std::conditional<w_bits<=8, uint8_t, typename std::conditional<w_bits<=16, uint16_t, typename std::conditional<w_bits<=32, uint32_t, uint64_t>::type>::type>::type>::type v_type;
    typedef bitstream<word_t, w_bits> stream_t;

    size_t m_cap = 0;
    size_t m_size = 0;
    uint8_t m_width = w_bits;
    stream_t bits;

    class proxy {
    private:
        int_array & m_arr;
        size_t m_idx;

    public:
        proxy(int_array & _arr, size_t _idx) : m_arr(_arr), m_idx(_idx) {
            assert(m_idx<m_arr.m_cap);
        }

        operator value_type() const {
            return m_arr.bits.read(m_idx * m_arr.m_width, (m_idx + 1) * m_arr.m_width - 1);
        }

        proxy& operator=(value_type value) {
            m_arr.write(m_idx, value);
            return *this;
        }

        proxy& operator=(const proxy& other) {
            m_arr.write(m_idx, (value_type)other);
            return *this;
        }

        proxy& operator&=(size_t val) {
            size_t elm = m_arr.bits.read(m_idx * m_arr.m_width, (m_idx + 1) * m_arr.m_width - 1);
            m_arr.write(m_idx, elm & val);
            return *this;
        }

        proxy& operator++() {
            auto val = (value_type)*this;
            m_arr.write(m_idx, val+1);
            return *this;
        }

        value_type operator++(int) {
            auto val = (value_type)*this;
            ++(*this);
            return val;
        }

        proxy& operator--() {
            auto val = (value_type)*this;
            m_arr.write(m_idx, val-1);
            return *this;
        }

        value_type operator--(int) {
            auto val = (value_type)*this;
            --(*this);
            return val;
        }

        proxy& operator+=(const value_type x) {
            auto val = (value_type)*this;
            m_arr.write(m_idx, val+x);
            return *this;
        }

        proxy& operator-=(const value_type x) {
            auto val = (value_type)*this;
            m_arr.write(m_idx, val-x);
            return *this;
        }

        bool operator==(const proxy& x)const {
            return value_type(*this) == value_type(x);
        }

        bool operator<(const proxy& x)const {
            return value_type(*this) < value_type(x);
        }
    };

    //simple constructor (cap is the capacity of the vector)
    int_array(size_t cap_, uint8_t width_) : m_cap(cap_), m_width(width_){
        assert(m_width <= w_bits);
        bits.stream_size = INT_CEIL(m_cap * width_, stream_t::word_bits);
        bits.stream = allocator::allocate<word_t>(bits.stream_size, true, 0);
    }

    //initialize from allocated memory
    int_array(word_t *ptr, size_t size_, uint8_t width_): m_cap(size_), m_size(size_), m_width(width_){
        assert(m_width <= w_bits);
        bits.stream_size = INT_CEIL(m_size*m_width, stream_t::word_bits);
        bits.stream = ptr;
    }

    //initialize with default value and size
    int_array(size_t size_, value_type def_val, uint8_t width_) : m_cap(size_), m_size(size_), m_width(width_){
        assert(m_width <= w_bits);
        bits.stream_size = INT_CEIL(m_size * width_, stream_t::word_bits);
        bits.stream = allocator::allocate<word_t>(bits.stream_size, true, 0);
        initialize(def_val, m_size);
    }

    //initializer-list constructor
    int_array(std::initializer_list<size_t> list) {
        assert(list.size()>0);
        size_t max=0;
        for(auto const& sym : list) if(sym>max) max=sym;

        m_width = std::max(4UL, sizeof(unsigned int)*8 - __builtin_clz(max));
        assert(m_width<=w_bits);
        m_size = list.size();
        m_cap = m_size;

        bits.stream_size = n_words();
        bits.stream = allocator::allocate<word_t>(bits.stream_size, false, 0);
        size_t i=0;
        for(auto const& val : list){
            write(i++,val);
        }
        mask_tail();
    }

    ~int_array() {
        allocator::deallocate<word_t>(bits.stream);
    }

    //default constructor
    int_array()=default;

    //move constructor
    explicit int_array(int_array<word_t>&& other) noexcept {
        bits.stream=nullptr;
        bits.stream_size=0;
        move(std::forward<int_array<word_t>>(other));
    };

    //copy constructor
    int_array(const int_array& other) noexcept {
        copy(other);
    };

    void move(int_array<word_t>&& other) noexcept {
        m_size = std::exchange(other.m_size, 0);
        m_cap = std::exchange(other.m_cap, 0);
        m_width = std::exchange(other.m_width, 0);
        bits.stream_size = std::exchange(other.bits.stream_size, 0);
        if(bits.stream!= nullptr){
            allocator::deallocate(bits.stream);
        }
        bits.stream = std::exchange(other.bits.stream, nullptr);
    }

    inline void copy(const int_array<word_t> &other) {
        m_size = other.m_size;
        m_cap = other.cap;
        m_width = other.m_width;
        bits.stream_size = n_words();

        if(bits.stream!= nullptr){
            allocator::deallocate(bits.stream);
        }

        bits.stream = allocator::allocate<word_t>(bits.stream_size, false, 0);
        memcpy(bits.stream, other.stream, bits.stream * sizeof(word_t));
    }

    inline void clear() {
        m_size=0;
    }

    inline void erase() {
        m_size=0;
        m_cap = 0;
        allocator::deallocate(bits.stream);
        bits.stream = nullptr;
        bits.stream_size = 0;
    }

    //copy assignment operator
    int_array& operator=(const int_array<word_t> & other){
        if(this!=&other){
            copy(other);
        }
        return *this;
    }

    //move assignment operator
    int_array& operator=(int_array<word_t> && other) noexcept{
        if(this!=&other){
            move(std::forward<int_array<word_t>>(other));
        }
        return *this;
    }

    inline bool operator==(const int_array<word_t> &other) const {
        if(m_size!=other.m_size) return false;
        size_t tot_bits = m_size * m_width;
        return bits.compare_chunk(other.bits.stream, 0, tot_bits);
    }

    [[nodiscard]] inline void* data() const{
        return bits.stream;
    }

    [[nodiscard]] inline word_t* stream() const{
        return bits.stream;
    }

    inline void set_data(const word_t* new_data, size_t size) {
        bits.stream = new_data;
        bits.stream_size = size;
    }

    [[nodiscard]] inline size_t width() const{
        return m_width;
    }

    [[nodiscard]] inline size_t n_words() const{
        return m_size==0? 0: INT_CEIL(n_bits(), stream_t::word_bits);
    }

    [[nodiscard]] inline size_t n_bytes() const{//number of bytes used by the data (ceil)
        return INT_CEIL(n_bits(), 8);
    }

    [[nodiscard]] inline size_t n_bits() const{//number of bits used by the data (exact)
        return m_size * m_width;
    }

    inline void mask_tail() {
        size_t tot_bits = m_size * m_width;
        size_t n_cells = tot_bits / stream_t::word_bits;//floor

        if (tot_bits > n_cells * stream_t::word_bits) {
            bits.write(tot_bits, (n_cells + 1) * stream_t::word_bits - 1, 0);
        }
    }

    inline void push_back(value_type value) {
        if((m_size+1)>m_cap){
            m_cap = m_cap==0? 2: m_cap*2;
            reserve(m_cap);
        }
        size_t i = m_size * m_width;
        size_t j = (m_size+1) * m_width - 1;
        bits.write(i, j , value);
        m_size++;
    }

    inline void pop_back() {
        if(m_size>=1) m_size--;
    }

    void reserve(size_t new_cap) {
        //reserve memory for new_size number of elements
        size_t new_buffer_size = 0;
        if(new_cap>0){
            new_buffer_size = INT_CEIL(new_cap * m_width, stream_t::word_bits);
        }

        if(new_buffer_size>bits.stream_size){
            bits.stream = allocator::reallocate<word_t>(bits.stream, bits.stream_size, new_buffer_size, false, 0);
            bits.stream_size = new_buffer_size;
        }
        m_cap = new_cap;
    }

    void resize(size_t new_size) {
        //reserve memory for new_size number of elements
        size_t new_buffer_size = INT_CEIL(new_size * m_width, stream_t::word_bits);
        bits.stream = allocator::reallocate<word_t>(bits.stream, bits.stream_size, new_buffer_size, false, 0);
        bits.stream_size = new_buffer_size;
        m_size = new_size;
        m_cap = m_size;
    }

    [[nodiscard]] inline size_t size() const{
        return m_size;
    }

    [[nodiscard]] inline size_t capacity() const{
        return m_cap;
    }

    [[nodiscard]] inline bool empty() const{
        return m_size==0;
    }

    inline value_type back() const{
        assert(m_size>0);
        return bits.read((m_size-1)*m_width, m_size*m_width-1);
    }

    inline void write(size_t idx, value_type value){
        assert(idx<m_cap);
        if(idx>=m_size) m_size = idx+1;
        bits.write(idx * m_width, (idx + 1) * m_width - 1, value);
    }

    inline value_type read(const size_t idx) const {
        assert(idx<m_size);
        return bits.read(idx * m_width, (idx + 1) * m_width - 1);
    }

    inline value_type operator[](const size_t idx) const {
        assert(idx<m_size);
        return bits.read(idx * m_width, (idx + 1) * m_width - 1);
    }

    inline proxy operator[](const size_t idx) {
        return proxy(*this, idx);
    }

    size_type serialize(std::ostream &out) const{
        size_t written_bytes = bits.serialize(out);
        written_bytes +=serialize_elm(out, m_size);
        written_bytes +=serialize_elm(out, m_cap);
        written_bytes +=serialize_elm(out, m_width);
        return written_bytes;
    }

    void load(std::istream &in){
        bits.load(in);
        load_elm(in, m_size);
        load_elm(in, m_cap);
        load_elm(in, m_width);
    }

    void initialize(value_type val, size_t n_elems){
        if(val==0){
            size_t n_cells = INT_CEIL(n_elems*m_width, stream_t::word_bits);
            assert(n_cells<=bits.stream_size);
            memset(bits.stream, 0, n_cells*sizeof(word_t));
        }else{
            std::cout<<"not implemented yet"<<std::endl;
            exit(1);
        }
    }

    void swap(int_array<word_t>& other){
        move(std::move(other));
    }

    void set_width(uint8_t new_width){
        assert(new_width<=w_bits);
        m_width = new_width;
    }
};

typedef int_array<size_t, 1> bit_array;

#endif //LPG_COMPRESSOR_INT_ARRAY_H
