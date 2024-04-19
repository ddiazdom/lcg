//
// Created by diego on 05-03-20.
//

#ifndef LMS_COMPRESSOR_SE_HASH_TABLE_H
#define LMS_COMPRESSOR_SE_HASH_TABLE_H
#include <iostream>
#include <cstring>
#include <cassert>
#include <fstream>
#include <utility>
#include "xxHash-dev/xxhash.h"
#include "bitstream.h"
#include "mxm3.h"
#include <algorithm>
#include <type_traits>
#include <utility>
#include <vector>
#ifdef __linux__
#include <malloc.h>
#endif

//#ifdef __ARM_NEON
//#include <arm_neon.h>
//#include <intrin.h>
//#endif
//#include <arm_neon.h>

#define HASH(stream, bytes) XXH3_64bits(stream, bytes);
//#define HASH(stream, bytes) mx3::hash((const uint8_t*)stream, bytes, 0);

template<class ht_type>
class hash_table_iterator{

    using buffer_t =                    typename ht_type::buff_t;
    using value_t =                     typename ht_type::val_type;

    struct data_t{
        void * pointer=nullptr;
        size_t n_bytes{};
        value_t value{};
    };

    typedef std::pair<size_t, size_t>   locus_t;

    typedef std::forward_iterator_tag   iterator_category;
    typedef const locus_t&              reference;

    const ht_type&                      ht;
    locus_t                             locus;
    data_t                              m_data;

    const data_t& data =                m_data;

public:

    hash_table_iterator(const ht_type& ht_, size_t bit_pos): ht(ht_) {
        if(bit_pos+1<ht.next_av_bit){
            size_t n_bits = ht.data.read(bit_pos, bit_pos+ht.description_bits()-1);
            locus.first = bit_pos+ht.description_bits();
            locus.second = locus.first+n_bits;
        }else{
            locus.first = ht.next_av_bit+ht.description_bits()-1;
        }
        m_data.pointer = nullptr;
    };

    inline size_t operator*() {
        return locus.first-ht.description_bits();
    }

    const data_t& pair(){
        if(locus.first<ht.next_av_bit){
            if(m_data.pointer==nullptr){
                m_data.pointer = malloc(INT_CEIL(ht.max_key_bits, ht_type::stream_t::word_bits)*sizeof(buffer_t));
            }
            ht.data.read_chunk(m_data.pointer, locus.first, locus.second-1);
            m_data.n_bytes = INT_CEIL(locus.second-locus.first, 8);
            ht.data.read_chunk(&m_data.value, locus.second, locus.second+ht.value_bits()-1);
        }
        return data;
    }

    const value_t& value(){
        if(locus.first<ht.next_av_bit){

            //clean the area
            char * tmp = reinterpret_cast<char *>(&m_data.value);
            memset(tmp, 0, sizeof(m_data.value));
            //

            //!!!TODO always clean the unused area if we use less bytes than sizeof(value_type)
            ht.data.read_chunk(&m_data.value, locus.second, locus.second+ht.value_bits()-1);
        }
        return m_data.value;
    }


    hash_table_iterator<ht_type> operator++(int) {
        hash_table_iterator<ht_type> tmp(*this);
        ++*this;
        return tmp;
    }

    inline bool operator==(const hash_table_iterator<ht_type>& other) const{
        return locus.first == other.locus.first;
    }

    inline bool operator!=(const hash_table_iterator<ht_type>& other) const{
        return locus.first != other.locus.first;
    }

    hash_table_iterator<ht_type>& operator=(const value_t &val) { //update the value associated with the current key
        if(locus.first<ht.next_av_bit){
            const void * tmp_ptr = &val;
            ht.data.write_chunk(tmp_ptr, locus.second, locus.second+ht.value_bits()-1);
        }
        return *this;
    }

    const hash_table_iterator<ht_type>& operator++() {
        if(locus.first<ht.next_av_bit){
            locus.first = locus.second+ht.value_bits()+ht.description_bits();
            if(locus.first<ht.next_av_bit){
                size_t key_bits = ht.data.read(locus.first-ht.description_bits(), locus.first-1);
                locus.second = locus.first+key_bits;
            }
        }
        return *this;
    }

    ~hash_table_iterator(){
        if(m_data.pointer!= nullptr){
            free(m_data.pointer);
        }
    }
};

template<class value_t,
        size_t val_bits=sizeof(value_t)*8,
        size_t desc_bits=32>
class hash_table{

private:

    typedef hash_table<value_t, val_bits, desc_bits> ht_type;
    typedef size_t                                   buffer_t;
    typedef bitstream<size_t>                        stream_t;
    typedef hash_table_iterator<ht_type>             iterator;
    friend                                           iterator;

    //limit bucket distance allowed for the data in the hash table
    const size_t limit_bck_dist = 0xFFFF;

    //64-bit pointer for the
    // hash table
    size_t* table = nullptr;

    // position of the next
    // available bit in the data buffer (1-based)
    size_t next_av_bit=0;

    //number of buckets in the
    // hash table always a power of two
    size_t n_buckets=0;

    //maximum bucket distance seen so far
    // for key (for the robin hood algorithm)
    size_t max_bck_dist=0;

    //number of elements
    // inserted in the hash table
    size_t n_elms=0;

    //current load factor
    float m_load_factor=0;

    //maximum allowed load factor
    float m_max_load_factor=0;

    //number of elements of the longest
    // string inserted so far in the hash table
    size_t max_key_bits=0;

    //stream of bits with the data
    stream_t data;

    //number of bits used in the key description
    size_t d_bits=0;

    //the variable query is a pointer containing the queried string
    inline bool equal(const void *query, size_t query_bits, size_t data_offset) const{
        size_t key_bits = data.read(data_offset, data_offset+d_bits-1);
        if(key_bits!=query_bits) return false;
        return data.compare_chunk(query, data_offset+d_bits, query_bits);
    }

    inline size_t insert_int(const void* key, size_t key_bits, const value_t& value){

        size_t pair_bits = key_bits+val_bits+d_bits;
        size_t idx = next_av_bit - 1 + pair_bits;

        if(idx >= data.n_bits()){
            resize_data_buffer(pair_bits);
        }

        idx = next_av_bit-1;

        //write metadata (length, alphabet width) of the key
        data.write(idx, idx+d_bits-1, key_bits);
        idx+=d_bits;

        //write the key
        data.write_chunk(key, idx, idx+key_bits-1);
        idx+=key_bits;

        //write the value
        if constexpr (val_bits<=64){
            data.write( idx, idx+val_bits-1, value);
        }else{
            data.write_chunk(&value, idx, idx+val_bits-1);
        }

        if(key_bits>max_key_bits) max_key_bits = key_bits;

        next_av_bit += pair_bits;
        assert((next_av_bit-pair_bits)>0);


        return next_av_bit-pair_bits;
    }

    inline void rehash(){

        max_bck_dist=0;
        size_t data_offset=0;
        void * tmp_key = malloc(INT_CEIL(max_key_bits, stream_t::word_bits)*sizeof(buffer_t));

        size_t dist, bck_dist, bck_offset, tmp_offset, idx, hash;

        for(size_t k=0;k<n_elms;k++){

            size_t key_bits = data.read(data_offset, data_offset+d_bits-1);
            size_t key_bytes = INT_CEIL(key_bits, 8);

            //this clean the tail of the buffer
            char * tmp = reinterpret_cast<char*>(tmp_key);
            tmp[key_bytes-1] = 0;
            //

            data.read_chunk(tmp_key, data_offset + d_bits, data_offset + d_bits + key_bits - 1);

            //hash = XXH3_64bits(tmp_key, key_bytes);
            hash = HASH(tmp_key, key_bytes);

            idx = hash & (n_buckets - 1);
            tmp_offset = data_offset + 1;

            if(table[idx]==0){
                table[idx] = tmp_offset;
            }else{
                dist = 0;
                while(table[idx]!=0){
                    bck_offset = (table[idx] & 0xFFFFFFFFFFFul);
                    bck_dist = table[idx] >> 44UL;

                    if(bck_dist<dist){ //steal to the rich
                        table[idx] = (dist<<44UL) | tmp_offset;
                        if(dist>max_bck_dist) max_bck_dist = dist;
                        tmp_offset = bck_offset;
                        dist = bck_dist+1;
                    }else{
                        dist++;
                    }
                    idx = (idx+1) & (n_buckets - 1);
                }

                table[idx] = (dist<<44UL) | tmp_offset ;
                if(dist>max_bck_dist) max_bck_dist = dist;
            }

            data_offset+=d_bits+key_bits+val_bits;
        }
        free(tmp_key);
    }

    [[nodiscard]] inline size_t new_table_size() const{
        //TODO if the hash table becomes too high, then I will change the growth policy to this:
        //size_t new_size = prime_generator::get_next_prime(size_t(double(n_buckets)/(m_max_load_factor-0.1)));
        //std::cout<<"Just testing "<<n_buckets<<" -> "<<new_size<<" -- "<<double(n_elms)/new_size<<", "<<(n_buckets<<1UL)<<std::endl;
        return n_buckets<<1UL;
    }

    //increase the data buffer to insert bits_to_fit bits
    void resize_data_buffer(size_t bits_to_fit) {

        //bytes we need to fit in the data buffer
        size_t bytes_to_fit = INT_CEIL(bits_to_fit, 8);

        //number of bytes used by the data buffer
        size_t used_data_bytes = INT_CEIL(next_av_bit, 8);

        //minimum number of bytes we require for inserting the data
        size_t min_data_bytes = used_data_bytes+bytes_to_fit;

        auto new_size = size_t(data.stream_size*sizeof(buffer_t)*1.5);

        //maximum number of bytes we can allocate
        size_t new_data_bytes = INT_CEIL(std::max(min_data_bytes, new_size), sizeof(buffer_t))*sizeof(buffer_t);

        data.stream = reinterpret_cast<buffer_t*>(realloc(data.stream, new_data_bytes));
        data.stream_size = new_data_bytes/sizeof(buffer_t);
    }

    void move(hash_table&& other){
        std::swap(n_buckets, other.n_buckets);
        std::swap(table, other.table);
        std::swap(next_av_bit, other.next_av_bit);
        std::swap(max_bck_dist, other.max_bck_dist);
        std::swap(n_elms, other.n_elms);
        std::swap(m_load_factor, other.m_load_factor);
        std::swap(m_max_load_factor, other.m_max_load_factor);
        std::swap(max_key_bits, other.max_key_bits);
        data.swap(other.data);
        std::swap(d_bits, other.d_bits);
    }

    void copy(hash_table& other){

        n_buckets = other.n_buckets;

        if(table==nullptr){
            table = (size_t*)malloc(sizeof(size_t)*n_buckets);
        }else{
            table = reinterpret_cast<size_t*>(realloc(table, n_buckets));
        }

        //TODO you missed copying the data buffer
        //

        next_av_bit = other.next_av_bit;
        max_bck_dist = other.max_bck_dist;
        n_elms = other.n_elms;
        m_load_factor = other.m_load_factor;
        m_max_load_factor = other.m_max_load_factor;
        max_key_bits = other.max_key_bits;
        d_bits = other.d_bits;
        data = other.data;
    }

    void init() {

        //size for the hash table
        n_buckets = 4;

        //number of bytes allocated for the hash table
        size_t table_bytes = (n_buckets * sizeof(size_t));
        table = reinterpret_cast<size_t*>(malloc(table_bytes));
        memset(table, 0, table_bytes);

        //number of bytes allocated for the data buffer
        data.stream_size = 8;
        size_t data_bytes = data.stream_size*sizeof(buffer_t);
        data.stream = reinterpret_cast<buffer_t *>(malloc(data_bytes));
        memset(data.stream, 0, data_bytes);

        next_av_bit = 1;
        m_load_factor = 0;
        max_bck_dist = 0;
        n_elms = 0;
        max_key_bits = 0;
    }

public:

    typedef value_t       val_type;
    typedef buffer_t      buff_t;

    explicit hash_table(float max_load_factor=0.8, size_t _d_bits=desc_bits){//size of the buffer
        d_bits = _d_bits;
        m_max_load_factor = max_load_factor;
        init();
    }

    hash_table(hash_table&& other) noexcept{
        move(std::forward<hash_table>(other));
    }

    //copy assignment operator
    hash_table& operator=(hash_table & other){
        if(this!=&other){
            copy(other);
        }
        return *this;
    }

    //move assignment operator
    hash_table& operator=(hash_table && other) noexcept{
        if(this!=&other){
            move(std::forward<hash_table>(other));
        }
        return *this;
    }

    //swap method
    void swap(hash_table& other){
        move(std::forward<hash_table>(other));
    }

    inline iterator begin()  {
        return iterator(*this, 0);
    }

    inline iterator end() const  {
        return iterator(*this, next_av_bit);
    }

    [[nodiscard]] inline size_t value_bits() const{
        return val_bits;
    }

    [[nodiscard]] inline size_t description_bits() const{
        return d_bits;
    }

    size_t serialize(std::ostream &out) const{
        size_t written_bytes = serialize_raw_vector(out, table, n_buckets);
        written_bytes += serialize_elm(out, next_av_bit);
        written_bytes += serialize_elm(out, max_bck_dist);
        written_bytes += serialize_elm(out, n_elms);
        written_bytes += serialize_elm(out, m_load_factor);
        written_bytes += serialize_elm(out, m_max_load_factor);
        written_bytes += serialize_elm(out, max_key_bits);
        data.serialize(out);
        written_bytes += serialize_elm(out, d_bits);
        return written_bytes;
    }

    void load(std::istream &in){
        load_raw_vector(in, table, n_buckets);
        load_elm(in, next_av_bit);
        load_elm(in, max_bck_dist);
        load_elm(in, n_elms);
        load_elm(in, m_load_factor);
        load_elm(in, m_max_load_factor);
        load_elm(in, max_key_bits);
        data.load(in);
        load_elm(in, d_bits);
    }

    //returns the update value in the hash table
    template<class val_t = value_t>
    typename std::enable_if<std::is_unsigned<val_t>::value, size_t>::type
    increment_value(const void* key, const size_t key_bits, value_t new_value) {
        auto res = insert(key, key_bits, new_value);
        if(res.second) return new_value;

        value_t old_value = data.read(res.first+d_bits+key_bits, res.first+d_bits+key_bits+val_bits-1);
        new_value += old_value;
        data.write(res.first+d_bits+key_bits, res.first+d_bits+key_bits+val_bits-1, new_value);
        return new_value;
    }

    void resize_table(size_t new_n_buckets) {
        assert(new_n_buckets>0 && (new_n_buckets & (new_n_buckets-1))==0);
        size_t new_table_bytes = new_n_buckets*sizeof(size_t);

        table = reinterpret_cast<size_t*>(realloc(table, new_table_bytes));
        n_buckets = new_n_buckets;
        memset(table, 0, new_n_buckets*sizeof(size_t));

        if(n_elms>0){
            rehash();
            m_load_factor = float(n_elms) / n_buckets;
        }
    };

    std::pair<size_t, bool> insert(const void* key, const size_t& key_bits, const value_t& val){

        size_t hash = HASH(key, INT_CEIL(key_bits, 8));

        size_t idx = hash & (n_buckets - 1);
        size_t in_offset=0;//locus where the key is inserted

        if(table[idx]==0) {
            in_offset  = insert_int(key, key_bits, val);
            table[idx] =  in_offset;
        }else{//resolve collision

            size_t dist=0, offset;
            size_t bck_dist, bck_offset;
            bool inserted = false;

            while(table[idx]!=0){

                bck_offset = (table[idx] & 0xFFFFFFFFFFFul);
                bck_dist = table[idx] >> 44UL;

                if(!inserted && bck_dist==dist && equal(key, key_bits , bck_offset-1)){
                    return {bck_offset-1, false};
                }else if(bck_dist<dist){ //steal to the rich

                    if(!inserted){//swap the keys
                        in_offset = insert_int(key, key_bits, val);
                        offset = in_offset;
                        inserted = true;
                    }

                    assert(dist<=limit_bck_dist);
                    if(dist>max_bck_dist) max_bck_dist = dist;

                    table[idx] = (dist<<44UL) | offset;
                    offset = bck_offset;
                    dist = bck_dist;
                }
                dist++;
                idx = (idx+1) & (n_buckets - 1);
            }

            if(!inserted){
                in_offset = insert_int(key, key_bits, val);
                offset = in_offset;
            }

            assert(dist<=limit_bck_dist);
            if(dist>max_bck_dist) max_bck_dist = dist;

            table[idx] = (dist<<44UL) | offset ;
        }

        n_elms++;
        m_load_factor = float(n_elms) / n_buckets;

        if(m_load_factor>=m_max_load_factor){
            resize_table(new_table_size());
        }

        assert(in_offset>0);
        return {in_offset-1, true};
    }

    [[nodiscard]] inline float load_factor() const{
        return m_load_factor;
    }

    [[nodiscard]] inline float max_load_factor() const{
        return m_max_load_factor;
    }

    [[nodiscard]] inline size_t size() const {
        return n_elms;
    }

    inline void shrink_databuff() {
        //next_av_bit is zero-based
        size_t written_bits = next_av_bit-1;
        size_t new_size = INT_CEIL(written_bits, stream_t::word_bits);

        if(data.stream_size!=new_size) {
            data.stream = reinterpret_cast<buffer_t *>(realloc(data.stream, new_size * sizeof(buffer_t)));
            data.stream_size = new_size;
            //fill the tail to avoid valgrind warnings
            size_t buffer_bits = new_size*stream_t::word_bits;
            assert(buffer_bits>=written_bits);
            if(buffer_bits>written_bits){//clean the buffer tail
                data.write(written_bits, buffer_bits-1, 0);
            }
        }
    }

    inline void reset(){
        n_buckets = 4;
        n_elms = 0;
        max_bck_dist = 0;
        m_load_factor = 0;
        next_av_bit = 1;

        size_t table_bytes = n_buckets*sizeof(size_t);
        table = reinterpret_cast<size_t*>(realloc(table, table_bytes));

        size_t data_bytes = table_bytes*2;
        data.stream = reinterpret_cast<buffer_t*>(realloc(data.stream, data_bytes));
        data.stream_size = data_bytes/sizeof(buffer_t);

        memset(table, 0, table_bytes);
        memset(data.stream, 0, data_bytes);
    }

    [[nodiscard]] inline size_t max_bucket_dist() const {
        return max_bck_dist;
    }

    inline std::pair<size_t, bool> key2value_test(const void* key, size_t key_bits) const {

        size_t hash = HASH(key, INT_CEIL(key_bits, 8));
        size_t idx = hash & (n_buckets - 1);

        size_t offset, bck_dist;

        bck_dist = (table[idx] >> 44U);
        if(bck_dist<1 || bck_dist==max_bck_dist) return {0, false};
        offset = (table[idx] & 0xFFFFFFFFFFFul);
        if(1==bck_dist && equal(key, key_bits, offset-1)) return {offset - 1, true};
        idx = (idx+1) & (n_buckets - 1);

        bck_dist = (table[idx] >> 44U);
        if(bck_dist<2 || bck_dist==max_bck_dist) return {0, false};
        offset = (table[idx] & 0xFFFFFFFFFFFul);
        if(2==bck_dist && equal(key, key_bits, offset-1)) return {offset - 1, true};
        idx = (idx+1) & (n_buckets - 1);

        bck_dist = (table[idx] >> 44U);
        if(bck_dist<3 || bck_dist==max_bck_dist) return {0, false};
        offset = (table[idx] & 0xFFFFFFFFFFFul);
        if(3==bck_dist && equal(key, key_bits, offset-1)) return {offset - 1, true};
        idx = (idx+1) & (n_buckets - 1);

        bck_dist = (table[idx] >> 44U);
        if(bck_dist<4 || bck_dist==max_bck_dist) return {0, false};
        offset = (table[idx] & 0xFFFFFFFFFFFul);
        if(4==bck_dist && equal(key, key_bits, offset-1)) return {offset - 1, true};
        idx = (idx+1) & (n_buckets - 1);

        size_t i = 6;
        while(i<=max_bck_dist && table[idx]!=0){
            bck_dist = (table[idx] >> 44U);
            if(bck_dist<i){
                return {0, false};
            }
            offset = (table[idx] & 0xFFFFFFFFFFFul);
            if(i==bck_dist && equal(key, key_bits, offset-1)){
                return {offset - 1, true};
            }
            idx = (idx+1) & (n_buckets - 1);
            i++;
        }
        return {0, false};
    }

    inline std::pair<size_t, bool> find(const void* key, size_t key_bits) const {

        size_t hash = HASH(key, INT_CEIL(key_bits, 8));
        size_t idx = hash & (n_buckets - 1);

        if(table[idx]==0){
            return {0, false};
        }else{
            size_t offset, i=0, bck_dist;
            while(i<=max_bck_dist && table[idx]!=0){
                bck_dist = (table[idx] >> 44U);
                if(bck_dist<i){
                    return {0, false};
                }
                offset = (table[idx] & 0xFFFFFFFFFFFul);
                if(i==bck_dist && equal(key, key_bits, offset-1)){
                    return {offset - 1, true};
                }
                idx = (idx+1) & (n_buckets - 1);
                i++;
            }
            return {0, false};
        }
    }

    inline bool key2value(const void* key, size_t key_bits, value_t& value) const {

        size_t hash = HASH(key, INT_CEIL(key_bits, 8));

        size_t idx = hash & (n_buckets - 1);
        size_t i=0;
        size_t bck_dist = (table[idx] >> 44U);
        size_t offset = (table[idx] & 0xFFFFFFFFFFFul);
        while(table[idx]!=0 && i<=bck_dist){
            if(i==bck_dist && equal(key, key_bits, offset-1)){
                size_t value_start = offset-1+d_bits+key_bits;
                if constexpr (val_bits<=64){
                    value = data.read(value_start, value_start+val_bits-1);
                }else{
                    data.read_chunk(value, value_start, value_start+val_bits-1);
                }
                return true;
            }
            idx = (idx+1) & (n_buckets - 1);
            i++;
            bck_dist = (table[idx] >> 44U);
            offset = (table[idx] & 0xFFFFFFFFFFFul);
        }
        return false;
    }

    //offset is the bit where the pair description starts
    inline void insert_value_at(size_t offset, value_t val) {
        size_t value_start = offset + d_bits + data.read(offset, offset + d_bits - 1);
        if constexpr (val_bits<=64){
            data.write(value_start, value_start + val_bits - 1, val);
        }else{
            data.write_chunk(&val, value_start, value_start + val_bits - 1);
        }
    }

    //offset is the bit where the pair description starts
    inline void get_value_from(size_t offset, value_t& dest) const {
        size_t value_start = offset + d_bits + data.read(offset, offset + d_bits - 1);
        if constexpr (val_bits<=64){
            dest = data.read(value_start, value_start + val_bits - 1);
        }else{
            data.read_chunk(&dest, value_start, value_start + val_bits - 1);
        }
    }

    //offset is the bit where the pair description starts
    inline size_t get_key_from(size_t offset, void* dest) const {
        size_t key_bits = data.read(offset, offset + d_bits - 1);;
        data.read_chunk(dest, offset+d_bits, offset+d_bits+key_bits-1);
        return key_bits;
    }

    ~hash_table(){
        if(table!= nullptr){
            free(table);
        }
        if(data.stream!= nullptr){
            free(data.stream);
        }
    }

    [[nodiscard]] inline size_t tot_buckets() const{
        return n_buckets;
    }

    inline const stream_t& get_data(){
        return data;
    }

    inline void clean(){
        memset(table, 0, n_buckets*sizeof(size_t));
        memset(data.stream, 0, data.stream_size*sizeof(buffer_t));
        n_elms = 0;
        max_bck_dist = 0;
        m_load_factor = 0;
        next_av_bit = 1;
    }

    void store_data_to_file(const std::string& output){
        std::filebuf fb;
        fb.open(output, std::ios::out | std::ios::binary);
        std::ostream ofs(&fb);
        ofs.write(reinterpret_cast<char *>(table), n_buckets*sizeof(size_t));
        ofs.tellp();
        shrink_databuff();
        data.serialize(ofs);
        fb.close();
    }

    size_t hash_table_tot_bytes(){
        size_t tot_bytes = n_buckets*sizeof(size_t);
        tot_bytes+= sizeof(data.stream_size);
        tot_bytes+= data.stream_size * sizeof(buffer_t);
        return  tot_bytes;
    }

    void ht_stats(size_t top){
        std::vector<size_t> freqs(max_bck_dist+1, 0);
        size_t dist;
        for(size_t i=0;i<n_buckets;i++){
            if(table[i]!=0){
                dist = table[i]>>44UL;
                assert(dist<=max_bck_dist);
                freqs[dist]++;
            }
        }

        std::vector<std::pair<size_t, size_t>> pairs(freqs.size());
        for(size_t i=0;i<pairs.size();i++){
            pairs[i] = {i, freqs[i]};
        }
        std::sort(pairs.begin(), pairs.end(), [&](auto &a, auto &b){
            return a.second>b.second;
        });

        std::cout<<"\nsize "<<size()<<std::endl;
        std::cout<<"n_buckets "<<tot_buckets()<<std::endl;
        std::cout<<"load factor: "<<load_factor()<<std::endl;
        std::cout<<"max_bck_dist: "<<max_bck_dist<<std::endl;
        for(size_t i=0;i<=std::min(top, max_bck_dist);i++){
            std::cout<<"Bck dist:"<<pairs[i].first<<" freq:"<<pairs[i].second<<" ( "<<(double(pairs[i].second)/n_elms)*100<<"% )"<<std::endl;
        }
    }

    void load_data_from_file(const std::string& input){
        assert(table==nullptr && data.stream== nullptr);
        std::ifstream ifs(input, std::ios_base::binary);
        table = reinterpret_cast<size_t*>(malloc(n_buckets*sizeof(size_t)));
        ifs.read(reinterpret_cast<char *>(table), n_buckets*sizeof(size_t));
        data.load(ifs);
        ifs.close();
    }

    [[nodiscard]] size_t data_bytes() const {
        return (data.stream_size*stream_t::word_bits)/8;
    }

    void destroy_table() {
        free(table);
        table = nullptr;
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    void destroy_data() {
        free(data.stream);
        data.stream = nullptr;
#ifdef __linux__
        malloc_trim(0);
#endif
    }
};

template<class value_t,
         size_t val_bits=sizeof(value_t)*8,
         size_t desc_bits=32>
class buffered_hash_table{

private:

    typedef buffered_hash_table<value_t, val_bits, desc_bits> ht_type;
    typedef size_t                                            buffer_t;
    typedef bitstream<size_t>                                 stream_t;
    typedef hash_table_iterator<ht_type>                      iterator;
    friend                                                    iterator;

    //limit bucket distance allowed for the data in the hash table
    const size_t limit_bck_dist = 0xFFFF;

    bool static_buffer = false;

    //64-bit pointer for the
    // hash table
    size_t* table = nullptr;

    //maximum buffer size in bytes
    // used for the data structure
    size_t max_buffer_bytes=0;

    // position of the next
    // available bit in the data buffer (1-based)
    size_t next_av_bit=0;

    //number of buckets in the
    // hash table always a power of two
    size_t n_buckets=0;

    //maximum bucket distance seen so far
    // for key (for the robin hood algorithm)
    size_t max_bck_dist=0;

    //number of elements
    // inserted in the hash table
    size_t n_elms=0;

    //current load factor
    float m_load_factor=0;

    //maximum allowed load factor
    float m_max_load_factor=0;

    //number of elements of the longest
    // string inserted so far in the hash table
    size_t max_key_bits=0;

    //file where the data is dumped
    std::ofstream dump_fs;
    std::string file;

    //stream of bits with the data
    stream_t data;

    //number of bits used in the key description
    size_t d_bits=0;

    //maximum number of bits you can fit in the data buffer
    size_t buff_limit{};

    //the variable query is a pointer containing the queried string
    inline bool equal(const void *query, size_t query_bits, size_t data_offset) const{
        size_t key_bits = data.read(data_offset, data_offset+d_bits-1);
        if(key_bits!=query_bits) return false;
        return data.compare_chunk(query, data_offset+d_bits, query_bits);
    }

    inline std::pair<size_t, bool> insert_int(const void* key, size_t key_bits, const value_t& value){

        size_t pair_bits = key_bits+val_bits+d_bits;
        assert(pair_bits<buff_limit);

        size_t idx = next_av_bit - 1 + pair_bits;
        bool dump=false;//flag to indicate that the data was dumped
        if(idx >= data.n_bits()){
            if(static_buffer){
                dump = true;
                dump_hash();
            }else{
                //we exceed the size of the data buffer
                dump = resize_data_buffer(pair_bits);
            }
        }

        idx = next_av_bit-1;

        //write metadata (length, alphabet width) of the key
        data.write(idx, idx+d_bits-1, key_bits);
        idx+=d_bits;

        //write the key
        data.write_chunk(key, idx, idx+key_bits-1);
        idx+=key_bits;

        //write the value
        if constexpr (val_bits<=64){
            data.write( idx, idx+val_bits-1, value);
        }else{
            data.write_chunk(&value, idx, idx+val_bits-1);
        }

        if(key_bits>max_key_bits) max_key_bits = key_bits;

        next_av_bit += pair_bits;
        assert((next_av_bit-pair_bits)>0);

        return {next_av_bit-pair_bits, dump};
    }

    inline void rehash(){

        assert(!static_buffer);

        max_bck_dist=0;
        size_t data_offset=0;
        void * tmp_key = malloc(INT_CEIL(max_key_bits, stream_t::word_bits)*sizeof(buffer_t));

        size_t dist, bck_dist, bck_offset, tmp_offset, idx, hash;

        for(size_t k=0;k<n_elms;k++){

            size_t key_bits = data.read(data_offset, data_offset+d_bits-1);
            size_t key_bytes = INT_CEIL(key_bits, 8);

            //this clean the tail of the buffer
            char * tmp = reinterpret_cast<char*>(tmp_key);
            tmp[key_bytes-1] = 0;
            //

            data.read_chunk(tmp_key, data_offset + d_bits, data_offset + d_bits + key_bits - 1);

            //hash = XXH3_64bits(tmp_key, key_bytes);
            hash = HASH(tmp_key, key_bytes);

            idx = hash & (n_buckets - 1);
            tmp_offset = data_offset + 1;

            if(table[idx]==0){
                table[idx] = tmp_offset;
            }else{
                dist = 0;
                while(table[idx]!=0){
                    bck_offset = (table[idx] & 0xFFFFFFFFFFFul);
                    bck_dist = table[idx] >> 44UL;

                    if(bck_dist<dist){ //steal to the rich
                        table[idx] = (dist<<44UL) | tmp_offset;
                        if(dist>max_bck_dist) max_bck_dist = dist;
                        tmp_offset = bck_offset;
                        dist = bck_dist+1;
                    }else{
                        dist++;
                    }
                    idx = (idx+1) & (n_buckets - 1);
                }

                table[idx] = (dist<<44UL) | tmp_offset ;
                if(dist>max_bck_dist) max_bck_dist = dist;
            }

            data_offset+=d_bits+key_bits+val_bits;
        }
        free(tmp_key);
    }

    inline bool exceed_buffer(size_t key_bits) const {

        size_t pair_bits = key_bits+val_bits+d_bits;
        size_t idx = next_av_bit - 1 + pair_bits;
        size_t table_bytes = n_buckets*sizeof(size_t);

        if((float(n_elms+1)/n_buckets) > m_max_load_factor) {
            if (static_buffer) return true;
            size_t new_n_buckets = new_table_size();
            table_bytes = new_n_buckets*sizeof(size_t);
        }

        if(idx>=data.n_bits() && static_buffer) return true;

        //number of bytes used by the data buffer, including the new element
        size_t data_bytes = INT_CEIL(idx, stream_t::word_bits)*sizeof(buffer_t);

        return max_buffer_bytes<(data_bytes+table_bytes);
    }

    inline size_t new_table_size() const{
        //TODO if the hash table becomes too high, then I will change the growth policy to this:
        //size_t new_size = prime_generator::get_next_prime(size_t(double(n_buckets)/(m_max_load_factor-0.1)));
        //std::cout<<"Just testing "<<n_buckets<<" -> "<<new_size<<" -- "<<double(n_elms)/new_size<<", "<<(n_buckets<<1UL)<<std::endl;
        return n_buckets<<1UL;
    }

    //number of bytes available for the hash table
    inline size_t av_tb_bytes() const{
        size_t data_bytes = data.stream_size*sizeof(buffer_t);
        size_t av_bytes = max_buffer_bytes-data_bytes;
        return av_bytes;
    }

    //increase the data buffer to insert bits_to_fit bits
    bool resize_data_buffer(size_t bits_to_fit) {
        assert(!static_buffer);
        bool data_dumped = false;
        size_t bytes_to_fit = INT_CEIL(bits_to_fit, stream_t::word_bits)*sizeof(buffer_t);

        //number of bytes allocated for the hash table
        size_t table_bytes = n_buckets * sizeof(size_t);

        //number of bytes used by the data buffer
        size_t used_data_bytes = INT_CEIL(next_av_bit, stream_t::word_bits)*sizeof(buffer_t);

        //number of available bytes
        size_t av_bytes = max_buffer_bytes-(used_data_bytes + table_bytes);

        if(bytes_to_fit>av_bytes){
            if(bytes_to_fit>(max_buffer_bytes/2)){
                std::cout<<"One of the keys requires more than half the buffer .. consider increasing the buffer size"<<std::endl;
            }
            data_dumped = true;
            dump_hash();//dump the hash table to disk
        }else{
            //minimum number of bytes we require for inserting the data
            size_t min_data_bytes = used_data_bytes+bytes_to_fit;

            auto new_size = size_t(data.stream_size*sizeof(buffer_t)*1.5);
            //maximum number of bytes we can allocate
            size_t new_data_bytes = std::min(std::max(min_data_bytes, new_size), max_buffer_bytes-table_bytes);
            data.stream = reinterpret_cast<buffer_t*>(realloc(data.stream, new_data_bytes));
            data.stream_size = new_data_bytes/sizeof(buffer_t);
        }
        return data_dumped;
    }

    void move(buffered_hash_table&& other){
        assert(!static_buffer);
        std::swap(n_buckets, other.n_buckets);
        std::swap(table, other.table);
        std::swap(max_buffer_bytes, other.max_buffer_bytes);
        std::swap(next_av_bit, other.next_av_bit);
        std::swap(max_bck_dist, other.max_bck_dist);
        std::swap(n_elms, other.n_elms);
        std::swap(m_load_factor, other.m_load_factor);
        std::swap(m_max_load_factor, other.m_max_load_factor);
        std::swap(max_key_bits, other.max_key_bits);
        std::swap(dump_fs, other.dump_fs);
        std::swap(file, other.file);
        data.swap(other.data);
        std::swap(d_bits, other.d_bits);
        std::swap(static_buffer, other.static_buffer);
    }

    void copy(buffered_hash_table& other){
        assert(!static_buffer);
        n_buckets = other.n_buckets;
        if(table==nullptr){
            table = (size_t*)malloc(sizeof(size_t)*n_buckets);
        }else{
            table = reinterpret_cast<size_t*>(realloc(table, n_buckets));
        }
        max_buffer_bytes = other.max_buffer_bytes;
        next_av_bit = other.next_av_bit;
        max_bck_dist = other.max_bck_dist;
        n_elms = other.n_elms;
        m_load_factor = other.m_load_factor;
        m_max_load_factor = other.m_max_load_factor;
        max_key_bits = other.max_key_bits;
        file = other.file;
        dump_fs.open(file, std::ios::app | std::ios::binary | std::ios::out);
        d_bits = other.d_bits;
        static_buffer = other.static_buffer;
        data = other.data;
    }

    void static_init(size_t buffer_size, void* buff_addr){

        //check the address is aligned
        assert(((uintptr_t)buff_addr % sizeof(size_t))==0);
        memset(buff_addr, 0, buffer_size);

        if(!file.empty()){
            dump_fs.open(file, std::ios::out | std::ios::binary);
            assert(dump_fs.good());
        }

        //define the size of the hash table
        n_buckets = (buffer_size/4)/sizeof(size_t); //25% buffer for the hash table
        n_buckets = next_power_of_two(n_buckets); //resize the table to the next power of two (50% of the buffer for the table, at most)
        assert(n_buckets>=4);
        table = reinterpret_cast<size_t*>(buff_addr);

        //define the size of the data
        size_t rem_buff = buffer_size - n_buckets*sizeof(size_t);//the other half is for the data

        //the address should be aligned!!
        assert(((uintptr_t)(table+n_buckets) % sizeof(buff_t))==0);
        data.stream_size = rem_buff/sizeof(buff_t);//the number of elements is floored
        data.stream = reinterpret_cast<buffer_t*>(table+n_buckets);

        buff_limit = (data.stream_size-1)*sizeof(buff_t)*8;
        next_av_bit = 1;
        m_load_factor = 0;
        max_bck_dist = 0;
        n_elms = 0;
        max_key_bits = 0;
        max_buffer_bytes = data.stream_size*sizeof(buff_t)+n_buckets*sizeof(size_t);
        assert(max_buffer_bytes<=buffer_size);
    }


    void dynamic_init(size_t buffer_size){
        assert(!static_buffer);

        n_buckets = 4;//minimum size for the hash table

        size_t table_bytes = (n_buckets * sizeof(size_t)); //number of bytes allocated for the hash table

        //number of bytes allocated for the data buffer
        size_t data_bytes = table_bytes*2;

        if(buffer_size!=0){
            //approximate to a number divisible by sizeof(buff_t)
            max_buffer_bytes = INT_CEIL(buffer_size, sizeof(buff_t))*sizeof(buff_t);
            data_bytes = std::min(max_buffer_bytes-table_bytes, data_bytes);

            if(!file.empty()){
                dump_fs.open(file, std::ios::out | std::ios::binary);
                assert(dump_fs.good());
            }
        }else{
            assert(file.empty());
            max_buffer_bytes = std::numeric_limits<size_t>::max();
        }

        table = reinterpret_cast<size_t*>(malloc(table_bytes));
        memset(table, 0, table_bytes);

        data.stream = reinterpret_cast<buffer_t *>(malloc(data_bytes));
        data.stream_size = data_bytes/sizeof(buffer_t);
        memset(data.stream, 0, data_bytes);

        next_av_bit = 1;
        m_load_factor = 0;
        max_bck_dist = 0;
        n_elms = 0;
        max_key_bits = 0;
    }

public:

    typedef value_t       val_type;
    typedef buffer_t      buff_t;

    explicit buffered_hash_table(size_t buffer_size=0, std::string file_="", float max_load_factor=0.8, void* buff_addr=nullptr, size_t _d_bits=desc_bits){//size of the buffer
        //buffer_size==0 means there is no limit in memory consumption
        d_bits = _d_bits;
        file = std::move(file_);
        m_max_load_factor = max_load_factor;
        if(buff_addr== nullptr){
            static_buffer = false;
            dynamic_init(buffer_size);
        }else{
            assert(buffer_size>=64);//must be at least 32 bytes
            static_buffer = true;
            static_init(buffer_size, buff_addr);
        }
        assert(static_buffer);
    }

    buffered_hash_table(buffered_hash_table&& other) noexcept{
        move(std::forward<buffered_hash_table>(other));
    }

    //copy assignment operator
    buffered_hash_table& operator=(buffered_hash_table & other){
        if(this!=&other){
            copy(other);
        }
        return *this;
    }

    //move assignment operator
    buffered_hash_table& operator=(buffered_hash_table && other) noexcept{
        if(this!=&other){
            move(std::forward<buffered_hash_table>(other));
        }
        return *this;
    }

    //swap method
    void swap(buffered_hash_table& other){
        move(std::forward<buffered_hash_table>(other));
    }

    inline iterator begin()  {
        return iterator(*this, 0);
    }

    inline iterator end() const  {
        return iterator(*this, next_av_bit);
    }

    inline size_t value_bits() const{
        return val_bits;
    }

    inline size_t description_bits() const{
        return d_bits;
    }

    size_t serialize(std::ostream &out) const{
        assert(!static_buffer);
        size_t written_bytes = serialize_elm(out, static_buffer);
        written_bytes += serialize_raw_vector(out, table, n_buckets);
        written_bytes += serialize_elm(out, max_buffer_bytes);
        written_bytes += serialize_elm(out, next_av_bit);
        written_bytes += serialize_elm(out, max_bck_dist);
        written_bytes += serialize_elm(out, n_elms);
        written_bytes += serialize_elm(out, m_load_factor);
        written_bytes += serialize_elm(out, m_max_load_factor);
        written_bytes += serialize_elm(out, max_key_bits);
        written_bytes += serialize_plain_vector(out, file);
        data.serialize(out);
        written_bytes += serialize_elm(out, d_bits);

        return written_bytes;
    }

    void load(std::istream &in){
        assert(!static_buffer);
        load_elm(in, static_buffer);
        load_raw_vector(in, table, n_buckets);
        load_elm(in, max_buffer_bytes);
        load_elm(in, next_av_bit);
        load_elm(in, max_bck_dist);
        load_elm(in, n_elms);
        load_elm(in, m_load_factor);
        load_elm(in, m_max_load_factor);
        load_elm(in, max_key_bits);
        load_plain_vector(in, file);
        data.load(in);
        load_elm(in, d_bits);
        if(!file.empty()){
            dump_fs.open(file, std::ios::app | std::ios::out | std::ios::binary);
        }
    }

    template<class val_t = value_t>
    typename std::enable_if<std::is_unsigned<val_t>::value, void>::type
    increment_value(const void* key, const size_t key_bits, value_t new_value) {
        auto res = insert(key, key_bits, new_value);
        if(res.second) return;

        value_t old_value = data.read(res.first+d_bits+key_bits, res.first+d_bits+key_bits+val_bits-1);
        new_value += old_value;
        data.write(res.first+d_bits+key_bits, res.first+d_bits+key_bits+val_bits-1, new_value);
    }

    /*std::pair<size_t, bool> insert_new(const void* key, const size_t& key_bits, const value_t& val){

        //size_t hash = XXH3_64bits(key, INT_CEIL(key_bits, 8));
        size_t hash = HASH(key, INT_CEIL(key_bits, 8));

        size_t idx = hash & (n_buckets - 1), in_offset;

        if(table[idx]==0){
            auto res  = insert_int(key, key_bits, val);
            in_offset = res.first;
            table[idx] =  in_offset;
        }else{

            size_t idx1 = (idx+1) & (n_buckets-1);
            size_t idx2 = (idx+2) & (n_buckets-1);
            size_t idx3 = (idx+3) & (n_buckets-1);
            unsigned int bck_close = (table[idx1]==0) | ((table[idx2]==0)<<1UL) | ((table[idx3]==0)<<2);
            if(bck_close!=3){
                size_t next_av = sizeof(unsigned int) - __builtin_clz(bck_close | 1);
                size_t dist1 = table[idx1]>>44UL;
                size_t dist2 = table[idx2]>>44UL;
                size_t dist3 = table[idx3]>>44UL;
                bool check_col = dist1==1U | dist2==2UL | dist3==3UL;
            }else{
                size_t bck_dist = table[idx]>>44UL, bck_offset, dist=0;
                while(table[idx]!=0 && bck_dist>dist){
                    idx = (idx+1) & (n_buckets-1);
                    bck_dist = table[idx] >> 44UL;
                    dist++;
                }
                while(table[idx]!=0 && bck_dist==dist){
                    bck_offset = table[idx] & 0xFFFFFFFFFFFul;
                    if(equal(key, key_bits , bck_offset-1)){
                        return {bck_offset-1, false};
                    }
                    idx = (idx+1) & (n_buckets-1);
                    bck_dist = table[idx] >> 44UL;
                    dist++;
                }

                assert(dist>0 && bck_dist<dist);
                auto res = insert_int(key, key_bits, val);
                size_t tmp_offset = res.first;
                in_offset = tmp_offset;

                while(table[idx]!=0){
                    bck_dist = table[idx] >> 44UL;
                    if(bck_dist<dist){
                        bck_offset = table[idx] & 0xFFFFFFFFFFFul;
                        table[idx] = (dist<<44UL) | tmp_offset;
                        if(dist>max_bck_dist) max_bck_dist = dist;
                        dist = bck_dist;
                        tmp_offset = bck_offset;
                    }
                    dist++;
                    idx = (idx+1) & (n_buckets-1);
                }

                assert(table[idx]==0);
                table[idx] = (dist<<44UL) | tmp_offset;
                if(dist>max_bck_dist) max_bck_dist = dist;
            }
        }

        n_elms++;
        m_load_factor = float(n_elms) / n_buckets;

        if(m_load_factor>=m_max_load_factor){
            size_t new_size = new_table_size();
            assert(new_size<=av_tb_bytes());
            resize_table(new_size);
        }
        return {in_offset-1, true};
    }*/

    void resize_table(size_t new_n_buckets) {
        assert(!static_buffer && new_n_buckets>0 && (new_n_buckets & (new_n_buckets-1))==0);
        size_t new_table_bytes = new_n_buckets*sizeof(size_t);

        assert(new_table_bytes <= av_tb_bytes());
        table = reinterpret_cast<size_t*>(realloc(table, new_table_bytes));
        n_buckets = new_n_buckets;
        memset(table, 0, new_n_buckets * sizeof(size_t));
        if(n_elms>0){
            rehash();
            m_load_factor = float(n_elms) / n_buckets;
        }
    };

    size_t longest_key() const {
        return max_key_bits;
    }

    std::pair<size_t, bool> insert(const void* key, const size_t& key_bits, const value_t& val){

        //size_t hash = XXH3_64bits(key, INT_CEIL(key_bits, 8));
        size_t hash = HASH(key, INT_CEIL(key_bits, 8));

        size_t idx = hash & (n_buckets - 1);
        size_t in_offset=0;//locus where the key is inserted

        if(table[idx]==0) {
            auto res  = insert_int(key, key_bits, val);
            //if(res.second) idx = hash & (n_buckets-1);
            in_offset = res.first;
            table[idx] =  in_offset;
        }else{ //resolve collision

            size_t dist=0, offset;
            size_t bck_dist, bck_offset;
            bool inserted = false;

            while(table[idx]!=0){

                bck_offset = (table[idx] & 0xFFFFFFFFFFFul);
                bck_dist = table[idx] >> 44UL;

                if(!inserted && bck_dist==dist && equal(key, key_bits , bck_offset-1)){
                    return {bck_offset-1, false};
                }else if(bck_dist<dist){ //steal to the rich

                    if(!inserted){//swap the keys
                        auto res = insert_int(key, key_bits, val);
                        in_offset = res.first;
                        if(!res.second){//was the data dumped?
                            offset = res.first;
                            inserted = true;
                        }else{
                            //data was dumped, it is not necessary to do the swap
                            table[hash & (n_buckets-1)] =  res.first;
                            goto finish;
                        }
                    }

                    assert(dist<=limit_bck_dist);
                    if(dist>max_bck_dist) max_bck_dist = dist;

                    table[idx] = (dist<<44UL) | offset;
                    offset = bck_offset;
                    dist = bck_dist;
                }
                dist++;
                idx = (idx+1) & (n_buckets - 1);
            }

            if(!inserted){
                auto res = insert_int(key, key_bits, val);
                offset = res.first;
                in_offset = offset;
                if(res.second){//data was dumped
                    idx = hash & (n_buckets-1);
                    dist = 0;
                }
            }

            assert(dist<=limit_bck_dist);
            if(dist>max_bck_dist) max_bck_dist = dist;

            table[idx] = (dist<<44UL) | offset ;
        }

        finish:
        n_elms++;
        m_load_factor = float(n_elms) / n_buckets;

        bool dump;
        if(m_load_factor>=m_max_load_factor){

            if(static_buffer){
                dump=true;
            }else{
                size_t new_size = new_table_size();
                if((new_size*sizeof(size_t)) > av_tb_bytes()){
                    //the new size doesn't fit the buffer
                    dump = true;
                }else{
                    //rehashing does not affect the position of the pair
                    dump = false;
                    resize_table(new_size);
                }
            }

            if(dump){
                //invalidate the last insertion and dump the hash table
                assert(next_av_bit>1);
                next_av_bit -= (key_bits+val_bits+d_bits);//this step crops the data buffer, so the last insertion is not included in the dump
                dump_hash();
                assert(n_elms==0 && m_load_factor==0);

                //insert again the new element
                auto res = insert_int(key, key_bits, val);
                assert(!res.second);
                in_offset = res.first;
                assert(table[hash & (n_buckets-1)]==0);
                table[hash & (n_buckets-1)] =  res.first;

                n_elms++;
                m_load_factor = float(n_elms) / n_buckets;
            }
        }

        assert(in_offset>0);
        return {in_offset-1, true};
    }

    inline float load_factor() const{
        return m_load_factor;
    }

    inline float max_load_factor() const{
        return m_max_load_factor;
    }

    inline size_t size() const {
        return n_elms;
    }

    inline void shrink_databuff() {
        assert(!static_buffer);
        size_t new_size = INT_CEIL(next_av_bit, stream_t::word_bits)+1;
        data.stream = reinterpret_cast<buffer_t*>(realloc(data.stream, new_size*sizeof(buffer_t)));
        data.stream_size = new_size;
    }

    inline void reset(){
        assert(!static_buffer);
        n_buckets = 4;
        n_elms = 0;
        max_bck_dist = 0;
        m_load_factor = 0;
        next_av_bit = 1;

        size_t table_bytes = n_buckets*sizeof(size_t);
        table = reinterpret_cast<size_t*>(realloc(table, table_bytes));

        size_t data_bytes = table_bytes*2;
        data.stream = reinterpret_cast<buffer_t*>(realloc(data.stream, data_bytes));
        data.stream_size = data_bytes/sizeof(buffer_t);
        memset(table, 0, table_bytes);
        memset(data.stream, 0, data_bytes);
    }

    inline size_t max_bucket_dist() const {
        return max_bck_dist;
    }

    inline std::pair<size_t, bool> find(const void* key, size_t key_bits) const {

        //size_t hash = XXH3_64bits(key, INT_CEIL(key_bits, 8));
        size_t hash = HASH(key, INT_CEIL(key_bits, 8));

        size_t idx = hash & (n_buckets - 1);

        if(table[idx]==0){
            return {0, false};
        }else{
            size_t offset, i=0, bck_dist;
            while(i<=max_bck_dist && table[idx]!=0){
                bck_dist = (table[idx] >> 44U);
                if(bck_dist<i){
                    return {0, false};
                }
                offset = (table[idx] & 0xFFFFFFFFFFFul);
                if(i==bck_dist && equal(key, key_bits, offset-1)){
                    return {offset - 1, true};
                }
                idx = (idx+1) & (n_buckets - 1);
                i++;
            }
            return {0, false};
        }
    }

    inline bool key2value(const void* key, size_t key_bits, value_t& value) const {

        //size_t hash = XXH3_64bits(key, INT_CEIL(key_bits, 8));
        size_t hash = HASH(key, INT_CEIL(key_bits, 8));

        size_t idx = hash & (n_buckets - 1);
        size_t i=0;
        size_t bck_dist = (table[idx] >> 44U);
        size_t offset = (table[idx] & 0xFFFFFFFFFFFul);
        while(table[idx]!=0 && i<=bck_dist){
            if(i==bck_dist && equal(key, key_bits, offset-1)){
                size_t value_start = offset-1+d_bits+key_bits;
                if constexpr (val_bits<=64){
                    value = data.read(value_start, value_start+val_bits-1);
                }else{
                    data.read_chunk(value, value_start, value_start+val_bits-1);
                }
                return true;
            }
            idx = (idx+1) & (n_buckets - 1);
            i++;
            bck_dist = (table[idx] >> 44U);
            offset = (table[idx] & 0xFFFFFFFFFFFul);
        }
        return false;
    }

    inline bool get_value(const void* key, size_t key_bits, value_t& value) const {

        //size_t hash = XXH3_64bits(key, INT_CEIL(key_bits, 8));
        size_t hash = HASH(key, INT_CEIL(key_bits, 8));

        size_t idx = hash & (n_buckets - 1);

        if(table[idx]==0){
            return false;
        }else{
            size_t offset, i=0, bck_dist;
            while(i<=max_bck_dist && table[idx]!=0){
                bck_dist = (table[idx] >> 44U);
                if(bck_dist<i){
                    return false;
                }
                offset = (table[idx] & 0xFFFFFFFFFFFul);
                if(i==bck_dist && equal(key, key_bits, offset-1)){
                    if constexpr (val_bits<=64){
                        value = data.read(offset-1+d_bits+key_bits, offset-1+d_bits+key_bits+val_bits-1);
                    }else{
                        data.read_chunk(value, offset-1+d_bits+key_bits, offset-1+d_bits+key_bits+val_bits-1);
                    }
                    return true;
                }
                idx = (idx+1) & (n_buckets - 1);
                i++;
            }
            return false;
        }
    }

    //offset is the bit where the pair description starts
    inline void insert_value_at(size_t offset, value_t val){
        size_t value_start = offset + d_bits + data.read(offset, offset + d_bits - 1);
        if constexpr (val_bits<=64){
            data.write(value_start, value_start + val_bits - 1, val);
        }else{
            data.write_chunk(&val, value_start, value_start + val_bits - 1);
        }
    }

    //offset is the bit where the pair description starts
    inline void get_value_from(size_t offset, value_t& dest) const {
        size_t value_start = offset + d_bits + data.read(offset, offset + d_bits - 1);
        if constexpr (val_bits<=64){
            dest = data.read(value_start, value_start + val_bits - 1);
        }else{
            data.read_chunk(&dest, value_start, value_start + val_bits - 1);
        }
    }

    //offset is the bit where the pair description starts
    inline size_t get_key(size_t offset, void* dest) const {
        size_t key_bits = data.read(offset, offset + d_bits - 1);;
        data.read_chunk(dest, offset+d_bits, offset+d_bits+key_bits-1);
        return key_bits;
    }

    ~buffered_hash_table(){
        if(!static_buffer){
            flush();
            if(table!= nullptr){
                free(table);
            }
            if(data.stream!= nullptr){
                free(data.stream);
            }
        }
        dump_fs.close();
    }

    void flush(){
        if(next_av_bit>1 && max_buffer_bytes!=std::numeric_limits<size_t>::max()){
            dump_hash();
            //after the dump there might remain some bits of the buffer tail,
            // which are moved to the buffer head
            if(next_av_bit>1 && !file.empty()){
                dump_fs.write((char *) data.stream, INT_CEIL(next_av_bit-1, 8));
                assert(dump_fs.good());
            }
            dump_fs.flush();
        }
    }

    void dump_hash(){

        //next_av_bit is one-based!!
        assert(next_av_bit>=1);
        size_t written_bits = next_av_bit-1;
        size_t tot_bytes = INT_CEIL(written_bits, 8);
        size_t bytes_to_write;

        buffer_t tail=0;
        if(tot_bytes>0) {
            size_t rem = written_bits % bitstream<buffer_t>::word_bits;
            if (rem==0) {
                bytes_to_write = tot_bytes;
                next_av_bit=1;
            } else {
                //leave the tail in the buffer
                size_t n_cells = INT_CEIL(written_bits, bitstream<buffer_t>::word_bits);
                bytes_to_write = (n_cells-1)*sizeof(buffer_t);
                tail = data.stream[n_cells-1] & bitstream<buffer_t>::masks[rem];
                next_av_bit = rem+1;
            }

            if(!file.empty()){
                dump_fs.write((char *) data.stream, bytes_to_write);
                assert(dump_fs.good());
            }
        }
        n_elms = 0;
        max_bck_dist = 0;
        m_load_factor = 0;

        // a buffer size of (1<<64)-1 bytes indicates
        // that there is no limit for the size
        // of the buffer
        if(max_buffer_bytes!=std::numeric_limits<buffer_t>::max() && !static_buffer){
            //half of the buffer for the hash table
            n_buckets = INT_CEIL( (max_buffer_bytes/2), sizeof(size_t));
            n_buckets = 1UL << sym_width(n_buckets);//next power of two
            table = reinterpret_cast<size_t*>(realloc(table, n_buckets*sizeof(size_t)));

            //half of the buffer for the data
            data.stream_size = INT_CEIL((max_buffer_bytes/2), sizeof(buffer_t));
            data.stream = reinterpret_cast<buffer_t *>(realloc(data.stream, data.stream_size*sizeof(buffer_t)));
        }

        memset(data.stream, 0, data.stream_size*sizeof(size_t));
        data.stream[0] = tail;

        memset(table, 0, n_buckets * sizeof(size_t));
    };

    inline size_t tot_buckets() const{
        return n_buckets;
    }

    inline const stream_t& get_data(){
        return data;
    }

    inline void clean(){
        memset(table, 0, n_buckets*sizeof(size_t));
        memset(data.stream, 0, data.stream_size*sizeof(buffer_t));
        n_elms = 0;
        max_bck_dist = 0;
        m_load_factor = 0;
        next_av_bit = 1;
    }

    const std::string& dump_file() const{
        return file;
    }

    void store_data_to_file(const std::string& output){
        assert(!static_buffer);
        std::filebuf fb;
        fb.open(output, std::ios::out | std::ios::binary);
        std::ostream ofs(&fb);
        ofs.write(reinterpret_cast<char *>(table), n_buckets*sizeof(size_t));
        ofs.tellp();
        data.serialize(ofs);
        fb.close();
    }

    size_t hash_table_tot_bytes(){
        size_t tot_bytes = n_buckets*sizeof(size_t);
        tot_bytes+= sizeof(data.stream_size);
        tot_bytes+= data.stream_size * sizeof(buffer_t);
        return  tot_bytes;
    }

    void ht_stats(size_t top){
        std::vector<size_t> freqs(max_bck_dist+1, 0);
        size_t dist;
        for(size_t i=0;i<n_buckets;i++){
            if(table[i]!=0){
                dist = table[i]>>44UL;
                assert(dist<=max_bck_dist);
                freqs[dist]++;
            }
        }

        std::vector<std::pair<size_t, size_t>> pairs(freqs.size());
        for(size_t i=0;i<pairs.size();i++){
            pairs[i] = {i, freqs[i]};
        }
        std::sort(pairs.begin(), pairs.end(), [&](auto &a, auto &b){
            return a.second>b.second;
        });

        std::cout<<"\nsize "<<size()<<std::endl;
        std::cout<<"n_buckets "<<tot_buckets()<<std::endl;
        std::cout<<"load factor: "<<load_factor()<<std::endl;
        std::cout<<"max_bck_dist: "<<max_bck_dist<<std::endl;
        for(size_t i=0;i<=std::min(top, max_bck_dist);i++){
            std::cout<<"Bck dist:"<<pairs[i].first<<" freq:"<<pairs[i].second<<" ( "<<(double(pairs[i].second)/n_elms)*100<<"% )"<<std::endl;
        }
    }

    void load_data_from_file(const std::string& input){
        assert(table==nullptr && data.stream== nullptr && !static_buffer);
        std::ifstream ifs(input, std::ios_base::binary);
        table = reinterpret_cast<size_t*>(malloc(n_buckets*sizeof(size_t)));
        ifs.read(reinterpret_cast<char *>(table), n_buckets*sizeof(size_t));
        data.load(ifs);
        ifs.close();
    }

    size_t data_bytes() const {
        return (data.stream_size*stream_t::word_bits)/8;
    }

    void destroy_table() {
        free(table);
        table = nullptr;
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    void destroy_data() {
        free(data.stream);
        data.stream = nullptr;
#ifdef __linux__
        malloc_trim(0);
#endif
    }
};
#endif //LMS_COMPRESSOR_SE_HASH_TABLE_H
