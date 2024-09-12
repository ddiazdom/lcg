//
// Created by Diaz, Diego on 3.4.2023.
//

#ifndef SIMPLE_EXAMPLE_TS_STRING_MAP_H
#define SIMPLE_EXAMPLE_TS_STRING_MAP_H
#include <atomic>
#include <vector>
#include <thread>
#include <cassert>
#include <cstring>
#include <fstream>
#include "cdt_common.hpp"
#include "../external/xxHash-dev/xxhash.h"

thread_local size_t tr_id;

#ifdef __linux__
#include <malloc.h>
#endif

#if defined(__x86_64__) || defined(_M_X64)
#define YIELD_PROCESSOR __builtin_ia32_pause()
#else
#define YIELD_PROCESSOR __asm__ __volatile__("yield")
#endif

struct ttas_lock{
    //std::atomic_flag flag=ATOMIC_FLAG_INIT;
    std::atomic<bool> flag{false};
    inline void lock() {
        for (;;) {
            if(!flag.exchange(true, std::memory_order_acquire)) {
                break;
            }
            while(flag.load(std::memory_order_relaxed)) {
                YIELD_PROCESSOR;
            }
        }
        //while(flag.test_and_set(std::memory_order_acquire));
    }

    inline void unlock(){
        flag.store(false, std::memory_order_release);
    }
};

template<class value_t>
class proxy {
private:
    char * addr= nullptr;
    uint8_t bytes = 0;

public:
    proxy()= default;

    proxy(char * _addr, uint8_t _bytes) : addr(_addr), bytes(_bytes) {}

    operator value_t() const {
        value_t val = 0;
        memcpy(&val, addr, bytes);
        return val;
    }

    proxy& operator=(value_t value) {
        memcpy(addr, &value, bytes);
        return *this;
    }

    /*proxy& operator=(const proxy& other) {
        auto val = (value_t)other;
        memcpy(addr, &val, bytes);
        return *this;
    }*/

    bool operator==(const proxy& x)const {
        return value_type(*this) == value_type(x);
    }

    bool operator<(const proxy& x)const {
        return value_type(*this) < value_type(x);
    }
};

template<class map_type>
class string_map_iterator{

    using value_t =                              typename map_type::value_type;

    const map_type&                              map;
    char * sm_buff=nullptr;
    char * sm_buff_limit= nullptr;
    size_t sm_idx = 0;
    size_t entry_idx{};
    size_t tot_entries;
    uint8_t key_len_bytes=0;
    uint8_t value_bytes=0;

public:


    typedef std::forward_iterator_tag                                  iterator_category;
    typedef std::pair<std::pair<const char *, size_t>, proxy<value_t>> entry_type;
    typedef entry_type                                                 reference;
    entry_type                                                         entry;

    string_map_iterator(const map_type& map_, size_t entry_idx): map(map_),
                                                                 entry_idx(entry_idx),
                                                                 tot_entries(map.size()) {
        if(entry_idx<tot_entries){

            size_t entry_tmp=0;
            sm_idx=0;
            while(map.sub_map_array[sm_idx]->empty() ||
                  (entry_tmp+map.sub_map_array[sm_idx]->size())<entry_idx){
                entry_tmp+=map.sub_map_array[sm_idx++]->size();
            }

            if(sm_idx<map.sub_map_array.size()){

                sm_buff = map.sub_map_array[sm_idx]->get_buffer();
                sm_buff_limit = sm_buff+map.sub_map_array[sm_idx]->used_bytes_in_buffer();
                key_len_bytes = map.sub_map_array[sm_idx]->bytes_for_key_len();
                value_bytes = map.sub_map_array[sm_idx]->bytes_for_value();

                size_t len;
                while(entry_tmp<entry_idx){
                    len = 0;
                    memcpy(&len, sm_buff, key_len_bytes);
                    sm_buff+=key_len_bytes+len+value_bytes;
                    entry_tmp++;
                }
                assert(entry_idx==entry_tmp);

                len = 0;
                memcpy(&len, sm_buff, key_len_bytes);
                sm_buff+=key_len_bytes;
                entry.first = {sm_buff, len};
                sm_buff+=len;
                entry.second = proxy<value_t>(sm_buff, value_bytes);
                sm_buff+=value_bytes;

                if(sm_buff==sm_buff_limit){
                    sm_idx++;
                    while(sm_idx<map.sub_map_array.size() && map.sub_map_array[sm_idx]->empty()){
                        sm_idx++;
                    }
                    if(sm_idx<map.sub_map_array.size()){
                        sm_buff = map.sub_map_array[sm_idx]->get_buffer();
                        sm_buff_limit = sm_buff+map.sub_map_array[sm_idx]->used_bytes_in_buffer();
                        key_len_bytes = map.sub_map_array[sm_idx]->bytes_for_key_len();
                        value_bytes = map.sub_map_array[sm_idx]->bytes_for_value();
                    }
                }
            } else{
                entry = {{nullptr, 0}, proxy<value_t>(nullptr, 0)};
            }
        }
    };

    string_map_iterator(const string_map_iterator& other): map(other.map),
                                                           entry(other.entry),
                                                           sm_idx(other.sm_idx),
                                                           entry_idx(other.entry_idx),
                                                           tot_entries(other.tot_entries){};

    inline reference operator*() {
        return entry;
    }

    const string_map_iterator operator++(int) {
        string_map_iterator<map_type> tmp(*this);
        ++*this;
        return tmp;
    }

    inline bool operator==(const string_map_iterator& other) const{
        return entry_idx == other.entry_idx;
    }

    inline bool operator!=(const string_map_iterator& other) const{
        return entry_idx != other.entry_idx;
    }

    string_map_iterator& operator=(string_map_iterator const& other) {
        if(other==this) return *this;
        map = other.map;
        entry = other.entry;
        sm_buff = other.sm_buff;
        sm_buff_limit = other.sm_buff_limit;
        sm_idx = other.sm_idx;
        entry_idx = other.entry_idx;
        tot_entries = other.tot_entries;
        key_len_bytes = other.key_len_bytes;
        value_bytes = other.value_bytes;
        return *this;
    }

    const string_map_iterator& operator++() {
        entry_idx++;
        if(entry_idx < tot_entries){
            size_t len = 0;
            memcpy(&len, sm_buff, key_len_bytes);
            sm_buff+=key_len_bytes;
            entry.first = {sm_buff, len};
            sm_buff+=len;
            entry.second = proxy<value_t>(sm_buff, value_bytes);
            sm_buff+=value_bytes;

            if(sm_buff==sm_buff_limit){
                sm_idx++;
                while(sm_idx<map.sub_map_array.size() && map.sub_map_array[sm_idx]->empty()){
                    sm_idx++;
                }
                if(sm_idx<map.sub_map_array.size()){
                    sm_buff = map.sub_map_array[sm_idx]->get_buffer();
                    sm_buff_limit = sm_buff+map.sub_map_array[sm_idx]->used_bytes_in_buffer();
                    key_len_bytes = map.sub_map_array[sm_idx]->bytes_for_key_len();
                    value_bytes = map.sub_map_array[sm_idx]->bytes_for_value();
                }
            }
        } else {
            entry_idx = tot_entries;
            entry = {{nullptr, 0}, proxy<value_t>(nullptr, 0)};
        }
        return *this;
    }
};

template<typename val_type>
class nts_string_submap {
public:
    static constexpr size_t empty_entry = std::numeric_limits<size_t>::max();
    //static constexpr size_t empty_entry = 1099511627775; //The hash table is at most 1TB in size
    static constexpr uint8_t key_len_bytes=4;
    static constexpr uint8_t value_bytes=5;

    struct bucket_t{
        size_t key_offset=empty_entry;
        //val_type value=0;
        [[nodiscard]] inline bool empty() const {
            return key_offset==empty_entry;
        }
    };

private:

    typedef std::vector<bucket_t> table_t;
    table_t* m_table{nullptr};
    size_t n_elms = 0;
    float m_max_load_factor = 0.6;

    //buffer for the keys
    char * buffer = nullptr;
    size_t next_av_offset = 0;
    size_t buff_size = 0;
    //

    void rehash(size_t new_tab_size) {

        assert(new_tab_size>m_table->size());
        auto * new_table = new table_t(new_tab_size);

        //rehash the values
        for(auto & bucket : *m_table) {
            if(bucket.key_offset!= empty_entry){
                insert_entry_in_table_bucket(*new_table, new_tab_size, bucket);
            }
        }
        std::swap(m_table, new_table);
        delete new_table;
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    inline void insert_entry_in_table_bucket(table_t& new_table, size_t new_tab_size, bucket_t& bucket) const {

        size_t key_len=0;
        char * addr = buffer+bucket.key_offset;
        memcpy(&key_len, addr, key_len_bytes);
        size_t hash = XXH3_64bits(addr+key_len_bytes, key_len);

        size_t idx = hash & (new_tab_size-1);

        if(new_table[idx].key_offset==empty_entry){
            new_table[idx].key_offset = bucket.key_offset;
            //new_table[idx].value = bucket.value;
            return;
        }

        //find a new empty sm_bucket
        size_t j=1;
        while(true){
            idx = (hash + ((j*j+j)>>1UL)) & (new_tab_size-1);
            if(new_table[idx].key_offset==empty_entry){
                new_table[idx].key_offset = bucket.key_offset;
                //new_table[idx].value = bucket.value;
                break;
            }
            j++;
        }
    }

    char * create_new_buffer(size_t min_bytes) {
        assert(min_bytes>=buff_size);

        size_t new_buff_size = std::max<size_t>(buff_size, 1);
        while(new_buff_size<min_bytes){
            new_buff_size<<=1UL;
        }

        //buffer = (char *)realloc(buffer, new_buff_size);
        buffer = mem<char>::reallocate(buffer, new_buff_size);
        memset((void *)(buffer+next_av_offset), 0, (new_buff_size-next_av_offset));

        buff_size = new_buff_size;
        return buffer;
    }

public:

    /*
    [[nodiscard]] inline std::pair<const char *, size_t> get_key_pointer(const size_t idx) const {
        assert((*m_table)[idx].key_offset!=empty_entry);
        char * addr = buffer+(*m_table)[idx].key_offset;
        size_t len = 0;
        memcpy(&len, addr, key_len_bytes);
        return {addr+key_len_bytes, len};
    }

    [[nodiscard]] inline char * get_value_pointer(const size_t idx) const {
        assert((*m_table)[idx].key_offset!=empty_entry);
        char * addr = buffer+(*m_table)[idx].key_offset;
        size_t len = 0;
        memcpy(&len, addr, key_len_bytes);
        addr+=key_len_bytes+len;
        return addr;
    }*/

    table_t *const& table = m_table;

    explicit nts_string_submap(size_t min_cap=4, float max_lf=0.6, size_t m_buff_size=0) : m_max_load_factor(max_lf),
                                                                                           buff_size(m_buff_size) {

        m_table = new table_t(round_to_power_of_two(min_cap));
        if(buff_size!=0){
            //buffer = (char *) malloc(buff_size);
            buffer = mem<char>::allocate(buff_size);
        }
    }

    ~nts_string_submap(){
        delete m_table;
        //free(buffer);
        mem<char>::deallocate(buffer);
    }

    bool value_add(const uint8_t* key, size_t len, val_type val, size_t hash) {

        bool inserted = false;
        bool success = false;

        size_t j=0;
        size_t idx = hash & (m_table->size()-1);

        while(!success) {
            bucket_t &bucket = (*m_table)[idx];

            if(bucket.key_offset==empty_entry) {

                size_t key_bytes = value_bytes+key_len_bytes+len;
                size_t offset = next_av_offset;
                size_t min_bytes =  next_av_offset+key_bytes;

                if(min_bytes>=buff_size){//the key does not fit the buffer, allocate a new buffer
                    create_new_buffer(min_bytes);
                }
                next_av_offset += key_bytes;
                n_elms++;
                bucket.key_offset = offset;
                char * addr = buffer+offset;
                memcpy(addr, &len, key_len_bytes);
                addr+=key_len_bytes;
                memcpy(addr, key, len);
                addr+=len;
                memcpy(addr, &val, value_bytes);

                //the key insertion exceeds the max. load factor (i.e., rehash)
                if((float(n_elms)/float(m_table->size()))>=m_max_load_factor) {
                    rehash(next_power_of_two(m_table->size()));
                }
                success = true;
                inserted = true;

            } else if(memcmp(&len, buffer+bucket.key_offset, key_len_bytes)==0 &&
                      memcmp(key, buffer+bucket.key_offset+key_len_bytes, len)==0){
                val_type new_val = 0;
                memcpy(&new_val, buffer+bucket.key_offset+key_len_bytes+len, value_bytes);
                new_val+=val;
                memcpy(buffer+bucket.key_offset+key_len_bytes+len, &new_val, value_bytes);
                success = true;
            } else {
                j++;
                idx = (hash + ((j*j + j)>>1UL)) & (m_table->size()-1);
            }
        }
        return inserted;
    }

    std::pair<proxy<val_type>, bool> insert(const uint8_t* key, size_t len, val_type val, size_t hash) {

        bool inserted = false;
        bool success = false;
        char *val_ptr = nullptr;

        size_t j=0;
        size_t idx = hash & (m_table->size()-1);

        while(!success) {
            bucket_t &bucket = (*m_table)[idx];

            if(bucket.key_offset==empty_entry) {

                size_t key_bytes = value_bytes+key_len_bytes+len;
                size_t offset = next_av_offset;
                size_t min_bytes =  next_av_offset+key_bytes;

                if(min_bytes>=buff_size){//the key does not fit the buffer, allocate a new buffer
                    create_new_buffer(min_bytes);
                }
                next_av_offset += key_bytes;
                n_elms++;
                bucket.key_offset = offset;
                char * addr = buffer+offset;
                memcpy(addr, &len, key_len_bytes);
                addr+=key_len_bytes;
                memcpy(addr, key, len);
                addr+=len;
                memcpy(addr, &val, value_bytes);
                val_ptr = addr;

                //the key insertion exceeds the max. load factor (i.e., rehash)
                if((float(n_elms)/float(m_table->size()))>=m_max_load_factor) {
                    rehash(next_power_of_two(m_table->size()));
                }
                success = true;
                inserted = true;

            } else if(memcmp(&len, buffer+bucket.key_offset, key_len_bytes)==0 &&
                      memcmp(key, buffer+bucket.key_offset+key_len_bytes, len)==0){
                val_ptr = buffer+bucket.key_offset+key_len_bytes+len;
                success = true;
            } else {
                j++;
                idx = (hash + ((j*j + j)>>1UL)) & (m_table->size()-1);
            }
        }

        return {{val_ptr, value_bytes}, inserted};
    }

    [[nodiscard]] inline char * get_buffer() const {
        return buffer;
    }

    [[nodiscard]] inline size_t used_bytes_in_buffer() const {
        return next_av_offset;
    }

    [[nodiscard]] inline uint8_t bytes_for_key_len() const {
        return key_len_bytes;
    }

    [[nodiscard]] inline uint8_t bytes_for_value() const {
        return value_bytes;
    }

    void set_min_capacity(size_t new_cap){
        new_cap = round_to_power_of_two(new_cap);
        if(new_cap>m_table->size()){
            rehash(new_cap);
        }
    }

    bool find(const uint8_t * key, size_t len, size_t hash, val_type& value) const {

        size_t j=0;
        size_t idx = hash & (m_table->size()-1);
        while(true){
            if((*m_table)[idx].key_offset==empty_entry){
                return false;
            }else if(memcmp(&len, buffer+((*m_table)[idx].key_offset), key_len_bytes)==0 &&
                     memcmp(key, buffer+((*m_table)[idx].key_offset) + key_len_bytes, len)==0) {
                value = 0;
                memcpy(&value, buffer+((*m_table)[idx].key_offset)+key_len_bytes+len, value_bytes);
                return true;
            }
            j++;
            idx = (hash + ((j*j + j)>>1UL)) & (m_table->size()-1);
        }
    }

    [[nodiscard]] inline float load_factor() const {
        return float(n_elms)/float(m_table->size());
    }

    [[nodiscard]] inline float max_load_factor() const {
        return m_max_load_factor;
    };

    inline size_t size() {
        return n_elms;
    }

    inline bool empty()  {
        return n_elms==0;
    }

    [[nodiscard]] inline size_t capacity() const{
        return m_table->size();
    };

    size_t buffer_usage(){
        return next_av_offset;
    }

    size_t table_mem_usage(){
        return m_table->size()*sizeof(bucket_t);
    }

    size_t buffer_capacity(){
        return buff_size;
    }

    void destroy_table(){
        delete m_table;
        m_table = nullptr;
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    void destroy_buffer(){
        //free(buffer);
        mem<char>::deallocate(buffer);
        buffer = nullptr;
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    void shrink_to_fit(){
        //buffer = (char *)realloc(buffer, next_av_offset);
        buffer = mem<char>::reallocate(buffer, next_av_offset);
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    void load_table(std::ifstream & ifs){
        //read the table
        size_t n_buckets;
        ifs.read((char *)&n_buckets, sizeof(n_buckets));
        delete m_table;

        m_table = new table_t(n_buckets);
        ifs.read((char*)m_table->data(), sizeof(bucket_t)*n_buckets);
    }

    void load_buffer(std::ifstream & ifs){
        //read the variables
        ifs.read((char *)&n_elms, sizeof(n_elms));
        ifs.read((char *)&m_max_load_factor, sizeof(m_max_load_factor));

        //load the buffer with the keys
        ifs.read((char *)&next_av_offset, sizeof(next_av_offset));
        if(buffer!= nullptr){
            buffer = mem<char>::reallocate(buffer, next_av_offset);
        }else{
            buffer = mem<char>::allocate(next_av_offset);
        }
        buff_size = next_av_offset;
        ifs.read((char *)buffer, buff_size);
    }

    void load_submap(std::ifstream& ifs){
        load_table(ifs);
        load_buffer(ifs);
    }

    void store_table(std::ofstream & ofs){
        //serialize the table
        size_t n_buckets = table->size();
        ofs.write((char* )&n_buckets, sizeof(n_buckets));
        ofs.write((char *)table->data(), sizeof(bucket_t)*table->size());
    }

    void store_buffer(std::ofstream &ofs){
        //serialize the variables
        ofs.write((char *)&n_elms, sizeof(n_elms));
        ofs.write((char*)&m_max_load_factor, sizeof(m_max_load_factor));

        //serialize the buffer
        ofs.write((char*)&next_av_offset, sizeof(next_av_offset));
        ofs.write(buffer, next_av_offset);
        //NOTE we don't store buff_size because we set it equal to next_av_offset if we load it
    }

    void store_submap(std::ofstream& ofs){
        store_table(ofs);
        store_buffer(ofs);
    }

    /*inline void print_table() {
        for(auto & bucket : *m_table){
            if(bucket.key_offset!=empty_entry){
                char * tmp = (buffer+bucket.key_offset);
                for(size_t j=0;j<bucket.key_len;j++){
                    std::cout<<tmp[j]<<"";
                }
                std::cout<<" -> "<<bucket.value<<std::endl;
                //printf("%lu-%lu) %s -> %lu\n", i++, n_elms, , bucket.value);
            }
        }
    }*/
};

//TODO I have a problem with the bucket information. I need to add the string len (in bytes) to the bucket
template<typename val_type>
class ts_string_submap {

public:

    struct bucket_t {
        //std::atomic_flag access = ATOMIC_FLAG_INIT;
        ttas_lock access;
        char * key_ptr = nullptr;
        uint32_t key_len=0;
        val_type value=0;

        inline bool empty() const {
            return key_ptr== nullptr;
        }
    };

private:

    typedef std::vector<bucket_t> table_t;
    char * invalid_addr = reinterpret_cast<char *>(std::numeric_limits<uintptr_t>::max());

    struct protected_pointers{
        std::atomic_flag access = ATOMIC_FLAG_INIT;
        table_t * tab_ptr = nullptr;
    };

    table_t* m_table{nullptr};
    protected_pointers prot_ptrs[100]{};

    size_t n_elms = 0;
    float m_max_load_factor = 0.6;
    std::atomic_flag doing_rehash = ATOMIC_FLAG_INIT;//flag to rehash

    //buffer for the keys
    char * buffer = nullptr;
    //std::atomic_flag buffer_lock = ATOMIC_FLAG_INIT; // flag to increase the number of elements in the hash table
    ttas_lock dat_buff_lock;
    std::atomic<bool> moving_buffer = false;
    size_t next_av_offset = 0;
    size_t buff_size = 0;
    size_t prev_buff_size=0;
    char * prev_buff_ptr=nullptr;
    //

    std::atomic_flag mod_ht = ATOMIC_FLAG_INIT;

    void rehash(size_t new_tab_size, table_t* tab_ptr) {

        while(doing_rehash.test_and_set(std::memory_order_acquire));

        table_t * table_ptr = m_table;
        if(table_ptr!=tab_ptr || table_ptr->size()>=new_tab_size){
            doing_rehash.clear(std::memory_order_release);
            return;// a thread trying to rehash an old reference
        }

        while(mod_ht.test_and_set(std::memory_order_acquire));

        assert(new_tab_size>table_ptr->size());
        auto * new_table_ptr = new table_t(new_tab_size);

        //rehash the values
        for(auto & bucket : *table_ptr) {
            //while(bucket.access.test_and_set(std::memory_order_acquire));
            bucket.access.lock();
            if(bucket.key_ptr!= nullptr){
                insert_entry_in_table_bucket(*new_table_ptr, new_tab_size, bucket);
            }
            bucket.key_ptr=invalid_addr;
            bucket.access.unlock();
            //bucket.access.clear(std::memory_order_release);
        }

        //replace the old hash table
        m_table = new_table_ptr;
        for(auto & prot_tab_ptr : prot_ptrs){
            while(prot_tab_ptr.access.test_and_set(std::memory_order_acquire));
            prot_tab_ptr.tab_ptr = new_table_ptr;
            prot_tab_ptr.access.clear(std::memory_order_release);
        }
        //std::cout<<"I removed the old hash table "<<table_ptr->size()<<std::endl;
        //free(table_ptr);
        mem<table_t>::deallocate(table_ptr);
        mod_ht.clear(std::memory_order_release);
        doing_rehash.clear(std::memory_order_release);
    }

    inline void insert_entry_in_table_bucket(table_t& new_table, size_t new_tab_size, bucket_t& bucket) const {

        size_t hash = XXH3_64bits(bucket.key_ptr, strlen(bucket.key_ptr));

        size_t idx = hash & (new_tab_size-1);

        if(new_table[idx].key_ptr==nullptr){
            new_table[idx].key_ptr = bucket.key_ptr;
            new_table[idx].value = bucket.value;
            return;
        }

        //find a new empty sm_bucket
        size_t j=1;
        while(true){
            idx = (hash + ((j*j+j)>>1UL)) & (new_tab_size-1);
            if(new_table[idx].key_ptr==nullptr){
                new_table[idx].key_ptr = bucket.key_ptr;
                new_table[idx].value = bucket.value;
                break;
            }
            j++;
        }
    }

    char * create_new_buffer(size_t min_bytes) {

        //do not allocate a new address for the buffer while there is an active moving in progress
        bool expected = false;
        moving_buffer.compare_exchange_strong(expected, true, std::memory_order_acq_rel);
        if(expected){
            return nullptr;
        }

        assert(min_bytes>=buff_size);

        size_t new_buff_size = std::max<size_t>(buff_size, 1);
        while(new_buff_size<min_bytes){
            new_buff_size<<=1UL;
        }

        //auto *new_buff_ptr = (char *)malloc(new_buff_size);
        auto *new_buff_ptr = mem<char>::allocate(new_buff_size);
        memset(new_buff_ptr, 0, new_buff_size);

        prev_buff_ptr = buffer;
        buffer = new_buff_ptr;

        prev_buff_size = buff_size;
        buff_size = new_buff_size;

        return new_buff_ptr;
    }

    inline void move_old_buffer() {

        char *old_buff_ptr = prev_buff_ptr;
        char *new_buff_ptr = buffer;
        size_t old_size = prev_buff_size;

        //this will synchronize all the threads that were writing in the old buffer
        for(auto & prot_ptr : prot_ptrs) {
            while(prot_ptr.access.test_and_set(std::memory_order_acquire));
            prot_ptr.access.clear(std::memory_order_release);
        }

        // wait until no thread is writing into the old buffer
        if(old_buff_ptr!= nullptr) {

            auto reg_start = reinterpret_cast<uintptr_t>(old_buff_ptr);
            auto reg_end = reinterpret_cast<uintptr_t>(old_buff_ptr+old_size);

            //the hash table cannot change (i.e., rehash) while moving the keys from one buffer to the other
            while(mod_ht.test_and_set(std::memory_order_acquire));
            auto table_ptr = m_table;
            size_t cont=0;
            for(auto & bucket : *table_ptr) {
                //while(bucket.access.test_and_set(std::memory_order_acquire));
                bucket.access.lock();
                auto bck_int_ptr = reinterpret_cast<uintptr_t>(bucket.key_ptr);
                if(bck_int_ptr>=reg_start && bck_int_ptr<reg_end) {
                    ptrdiff_t diff = bucket.key_ptr-old_buff_ptr;
                    strcpy(new_buff_ptr+diff, bucket.key_ptr);
                    bucket.key_ptr =  new_buff_ptr + diff;
                    cont++;
                }
                //bucket.access.clear(std::memory_order_release);
                bucket.access.unlock();
            }
            assert(cont>0);
            mod_ht.clear(std::memory_order_release);

            //free(old_buff_ptr);
            mem<char>::deallocate(old_buff_ptr);
            //std::cout<<"Deleting the old key buffer"<<std::endl;
        }
        moving_buffer.store(false, std::memory_order_release);
    }

public:

    std::pair<char *, uint32_t> get_key_pointer(size_t idx) const {
        return {(*m_table)[idx].key_ptr, (*m_table)[idx].key_len};
    }

    [[nodiscard]] val_type * get_value_pointer(const size_t idx) const {
        return &((*m_table)[idx].value);
    }

    table_t * const& table = m_table;

    explicit ts_string_submap(size_t min_cap=4, float max_lf=0.6, size_t m_buff_size=0) : m_max_load_factor(max_lf),
                                                                                          buff_size(m_buff_size) {

        m_table = new table_t(round_to_power_of_two(min_cap));

        if(buff_size!=0){
            //buffer = (char *) malloc(buff_size);
            buffer = mem<char>::allocate(buff_size);
        }

        for(auto & prot_tab_ptr : prot_ptrs){
            prot_tab_ptr.tab_ptr = m_table;
        }
    }

    ~ts_string_submap(){
        delete m_table;
        //free(buffer);
        mem<char>::deallocate(buffer);
    }

    //thread-safe increment by "val" of the value associated with "key" in the hash table
    bool value_add(const uint8_t* key, size_t len, val_type val, size_t hash) {

        bool inserted = false;
        bool success = false;

        restart:
        size_t j=0;

        //protect the hash pointer from being freed
        while(prot_ptrs[tr_id].access.test_and_set(std::memory_order_acquire)){
            YIELD_PROCESSOR;
        };

        table_t * table_ptr = prot_ptrs[tr_id].tab_ptr;
        size_t idx = hash & (table_ptr->size()-1);

        while(!success) {
            //wait for the sm_bucket to be available and lock it
            bucket_t &bucket = (*table_ptr)[idx];
            //while(bucket.access.test_and_set(std::memory_order_acquire));
            bucket.access.lock();

            if(bucket.key_ptr==nullptr) {

                size_t key_bytes = len+1;
                char * buff_ptr;
                bool trigger_buff_move = false;
                size_t t_n_elms;

                //find a position in the buffer to insert the new key
                // only one thread at a time can execute this segment of code
                //while(buffer_lock.test_and_set(std::memory_order_acquire)){};
                dat_buff_lock.lock();

                size_t offset = next_av_offset;
                size_t min_bytes =  next_av_offset+key_bytes;

                if(min_bytes>=buff_size){//the key does not fit the buffer, allocate a new buffer
                    buff_ptr = create_new_buffer(min_bytes);
                    trigger_buff_move = true;
                    if(buff_ptr == nullptr){//there is a move of an old buffer in progress, so we couldn't allocate a new buffer
                        //buffer_lock.clear(std::memory_order_release);
                        dat_buff_lock.unlock();
                        //bucket.access.clear(std::memory_order_release);
                        bucket.access.unlock();
                        prot_ptrs[tr_id].access.clear(std::memory_order_release);
                        //TODO sleep for a moment
                        goto restart;
                    }
                }else{
                    buff_ptr = buffer;
                }

                next_av_offset += key_bytes;
                t_n_elms = ++n_elms;
                //buffer_lock.clear(std::memory_order_release);
                dat_buff_lock.unlock();
                //

                bucket.key_ptr = buff_ptr+offset;
                memcpy(bucket.key_ptr, key, len);

                bucket.value = val;
                //bucket.access.clear(std::memory_order_release);
                bucket.access.unlock();
                prot_ptrs[tr_id].access.clear(std::memory_order_release);

                //the key insertion triggers the creation of a new buffer.
                //move the keys from the old buffer to the new one
                if(trigger_buff_move){
                    move_old_buffer();
                }

                //the key insertion exceeds the max. load factor (i.e., rehash)
                if((float(t_n_elms)/float(table_ptr->size()))>=m_max_load_factor) {
                    rehash(next_power_of_two(table_ptr->size()), table_ptr);
                }

                success = true;
                inserted = true;

            } else if(bucket.key_ptr==invalid_addr) { // rehashing in progress (try again)
                //bucket.access.clear(std::memory_order_release);
                bucket.access.unlock();
                prot_ptrs[tr_id].access.clear(std::memory_order_release);
                //TODO sleep for a moment
                goto restart;
            } else {
                if(memcmp(key, bucket.key_ptr, len)==0 && *(bucket.key_ptr+len) == 0){
                    bucket.value += val;
                    //bucket.access.clear(std::memory_order_release);
                    bucket.access.unlock();
                    prot_ptrs[tr_id].access.clear(std::memory_order_release);
                    success = true;
                } else {
                    //bucket.access.clear(std::memory_order_release);
                    bucket.access.unlock();
                    j++;
                    idx = (hash + ((j*j + j)>>1UL)) & (table_ptr->size()-1);
                }
            }
        }
        return inserted;
    }

    void set_min_capacity(size_t new_cap) {
        new_cap = round_to_power_of_two(new_cap);
        auto table_ptr = m_table;
        if(new_cap>table_ptr->size()){
            rehash(new_cap, table_ptr);
        }
    }

    //warning: this function is not thread safe
    bool find(const uint8_t * key, size_t len, size_t hash, val_type& value) const {
        auto table_ptr = m_table;

        size_t j=0;
        size_t idx = hash & (table_ptr->size()-1);
        while(true){
            if((*table_ptr)[idx].key_ptr== nullptr){
                return false;
            }else if(memcmp((*table_ptr)[idx].key_ptr, key, len)==0 &&
                     *((*table_ptr)[idx].key_ptr + len)==0) {
                value = (*table_ptr)[idx].value;
                return true;
            }
            j++;
            idx = (hash + ((j*j + j)>>1UL)) & (table_ptr->size()-1);
        }
    }

    inline float load_factor() const {
        return float(n_elms)/float(m_table->size());
    }

    inline float max_load_factor() const {
        return m_max_load_factor;
    };

    inline size_t size() {
        dat_buff_lock.lock();
        //while(buffer_lock.test_and_set(std::memory_order_acquire));
        size_t res = n_elms;
        //buffer_lock.clear(std::memory_order_release);
        dat_buff_lock.unlock();

        return res;
    }

    inline bool empty()  {
        dat_buff_lock.lock();
        //while(buffer_lock.test_and_set(std::memory_order_acquire));
        bool res = n_elms==0;
        //buffer_lock.clear(std::memory_order_release);
        dat_buff_lock.unlock();
        return res;
    }

    inline size_t capacity() const{
        return m_table->size();
    };

    size_t buffer_usage(){
        dat_buff_lock.lock();
        //while(buffer_lock.test_and_set(std::memory_order_acquire));
        size_t usg = next_av_offset;
        //buffer_lock.clear(std::memory_order_release);
        dat_buff_lock.unlock();
        return usg;
    }

    size_t table_mem_usage(){
        size_t usg = m_table->size();
        return usg*sizeof(bucket_t);
    }

    size_t buffer_capacity(){
        dat_buff_lock.lock();
        //while(buffer_lock.test_and_set(std::memory_order_acquire));
        size_t bz = buff_size;
        //buffer_lock.clear(std::memory_order_release);
        dat_buff_lock.unlock();
        return bz;
    }

    void destroy_table(){
    }

    void destroy_buffer(){
    }

    void shrink_to_fit(){
    }

    inline void print_table() {
        dat_buff_lock.lock();
        //while(buffer_lock.test_and_set(std::memory_order_acquire));
        while(doing_rehash.test_and_set(std::memory_order_acquire));

        size_t i=0;
        for(auto & bucket : *m_table){
            if(bucket.key_ptr!=nullptr){
                printf("%lu-%lu) %s -> %lu\n", i++, n_elms, bucket.key_ptr, bucket.value);
            }
        }

        //buffer_lock.clear(std::memory_order_release);
        doing_rehash.clear(std::memory_order_release);
        dat_buff_lock.unlock();
    }
};

template<typename val_type,
         class sub_map_t,
         bool partition_key>
class string_map {

    typedef sub_map_t sub_map_type;
    typedef string_map_iterator<string_map<val_type, sub_map_t, partition_key>> iterator_type;
    friend class string_map_iterator<string_map<val_type, sub_map_t, partition_key>>;

    std::atomic<size_t> id_counter{0};
    std::vector<sub_map_type *> sub_map_array;
    float max_lf=0.6;
    size_t n_sub_maps{};
    size_t sm_mod{};
    size_t sm_shift{};

public:

    typedef typename sub_map_t::bucket_t bucket_t;
    typedef val_type                     value_type;

    explicit string_map(size_t min_sm_cap=4, float _max_lf=0.6, size_t _n_submaps=4, size_t min_sm_bz=0): max_lf(_max_lf),
                                                                                                          n_sub_maps(_n_submaps),
                                                                                                          sm_mod(n_sub_maps-1),
                                                                                                          sm_shift(__builtin_ctzll(n_sub_maps)){

        assert(!(n_sub_maps & (n_sub_maps-1)));
        sub_map_array.resize(n_sub_maps);
        if constexpr (!partition_key){
            static_assert(std::is_same<sub_map_t, ts_string_submap<value_type>>::value);
        }

        for(size_t i=0;i<n_sub_maps;i++){
            sub_map_array[i] = new sub_map_type(min_sm_cap, _max_lf, min_sm_bz);
        }
    }

    ~string_map(){
        for(size_t i=0;i<n_sub_maps;i++){
            delete sub_map_array[i];
        }
    }

    [[nodiscard]] inline float load_factor() const {
        return float(size())/float(capacity());
    }

    void register_thread(){
        size_t id = id_counter.fetch_add(1, std::memory_order_acq_rel);
        tr_id = id;
    }

    [[nodiscard]] inline float max_load_factor() const {
        return max_lf;
    };

    inline bool value_add(const uint8_t * key, size_t len, size_t new_val){
        size_t hash = XXH3_64bits(key, len);
        size_t j = (hash ^ (hash >> sm_shift)) & sm_mod;
        if constexpr (partition_key){
            if(j==tr_id){
                sub_map_array[j]->value_add(key, len, new_val, hash);
                return true;
            }else{
                return false;
            }
        } else{
            return sub_map_array[j]->value_add(key, len, new_val, hash);
        }
    }

    inline std::pair<proxy<value_type>, bool> insert(const uint8_t * key, size_t len, size_t new_val){
        size_t hash = XXH3_64bits(key, len);
        size_t j = (hash ^ (hash >> sm_shift)) & sm_mod;
        if constexpr (partition_key){
            if(j==tr_id){
                return sub_map_array[j]->insert(key, len, new_val, hash);
            }else{
                return {{nullptr, 0}, false};
            }
        } else{
            return sub_map_array[j]->insert(key, len, new_val, hash);
        }
    }

    inline void value_add_submap(const uint8_t * key, size_t len, size_t new_val){
        size_t hash = XXH3_64bits(key, len);
        size_t j = (hash ^ (hash >> sm_shift)) & sm_mod;
        if(j!=tr_id) return;
        sub_map_array[j]->value_add(key, len, new_val, hash);
    }

    bool find(const uint8_t * key, size_t len, val_type& val){
        size_t hash = XXH3_64bits(key, len);
        size_t j = (hash ^ (hash >> sm_shift)) & sm_mod;
        return sub_map_array[j]->find(key, len, hash, val);
    }

    [[nodiscard]] inline size_t size() const {
        size_t sz=0;
        for(size_t i=0;i<n_sub_maps;i++){
            sz+=sub_map_array[i]->size();
        }
        return sz;
    }

    [[nodiscard]] inline bool empty() const {
        bool em = false;
        for(size_t i=0;i<n_sub_maps;i++){
            em |= sub_map_array[i]->empty();
        }
        return !em;
    }

    [[nodiscard]] inline size_t capacity() const {
        size_t cap=0;
        for(size_t i=0;i<n_sub_maps;i++){
            cap+=sub_map_array[i]->capacity();
        }
        return cap;
    };

    size_t buffer_usage() {
        size_t bu=0;
        for(size_t i=0;i<n_sub_maps;i++){
            bu+=sub_map_array[i]->buffer_usage();
        }
        return bu;
    }

    size_t table_mem_usage(){
        size_t mu=0;
        for(size_t i=0;i<n_sub_maps;i++){
            mu+=sub_map_array[i]->table_mem_usage();
        }
        return mu;
    }

    size_t buffer_capacity(){
        size_t bc=0;
        for(size_t i=0;i<n_sub_maps;i++){
            bc+=sub_map_array[i]->buffer_capacity();
        }
        return bc;
    }

    iterator_type begin(){
        return iterator_type(*this, 0);
    }

    iterator_type end(){
        return iterator_type(*this, size());
    }

    void destroy_tables(){
        for(size_t i=0;i<n_sub_maps;i++){
            sub_map_array[i]->destroy_table();
        }
    }

    void destroy_buffers(){
        for(size_t i=0;i<n_sub_maps;i++){
            sub_map_array[i]->destroy_buffer();
        }
    }

    void destroy_map() {
        destroy_buffers();
        destroy_tables();
    }

    void shrink_to_fit(){
        for(size_t i=0;i<n_sub_maps;i++){
            sub_map_array[i]->shrink_to_fit();
        }
    }

    void load_tables(const std::string& input_file){
        std::ifstream ifs(input_file, std::ios::in | std::ios::binary);
        ifs.read((char *)&n_sub_maps, sizeof(n_sub_maps));
        sub_map_array.resize(n_sub_maps);
        for(size_t i=0;i<n_sub_maps;i++){
            sub_map_array[i]->load_table(ifs);
        }
        ifs.close();
    }

    void load_buffers(const std::string& input_file){
        std::ifstream ifs(input_file, std::ios::in | std::ios::binary);
        ifs.read((char *)&n_sub_maps, sizeof(n_sub_maps));
        sub_map_array.resize(n_sub_maps);
        for(size_t i=0;i<n_sub_maps;i++){
            sub_map_array[i]->load_buffer(ifs);
        }
        ifs.close();
    }

    void load_map(std::string& input_file){
        std::ifstream ifs(input_file, std::ios::in | std::ios::binary);
        ifs.read((char *)&n_sub_maps, sizeof(n_sub_maps));
        sub_map_array.resize(n_sub_maps);
        for(size_t i=0;i<n_sub_maps;i++){
            sub_map_array[i]->load_submap(ifs);
        }
        ifs.read((char *)&max_lf, sizeof(max_lf));
        sm_mod = n_sub_maps-1;
        sm_shift = __builtin_ctzll(n_sub_maps);
        ifs.close();
    }

    void store_tables(const std::string& output_file){
        std::ofstream ofs(output_file, std::ios::out | std::ios::binary);
        ofs.write((char *)&n_sub_maps, sizeof(n_sub_maps));
        for(size_t i=0;i<n_sub_maps;i++){
            sub_map_array[i]->store_table(ofs);
        }
        ofs.close();
    }

    void store_buffers(const std::string& output_file){
        std::ofstream ofs(output_file, std::ios::out | std::ios::binary);
        ofs.write((char *)&n_sub_maps, sizeof(n_sub_maps));
        for(size_t i=0;i<n_sub_maps;i++){
            sub_map_array[i]->store_buffer(ofs);
        }
        ofs.close();
    }

    void store_map(const std::string& output_file){
        std::ofstream ofs(output_file, std::ios::out | std::ios::binary);
        ofs.write((char *)&n_sub_maps, sizeof(n_sub_maps));
        for(size_t i=0;i<n_sub_maps;i++){
            sub_map_array[i]->store_submap(ofs);
        }
        ofs.write((char *)&max_lf, sizeof(max_lf));
        ofs.close();
    }

    void print_table(){
        for(size_t i=0;i<n_sub_maps;i++){
            sub_map_array[i]->print_table();
        }
    }

    void print_table_stats(){
        std::cout<<"Map stats"<<std::endl;
        std::cout<<"  Number of entries : "<<size()<<std::endl;
        std::cout<<"  Load factor : "<<load_factor()<<std::endl;
        std::cout<<"  Table capacity : "<<capacity()<<std::endl;
        std::cout<<"  Table usage (bytes) : "<<table_mem_usage()<<std::endl;
        std::cout<<"  Key buffer usage (bytes) : "<<buffer_usage()<<std::endl;
        std::cout<<"  Key buffer capacity (bytes) : "<<buffer_capacity()<<std::endl;
        report_submap_dist();
    }

    void report_submap_dist(){
        for(size_t i=0;i<n_sub_maps;i++){
            std::cout<<"Submap "<<i<<" prop="<<float(sub_map_array[i]->size())/float(size())<<std::endl;
        }
    }
};

//template<typename val>
//using ts_string_map = string_map<val, ts_string_submap<val>, false>;

template<typename val>
using par_string_map = string_map<val, nts_string_submap<val>, true>;

#endif //SIMPLE_EXAMPLE_TS_STRING_MAP_H
