//
// Created by Diaz, Diego on 9.1.2025.
//

#ifndef LCG_PHRASE_SET_VBYTE_H
#define LCG_PHRASE_SET_VBYTE_H

#include "xxhash.h"
#include "cds/cdt_common.hpp"
#include "cds/vbyte.h"
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cassert>
#ifdef __linux__
#include <malloc.h>
#endif

//20 MB
#define STR_THRESHOLD1 20971520

//200 MB
#define STR_THRESHOLD2 209715200

//1 GB
#define STR_THRESHOLD3 1073741824

class phrase_set_vbyte {
private:

    typedef std::vector<uint64_t> table_t;
    static constexpr uint64_t null_addr = std::numeric_limits<uint64_t>::max();
    static constexpr uint64_t d_one = 1UL<<44UL;

    uint8_t *phrase_stream = nullptr;
    uint64_t stream_size=0;
    uint64_t stream_cap=0;
    table_t m_table;
    float m_max_load_factor = 0.6;
    uint64_t elm_threshold=0;
    uint64_t frac_lf = 60;
    uint64_t n_phrases=0;
    uint64_t last_mt=0;
    uint64_t last_fp_pos=0;

    void rehash(uint64_t new_tab_size) {

        assert(new_tab_size>m_table.size());
        m_table.resize(new_tab_size);
        memset(m_table.data(), (int)null_addr, m_table.size()*sizeof(table_t::value_type));

        //rehash the values
        uint64_t phr_addr=0, phr_bytes, phr_len_bytes;
        uint64_t proc_phrases=0;
        while(phr_addr<stream_size){
            phr_bytes=0;
            phr_len_bytes = vbyte::read(&phrase_stream[phr_addr], phr_bytes);
            rehash_entry_rb(XXH3_64bits(&phrase_stream[phr_addr+phr_len_bytes], phr_bytes), phr_addr);
            phr_addr+=phr_len_bytes+phr_bytes;
            //read the metasymbol and move forward
            phr_addr+=vbyte::read(&phrase_stream[phr_addr], phr_bytes);
            proc_phrases++;
            assert(proc_phrases<=n_phrases);
        }
        assert(proc_phrases==n_phrases);
        assert(phr_addr==stream_size);
        elm_threshold = (m_table.size()*frac_lf)/100;
    }

    void increase_stream_cap(uint64_t min_cap){
        if(min_cap <= STR_THRESHOLD1){
            stream_cap = stream_size*2;//100% = 2x increase
        } else if(min_cap <= STR_THRESHOLD2){
            stream_cap = INT_CEIL((stream_size*15), 10);//%50=1.5x increase capacity
        } else if(min_cap <= STR_THRESHOLD3){
            stream_cap = INT_CEIL((stream_size*12), 10);//20%=1.2x increase capacity
        } else {
            stream_cap = INT_CEIL((stream_size*105), 100);//%5=1.05x increase capacity
        }
        stream_cap = std::max(min_cap, stream_cap);

        if(phrase_stream== nullptr){
            phrase_stream = mem<uint8_t>::allocate(stream_cap);
        }else{
            phrase_stream = mem<uint8_t>::reallocate(phrase_stream, stream_cap);
        }
    }

public:

    struct phrase_t{
        uint8_t* phrase= nullptr;
        uint64_t len=0;
        uint64_t mt=0;

        inline uint64_t parse(uint8_t* new_addr) {
            uint64_t read_bytes=0;
            read_bytes+= vbyte::read(new_addr, len);
            phrase=&new_addr[read_bytes];
            read_bytes+=len;
            read_bytes+= vbyte::read(&new_addr[read_bytes], mt);
            return read_bytes;
        }

        struct phrase_iterator{
            const uint8_t *addr;
            uint64_t sym = 0;
            uint64_t pos = 0;
            uint64_t size = 0;

            phrase_iterator(const uint8_t* _addr, uint64_t start_pos, uint64_t _stream_size): addr(_addr),
                                                                                              pos(start_pos),
                                                                                              size(_stream_size) {
                if(pos<size){
                    pos+=vbyte::read(&addr[pos], sym);
                }else{
                    pos = size+1;
                }
            }

            // Move to the next tuple
            inline void operator++() {
                sym = 0;
                if(pos<size){
                    pos+=vbyte::read(&addr[pos], sym);
                }else{
                    pos = size+1;
                }
            }

            inline bool operator==(const phrase_iterator& other) const {
                return other.addr==addr && other.pos==pos;
            }

            inline bool operator!=(const phrase_iterator& other) const {
                return other.pos!=pos || other.addr!=addr;
            }
        };

        [[nodiscard]] phrase_iterator begin() const {
            return {phrase, 0, len};
        }

        [[nodiscard]] phrase_iterator begin(uint64_t stream_pos) const {
            assert(stream_pos>=0 && stream_pos<len);
            return {phrase, stream_pos, len};
        }

        [[nodiscard]] phrase_iterator end() const {
            return {phrase, len, len};
        }
    };

    //iterate over the phrases of the set
    struct set_iterator {

    private:
        uint8_t* stream;
        uint64_t stream_pos=0;
        uint64_t stream_size;
        uint64_t curr_phr_pos = 0;
        phrase_t curr_phrase;

        void decode_phrase(){
            curr_phr_pos = stream_pos;
            if (stream_pos<stream_size){
                stream_pos+=curr_phrase.parse(&stream[stream_pos]);
                assert(stream_pos<=stream_size);
            }else{
                stream_pos = stream_size+1;
                curr_phrase.len = 0;
                curr_phrase.mt = 0;
                curr_phrase.phrase = nullptr;
            }
        };

    public:

        explicit set_iterator(uint8_t* _stream, uint64_t _start_pos, uint64_t _stream_size) : stream(_stream),
                                                                                              stream_pos(_start_pos),
                                                                                              stream_size(_stream_size){
            decode_phrase();
        }

        // Move to the next tuple
        inline void operator++() {
            decode_phrase();
        }

        inline phrase_t& operator*(){
            return curr_phrase;
        }

        [[nodiscard]] inline uint64_t pos() const {
            return curr_phr_pos;
        }

        inline bool operator==(const set_iterator& other) const {
            return other.stream==stream && other.stream_pos==stream_pos;
        }

        inline bool operator!=(const set_iterator& other) const {
            return other.stream_pos!=stream_pos || other.stream!=stream;
        }
    };

    explicit phrase_set_vbyte(uint64_t min_cap=4, float max_lf=0.85) {
        assert(min_cap>0);
        m_max_load_factor = max_lf;
        m_table = table_t(round_to_power_of_two(min_cap), null_addr);
        frac_lf = uint64_t(m_max_load_factor*100);
        elm_threshold = (m_table.size()*frac_lf)/100;
    }

    phrase_set_vbyte(phrase_set_vbyte&& other) noexcept {
        std::swap(phrase_stream, other.phrase_stream);
        std::swap(stream_size, other.stream_size);
        std::swap(stream_cap, other.stream_cap);
        m_table.swap(other.m_table);
        std::swap(m_max_load_factor, other.m_max_load_factor);
        std::swap(elm_threshold, other.elm_threshold);
        std::swap(frac_lf, other.frac_lf);
        std::swap(n_phrases, other.n_phrases);
        std::swap(last_mt, other.last_mt);
        std::swap(last_fp_pos, other.last_fp_pos);
    }

    phrase_set_vbyte(const phrase_set_vbyte& other) noexcept {
        copy(other);
    }

    void set_load_factor(float new_max_lf){
        m_max_load_factor = new_max_lf;
        frac_lf = uint64_t(m_max_load_factor*100);
        elm_threshold = (m_table.size()*frac_lf)/100;
        if(n_phrases>=elm_threshold) {
            rehash(next_power_of_two(m_table.size()));
        }
    }

    void copy(const phrase_set_vbyte& other){
        stream_size = other.stream_size;
        stream_cap = other.stream_cap;
        m_table = other.m_table;
        m_max_load_factor = other.m_max_load_factor;
        elm_threshold = other.elm_threshold;
        frac_lf = other.frac_lf;
        n_phrases = other.n_phrases;
        last_mt = other.last_mt;
        last_fp_pos = other.last_fp_pos;

        phrase_stream = mem<uint8_t>::allocate(stream_cap);
        memcpy(phrase_stream, other.phrase_stream, stream_size);
    }

    void swap(phrase_set_vbyte& other){
        std::swap(phrase_stream, other.phrase_stream);
        std::swap(stream_size, other.stream_size);
        std::swap(stream_cap, other.stream_cap);
        m_table.swap(other.m_table);
        std::swap(m_max_load_factor, other.m_max_load_factor);
        std::swap(elm_threshold, other.elm_threshold);
        std::swap(frac_lf, other.frac_lf);
        std::swap(n_phrases, other.n_phrases);
        std::swap(last_mt, other.last_mt);
        std::swap(last_fp_pos, other.last_fp_pos);
    }

    phrase_set_vbyte& operator=(const phrase_set_vbyte& other){
        if(this!=&other){
            copy(other);
        }
        return *this;
    }

    inline void rehash_entry_rb(uint64_t hash, uint64_t phr_addr){
        uint64_t idx = hash & (m_table.size() - 1);
        if(m_table[idx]!=null_addr) {
            uint64_t dist=0, bck_dist;
            bck_dist = m_table[idx] >> 44UL;
            while(bck_dist>=dist && m_table[idx]!=null_addr){
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }
            assert(dist<=65535);
            //assert(bck_dist<dist || m_table[idx]==null_addr);
            phr_addr |= dist<<44UL;
            uint64_t tmp;
            while(m_table[idx]!=null_addr){
                tmp = m_table[idx];
                m_table[idx] = phr_addr;
                phr_addr = tmp+d_one;
                idx = (idx+1) & (m_table.size()-1);
            }
        }
        m_table[idx] = phr_addr ;
    }

    //this function is for when we already know the phrase is *not* in the set
    inline uint64_t add_phrase(const uint8_t* phrase, uint64_t phr_bytes){

        uint64_t hash = XXH3_64bits(phrase, phr_bytes);
        uint64_t mt;

        uint64_t idx = hash & (m_table.size()-1);
        uint64_t bck_val = stream_size;
        if(m_table[idx]!=null_addr) {
            uint64_t dist=0, bck_dist;
            bck_dist = m_table[idx] >> 44UL;
            while(bck_dist>=dist && m_table[idx]!=null_addr){
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }
            assert(dist<=65535);
            //assert(bck_dist<dist || m_table[idx]==null_addr);
            bck_val |= dist<<44UL;
            uint64_t tmp;
            while(m_table[idx]!=null_addr){
                tmp = m_table[idx];
                m_table[idx] = bck_val;
                bck_val = tmp+d_one;
                idx = (idx+1) & (m_table.size()-1);
            }
        }

        m_table[idx] = bck_val;
        mt = n_phrases++;

        uint64_t mt_vb_sz = vbyte_len(mt);
        uint64_t len_vb_sz = vbyte_len(phr_bytes);
        uint64_t new_stream_size = stream_size+mt_vb_sz+len_vb_sz+phr_bytes;
        if(new_stream_size>=stream_cap){
            increase_stream_cap(new_stream_size);
        }

        vbyte::write(&phrase_stream[stream_size], phr_bytes, len_vb_sz);
        stream_size+=len_vb_sz;
        memcpy(&phrase_stream[stream_size], phrase, phr_bytes);
        stream_size+=phr_bytes;
        vbyte::write(&phrase_stream[stream_size], mt, mt_vb_sz);
        stream_size+=mt_vb_sz;

        //the insertion exceeds the max. load factor (i.e., rehash)
        if(n_phrases>=elm_threshold) {
            rehash(next_power_of_two(m_table.size()));
        }

        return mt;
    }

    inline uint64_t insert(const uint8_t* phrase, uint64_t phr_bytes) {

        uint64_t hash = XXH3_64bits(phrase, phr_bytes);
        uint64_t mt;

        uint64_t idx = hash & (m_table.size()-1);
        uint64_t bck_val = stream_size;
        if(m_table[idx]!=null_addr) {
            uint64_t bck_dist = m_table[idx] >> 44UL, dist=0;
            while(bck_dist>dist && m_table[idx]!=null_addr){
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }

            while(dist==bck_dist){
                uint64_t bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
                uint64_t ht_phr_bytes;
                uint64_t n_bytes = vbyte::read(&phrase_stream[bck_addr], ht_phr_bytes);
                if(phr_bytes==ht_phr_bytes && memcmp(&phrase_stream[bck_addr+n_bytes], phrase, phr_bytes)==0){
                    vbyte::read(&phrase_stream[bck_addr+n_bytes+phr_bytes], mt);
                    return mt;
                }
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }

            assert(dist<=65535);
            //assert(bck_dist<dist || m_table[idx]==null_addr);

            bck_val |= (dist<<44UL);
            uint64_t tmp;
            while(m_table[idx]!=null_addr){
                tmp = m_table[idx];
                m_table[idx] = bck_val;
                bck_val = tmp+d_one;
                idx = (idx+1) & (m_table.size()-1);
            }
        }
        m_table[idx] = bck_val;

        mt = n_phrases++;

        uint64_t mt_vb_sz = vbyte_len(mt);
        uint64_t len_vb_sz = vbyte_len(phr_bytes);
        uint64_t new_stream_size = stream_size+mt_vb_sz+len_vb_sz+phr_bytes;
        if(new_stream_size>=stream_cap){
            increase_stream_cap(new_stream_size);
        }

        vbyte::write(&phrase_stream[stream_size], phr_bytes, len_vb_sz);
        stream_size+=len_vb_sz;
        memcpy(&phrase_stream[stream_size], phrase, phr_bytes);
        stream_size+=phr_bytes;
        vbyte::write(&phrase_stream[stream_size], mt, mt_vb_sz);
        stream_size+=mt_vb_sz;

        //the insertion exceeds the max. load factor (i.e., rehash)
        if(n_phrases>=elm_threshold) {
            rehash(next_power_of_two(m_table.size()));
        }
        return mt;
    }

    inline uint32_t insert(const uint8_t* phrase, const uint64_t phr_bytes, const uint64_t hash) {

        uint64_t mt;
        uint64_t idx = hash & (m_table.size()-1);
        uint64_t bck_val = stream_size;
        if(m_table[idx]!=null_addr) {
            uint64_t bck_dist = m_table[idx] >> 44UL, dist=0;
            while(bck_dist>dist && m_table[idx]!=null_addr){
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }

            while(dist==bck_dist){
                uint64_t bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
                uint64_t ht_phr_bytes;
                uint64_t n_bytes = vbyte::read(&phrase_stream[bck_addr], ht_phr_bytes);
                if(phr_bytes==ht_phr_bytes && memcmp(&phrase_stream[bck_addr+n_bytes], phrase, phr_bytes)==0){
                    vbyte::read(&phrase_stream[bck_addr+n_bytes+phr_bytes], mt);
                    return mt;
                }
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }

            assert(dist<=65535);
            //assert(bck_dist<dist || m_table[idx]==null_addr);

            bck_val |= (dist<<44UL);
            uint64_t tmp;
            while(m_table[idx]!=null_addr){
                tmp = m_table[idx];
                m_table[idx] = bck_val;
                bck_val = tmp+d_one;
                idx = (idx+1) & (m_table.size()-1);
            }
        }
        m_table[idx] = bck_val;

        mt = n_phrases++;

        uint64_t mt_vb_sz = vbyte_len(mt);
        uint64_t len_vb_sz = vbyte_len(phr_bytes);
        uint64_t new_stream_size = stream_size+mt_vb_sz+len_vb_sz+phr_bytes;
        if(new_stream_size>=stream_cap){
            increase_stream_cap(new_stream_size);
        }

        vbyte::write(&phrase_stream[stream_size], phr_bytes, len_vb_sz);
        stream_size+=len_vb_sz;
        memcpy(&phrase_stream[stream_size], phrase, phr_bytes);
        stream_size+=phr_bytes;
        vbyte::write(&phrase_stream[stream_size], mt, mt_vb_sz);
        stream_size+=mt_vb_sz;

        //the insertion exceeds the max. load factor (i.e., rehash)
        if(n_phrases>=elm_threshold) {
            rehash(next_power_of_two(m_table.size()));
        }
        return mt;
    }

    inline bool find(const uint8_t* phrase, uint64_t phr_bytes, uint64_t& mt) const {

        uint64_t hash = XXH3_64bits(phrase, phr_bytes);
        uint64_t idx = hash & (m_table.size()-1);

        if(m_table[idx]!=null_addr) {
            uint64_t bck_dist = m_table[idx] >> 44UL, dist=0;
            while(bck_dist>dist && m_table[idx]!=null_addr){
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }

            while(dist==bck_dist){
                uint64_t bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
                uint64_t ht_phr_bytes;
                uint64_t n_bytes = vbyte::read(&phrase_stream[bck_addr], ht_phr_bytes);
                if(phr_bytes==ht_phr_bytes && memcmp(&phrase_stream[bck_addr+n_bytes], phrase, phr_bytes)==0){
                    vbyte::read(&phrase_stream[bck_addr+n_bytes+phr_bytes], mt);
                    return true;
                }

                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }
        }
        return false;
    }

    inline bool find(const uint8_t* phrase, const uint64_t phr_bytes, uint64_t& mt, const uint64_t hash) const {

        uint64_t idx = hash & (m_table.size()-1);
        if(m_table[idx]!=null_addr) {

            uint64_t bck_dist = m_table[idx] >> 44UL, dist=0;
            while(bck_dist>dist && m_table[idx]!=null_addr){
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }

            while(dist==bck_dist){
                uint64_t bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
                uint64_t ht_phr_bytes;
                uint64_t n_bytes = vbyte::read(&phrase_stream[bck_addr], ht_phr_bytes);
                if(phr_bytes==ht_phr_bytes && memcmp(&phrase_stream[bck_addr+n_bytes], phrase, phr_bytes)==0){
                    vbyte::read(&phrase_stream[bck_addr+n_bytes+phr_bytes], mt);
                    return true;
                }
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }
        }
        return false;
    }

    [[nodiscard]] inline float load_factor() const {
        return float(n_phrases)/float(m_table.size());
    }

    [[nodiscard]] inline float max_load_factor() const {
        return m_max_load_factor;
    };

    [[nodiscard]] inline uint64_t size() const {
        return n_phrases;
    }

    [[nodiscard]] inline bool empty()  const {
        return n_phrases==0;
    }

    [[nodiscard]] inline uint64_t capacity() const{
        return m_table.size();
    };

    [[nodiscard]] inline uint64_t stream_len() const {
        return stream_size;
    }

    [[nodiscard]] inline uint64_t table_mem_usage() const {
        return m_table.size()*sizeof(table_t::value_type);
    }

    [[nodiscard]] inline uint64_t phrases_mem_usage() const {
        return stream_cap;
    }

    [[nodiscard]] inline uint64_t mem_usage() const {
        return table_mem_usage()+phrases_mem_usage();
    }

    [[nodiscard]] inline uint64_t eff_mem_usage() const {
        return (n_phrases*sizeof(table_t::value_type)) + stream_size;
    }

    [[nodiscard]] inline const uint8_t* phr_stream() const {
        return phrase_stream;
    }

    inline void set_stream_capacity(uint64_t new_capacity){
        if(new_capacity>=stream_size){
            stream_cap = new_capacity;
            phrase_stream = mem<uint8_t>::reallocate(phrase_stream, stream_cap);
        }
    }

    void shrink_to_fit() {
        stream_cap = stream_size;
        phrase_stream = mem<uint8_t>::reallocate(phrase_stream, stream_cap);
    }

    [[nodiscard]] set_iterator begin() const {
        return set_iterator(phrase_stream, 0, stream_size);
    }

    [[nodiscard]] set_iterator begin(uint64_t stream_pos) const {
        assert(stream_pos>=0 && stream_pos<stream_size);
        return set_iterator(phrase_stream, stream_pos, stream_size);
    }

    [[nodiscard]] set_iterator end() const {
        return set_iterator(phrase_stream, stream_size, stream_size);
    }

    void destroy_table(){
        table_t().swap(m_table);
    }

    void destroy_stream(){
        if(phrase_stream!= nullptr){
            mem<uint8_t>::deallocate(phrase_stream);
        }
        phrase_stream= nullptr;
        stream_size = 0;
        stream_cap = 0;
        n_phrases = 0;
        last_fp_pos = 0;
        last_mt = 0;
    }

    void destroy(){
        destroy_table();
        destroy_stream();
    }

    [[nodiscard]] uint64_t buff_bytes_available() const {
        return (stream_cap-stream_size);
    }

    [[nodiscard]] uint64_t stream_capacity() const {
        return stream_cap;
    }

    ~phrase_set_vbyte(){
        destroy();
    }

    void psl_dist(){
        if(n_phrases>0){
            uint64_t freqs[2000]={0};
            for(unsigned long long entry : m_table){
                if(entry!=null_addr){
                    uint64_t bck_dist = entry >> 44UL;
                    if(bck_dist<2000){
                        freqs[bck_dist]++;
                    }
                }
            }
            uint64_t avg=0;
            float acc=0;
            for(uint64_t i=0;i<2000;i++){
                if(freqs[i]>0){
                    float frac = float(freqs[i])/float(n_phrases);
                    acc +=frac;
                    std::cout<<"psl: "<<i<<" --> "<<freqs[i]<<" "<<frac<<" "<<acc<<std::endl;
                    avg+=freqs[i]*i;
                }
            }
            std::cout<<"The average is "<<float(avg)/float(n_phrases)<<std::endl;
        }
    }

    /*[[nodiscard]] inline uint64_t tot_symbols() const {
        if constexpr (std::is_same<seq_type, uint8_t>::value){
            return stream_size - (n_phrases*sizeof(uint32_t)*2);
        }else{
            return stream_size - (n_phrases*2);
        }
    }*/

    uint64_t absorb_set(phrase_set_vbyte& other, const uint64_t* i_map,
                        uint64_t i_map_len, uint64_t* o_map, uint64_t o_map_len){

        if(other.empty()) return size();
        assert(o_map_len==(other.size()+1));

        auto set_it = other.begin();
        auto set_it_end = other.end();

        uint64_t new_seq_buff_sz = 2048;
        uint8_t *new_seq_buff = mem<uint8_t>::allocate(new_seq_buff_sz);
        uint64_t new_seq_bytes, proc_phrases=0, mt, sym_vb_sz;

        while(set_it!=set_it_end){
            auto phr = *set_it;
            auto phr_it = phr.begin();
            auto phr_it_end = phr.end();
            new_seq_bytes=0;
            while(phr_it!=phr_it_end){
                assert(phr_it.sym>0 && phr_it.sym<i_map_len);
                phr_it.sym = i_map[phr_it.sym]+1;
                sym_vb_sz = vbyte_len(phr_it.sym);
                if((new_seq_bytes+sym_vb_sz) > new_seq_buff_sz){
                    new_seq_buff_sz*=2;
                    new_seq_buff = mem<uint8_t>::reallocate(new_seq_buff, new_seq_buff_sz);
                }
                vbyte::write(&new_seq_buff[new_seq_bytes], phr_it.sym, sym_vb_sz);
                new_seq_bytes+=sym_vb_sz;
                ++phr_it;
            }
            mt = insert(new_seq_buff, new_seq_bytes);
            o_map[proc_phrases+1] = mt;
            proc_phrases++;
            assert(proc_phrases<=other.n_phrases);
            ++set_it;
        }

        assert(proc_phrases==other.n_phrases);
        mem<uint8_t>::deallocate(new_seq_buff);

        return size();
        /*uint64_t len, mt;
        uint64_t pos=0, proc_phrases=0;
        while(pos<other.stream_size){
            if constexpr (std::is_same<seq_type, uint8_t>::value){
                memcpy(&len, &other.phrase_stream[pos], sizeof(uint32_t));//read the length
                pos+=sizeof(uint32_t);//skip length
                mt = insert(&other.phrase_stream[pos], len);
                pos+=len+sizeof(uint32_t);//skip the phrase and mt
            }else{
                len = other.phrase_stream[pos];//read the length
                pos++;//skip length
                uint64_t end = pos+len;
                for(uint64_t j=pos;j<end;j++){
                    assert(other.phrase_stream[j]>0 && other.phrase_stream[j]<i_map_len);
                    other.phrase_stream[j] = i_map[other.phrase_stream[j]]+1;
                }
                mt = insert(&other.phrase_stream[pos], len);
                pos+=len+1;//skip the phrase and mt
            }
            o_map[proc_phrases+1] = mt;
            proc_phrases++;
            assert(proc_phrases<=other.n_phrases);
        }
        assert(proc_phrases==other.n_phrases);
        assert(pos==other.stream_size);
        return size();*/
    }

    void clear(){
        stream_size=0;
        stream_cap = std::min<uint64_t>(4, stream_cap>>3);//keep 1/8 of the buffer
        phrase_stream = mem<uint8_t>::reallocate(phrase_stream, stream_cap);
        uint64_t new_tab_size = std::max<uint64_t>(4, (m_table.size()>>3));// keep 1/8 of the table
        assert(is_power_of_two(new_tab_size));//dummy check
        m_table.resize(new_tab_size);
        memset(m_table.data(), 0xFF, m_table.size()*sizeof(table_t::value_type));
        elm_threshold = (m_table.size()*frac_lf)/100;
        n_phrases = 0;
        last_mt = 0;
        last_fp_pos =0;
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    inline void update_fps(const uint64_t* prev_fps, uint64_t& len_prev_fps, uint64_t*& fps, uint64_t& len_fps) {
        if(last_mt<n_phrases){
            assert(last_fp_pos<=stream_size);

            fps = mem<uint64_t>::reallocate(fps, n_phrases+1);
            len_fps = n_phrases+1;

            auto set_it = set_iterator(phrase_stream, last_fp_pos, stream_size);
            auto set_it_end = end();
            std::vector<uint64_t> fp_sequence;
            while(set_it!=set_it_end){
                auto phr = *set_it;

                /*for(uint64_t j=0;j<phr.len;j++){
                    assert(phr.phrase[j]>0 && phr.phrase[j]<len_prev_fps);
                    fp_sequence.push_back(prev_fps[phr.phrase[j]]);
                }*/
                auto phr_it = phr.begin();
                auto phr_it_end = phr.end();
                while(phr_it!=phr_it_end){
                    assert(phr_it.sym>0 && phr_it.sym<len_prev_fps);
                    fp_sequence.push_back(prev_fps[phr_it.sym]);
                    ++phr_it;
                }
                fps[last_mt+1] = XXH3_64bits(fp_sequence.data(), fp_sequence.size()*sizeof(uint64_t));
                last_mt++;
                ++set_it;
                fp_sequence.clear();
            }
            last_fp_pos=stream_size;
        }
    }

    inline void update_fps_with_sink(const uint64_t* prev_fps_sink, const uint64_t& len_prev_fps_sink,
                                     uint64_t alpha_sink, const uint64_t* prev_fps, uint64_t& len_prev_fps,
                                     uint64_t*& fps, uint64_t& len_fps) {

        const uint64_t *p_fps[2] = {prev_fps, prev_fps_sink};
        const uint64_t p_fps_len[2] = {len_prev_fps, len_prev_fps_sink};

        if(last_mt<n_phrases) {
            assert(last_fp_pos<=stream_size);

            fps = mem<uint64_t>::reallocate(fps, n_phrases+1);
            len_fps = n_phrases+1;

            auto set_it = set_iterator(phrase_stream, last_fp_pos, stream_size);
            auto set_it_end = end();
            std::vector<uint64_t> fp_sequence;

            while(set_it!=set_it_end){

                auto phr = *set_it;
                auto phr_it = phr.begin();
                auto phr_it_end = phr.end();

                while(phr_it!=phr_it_end){
                    if(phr_it.sym<=alpha_sink){
                        assert(phr_it.sym>0 && phr_it.sym<p_fps_len[1]);
                        fp_sequence.push_back(p_fps[1][phr_it.sym]);
                    }else{
                        uint64_t rank = phr_it.sym-alpha_sink;
                        assert(rank>0 && rank<p_fps_len[0]);
                        fp_sequence.push_back(p_fps[0][rank]);
                    }
                    ++phr_it;
                }

                /*for(uint64_t j=0;j<phr.len;j++){
                    if(phr.phrase[j]<=alpha_sink){
                        assert(phr.phrase[j]>0 && phr.phrase[j]<p_fps_len[1]);
                        fp_sequence.push_back(p_fps[1][phr.phrase[j]]);
                    }else{
                        uint32_t rank = phr.phrase[j]-alpha_sink;
                        assert(rank>0 && rank<p_fps_len[0]);
                        fp_sequence.push_back(p_fps[0][rank]);
                    }
                }*/

                fps[last_mt+1] = XXH3_64bits(fp_sequence.data(), fp_sequence.size()*sizeof(uint64_t));
                last_mt++;
                ++set_it;
                fp_sequence.clear();
            }
            last_fp_pos=stream_size;
        }
    }

    uint64_t serialize(std::ofstream &ofs){
        uint64_t written_bytes = 0;
        written_bytes+=serialize_elm(ofs, stream_size);
        written_bytes+=serialize_elm(ofs, stream_cap);
        written_bytes+=serialize_elm(ofs, m_max_load_factor);
        written_bytes+=serialize_elm(ofs, elm_threshold);
        written_bytes+=serialize_elm(ofs, frac_lf);
        written_bytes+=serialize_elm(ofs, n_phrases);
        written_bytes+=serialize_elm(ofs, last_mt);
        written_bytes+=serialize_elm(ofs, last_fp_pos);
        written_bytes+= serialize_raw_vector(ofs, phrase_stream, stream_size);
        written_bytes+= serialize_plain_vector(ofs, m_table);
        return written_bytes;
    }

    void load(std::ifstream & ifs){
        load_elm(ifs, stream_size);
        load_elm(ifs, stream_cap);
        load_elm(ifs, m_max_load_factor);
        load_elm(ifs, elm_threshold);
        load_elm(ifs, frac_lf);
        load_elm(ifs, n_phrases);
        load_elm(ifs, last_mt);
        load_elm(ifs, last_fp_pos);

        load_raw_vector(ifs, phrase_stream, stream_size);
        load_plain_vector(ifs, m_table);
        stream_cap = stream_size;
    }
};
#endif //LCG_PHRASE_SET_VBYTE_H
