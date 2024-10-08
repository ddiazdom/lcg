//
// Created by Diaz, Diego on 17.9.2023.
//

#ifndef LCG_LZL_MAP_H
#define LCG_LZL_MAP_H
#include "../../external/xxHash-dev/xxhash.h"
#include "cds/cdt_common.hpp"
#include "buff_vector.h"
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

template<class seq_type>
class phrase_set {

public:

    seq_type *phrase_stream = nullptr;
    size_t stream_size=0;
    size_t stream_cap=0;
    static constexpr uint8_t seq_bytes=sizeof(seq_type);
    static constexpr uint64_t null_addr = std::numeric_limits<uint64_t>::max();
    typedef std::vector<uint64_t> table_t;

    struct phrase_t{
        seq_type* phrase;
        uint32_t len;
        uint32_t mt;
        phrase_t(seq_type* _phrase, uint32_t _len, uint32_t _mt): phrase(_phrase), len(_len), mt(_mt){}
    };

    struct iterator {

    private:
        seq_type* stream;
        size_t stream_pos=0;
        size_t stream_size;
        phrase_t curr_phrase;

        void decode_phrase(){
            if (stream_pos<stream_size){
                if constexpr (std::is_same<seq_type, uint8_t>::value){
                    memcpy(&curr_phrase.len, &stream[stream_pos], sizeof(uint32_t));
                    stream_pos+=sizeof(uint32_t);
                    curr_phrase.phrase = &stream[stream_pos];
                    stream_pos+=curr_phrase.len;
                    memcpy(&curr_phrase.mt, &stream[stream_pos], sizeof(uint32_t));
                    stream_pos+=sizeof(uint32_t);
                } else {
                    curr_phrase.len = stream[stream_pos];
                    stream_pos++;
                    curr_phrase.phrase = &stream[stream_pos];
                    stream_pos+=curr_phrase.len;
                    curr_phrase.mt = stream[stream_pos];
                    stream_pos++;
                }
                assert(stream_pos<=stream_size);
            }else{
                stream_pos = stream_size+1;
                curr_phrase.len = 0;
                curr_phrase.mt = 0;
                curr_phrase.phrase = nullptr;
            }
        };

    public:

        explicit iterator(seq_type* _stream, size_t _start_pos, size_t _stream_size) : stream(_stream),
                                                                                       stream_pos(_start_pos),
                                                                                       stream_size(_stream_size),
                                                                                       curr_phrase(nullptr, 0, 0) {
            decode_phrase();
        }

        // Move to the next tuple
        inline void operator++() {
            decode_phrase();
        }

        inline phrase_t& operator*(){
            return curr_phrase;
        }

        inline bool operator==(const iterator& other) {
            return other.stream==stream && other.stream_pos==stream_pos;
        }

        inline bool operator!=(const iterator& other) {
            return other.stream!=stream || other.stream_pos!=stream_pos;
        }
    };

private:

    table_t m_table;
    float m_max_load_factor = 0.6;
    size_t elm_threshold=0;
    size_t frac_lf = 60;
    size_t n_phrases=0;
    size_t last_mt=0;
    size_t last_fp_pos=0;
    size_t n_rehashes=0;

    void rehash(size_t new_tab_size) {
        //std::cout<<"rehashing "<<n_phrases<<" phrases to a hash table of size "<<new_tab_size<<std::endl;
        assert(new_tab_size>m_table.size());
        m_table.resize(new_tab_size);
        memset(m_table.data(), (int)null_addr, m_table.size()*sizeof(table_t::value_type));

        //rehash the values
        uint32_t len, phr_addr, pos=0;
        size_t proc_phrases=0;
        while(pos<stream_size){
            phr_addr=pos;
            if constexpr (std::is_same<seq_type, uint8_t>::value){
                memcpy(&len, &phrase_stream[pos], sizeof(uint32_t));//read the length
                pos+=sizeof(uint32_t);//skip length
                //insert_entry_in_table_bucket(XXH3_64bits(&phrase_stream[pos], len), phr_addr);
                rehash_entry_rb(XXH3_64bits(&phrase_stream[pos], len), phr_addr);
                pos+=len+sizeof(uint32_t);//skip the phrase and mt
            }else{
                len = phrase_stream[pos];//read the length
                pos++;//skip length
                //insert_entry_in_table_bucket(XXH3_64bits(&phrase_stream[pos], len*seq_bytes), phr_addr);
                rehash_entry_rb(XXH3_64bits(&phrase_stream[pos], len*seq_bytes), phr_addr);
                pos+=len+1;//skip the phrase and mt
            }
            proc_phrases++;
            assert(proc_phrases<=n_phrases);
        }
        assert(proc_phrases==n_phrases);
        assert(pos==stream_size);
        elm_threshold = (m_table.size()*frac_lf)/100;
        n_rehashes++;
    }

    inline void insert_entry_in_table_bucket(uint64_t hash, uint64_t phr_addr) {

        size_t idx = hash & (m_table.size()-1);
        if(m_table[idx]==null_addr){
            m_table[idx] = phr_addr;
            return;
        }

        //find a new empty sm_bucket
        size_t j=1;
        while(true){
            idx = (hash + ((j*j+j)>>1UL)) & (m_table.size()-1);
            if(m_table[idx]==null_addr){
                m_table[idx] = phr_addr;
                break;
            }
            j++;
        }
    }

    void increase_stream_cap(size_t min_cap){
        size_t bytes = stream_cap*sizeof(seq_type);
        if(bytes <= STR_THRESHOLD1){
            stream_cap = stream_size*2;
        } else if(bytes <= STR_THRESHOLD2){
            stream_cap = static_cast<size_t>(stream_size*1.5);
        } else if(bytes <= STR_THRESHOLD3){
            stream_cap = static_cast<size_t>(stream_size*1.2);
        } else {
            stream_cap = static_cast<size_t>(stream_size*1.05);
        }
        stream_cap = std::max(min_cap, stream_cap);

        if(phrase_stream== nullptr){
            phrase_stream = mem<seq_type>::allocate(stream_cap);
        }else{
            phrase_stream = mem<seq_type>::reallocate(phrase_stream, stream_cap);
        }
    }

public:

    explicit phrase_set(size_t min_cap=4, float max_lf=0.85){
        assert(min_cap>0);
        m_max_load_factor = max_lf;
        m_table = table_t(round_to_power_of_two(min_cap), null_addr);
        frac_lf = size_t(m_max_load_factor*100);
        elm_threshold = (m_table.size()*frac_lf)/100;
    }

    /*inline uint32_t insert(seq_type* q_phrase, size_t q_len, bool& inserted) {
        size_t q_bytes = q_len*seq_bytes;
        size_t hash = XXH3_64bits(q_phrase, q_bytes);

        size_t j = 0;
        size_t idx = hash & (m_table.size()-1);

        while(m_table[idx]!=null_addr) {
            uint32_t len;
            uint64_t phr_addr = m_table[idx];
            if constexpr (std::is_same<seq_type, uint8_t>::value){
                memcpy(&len, &phrase_stream[phr_addr], sizeof(uint32_t));
                phr_addr+=sizeof(uint32_t);
                if(q_len == len &&
                   memcmp(q_phrase, &phrase_stream[phr_addr], q_bytes)==0){
                    inserted = false;
                    uint32_t mt;
                    memcpy(&mt, &phrase_stream[phr_addr+len], sizeof(uint32_t));
                    return mt;
                }
            }else{
                len = phrase_stream[phr_addr];
                phr_addr++;
                if(q_len == len &&
                   memcmp(q_phrase, &phrase_stream[phr_addr], q_bytes)==0){
                    inserted = false;
                    return phrase_stream[phr_addr+len];
                }
            }
            j++;
            idx = (hash + ((j*j + j)>>1UL)) & (m_table.size()-1);
        }

        m_table[idx] = stream_size;
        uint32_t mt = n_phrases++;

        if constexpr (std::is_same<seq_type, uint8_t>::value){
            size_t n_words = (2*sizeof(uint32_t)) + q_len;
            if((stream_size+n_words)>=stream_cap){
                increase_stream_cap(stream_size+n_words);
            }
            memcpy(&phrase_stream[stream_size], &q_len, sizeof(uint32_t));
            stream_size+=sizeof(uint32_t);
            memcpy(&phrase_stream[stream_size], q_phrase, q_bytes);
            stream_size+=q_len;
            memcpy(&phrase_stream[stream_size], &mt, sizeof(uint32_t));
            stream_size+=sizeof(uint32_t);
        } else {
            size_t n_words = 2+q_len;
            if((stream_size+n_words)>=stream_cap){
                increase_stream_cap(stream_size+n_words);
            }
            phrase_stream[stream_size] = q_len;
            stream_size++;
            memcpy(&phrase_stream[stream_size], q_phrase, q_bytes);
            stream_size+=q_len;
            phrase_stream[stream_size] = mt;
            stream_size++;
        }

        //the insertion exceeds the max. load factor (i.e., rehash)
        if(n_phrases>=elm_threshold) {
            rehash(next_power_of_two(m_table.size()));
        }

        inserted = true;
        return mt;
    }*/

    inline static bool compare(seq_type * q_phrase, size_t q_len, seq_type* phr_addr){
        if constexpr (std::is_same<seq_type, uint8_t>::value){
            uint32_t len;
            memcpy(&len, phr_addr, sizeof(uint32_t));
            phr_addr+=sizeof(uint32_t);
            return (q_len == len && memcmp(q_phrase, phr_addr, q_len)==0);
        } else {
            return (q_len == *phr_addr && memcmp(q_phrase, (phr_addr+1), q_len*seq_bytes)==0);
        }
    }

    inline void rehash_entry_rb(uint64_t hash, size_t phr_addr){

        size_t idx = hash & (m_table.size() - 1);
        size_t dist=0, bck_addr, bck_dist;

        if(m_table[idx]!=null_addr) {

            bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
            bck_dist = m_table[idx] >> 44UL;
            while(bck_dist>=dist && m_table[idx]!=null_addr){
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
                bck_dist = m_table[idx] >> 44UL;
            }

            assert(bck_dist<dist || m_table[idx]==null_addr);

            while(m_table[idx]!=null_addr){
                m_table[idx] = (dist<<44UL) | phr_addr;
                phr_addr = bck_addr;
                dist = bck_dist+1;

                idx = (idx+1) & (m_table.size()-1);
                bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
                bck_dist = m_table[idx] >> 44UL;
            }
        }
        assert(dist<=65535);
        m_table[idx] = (dist<<44UL) | phr_addr ;
    }

    uint32_t insert(seq_type* q_phrase, size_t q_len) {

        size_t q_bytes = q_len*seq_bytes;
        uint64_t hash = XXH3_64bits(q_phrase, q_bytes);
        size_t idx = hash & (m_table.size()-1), bck_dist, bck_addr, dist=0, addr=stream_size;

        if(m_table[idx]!=null_addr) {

            bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
            bck_dist = m_table[idx] >> 44UL;
            while(bck_dist>dist && m_table[idx]!=null_addr){
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
                bck_dist = m_table[idx] >> 44UL;
            }

            while(dist==bck_dist){
                if(compare(q_phrase, q_len, &phrase_stream[bck_addr])){
                    if(std::is_same<seq_type, uint8_t>::value){
                        uint32_t mt;
                        memcpy(&mt, &phrase_stream[bck_addr+sizeof(uint32_t)+q_len], sizeof(uint32_t));
                        return mt;
                    } else {
                        return phrase_stream[bck_addr+1+q_len];
                    }
                }
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
                bck_dist = m_table[idx] >> 44UL;
            }

            assert(bck_dist<dist || m_table[idx]==null_addr);

            while(m_table[idx]!=null_addr){
                m_table[idx] = (dist<<44UL) | addr;
                addr = bck_addr;
                dist = bck_dist+1;

                idx = (idx+1) & (m_table.size()-1);
                bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
                bck_dist = m_table[idx] >> 44UL;
            }
        }

        assert(dist<=65535);
        m_table[idx] = (dist<<44UL) | addr ;

        uint32_t mt = n_phrases++;
        if constexpr (std::is_same<seq_type, uint8_t>::value) {
            size_t n_words = 2*sizeof(uint32_t) + q_len;
            if((stream_size+n_words)>=stream_cap){
                increase_stream_cap(stream_size+n_words);
            }
            memcpy(&phrase_stream[stream_size], &q_len, sizeof(uint32_t));
            stream_size+=sizeof(uint32_t);
            memcpy(&phrase_stream[stream_size], q_phrase, q_len);
            stream_size+=q_len;
            memcpy(&phrase_stream[stream_size], &mt, sizeof(uint32_t));
            stream_size+=sizeof(uint32_t);
        } else {
            size_t n_words = 2+q_len;
            if((stream_size+n_words)>=stream_cap){
                increase_stream_cap(stream_size+n_words);
            }
            phrase_stream[stream_size] = q_len;
            stream_size++;
            memcpy(&phrase_stream[stream_size], q_phrase, q_len*seq_bytes);
            stream_size+=q_len;
            phrase_stream[stream_size] = mt;
            stream_size++;
        }

        //the insertion exceeds the max. load factor (i.e., rehash)
        if(n_phrases>=elm_threshold) {
            rehash(next_power_of_two(m_table.size()));
        }

        return mt;
    }

    void set_min_capacity(size_t new_cap){
        new_cap = round_to_power_of_two(new_cap);
        if(new_cap>m_table.size()){
            rehash(new_cap);
        }
    }

    /*bool find(off_t source, size_t len, size_t& mt) const {

        size_t hash = XXH3_64bits(&seq[source], len*seq_bytes);
        size_t j=0;
        size_t idx = hash & (m_table.size()-1);

        while(true){
            if(m_table[idx]==null_source){
                return false;
            }else{
                const phrase_t & phrase = phrases[m_table[idx]];
                if(len==phrase.len &&
                   memcmp(&seq[source], &seq[phrase.source], len*seq_bytes)==0) {
                    mt = m_table[idx];
                    return true;
                }
                j++;
                idx = (hash + ((j*j + j)>>1UL)) & (m_table.size()-1);
            }
        }
    }*/

    [[nodiscard]] inline float load_factor() const {
        return float(n_phrases)/float(m_table.size());
    }

    [[nodiscard]] inline float max_load_factor() const {
        return m_max_load_factor;
    };

    [[nodiscard]] inline size_t size() const {
        return n_phrases;
    }

    [[nodiscard]] inline bool empty()  const {
        return stream_size==0;
    }

    [[nodiscard]] inline size_t capacity() const{
        return m_table.size();
    };

    [[nodiscard]] inline size_t table_mem_usage() const {
        return m_table.size()*sizeof(table_t::value_type);
    }

    [[nodiscard]] inline size_t phrases_mem_usage() const {
        return stream_cap*sizeof(seq_type);
    }

    [[nodiscard]] inline size_t mem_usage() const {
        return table_mem_usage()+phrases_mem_usage();
    }

    [[nodiscard]] inline size_t eff_mem_usage() const {
        return (n_phrases*sizeof(table_t::value_type)) + (stream_size*sizeof(seq_type));
    }

    void shrink_to_fit() {
        stream_cap = stream_size;
        phrase_stream = mem<seq_type>::reallocate(phrase_stream, stream_cap);
    }

    [[nodiscard]] iterator begin() const {
        return iterator(phrase_stream, 0, stream_size);
    }

    [[nodiscard]] iterator begin(size_t stream_pos) const {
        assert(stream_pos>=0 && stream_pos<stream_size);
        return iterator(phrase_stream, stream_pos, stream_size);
    }

    [[nodiscard]] iterator end() const {
        return iterator(phrase_stream, stream_size, stream_size);
    }

    void destroy_table(){
        std::vector<uint64_t>().swap(m_table);
    }

    [[nodiscard]] size_t buff_bytes_available() const {
        return (stream_cap-stream_size)*sizeof(seq_type);
    }

    ~phrase_set(){
        if(phrase_stream!= nullptr){
            mem<seq_type>::deallocate(phrase_stream);
            phrase_stream= nullptr;
            destroy_table();
        }
    }

    void psl_dist(){
        size_t freqs[2000]={0};
        for(unsigned long long entry : m_table){
            if(entry!=null_addr){
                size_t bck_dist = entry >> 44UL;
                if(bck_dist<2000){
                    freqs[bck_dist]++;
                }
            }
        }
        for(size_t i=0;i<2000;i++){
            if(freqs[i]>0){
                std::cout<<"psl: "<<i<<" --> "<<freqs[i]<<std::endl;
            }
        }
    }

    inline void update_fps(const uint64_t* prev_fps, size_t& len_prev_fps, uint64_t*& fps, size_t& len_fps) {
        if(last_mt<n_phrases){
            assert(last_fp_pos<=stream_size);

            fps = mem<uint64_t>::reallocate(fps, n_phrases+1);
            len_fps = n_phrases+1;

            auto it = iterator(phrase_stream, last_fp_pos, stream_size);
            auto it_end = end();
            std::vector<uint64_t> fp_sequence;
            while(it!=it_end){
                auto phr = *it;
                for(size_t j=0;j<phr.len;j++){
                    assert(phr.phrase[j]>0 && phr.phrase[j]<len_prev_fps);
                    fp_sequence.push_back(prev_fps[phr.phrase[j]]);
                }
                fps[last_mt+1] = XXH3_64bits(fp_sequence.data(), fp_sequence.size()*sizeof(uint64_t));
                last_mt++;
                ++it;
                fp_sequence.clear();
            }
            last_fp_pos=stream_size;
        }
    }
};

#endif //LCG_LZL_MAP_H
