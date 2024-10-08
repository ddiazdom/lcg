//
// Created by Diaz, Diego on 17.9.2023.
//

#ifndef LCG_LZL_MAP_H
#define LCG_LZL_MAP_H
#include "../../external/xxHash-dev/xxhash.h"
#include "cds/cdt_common.hpp"
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

private:

    typedef std::vector<uint64_t> table_t;
    static constexpr uint8_t seq_bytes=sizeof(seq_type);
    static constexpr uint64_t null_addr = std::numeric_limits<uint64_t>::max();
    static constexpr uint64_t d_one = 1UL<<44UL;

    seq_type *phrase_stream = nullptr;
    size_t stream_size=0;
    size_t stream_cap=0;
    table_t m_table;
    float m_max_load_factor = 0.6;
    size_t elm_threshold=0;
    size_t frac_lf = 60;
    size_t n_phrases=0;
    size_t last_mt=0;
    size_t last_fp_pos=0;

    void rehash(size_t new_tab_size) {
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
                rehash_entry_rb(XXH3_64bits(&phrase_stream[pos], len), phr_addr);
                pos+=len+sizeof(uint32_t);//skip the phrase and mt
            }else{
                len = phrase_stream[pos];//read the length
                pos++;//skip length
                rehash_entry_rb(XXH3_64bits(&phrase_stream[pos], len*seq_bytes), phr_addr);
                pos+=len+1;//skip the phrase and mt
            }
            proc_phrases++;
            assert(proc_phrases<=n_phrases);
        }
        assert(proc_phrases==n_phrases);
        assert(pos==stream_size);
        elm_threshold = (m_table.size()*frac_lf)/100;
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
        size_t curr_phr_pos = 0;
        phrase_t curr_phrase;

        void decode_phrase(){
            curr_phr_pos = stream_pos;
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

        [[nodiscard]] inline size_t pos() const {
            return curr_phr_pos;
        }

        inline bool operator==(const iterator& other) {
            return other.stream==stream && other.stream_pos==stream_pos;
        }

        inline bool operator!=(const iterator& other) {
            return other.stream!=stream || other.stream_pos!=stream_pos;
        }
    };

    explicit phrase_set(size_t min_cap=4, float max_lf=0.85) {
        assert(min_cap>0);
        m_max_load_factor = max_lf;
        m_table = table_t(round_to_power_of_two(min_cap), null_addr);
        frac_lf = size_t(m_max_load_factor*100);
        elm_threshold = (m_table.size()*frac_lf)/100;
    }

    phrase_set(phrase_set&& other) noexcept {
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

    phrase_set(const phrase_set& other) noexcept {
        copy(other);
    }

    void copy(const phrase_set& other){
        stream_size = other.stream_size;
        stream_cap = other.stream_cap;
        m_table = other.m_table;
        m_max_load_factor = other.m_max_load_factor;
        elm_threshold = other.elm_threshold;
        frac_lf = other.frac_lf;
        n_phrases = other.n_phrases;
        last_mt = other.last_mt;
        last_fp_pos = other.last_fp_pos;

        if(phrase_stream!= nullptr){
            mem<seq_type>::deallocate(phrase_stream);
        }
        phrase_stream = mem<seq_type>::allocate(stream_cap);
        memcpy(phrase_stream, other.phrase_stream, sizeof(seq_type)*stream_size);
    }

    void swap(phrase_set& other){
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

    phrase_set& operator=(const phrase_set& other){
        if(this!=&other){
            copy(other);
        }
        return *this;
    }

    inline static bool compare(const seq_type * q_phrase, size_t q_len, seq_type* phr_addr){
        if constexpr (std::is_same<seq_type, uint8_t>::value){
            uint32_t len;
            memcpy(&len, phr_addr, sizeof(uint32_t));
            return (q_len == len && memcmp(q_phrase, (phr_addr+sizeof(uint32_t)), q_len)==0);
        } else {
            return (q_len == *phr_addr && memcmp(q_phrase, (phr_addr+1), q_len*seq_bytes)==0);
        }
    }

    inline void rehash_entry_rb(uint64_t hash, uint64_t phr_addr){

        size_t idx = hash & (m_table.size() - 1);
        if(m_table[idx]!=null_addr) {
            size_t dist=0, bck_dist;
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
    inline size_t add_phrase(const seq_type* phrase, size_t len){

        size_t q_bytes = len*seq_bytes;
        uint64_t hash = XXH3_64bits(phrase, q_bytes);
        uint32_t mt;

        size_t idx = hash & (m_table.size()-1);
        uint64_t bck_val = stream_size;
        if(m_table[idx]!=null_addr) {
            size_t dist=0, bck_dist;
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
        if constexpr (std::is_same<seq_type, uint8_t>::value) {
            size_t n_words = 2*sizeof(uint32_t) + len;
            if((stream_size+n_words)>=stream_cap){
                increase_stream_cap(stream_size+n_words);
            }
            memcpy(&phrase_stream[stream_size], &len, sizeof(uint32_t));
            stream_size+=sizeof(uint32_t);
            memcpy(&phrase_stream[stream_size], phrase, len);
            stream_size+=len;
            memcpy(&phrase_stream[stream_size], &mt, sizeof(uint32_t));
            stream_size+=sizeof(uint32_t);
        } else {
            size_t n_words = 2+len;
            if((stream_size+n_words)>=stream_cap){
                increase_stream_cap(stream_size+n_words);
            }
            phrase_stream[stream_size] = len;
            stream_size++;
            memcpy(&phrase_stream[stream_size], phrase, len*seq_bytes);
            stream_size+=len;
            phrase_stream[stream_size] = mt;
            stream_size++;
        }

        //the insertion exceeds the max. load factor (i.e., rehash)
        if(n_phrases>=elm_threshold) {
            rehash(next_power_of_two(m_table.size()));
        }

        return mt;
    }

    inline uint32_t insert(const seq_type* q_phrase, size_t q_len) {

        size_t q_bytes = q_len*seq_bytes;
        uint64_t hash = XXH3_64bits(q_phrase, q_bytes);
        uint32_t mt;

        size_t idx = hash & (m_table.size()-1);
        uint64_t bck_val = stream_size;
        if(m_table[idx]!=null_addr) {
            size_t bck_dist = m_table[idx] >> 44UL, dist=0;
            while(bck_dist>dist && m_table[idx]!=null_addr){
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }

            while(dist==bck_dist){
                uint64_t bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
                if(compare(q_phrase, q_len, &phrase_stream[bck_addr])){
                    if constexpr (std::is_same<seq_type, uint8_t>::value){
                        memcpy(&mt, &phrase_stream[bck_addr+sizeof(uint32_t)+q_len], sizeof(uint32_t));
                    } else {
                        mt = phrase_stream[bck_addr+1+q_len];
                    }
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

    inline bool find(const seq_type* q_phrase, size_t q_len, uint32_t& mt) const {

        size_t q_bytes = q_len*seq_bytes;
        uint64_t hash = XXH3_64bits(q_phrase, q_bytes);
        size_t idx = hash & (m_table.size()-1);

        if(m_table[idx]!=null_addr) {
            size_t bck_dist = m_table[idx] >> 44UL, dist=0;
            while(bck_dist>dist && m_table[idx]!=null_addr){
                dist++;
                idx = (idx+1) & (m_table.size()-1);
                bck_dist = m_table[idx] >> 44UL;
            }

            while(dist==bck_dist){
                uint64_t bck_addr = (m_table[idx] & 0xFFFFFFFFFFFul);
                if(compare(q_phrase, q_len, &phrase_stream[bck_addr])){
                    if constexpr (std::is_same<seq_type, uint8_t>::value){
                        memcpy(&mt, &phrase_stream[bck_addr+sizeof(uint32_t)+q_len], sizeof(uint32_t));
                    } else {
                        mt = phrase_stream[bck_addr+1+q_len];
                    }
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

    [[nodiscard]] inline size_t size() const {
        return n_phrases;
    }

    [[nodiscard]] inline bool empty()  const {
        return n_phrases==0;
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

    inline const seq_type* phr_stream() const {
        return phrase_stream;
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
        if(n_phrases>0){
            size_t freqs[2000]={0};
            for(unsigned long long entry : m_table){
                if(entry!=null_addr){
                    size_t bck_dist = entry >> 44UL;
                    if(bck_dist<2000){
                        freqs[bck_dist]++;
                    }
                }
            }
            size_t avg=0;
            float acc=0;
            for(size_t i=0;i<2000;i++){
                if(freqs[i]>0){
                    float frac = float(freqs[i])/float(n_phrases);
                    acc +=frac;
                    std::cout<<"psl: "<<i<<" --> "<<freqs[i]<<" "<<frac<<" "<<acc<<std::endl;
                    avg+=freqs[i]*i;
                }
            }
            std::cout<<"The average is "<<float(avg)/float(n_phrases)<<std::endl;
            //std::cout<<" A "<<float(a)/float(a+b+c)<<std::endl;
            //std::cout<<" B "<<float(b)/float(a+b+c)<<std::endl;
            //std::cout<<" C "<<float(c)/float(a+b+c)<<std::endl;
        }
    }

    [[nodiscard]] inline size_t tot_symbols() const {
        if constexpr (std::is_same<seq_type, uint8_t>::value){
            return stream_size- (n_phrases*sizeof(uint32_t)*2);
        }else{
            return stream_size- (n_phrases*2);
        }
    }

    size_t absorb_set(phrase_set<seq_type>& other,
                      const uint64_t* i_map, size_t i_map_len,
                      uint64_t* o_map, size_t o_map_len){

        if(other.empty()) return size();
        assert(o_map_len==(other.size()+1));
        uint32_t len, pos=0, mt;
        size_t proc_phrases=0;
        while(pos<other.stream_size){
            if constexpr (std::is_same<seq_type, uint8_t>::value){
                memcpy(&len, &other.phrase_stream[pos], sizeof(uint32_t));//read the length
                pos+=sizeof(uint32_t);//skip length
                mt = insert(&other.phrase_stream[pos], len);
                pos+=len+sizeof(uint32_t);//skip the phrase and mt
            }else{
                len = other.phrase_stream[pos];//read the length
                pos++;//skip length
                size_t end = pos+len;
                for(size_t j=pos;j<end;j++){
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


        return size();
    }

    void clear(){
        stream_size=0;
        stream_cap = std::min<size_t>(4, stream_cap>>3);//keep 1/8 of the buffer
        phrase_stream = mem<seq_type>::reallocate(phrase_stream, stream_cap);
        size_t new_tab_size = std::max<size_t>(4, (m_table.size()>>3));// keep 1/8 of the table
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
