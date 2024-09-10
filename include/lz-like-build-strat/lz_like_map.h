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

template<class seq_type>
class lz_like_map {

public:

    seq_type *seq;//pointer to the source sequence
    static constexpr uint8_t seq_bytes=sizeof(seq_type);

    static constexpr uint32_t null_source = std::numeric_limits<uint32_t>::max();

    struct phrase_t {
        uint32_t source = null_source;
        uint32_t len;
        phrase_t(uint32_t _source, uint32_t _len): source(_source), len(_len){}

        [[nodiscard]] std::string to_string() const{
            std::string str;
            str+="source: "+std::to_string(source)+", ";
            str+="len: "+std::to_string(len)+", ";
            return str;
        }
    };

    typedef std::vector<uint32_t> table_t;
    typedef buff_vector<phrase_t> phrase_list_t;
    //typedef std::vector<phrase_t> phrase_list_t;

private:

    table_t m_table;
    phrase_list_t phrases;
    float m_max_load_factor = 0.6;
    size_t elm_threshold=0;
    size_t frac_lf = 60;

    void rehash(size_t new_tab_size) {

        assert(new_tab_size>m_table.size());
        m_table.resize(new_tab_size);
        memset(m_table.data(), (int)null_source, m_table.size()*sizeof(uint32_t));

        //rehash the values
        uint32_t p_idx=0;
        for(auto const & phrase : phrases) {
            insert_entry_in_table_bucket(phrase, p_idx);
            p_idx++;
        }

        elm_threshold = (m_table.size()*frac_lf)/100;
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    inline void insert_entry_in_table_bucket(const phrase_t& phrase, const uint32_t p_idx) {

        size_t hash = XXH3_64bits(&seq[phrase.source], phrase.len*seq_bytes);
        size_t idx = hash & (m_table.size()-1);

        if(m_table[idx]==null_source){
            m_table[idx] = p_idx;
            return;
        }

        //find a new empty sm_bucket
        size_t j=1;
        while(true){
            idx = (hash + ((j*j+j)>>1UL)) & (m_table.size()-1);
            if(m_table[idx]==null_source){
                m_table[idx] = p_idx;
                break;
            }
            j++;
        }
    }

public:

    const phrase_list_t& phrase_set = phrases;

    explicit lz_like_map(seq_type* _seq, char * buffer=nullptr, size_t buff_bytes=0, size_t min_cap=4, float max_lf=0.6) : seq(_seq),
                                                                                                                           phrases(buffer, buff_bytes),
                                                                                                                           m_max_load_factor(max_lf){
        m_table = table_t(round_to_power_of_two(min_cap), null_source);
        frac_lf = size_t(m_max_load_factor*100);
        elm_threshold = (m_table.size()*frac_lf)/100;
    }

     inline uint32_t insert(off_t q_source, size_t q_len, bool& inserted) {
        size_t q_bytes = q_len*seq_bytes;
        size_t hash = XXH3_64bits(&seq[q_source], q_bytes);

        size_t j = 0;
        size_t idx = hash & (m_table.size()-1);
        while(m_table[idx]!=null_source) {
            if(q_len == phrases[m_table[idx]].len &&
               memcmp(&seq[q_source], &seq[phrases[m_table[idx]].source], q_bytes)==0){
                inserted = false;
                //phrases[m_table[idx]].repeated=true;
                // the reference is always the rightmost occurrence
                // in the text's scan
                // phrase.source = q_source;
                return m_table[idx];
            }
            j++;
            idx = (hash + ((j*j + j)>>1UL)) & (m_table.size()-1);
        }

        m_table[idx] = phrases.size();
        phrases.emplace_back(q_source, q_len);

        //the insertion exceeds the max. load factor (i.e., rehash)
        if(phrases.size()>=elm_threshold) {
            rehash(next_power_of_two(m_table.size()));
        }

        inserted = true;
        return phrases.size()-1;
    }

    /*inline uint32_t unhashed_insert(off_t q_source, size_t q_len) {
        assert(q_len<=1073741823);
        phrases.emplace_back(q_source, q_len, false, true);
        return phrases.size()-1;
    }*/

    void set_min_capacity(size_t new_cap){
        new_cap = round_to_power_of_two(new_cap);
        if(new_cap>m_table.size()){
            rehash(new_cap);
        }
    }

    bool find(off_t source, size_t len, size_t& mt) const {

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
    }

    [[nodiscard]] inline float load_factor() const {
        return float(phrases.size())/float(m_table.size());
    }

    [[nodiscard]] inline float max_load_factor() const {
        return m_max_load_factor;
    };

    [[nodiscard]] inline size_t size() const {
        return phrases.size();
    }

    [[nodiscard]] inline bool empty()  const {
        return phrases.empty();
    }

    [[nodiscard]] inline size_t capacity() const{
        return m_table.size();
    };

    size_t table_mem_usage(){
        return m_table.size()*sizeof(uint32_t);
    }

    size_t phrases_mem_usage(){
        return phrases.capacity()*sizeof(phrase_t);
    }

    void shrink_to_fit(){
        phrases.shrink_to_fit();
    }

    void insert_dummy_entry(phrase_t dummy){
        phrases.push_back(dummy);
    }

    size_t mem_usage(){
        return table_mem_usage()+phrases_mem_usage();
    }

    void destroy_table(){
        std::vector<uint32_t>().swap(m_table);
#ifdef __linux__
        malloc_trim(0);
#endif
    }
};

#endif //LCG_LZL_MAP_H
