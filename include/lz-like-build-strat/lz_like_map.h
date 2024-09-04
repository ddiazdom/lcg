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

template<class data_type>
class lz_like_map {

public:

    data_type *data;//pointer to the source sequence
    static constexpr uint8_t data_bytes=sizeof(data_type);

    static constexpr uint32_t null_source = std::numeric_limits<uint32_t>::max();

    struct phrase_t {
        uint32_t source = null_source;
        uint32_t len:30;
        bool repeated:1;
        bool bypassed:1;
        phrase_t(uint32_t _source, uint32_t _len, bool _repeated, bool _bypassed): source(_source), len(_len), repeated(_repeated), bypassed(_bypassed){}

        [[nodiscard]] std::string to_string() const{
            std::string str;
            str+="source: "+std::to_string(source)+", ";
            str+="len: "+std::to_string(len)+", ";
            str+="repeated: "+std::to_string(repeated)+", ";
            str+="bypassed: "+std::to_string(bypassed);



            return str;
        }
    };

    typedef std::vector<uint32_t> table_t;
    typedef std::vector<phrase_t> phrase_list_t;

private:

    table_t m_table;
    phrase_list_t phrases;
    float m_max_load_factor = 0.6;
    size_t elm_threshold=0;
    size_t frac_lf = 60;
    size_t n_hashed=0;

    void rehash(size_t new_tab_size) {

        assert(new_tab_size>m_table.size());
        m_table.resize(new_tab_size);
        memset(m_table.data(), (int)null_source, m_table.size()*sizeof(uint32_t));

        //rehash the values
        uint32_t p_idx=0;
        for(auto const & phrase : phrases) {
            if(!phrase.bypassed){
                insert_entry_in_table_bucket(phrase, p_idx);
            }
            p_idx++;
        }

        elm_threshold = (m_table.size()*frac_lf)/100;
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    inline void insert_entry_in_table_bucket(const phrase_t& phrase, const uint32_t p_idx) {

        size_t hash = XXH3_64bits(&data[phrase.source], phrase.len*data_bytes);
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

    explicit lz_like_map(data_type* _data, size_t min_cap=4, float max_lf=0.6) : data(_data), m_max_load_factor(max_lf) {
        m_table = table_t(round_to_power_of_two(min_cap), null_source);
        frac_lf = size_t(m_max_load_factor*100);
        elm_threshold = (m_table.size()*frac_lf)/100;
    }

     inline uint32_t insert(off_t q_source, size_t q_len, bool& inserted) {

        size_t q_bytes = q_len*data_bytes;
        size_t hash = XXH3_64bits(&data[q_source], q_bytes);

        size_t j = 0;
        size_t idx = hash & (m_table.size()-1);
        while(m_table[idx]!=null_source) {
            if(q_len == phrases[m_table[idx]].len &&
               memcmp(&data[q_source], &data[phrases[m_table[idx]].source], q_bytes)==0){
                inserted = false;
                phrases[m_table[idx]].repeated=true;
                // the reference is always the rightmost occurrence
                // in the text's scan
                // phrase.source = q_source;
                return m_table[idx];
            }
            j++;
            idx = (hash + ((j*j + j)>>1UL)) & (m_table.size()-1);
        }

        m_table[idx] = phrases.size();
        phrases.emplace_back(q_source, q_len, false, false);
        n_hashed++;

        //the insertion exceeds the max. load factor (i.e., rehash)
        if(n_hashed>=elm_threshold) {
            rehash(next_power_of_two(m_table.size()));
        }

        inserted = true;
        return phrases.size()-1;
    }

    inline uint32_t unhashed_insert(off_t q_source, size_t q_len) {
        phrases.emplace_back(q_source, q_len, false, true);
        return phrases.size()-1;
    }

    void set_min_capacity(size_t new_cap){
        new_cap = round_to_power_of_two(new_cap);
        if(new_cap>m_table.size()){
            rehash(new_cap);
        }
    }

    bool find(off_t source, size_t len, size_t& mt) const {

        size_t hash = XXH3_64bits(&data[source], len*data_bytes);
        size_t j=0;
        size_t idx = hash & (m_table.size()-1);

        while(true){
            if(m_table[idx]==null_source){
                return false;
            }else{
                const phrase_t & phrase = phrases[m_table[idx]];
                if(len==phrase.len &&
                   memcmp(&data[source], &data[phrase.source], len*data_bytes)==0) {
                    mt = m_table[idx];
                    return true;
                }
                j++;
                idx = (hash + ((j*j + j)>>1UL)) & (m_table.size()-1);
            }
        }
    }

    [[nodiscard]] inline float load_factor() const {
        return float(n_hashed)/float(m_table.size());
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

    float hashed_fraction(){
        return float(n_hashed)/float(phrases.size());
    }

    size_t phrases_mem_usage(){
        return phrases.capacity()*sizeof(phrase_t);
    }

    void shrink_to_fit(){
        phrases.shrink_to_fit();
#ifdef __linux__
        malloc_trim(0);
#endif
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
