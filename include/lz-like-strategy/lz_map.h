//
// Created by Diaz, Diego on 17.9.2023.
//

#ifndef LCG_LZL_MAP_H
#define LCG_LZL_MAP_H

class lz_map {

public:

    uint8_t *data;
    static constexpr uint32_t null_source = std::numeric_limits<uint32_t>::max();
    struct bucket_t{
        uint32_t source = null_source;
        uint32_t len{};
        uint32_t mt_sym{};
    };

private:

    typedef std::vector<bucket_t> table_t;
    table_t* m_table{nullptr};
    size_t n_elms = 0;
    uint32_t next_av_mt=1;
    float m_max_load_factor = 0.6;

    void rehash(size_t new_tab_size) {

        assert(new_tab_size>m_table->size());
        auto * new_table = new table_t(new_tab_size);

        //rehash the values
        for(auto & bucket : *m_table) {
            if(bucket.source!=null_source){
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

        size_t hash = XXH3_64bits(&data[bucket.source], bucket.len);
        size_t idx = hash & (new_tab_size-1);

        if(new_table[idx].source==null_source){
            new_table[idx] = bucket;
            return;
        }

        //find a new empty sm_bucket
        size_t j=1;
        while(true){
            idx = (hash + ((j*j+j)>>1UL)) & (new_tab_size-1);
            if(new_table[idx].source==null_source){
                new_table[idx] = bucket;
                break;
            }
            j++;
        }
    }

public:

    table_t *const& table = m_table;

    explicit lz_map(uint8_t* _data, size_t min_cap=4, float max_lf=0.6) : data(_data), m_max_load_factor(max_lf) {
        m_table = new table_t(round_to_power_of_two(min_cap));
    }

    ~lz_map(){
        delete m_table;
    }

     uint32_t insert(off_t source, size_t len, bool& inserted) {

        inserted = false;
        bool success = false;
        uint32_t mt;

        size_t hash = XXH3_64bits(&data[source], len);
        size_t j= 0;
        size_t idx = hash & (m_table->size()-1);

        while(!success) {
            bucket_t &bucket = (*m_table)[idx];

            if(bucket.source==null_source) {
                bucket.source = source;
                bucket.len = len;
                bucket.mt_sym = next_av_mt++;
                mt = bucket.mt_sym;

                //the key insertion exceeds the max. load factor (i.e., rehash)
                n_elms++;
                if((float(n_elms)/float(m_table->size()))>=m_max_load_factor) {
                    rehash(next_power_of_two(m_table->size()));
                }
                success = true;
                inserted = true;
            } else if(len == bucket.len &&
                      memcmp(&data[source], &data[bucket.source], len)==0){
                mt = bucket.mt_sym;
                success = true;
                inserted = false;
            } else {
                j++;
                idx = (hash + ((j*j + j)>>1UL)) & (m_table->size()-1);
            }
        }
        return mt;
    }

    void set_min_capacity(size_t new_cap){
        new_cap = round_to_power_of_two(new_cap);
        if(new_cap>m_table->size()){
            rehash(new_cap);
        }
    }

    bool find(const off_t source, size_t len, size_t& bucket) const {

        size_t hash = XXH3_64bits(&data[source], len);
        size_t j=0;
        size_t idx = hash & (m_table->size()-1);

        while(true){
            if((*m_table)[idx].source==null_source){
                return false;
            }else if(len==(*m_table)[idx].len &&
                     memcmp(&data[source], &data[(*m_table)[idx].source], len)==0) {
                bucket = idx;
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

    [[nodiscard]] inline size_t size() const {
        return n_elms;
    }

    [[nodiscard]] inline bool empty()  const {
        return n_elms==0;
    }

    [[nodiscard]] inline size_t capacity() const{
        return m_table->size();
    };

    size_t table_mem_usage(){
        return m_table->size()*sizeof(bucket_t);
    }

    void destroy_table(){
        delete m_table;
        m_table = nullptr;
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
        ifs.read((char*)m_table->data(), long(sizeof(bucket_t)*n_buckets));
    }

    void store_table(std::ofstream & ofs) const{
        //serialize the table
        size_t n_buckets = table->size();
        ofs.write((char* )&n_buckets, sizeof(n_buckets));
        ofs.write((char *)table->data(), long(sizeof(bucket_t)*table->size()));
    }
};

#endif //LCG_LZL_MAP_H
