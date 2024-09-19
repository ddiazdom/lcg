//
// Created by Diaz, Diego on 18.9.2024.
//

#ifndef LCG_EXT_GRAM_H
#define LCG_EXT_GRAM_H

struct i_gram_stream{

    bitstream<size_t> buffer;
    partial_gram<uint8_t> comp_gram;
    int fd_r;
    std::vector<lvl_metadata_type> metadata;
    size_t lvl=0;
    size_t read_bytes=0;
    size_t curr_lvl=0;

    explicit i_gram_stream(int _fd_r) : fd_r(_fd_r){
        comp_gram.load_metadata_fd(fd_r);
        read_bytes = comp_gram.bytes_metadata();
        lvl = comp_gram.lvl;
        metadata = comp_gram.metadata;
    }

    void load_and_transform_curr_rule_set(uint32_t* lvl_rules, size_t stream_size, buff_vector<uint32_t>& map,
                                          vector_uint64_t& partition){

        size_t sym, pos=0, str_pos=1;
        comp_gram.load_next_rule_set_fd(fd_r, curr_lvl, buffer);
        uint32_t len;

        size_t n_words = bitstream<size_t>::bits2words(metadata[curr_lvl+1].n_bits());;
        read_bytes += bitstream<size_t>::words2bytes(n_words);

        off_t par=-1;
        size_t elm_per_block = INT_CEIL(metadata[curr_lvl+1].n_rules, partition.size());

        for(size_t j=0;j<metadata[curr_lvl+1].n_rules;j++){
            len=0;
            do{
                sym = buffer.read(pos, pos+metadata[curr_lvl+1].sym_width-1);

                lvl_rules[str_pos+len] = map[(sym>>1UL)-1];
                pos+=metadata[curr_lvl+1].sym_width;
                len++;
            } while(!(sym & 1UL));

            par+= (j % elm_per_block)==0;
            partition[par] = str_pos;
            assert(len<=151);

            lvl_rules[str_pos-1] = len;
            str_pos+=len+1;
        }
        curr_lvl++;
        assert((str_pos-1)==stream_size);
    }

    // Load the phrases in the level with their lengths.
    // It also stores pointers to a partition of the phrases,
    // so we can operate in parallel over the set
    template<class sym_type>
    void load_curr_rule_set(sym_type* rules_stream, size_t stream_size, vector_uint64_t & partition){

        size_t len_words;
        if constexpr (std::is_same<sym_type, uint8_t>::value){
            len_words=sizeof(uint32_t);
        }else{
            len_words=1;
        }

        size_t sym, pos=0, str_pos=len_words;
        assert(!partition.empty());
        comp_gram.load_next_rule_set_fd(fd_r, curr_lvl, buffer);
        uint32_t len;

        size_t n_words = bitstream<size_t>::bits2words(metadata[curr_lvl+1].n_bits());;
        read_bytes += bitstream<size_t>::words2bytes(n_words);

        off_t par=-1;
        size_t elm_per_block = INT_CEIL(metadata[curr_lvl+1].n_rules, partition.size());

        for(size_t j=0;j<metadata[curr_lvl+1].n_rules;j++){
            len=0;
            do{
                sym = buffer.read(pos, pos+metadata[curr_lvl+1].sym_width-1);
                rules_stream[str_pos+len] = sym>>1UL;
                if constexpr (std::is_same<sym_type, uint32_t>::value) {
                    rules_stream[str_pos+len]--;//the symbols are one-based to differentiate them from the sep symbol
                }
                assert(rules_stream[str_pos+len]<metadata[curr_lvl].n_rules);
                pos+=metadata[curr_lvl+1].sym_width;
                len++;
            } while(!(sym & 1UL));
            assert(len<=151);

            par+= (j % elm_per_block)==0;
            partition[par] = str_pos;

            if constexpr (std::is_same<sym_type, uint8_t>::value){
                memcpy(&rules_stream[str_pos], &len, sizeof(uint32_t));
            }else{
                rules_stream[str_pos] = len;
            }
            str_pos+=len+len_words;
        }
        //std::cout<<str_pos<<" "<<stream_size<<" "<<metadata[i+1].n_rules<<std::endl;
        curr_lvl++;
        assert((str_pos-len_words)==stream_size);
    }

    template<class sym_type>
    void compute_level_fps(sym_type* lvl_stream,
                           size_t stream_size,
                           buff_vector<uint64_t> & prev_fps,
                           buff_vector<uint64_t> & new_fps,
                           vector_uint64_t & partition) const {

        vector_uint64_t fp_seq(comp_gram.longest_rule);
        size_t pos = 0, end, tmp_pos, rule=0;//, tmp;
        uint32_t len;
        while(pos<stream_size){
            //tmp = pos;

            if constexpr (std::is_same<sym_type, uint8_t>::value){
                memcpy(&len, &lvl_stream[pos], sizeof(uint32_t));
                pos+=sizeof(uint32_t);
            }else{
                len = lvl_stream[pos];
                pos+=1;
            }

            assert(len>0);
            end = pos+len;
            tmp_pos=0;
            while(pos<end){
                //std::cout<<lvl_stream[pos]<<" ";
                fp_seq[tmp_pos++] = prev_fps[lvl_stream[pos++]];
            }
            //std::cout<<" ---> "<<len<<":"<<tmp_pos<<":"<<rule<<std::endl;
            assert(rule<new_fps.size());
            new_fps[rule++] = XXH3_64bits(fp_seq.data(), sizeof(uint64_t)*len);
        }
    }

    void load_compressed_string(buff_vector<uint32_t>& rules_stream){
        assert((curr_lvl+1)==lvl);
        size_t sym, pos=0, str_pos=0;
        comp_gram.load_next_rule_set_fd(fd_r, curr_lvl, buffer);
        rules_stream.resize(metadata[curr_lvl+1].tot_symbols);

        size_t n_words = bitstream<size_t>::bits2words(metadata[curr_lvl+1].n_bits());;
        read_bytes += bitstream<size_t>::words2bytes(n_words);

        do{
            sym = buffer.read(pos, pos+metadata[curr_lvl+1].sym_width-1);
            rules_stream[str_pos] = (sym>>1UL)-1;
            assert(rules_stream[str_pos]<metadata[curr_lvl].n_rules);
            pos+=metadata[curr_lvl+1].sym_width;
            str_pos++;
        } while(!(sym & 1UL));

        assert(str_pos==rules_stream.size());
        curr_lvl++;
        assert(curr_lvl==(metadata.size()-1));
    }

    lvl_metadata_type* metadata_curr_lvl(){
        return &metadata[curr_lvl+1];
    }

    size_t bytes_curr_level(){
        assert(curr_lvl<(lvl-1));
        if(curr_lvl==0){
            return metadata[1].tot_symbols + (metadata[1].n_rules*sizeof(uint32_t));
        }else{
            return sizeof(uint32_t)*(metadata[curr_lvl+1].tot_symbols+metadata[curr_lvl+1].n_rules);
        }
    }

    [[nodiscard]] size_t disk_usage() const{
        return comp_gram.bytes();
    }

    ~i_gram_stream(){
        buffer.destroy();
    }
};

#define WRITE_B\
    assert(str_pos_b<rule_set_b.size());          \
    assert(str_pos_c<rule_set_c.size());          \
    if constexpr (std::is_same<sym_type, uint8_t>::value){\
        memcpy(&len_b, &rule_set_b[str_pos_b], sizeof(uint32_t));\
        memcpy(&rule_set_c[str_pos_c], &len_b, sizeof(uint32_t));\
        str_pos_b+=sizeof(uint32_t);\
        str_pos_c+=sizeof(uint32_t);\
        memcpy(&rule_set_c[str_pos_c], &rule_set_b[str_pos_b], len_b);\
        str_pos_c+=len_b;\
        str_pos_b+=len_b;\
    } else {\
        len_b = rule_set_b[str_pos_b];\
        rule_set_c[str_pos_c] = len_b;\
        str_pos_b++;\
        str_pos_c++;\
        memcpy(&rule_set_c[str_pos_c], &rule_set_b[str_pos_b], len_b*sizeof(uint32_t));\
        str_pos_c+=len_b;\
        str_pos_b+=len_b;\
    }          \
    assert(len_b<=151);\
    fps_c.push_back(fps_b[rule_b]);\
    met_c.tot_symbols+=len_b;\
    merge_marks.push_back(in_b);\
    rule_b++;\

#define WRITE_A \
    assert(str_pos_a<rule_set_a.size());          \
    assert(str_pos_c<rule_set_c.size());          \
    if constexpr (std::is_same<sym_type, uint8_t>::value){\
        memcpy(&len_a, &rule_set_a[str_pos_a], sizeof(uint32_t));\
        memcpy(&rule_set_c[str_pos_c], &len_a, sizeof(uint32_t));\
        str_pos_a+=sizeof(uint32_t);\
        str_pos_c+=sizeof(uint32_t);\
        memcpy(&rule_set_c[str_pos_c], &rule_set_a[str_pos_a], len_a);\
        str_pos_c+=len_a;\
        str_pos_a+=len_a;\
    } else{\
        len_a = rule_set_a[str_pos_a];\
        rule_set_c[str_pos_c] = len_a;\
        str_pos_a++;\
        str_pos_c++;\
        memcpy(&rule_set_c[str_pos_c], &rule_set_a[str_pos_a], len_a*sizeof(uint32_t)); \
        str_pos_c+=len_a;\
        str_pos_a+=len_a;\
    }\
    assert(len_a<=151);\
    fps_c.push_back(fps_a[rule_a]);\
    met_c.tot_symbols+=len_a;\
    merge_marks.push_back(in_a);\
    rule_a++;\

struct extensible_gram {

    buff_vector<uint8_t> ter_rules;
    std::vector<buff_vector<uint32_t>> nt_rules;
    std::vector<buff_vector<uint64_t>> gram_fps;

    std::vector<lvl_metadata_type> metadata;
    size_t lvl=0;
    size_t longest_rule=0;

    void transform_level(size_t i, buff_vector<uint32_t>& map){
        assert(i>=2 && (i-2)<nt_rules.size());
        size_t pos = 0, end, rule=0;
        uint32_t len;
        while(pos<nt_rules[i-2].size()){
            len = nt_rules[i-2][pos++];
            assert(len>0);
            end = pos+len;
            while(pos<end){
                nt_rules[i-2][pos] = map[nt_rules[i-2][pos]];
                pos++;
            }
            rule++;
        }
        assert(rule==metadata[i].n_rules);
    }

    // Load the phrases in the level with their lengths.
    // It also stores pointers to a partition of the phrases,
    // so we can operate in parallel over the set
    template<class sym_type>
    void load_curr_rule_set(sym_type* rules_stream, size_t stream_size, size_t curr_lvl,
                            bitstream<size_t>& buffer, vector_uint64_t & partition){

        size_t len_words;
        if constexpr (std::is_same<sym_type, uint8_t>::value){
            len_words=sizeof(uint32_t);
        }else{
            len_words=1;
        }

        size_t sym, pos=0, str_pos=len_words;
        assert(!partition.empty());
        uint32_t len;

        off_t par=-1;
        size_t elm_per_block = INT_CEIL(metadata[curr_lvl+1].n_rules, partition.size());

        for(size_t j=0;j<metadata[curr_lvl+1].n_rules;j++){
            len=0;
            do{
                sym = buffer.read(pos, pos+metadata[curr_lvl+1].sym_width-1);
                rules_stream[str_pos+len] = sym>>1UL;
                if constexpr (std::is_same<sym_type, uint32_t>::value){
                    rules_stream[str_pos+len]--;//the symbols are one-based to differentiate them from the sep symbol
                }
                assert(rules_stream[str_pos+len]<metadata[curr_lvl].n_rules);
                pos+=metadata[curr_lvl+1].sym_width;
                len++;
            } while(!(sym & 1UL));

            assert(len<=151);

            par+= (j % elm_per_block)==0;
            partition[par] = str_pos;

            if constexpr (std::is_same<sym_type, uint8_t>::value){
                memcpy(&rules_stream[str_pos-len_words], &len, sizeof(uint32_t));
            }else{
                rules_stream[str_pos-len_words] = len;
            }
            str_pos+=len+len_words;
        }
        assert((str_pos-len_words)==stream_size);
    }

    void load_compressed_string(buff_vector<uint32_t>& rules_stream, size_t curr_lvl, bitstream<size_t>& buffer){
        assert((curr_lvl+1)==lvl);
        size_t sym, pos=0, str_pos=0;
        rules_stream.resize(metadata[curr_lvl+1].tot_symbols);
        do{
            sym = buffer.read(pos, pos+metadata[curr_lvl+1].sym_width-1);
            rules_stream[str_pos] = (sym>>1UL)-1;
            assert(rules_stream[str_pos]<metadata[curr_lvl].n_rules);
            pos+=metadata[curr_lvl+1].sym_width;
            str_pos++;
        } while(!(sym & 1UL));
        assert(str_pos==rules_stream.size());
    }

    size_t load_from_file(int fd_r, size_t n_threads=1) {

        partial_gram<uint8_t> comp_gram;
        comp_gram.load_metadata_fd(fd_r);
        size_t read_bytes = comp_gram.bytes_metadata();
        longest_rule = comp_gram.longest_rule;

        metadata = comp_gram.metadata;
        lvl = comp_gram.lvl;
        size_t curr_lvl = 0;
        bitstream<size_t> buffer;

        nt_rules.resize(metadata.size()-2);//minus the alphabet and the first level
        gram_fps.resize(metadata.size()-1);//minus the concatenated string
        vector_uint64_t partition(n_threads, 0);

        gram_fps[0].resize(metadata[0].tot_symbols);
        for(size_t sym=0;sym<gram_fps[0].size();sym++){
            gram_fps[0][sym] = XXH3_64bits(&sym, sizeof(uint8_t));
        }

        comp_gram.load_next_rule_set_fd(fd_r, curr_lvl, buffer);
        read_bytes+= comp_bytes_level(curr_lvl);
        ter_rules.resize(bytes_level(curr_lvl));
        load_curr_rule_set(ter_rules.data, ter_rules.size(), curr_lvl++, buffer, partition);

        gram_fps[1].resize(metadata[1].n_rules);
        compute_level_fps(ter_rules.data, ter_rules.size(), gram_fps[0], gram_fps[1], partition);

        for(size_t i=2;i<metadata.size()-1;i++){
            comp_gram.load_next_rule_set_fd(fd_r, curr_lvl, buffer);
            read_bytes+= comp_bytes_level(curr_lvl);
            nt_rules[i-2].resize(bytes_level(curr_lvl)/sizeof(uint32_t));
            load_curr_rule_set(nt_rules[i-2].data, nt_rules[i-2].size(), curr_lvl++, buffer, partition);

            gram_fps[i].resize(metadata[i].n_rules);
            compute_level_fps(nt_rules[i-2].data, nt_rules[i-2].size(), gram_fps[i-1], gram_fps[i], partition);
        }

        //read the concatenated strings
        comp_gram.load_next_rule_set_fd(fd_r, curr_lvl, buffer);
        read_bytes+= comp_bytes_level(curr_lvl);
        load_compressed_string(nt_rules.back(), curr_lvl++, buffer);
        buffer.destroy();
        assert(read_bytes==comp_gram.bytes());
        return read_bytes;
    }

    template<class sym_type>
    void compute_level_fps(sym_type* lvl_stream,
                           size_t stream_size,
                           buff_vector<uint64_t> & prev_fps,
                           buff_vector<uint64_t> & new_fps,
                           vector_uint64_t & partition) const {

        vector_uint64_t fp_seq(longest_rule);
        size_t pos = 0, end, tmp_pos, rule=0;//, tmp;
        uint32_t len;
        while(pos<stream_size){
            if constexpr (std::is_same<sym_type, uint8_t>::value){
                memcpy(&len, &lvl_stream[pos], sizeof(uint32_t));
                pos+=sizeof(uint32_t);
            }else{
                len = lvl_stream[pos];
                pos+=1;
            }

            assert(len>0 && len<=151);
            end = pos+len;
            tmp_pos=0;
            while(pos<end){
                fp_seq[tmp_pos++] = prev_fps[lvl_stream[pos++]];
            }
            assert(rule<new_fps.size());
            new_fps[rule++] = XXH3_64bits(fp_seq.data(), sizeof(uint64_t)*len);
        }
    }

    size_t bytes_level(size_t curr_lvl){
        assert(curr_lvl<(lvl-1));
        if(curr_lvl==0){
            return metadata[1].tot_symbols + (metadata[1].n_rules*sizeof(uint32_t));
        }else{
            return sizeof(uint32_t)*(metadata[curr_lvl+1].tot_symbols+metadata[curr_lvl+1].n_rules);
        }
    }

    size_t comp_bytes_level(size_t curr_lvl){
        size_t n_words = bitstream<size_t>::bits2words(metadata[curr_lvl+1].n_bits());;
        size_t read_bytes = bitstream<size_t>::words2bytes(n_words);
        return read_bytes;
    }

    size_t ext_gram_size_in_bytes(){
        size_t bytes = ter_rules.eff_mem_usage();
        for(auto const& nt : nt_rules){
            bytes+=nt.eff_mem_usage();
        }
        for(auto const& fp_set : gram_fps){
            bytes+= fp_set.eff_mem_usage();
        }
        bytes+=sizeof(lvl_metadata_type)*metadata.size();
        return bytes;
    }

    template<class sym_type>
    void merge_level(buff_vector<sym_type>& rule_set_a, lvl_metadata_type& mt_a, buff_vector<uint64_t> & fps_a, buff_vector<uint32_t>& map_a,
                     buff_vector<sym_type>& rule_set_b, lvl_metadata_type& mt_b, buff_vector<uint64_t> & fps_b, buff_vector<uint32_t>& map_b) {

        size_t rule_a=0, rule_b=0, str_pos_a=0, str_pos_b=0, str_pos_c=0;
        uint32_t len_a=0, len_b=0;
        lvl_metadata_type met_c;

        uint8_t in_a=1, in_b=2, in_ab=3;

        buff_vector<uint8_t> merge_marks;
        merge_marks.increase_capacity(mt_a.n_rules+mt_b.n_rules);
        buff_vector<uint64_t> fps_c;

        buff_vector<sym_type> rule_set_c;
        size_t buff_size=0, words_len;
        if constexpr (std::is_same<sym_type, uint8_t>::value){
            buff_size+=(mt_a.n_rules+mt_b.n_rules)*sizeof(uint32_t);
            words_len = sizeof(uint32_t);
        }else{
            words_len = 1;
            buff_size+=(mt_a.n_rules+mt_b.n_rules);
        }
        buff_size+=mt_a.tot_symbols+mt_b.tot_symbols;
        rule_set_c.resize(buff_size);

        while(rule_a < mt_a.n_rules && rule_b < mt_b.n_rules) {

            //TODO checking
            uint32_t la, lb, tmpa, tmpb;
            if constexpr (std::is_same<sym_type, uint8_t>::value){
                memcpy(&la, &rule_set_a[str_pos_a], sizeof(uint32_t));
                memcpy(&lb, &rule_set_b[str_pos_b], sizeof(uint32_t));
            }else{
                la = rule_set_b[str_pos_a];
                lb = rule_set_b[str_pos_b];
            }
            tmpa=str_pos_a+words_len;
            tmpb=str_pos_b+words_len;

            std::cout<< fps_a[rule_a] << " " << fps_b[rule_b] << " " << rule_a << " " << rule_b << std::endl;
            for(size_t j=0;j<la;j++){
                std::cout<<rule_set_a[tmpa+j]<<" ";
            }
            std::cout<<" --> "<<la<<" ? "<<rule_set_a.size()<<" "<<str_pos_a+len_a<<std::endl;
            for(size_t j=0;j<lb;j++){
                std::cout<<rule_set_b[tmpb+j]<<" ";
            }
            std::cout<<" --> "<<lb<<" ? "<<rule_set_b.size()<<" "<<str_pos_b+len_b<<std::endl;
            //

            if(fps_a[rule_a]<fps_b[rule_b]){
                WRITE_A
            } else if(fps_b[rule_b]<fps_a[rule_a]) {
                WRITE_B
            } else {
                if constexpr (std::is_same<sym_type, uint8_t>::value){
                    memcpy(&len_a, &rule_set_a[str_pos_a], sizeof(uint32_t));
                    memcpy(&len_b, &rule_set_b[str_pos_b], sizeof(uint32_t));
                }else{
                    len_a = rule_set_a[str_pos_a];
                    len_b = rule_set_b[str_pos_b];
                }
                assert(len_a<=151 && len_b<=151);
                assert((str_pos_a+words_len)<rule_set_a.size());
                assert((str_pos_b+words_len)<rule_set_b.size());
                assert((str_pos_c+words_len)<rule_set_c.size());
                str_pos_a+=words_len;
                str_pos_b+=words_len;
                bool eq = (len_a == len_b) && (memcmp(&rule_set_a[str_pos_a], &rule_set_b[str_pos_b], len_a * sizeof(sym_type)) == 0);

                if(eq){
                    if constexpr (std::is_same<sym_type, uint8_t>::value){
                        memcpy(&rule_set_c[str_pos_c], &len_a, sizeof(uint32_t));
                    } else {
                        rule_set_c[str_pos_c] = len_a;
                    }
                    str_pos_c+=words_len;
                    memcpy(&rule_set_c[str_pos_c], &rule_set_a[str_pos_a], len_a*sizeof(sym_type));
                    str_pos_a+=len_a;
                    str_pos_b+=len_a;
                    str_pos_c+=len_a;

                    fps_c.push_back(fps_a[rule_a]);
                    met_c.tot_symbols += len_a;
                    merge_marks.push_back(in_ab);
                    rule_a++;
                    rule_b++;
                } else {
                    //a collision occurred (extremely unlikely, but not impossible)

                    std::cout<< "Collision warning:  "<< fps_a[rule_a] << " " << fps_b[rule_b] << " " << rule_a << " " << rule_b << std::endl;
                    for(size_t j=0;j<len_a;j++){
                        std::cout<<rule_set_a[str_pos_a+j]<<" ";
                    }
                    std::cout<<" --> "<<len_a<<" ? "<<rule_set_a.size()<<" "<<str_pos_a+len_a<<std::endl;
                    for(size_t j=0;j<len_b;j++){
                        std::cout<<rule_set_b[str_pos_b+j]<<" ";
                    }
                    std::cout<<" --> "<<len_b<<" ? "<<rule_set_b.size()<<" "<<str_pos_b+len_b<<std::endl;

                    auto n_comparisons = (off_t)std::min(len_a, len_b);
                    off_t diff_pos=-1;
                    for(off_t j=0;j<n_comparisons;j++){
                        if(rule_set_a[str_pos_a+j]!=rule_set_b[str_pos_b+j]){
                            diff_pos = j;
                            break;
                        }
                    }
                    bool a_is_smaller = (diff_pos!=-1 && rule_set_a[str_pos_a+diff_pos]<rule_set_b[str_pos_b+diff_pos]) ||
                                        (diff_pos==-1 && len_a<len_b);

                    str_pos_a-=words_len;
                    str_pos_b-=words_len;
                    if(a_is_smaller) {
                        WRITE_A
                    } else {
                        WRITE_B
                    }
                }
            }
        }

        while(rule_a<mt_a.n_rules){
            WRITE_A
        }

        while(rule_b<mt_b.n_rules) {
            WRITE_B
        }

        rule_set_c.len = str_pos_c;
        rule_set_c.shrink_to_fit();
        rule_set_a.swap(rule_set_c);

        fps_c.shrink_to_fit();
        fps_a.swap(fps_c);
        fps_c.destroy();

        size_t mt_sym_a=0, mt_sym_b=0, mg_mt_sym=0;
        map_a.resize(mt_a.n_rules);
        map_a.shrink_to_fit();
        map_b.resize(mt_b.n_rules);
        map_b.shrink_to_fit();
        for(const unsigned char& side : merge_marks){
            switch (side) {
                case 1:
                    map_a[mt_sym_a++] = mg_mt_sym;
                    break;
                case 2:
                    map_b[mt_sym_b++] = mg_mt_sym;
                    break;
                default:
                    assert(side==in_ab);
                    map_a[mt_sym_a++] = mg_mt_sym;
                    map_b[mt_sym_b++] = mg_mt_sym;
            }
            mg_mt_sym++;
        }
        met_c.n_rules = merge_marks.size();
        merge_marks.destroy();
        mt_a = met_c;
        mt_a.terminals = false;
    }

    size_t merge_grams(i_gram_stream& ng, size_t n_threads) {

        vector_uint64_t ng_part(n_threads);
        size_t ng_bytes = ng.disk_usage();

        buff_vector<uint32_t> map_a;
        buff_vector<uint32_t> map_b;
        buff_vector<uint64_t> ng_fps;

        if(lvl!=ng.lvl){
            if(lvl<ng.lvl){
                nt_rules.resize(ng.lvl-1);
                metadata.resize(ng.lvl+1);
                gram_fps.resize(ng.lvl);
                nt_rules[lvl-2].swap(nt_rules.back());
                std::swap(metadata[lvl], metadata.back());
            }else {
                ng.metadata.resize(lvl);
            }
        }

        //max levels that are not the sequences of compressed strings
        size_t max_lvl_a = lvl-1;
        size_t max_lvl_b = ng.lvl-1;
        size_t max_lvl = std::max(max_lvl_a, max_lvl_b);
        lvl_metadata_type *ng_mt_lvl;
        buff_vector<uint32_t> ng_cmp_string;

        {
            buff_vector<uint8_t> ng_buff;
            ng_buff.resize(ng.bytes_curr_level());

            //load the first level
            ng_mt_lvl = ng.metadata_curr_lvl();
            ng.load_curr_rule_set(ng_buff.data, ng_buff.size(), ng_part);
            ng_fps.resize(ng_mt_lvl->n_rules);
            ng.compute_level_fps(ng_buff.data, ng_buff.size(), gram_fps[0], ng_fps, ng_part);
            merge_level(ter_rules, metadata[1], gram_fps[1], map_a, ng_buff, *ng_mt_lvl, ng_fps, map_b);
        }

        size_t mg_lvl=2;
        buff_vector<uint32_t> ng_buff;

        while(mg_lvl<=max_lvl){

            if(max_lvl_a<mg_lvl){
                create_new_level(nt_rules[mg_lvl-2], mg_lvl, nt_rules.back(), gram_fps[mg_lvl-1], 0, map_a);
            }else{
                transform_level(mg_lvl, map_a);
            }

            if(max_lvl_b<mg_lvl){
                create_new_level(ng_buff, mg_lvl, ng_cmp_string, gram_fps[mg_lvl-1], 0, map_b);
            }else {

                ng_buff.resize(ng.bytes_curr_level()/sizeof(uint32_t));
                ng_mt_lvl = ng.metadata_curr_lvl();
                ng.load_and_transform_curr_rule_set(ng_buff.data, ng_buff.size(), map_b, ng_part);
                ng_fps.resize(ng_mt_lvl->n_rules);
                ng.compute_level_fps(ng_buff.data, ng_buff.size(), gram_fps[mg_lvl-1], ng_fps, ng_part);

                if(max_lvl_b==mg_lvl){
                    ng.load_compressed_string(ng_cmp_string);
                }
            }

            merge_level(nt_rules[mg_lvl-2], metadata[mg_lvl], gram_fps[mg_lvl], map_a, ng_buff, *ng_mt_lvl, ng_fps, map_b);
            metadata[mg_lvl].sym_width = sym_width(metadata[mg_lvl-1].n_rules);
            mg_lvl++;
        }

        concatenate_strings(map_a, ng_cmp_string, map_b);
        metadata[mg_lvl].sym_width = sym_width(metadata[mg_lvl].n_rules);
        lvl = mg_lvl;

        assert(ng.read_bytes==ng_bytes);
        return ng.read_bytes;
    }

    void create_new_level(buff_vector<uint32_t>& rule_stream, size_t new_lvl,
                          buff_vector<uint32_t>& comp_strings,
                          buff_vector<uint64_t>& prev_fps, uint64_t fp_seed,
                          buff_vector<uint32_t>& map){

        assert(metadata.back().n_rules == map.size());
        // new_level is to the previous level in the metadata because
        // the first element of the metadata vector has the terminal alphabet
        metadata[new_lvl].sym_width = sym_width(metadata.back().n_rules);
        metadata[new_lvl].n_rules = metadata.back().n_rules;
        metadata[new_lvl].tot_symbols = metadata.back().n_rules;
        metadata[new_lvl].terminals = false;

        gram_fps[new_lvl].resize(metadata[new_lvl].n_rules);
        rule_stream.resize(map.size()*2);

        std::vector<std::tuple<uint64_t, uint64_t, uint32_t>> perm(map.size());
        for(size_t i=0;i<map.size();i++){
            std::get<0>(perm[i]) = i;
            uint64_t fp = prev_fps[map[i]];
            std::get<1>(perm[i]) = XXH3_64bits(&fp, sizeof(uint64_t));
            std::get<2>(perm[i]) = map[i];
        }

        std::sort(perm.begin(), perm.end(), [&](auto const& a, auto const &b) -> bool{
            if(std::get<1>(a)!=std::get<1>(b)){
                return std::get<1>(a) < std::get<1>(b);//break ties using the level fingerprint
            }
            assert(std::get<0>(a)==std::get<0>(b) ||
                   prev_fps[std::get<2>(a)]!=prev_fps[std::get<2>(b)]);
            return prev_fps[std::get<2>(a)]<prev_fps[std::get<2>(b)];
        });

        size_t pos=0;
        std::vector<uint32_t> inv_perm(perm.size());
        for(size_t mt_sym=0;mt_sym<perm.size();mt_sym++){
            rule_stream[pos++] = 1;
            rule_stream[pos++] = std::get<0>(perm[mt_sym]);
            gram_fps[new_lvl][mt_sym] = std::get<1>(perm[mt_sym]);
            map[std::get<0>(perm[mt_sym])] = std::get<2>(perm[mt_sym]);
            inv_perm[std::get<0>(perm[mt_sym])] = mt_sym;
        }

        //update the compressed string
        size_t mt_sym;
        pos = 0;
        while(pos<comp_strings.size()){
            mt_sym = comp_strings[pos];
            mt_sym = inv_perm[mt_sym];
            comp_strings[pos]  = mt_sym;
            pos++;
        }
    }

    void concatenate_strings(buff_vector<uint32_t>& map_a, buff_vector<uint32_t>& stream_b, buff_vector<uint32_t>& map_b){

        buff_vector<uint32_t>& stream = nt_rules.back();
        size_t stream_a_size = stream.size();
        metadata.back().tot_symbols += stream_b.size();
        stream.resize(nt_rules.back().size()+stream_b.size());

        size_t a_pos =0;
        while(a_pos<stream_a_size){
            assert(stream[a_pos]<map_a.size());
            stream[a_pos] = map_a[stream[a_pos]];
            //std::cout<<stream[a_pos]<<" "<<map_a.size()<<" "<<stream_a_size<<" from A "<<" "<<a_pos<<std::endl;
            a_pos++;
        }

        size_t b_pos=0;
        while(a_pos<stream.size()){
            assert(stream_b[b_pos]<map_b.size());
            stream[a_pos] = map_b[stream_b[b_pos]];
            //std::cout<<stream[a_pos]<<" "<<map_b.size()<<" "<<stream_b.size()<<" from B "<<" "<<a_pos<<" "<<b_pos<<std::endl;
            a_pos++;
            b_pos++;
        }
    }
};



#endif //LCG_EXT_GRAM_H
