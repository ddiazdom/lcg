//
// Created by Diaz, Diego on 17.9.2023.
//

#ifndef LCG_LC_PARSING_H
#define LCG_LC_PARSING_H

#include <fcntl.h>

//threshold to collapse the local grammars
//10GB
#define COL_THRESHOLD_1 10737418240
//100GB
#define COL_THRESHOLD_2 107374182400
//1000GB = 1TB
#define COL_THRESHOLD_3 1099511627776

#ifdef DEBUG_MODE
#include "malloc_count.h"
#endif

#include "xxhash.h"
#include "collapse_gram.h"
#include "fastx_parser.h"

#include "cds/ts_queue.h"
#include "cds/ts_priority_queue.h"
#include "cds/utils.h"
#include "cds/vbyte_encoding.h"

#define PARSING_INFO \
do{\
std::string msg1 ="Parsed_input "+report_space((off_t)p_state.f_proc_syms)+\
", malloc_peak "+report_space((off_t)malloc_count_peak())+                 \
", byte_usage "+report_space((off_t)text_chunks[buff_id].gram.mem_usage())+\
", input_fraction "+std::to_string(input_frac)+                            \
", fraction_limit "+std::to_string(p_state.max_frac);                      \
logger<msg_lvl>::debug(msg1);                                              \
                                                                           \
std::ostringstream oss;                                                    \
oss<<"  Processed input: "<<report_space((off_t)p_state.f_proc_syms)\
<<" ("<< std::fixed<< std::setprecision(2)<<((float(p_state.f_proc_syms)/float(p_state.f_size))*100)<<"%)    ";  \
logger<msg_lvl, true, true>::info(oss.str());\
}while(0);


struct parsing_state {

    uint8_t sep_sym;
    size_t n_threads=1;
    off_t page_cache_limit=0;
    size_t chunk_size=0;//the chunk size computed for the
    size_t n_chunks=0;
    int fd_r=-1;
    size_t f_size=0;
    off_t f_rem_bytes=0;
    off_t f_read_bytes=0;
    size_t f_proc_syms=0; //number of symbols processed in the file (it includes extra format symbols)
    size_t f_eff_proc_syms=0;//number of effective symbols processed in the file
    off_t r_page_cache_bytes=0;
    size_t chunk_id=0;
    size_t sink_gram_mem_usage=0;
    size_t factor_inc = 1;
    size_t empty_entries = 0; //number of emtpy entries in the file (only valid for FASTX files)
    float max_frac=0;
    text_format txt_fmt=PLAIN;

    //for inputs up to 10GB, we collapse every 20 processed chunks
    //for inputs over 10GB, each thread grammar uses up to 2.5% of the input
    //for inputs over 100GB, each thread grammar uses up to 1.5% of the input
    //for inputs over 1000GB, each thread grammar uses up to 0.6% of the input
    const float i_fracs[4] = {0.1, 0.025, 0.015, 0.006};

    parsing_state(uint8_t s_sym, size_t _n_threads, off_t p_cache_lim): sep_sym(s_sym),
                                                                        n_threads(_n_threads),
                                                                        page_cache_limit(p_cache_lim) {}

    void new_file(std::string& i_file, text_format _txt_fmt, size_t ck_size, float i_frac){

        if(fd_is_valid(fd_r)){

#ifdef __linux__
            posix_fadvise(fd_r, 0, f_size, POSIX_FADV_DONTNEED);
#endif
            close(fd_r);
            fd_r = -1;
        }

        txt_fmt = _txt_fmt;
        fd_r = open(i_file.c_str(), O_RDONLY);
        f_size = file_size(i_file);

#ifdef __linux__
        posix_fadvise(fd_r, 0, f_size, POSIX_FADV_SEQUENTIAL);
#endif
        f_eff_proc_syms = 0;
        f_proc_syms = 0;
        f_read_bytes = 0;
        f_rem_bytes = (off_t)f_size;
        empty_entries = 0;

        //if the chunk size is undefined, we set it to 0.5% of the input file
        chunk_size = ck_size==0 ? size_t(ceil(0.005 * double(f_size))) : ck_size;
        // we put a cap on the chunk size of 200 MiB by design
        // (bigger chunks only use more RAM and do not speed up the compression process)
        chunk_size = std::min<size_t>(chunk_size, 1024*1024*200);

        //the cap for the number of threads is the number of chunks
        size_t tot_chunks = INT_CEIL(f_size, chunk_size);
        n_threads = std::min(n_threads, tot_chunks);

        n_chunks = n_threads+1;

        //we cannot set more chunks thant the total number of chunks in the input file
        n_chunks = std::min<unsigned long>(n_chunks, tot_chunks);

        //compute the approx. memory usage for the local grammar
        if(i_frac==0){
            //compute the rule to collapse local grammars into the sink grammar
            if(f_size<=COL_THRESHOLD_1){
                max_frac = i_fracs[0];
            } else if(f_size<=COL_THRESHOLD_2){
                max_frac = i_fracs[1];
            } else if(f_size<=COL_THRESHOLD_3){
                max_frac = i_fracs[2];
            } else{
                max_frac = i_fracs[3];
            }
            //The algorithm collapses the local grammars when:
            //the sum of the space usage of the local grammars exceed f_size*i_frac
            max_frac *=float(n_chunks);
        }
    }

    [[nodiscard]] inline size_t chunks_approx_mem_usage() const {
        return size_t(((chunk_size*115)/100)*n_chunks);
    }

    ~parsing_state(){
#ifdef __linux__
        posix_fadvise(fd_r, 0, f_size, POSIX_FADV_DONTNEED);
#endif
        close(fd_r);
    }

    void flush_page_cache(){
#ifdef __linux__
        std::cout<<"removing from page cache "<<r_page_cache_bytes<<" "<<f_read_bytes<<std::endl;
        posix_fadvise(fd_r, f_read_bytes-r_page_cache_bytes, r_page_cache_bytes, POSIX_FADV_DONTNEED);
        r_page_cache_bytes=0;
#endif
    }
};

struct phrase_overflow{
    uint32_t position;
    uint32_t length;
    uint32_t metasymbol;

    [[nodiscard]] inline uint32_t right_end() const {
        return position+length-1;
    }
};

void finish_byte_parse(text_chunk& chunk, off_t &parse_distance, std::vector<phrase_overflow>& phr_with_ovf){

    uint32_t mt_sym=1;
    off_t txt_size = chunk.e_bytes;

    off_t ovf_idx=0, next_ovf=-1;
    if(!phr_with_ovf.empty()){
        next_ovf = phr_with_ovf[ovf_idx].right_end();
    }

    off_t pos = txt_size-1;
    chunk.parse = (uint32_t *)(chunk.text+parse_distance);
    chunk.parse--;

    while(pos>=0){
        if(pos==next_ovf){
            pos -= phr_with_ovf[ovf_idx].length;
            *chunk.parse = phr_with_ovf[ovf_idx].metasymbol;
            next_ovf=-1;
            if(ovf_idx<off_t(phr_with_ovf.size()-1)){
                next_ovf = phr_with_ovf[++ovf_idx].right_end();
            }
        }  else {
            pos-=vbyte_decoder<uint32_t>::read_right2left(&chunk.text[pos], mt_sym);
            *chunk.parse = mt_sym;
        }
        chunk.parse--;
        while(pos>next_ovf && chunk.text[pos]==0) pos--;
    }

    chunk.parse++;
    assert(pos==-1);
    assert(chunk.parse[chunk.parse_size-1]==0);
    assert((uintptr_t) chunk.text<= (uintptr_t)chunk.parse && (uintptr_t)chunk.parse<= (uintptr_t)&chunk.text[txt_size-1]);
    assert((uintptr_t) chunk.text<= (uintptr_t)(chunk.parse+chunk.parse_size) &&
           (uintptr_t) (chunk.parse+chunk.parse_size) <= (uintptr_t)(chunk.text+chunk.buffer_bytes));
    assert(parse_distance==(off_t)chunk.dist(reinterpret_cast<uint8_t*>(chunk.parse+chunk.parse_size)));
}

template<bool query_sink>
void byte_par_r2l(text_chunk& chunk, off_t& n_strings, size_t sep_sym);

template<>
void byte_par_r2l<true>(text_chunk& chunk, off_t& n_strings, size_t sep_sym) {

    size_t next_av_in_sink = chunk.sink_gram.ter_dict.size();
    uint64_t * fps = chunk.gram.fps[chunk.round];
    uint8_t * text = chunk.text;
    uint64_t hash;
    off_t lb, rb = chunk.e_bytes-1, i=chunk.e_bytes-2, byte_offset, parse_size;
    uint8_t v_len;
    assert(text[i+1]==sep_sym && text[i]>text[i+1]);

    text[rb] = 128;//the vbyte code of the metasymbol mt=0 representing a separator in the next round of parsing
    n_strings=1;
    parse_size=4;//we will count in bytes of sizeof(uint32_t)
    off_t mbo[2]={4};

    std::vector<phrase_overflow> phr_with_ovf;
    bool r_cmp = true, l_cmp, new_str;
    uint32_t mt_sym, phrase_len;
    uint8_t mid_sym = text[i];
    while(--i>0 && text[i]==mid_sym);

    while(i>=0){
        l_cmp = fps[text[i]]>fps[mid_sym];
        if(l_cmp && !r_cmp){
            lb = i+1;
            new_str = mid_sym==sep_sym;
            lb += new_str;
            phrase_len=rb-lb;

            byte_offset = rb + parse_size;
            mbo[byte_offset > mbo[1]] = byte_offset;

            hash = XXH3_64bits(&text[lb], phrase_len);
            bool found = chunk.sink_gram.ter_dict.find(&text[lb], phrase_len, mt_sym, hash);
            if(!found){
                mt_sym = next_av_in_sink + chunk.gram.ter_dict.insert(&text[lb], phrase_len, hash);
            }
            mt_sym++;

            v_len = vbyte_len(mt_sym);
            if(__builtin_expect(v_len>phrase_len, 0)){
                //metasymbol does not fit its phrase
                phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
            } else {
                vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
                memset(&text[lb+v_len], 0, phrase_len-v_len);
            }

            parse_size+=4;
            rb = i+1;
            if(new_str){
                assert(text[rb]==sep_sym);
                text[rb] = 128;//vbyte code for 0 (the separator symbol in the next levels)
                n_strings++;
                parse_size+=4;
            }
        }
        r_cmp = l_cmp;
        mid_sym = text[i];
        while(--i>0 && text[i]==mid_sym);
    }

    lb = 0;
    phrase_len=rb-lb;
    byte_offset = rb + parse_size;
    mbo[byte_offset > mbo[1]] = byte_offset;

    hash = XXH3_64bits(&text[lb], phrase_len);
    bool found = chunk.sink_gram.ter_dict.find(&text[lb], phrase_len, mt_sym, hash);
    if(!found){
        mt_sym = next_av_in_sink + chunk.gram.ter_dict.insert(&text[lb], phrase_len, hash);
    }
    mt_sym++;

    v_len = vbyte_len(mt_sym);
    if(__builtin_expect(v_len>phrase_len, 0)){
        phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
    }else {
        vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
        memset(&text[lb+v_len], 0, phrase_len-v_len);
    }
    parse_size+=4;
    chunk.parse_size = parse_size/4;//divide by 4=sizeof(uint32_t) to avoid the <<

    //computes how many bytes we require for the parse
    off_t max_byte_offset = INT_CEIL(mbo[1], sizeof(uint32_t))*sizeof(uint32_t);
    chunk.increase_capacity(max_byte_offset);

    chunk.gram.update_fps_with_sink(chunk.round, chunk.sink_gram);
    finish_byte_parse(chunk, max_byte_offset, phr_with_ovf);
}

template<>
void byte_par_r2l<false>(text_chunk& chunk, off_t& n_strings, size_t sep_sym) {

    uint64_t * fps = chunk.gram.fps[chunk.round];
    uint8_t * text = chunk.text;
    off_t lb, rb = chunk.e_bytes-1, i=chunk.e_bytes-2, byte_offset, parse_size;

    uint8_t v_len;
    assert(text[i+1]==sep_sym && text[i]>text[i+1]);

    text[rb] = 128;//the vbyte code of the metasymbol mt=0 representing a separator in the next round of parsing
    n_strings=1;
    parse_size=4;//we will count in bytes of sizeof(uint32_t)
    off_t mbo[2]={4};

    std::vector<phrase_overflow> phr_with_ovf;
    bool r_cmp = true, l_cmp, new_str;
    uint32_t mt_sym, phrase_len;
    uint8_t mid_sym = text[i];
    while(--i>0 && text[i]==mid_sym);

    while(i>=0){
        l_cmp = fps[text[i]]>fps[mid_sym];
        if(l_cmp && !r_cmp){
            lb = i+1;
            new_str = mid_sym==sep_sym;
            lb += new_str;
            phrase_len=rb-lb;

            byte_offset = rb + parse_size;
            mbo[byte_offset > mbo[1]] = byte_offset;
            mt_sym = chunk.gram.ter_dict.insert(&text[lb], phrase_len)+1;

            v_len = vbyte_len(mt_sym);
            if(__builtin_expect(v_len>phrase_len, 0)){
                //metasymbol does not fit its phrase
                phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
            } else {
                vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
                memset(&text[lb+v_len], 0, phrase_len-v_len);
            }

            parse_size+=4;
            rb = i+1;
            if(new_str){
                assert(text[rb]==sep_sym);
                text[rb] = 128;//vbyte code for 0 (the separator symbol in the next levels)
                n_strings++;
                parse_size+=4;
            }
        }
        r_cmp = l_cmp;
        mid_sym = text[i];
        while(--i>0 && text[i]==mid_sym);
    }

    lb = 0;
    phrase_len=rb-lb;
    byte_offset = rb + parse_size;
    mbo[byte_offset > mbo[1]] = byte_offset;
    mt_sym = chunk.gram.ter_dict.insert(&text[lb], phrase_len) + 1;

    v_len = vbyte_len(mt_sym);
    if(__builtin_expect(v_len>phrase_len, 0)){
        phr_with_ovf.push_back({uint32_t(lb), phrase_len, mt_sym});
    }else {
        vbyte_decoder<uint32_t>::write_right2left(&text[lb], mt_sym, v_len);
        memset(&text[lb+v_len], 0, phrase_len-v_len);
    }
    parse_size+=4;
    chunk.parse_size = parse_size/4;//divide by 4=sizeof(uint32_t) to avoid the << shift during the scan

    //computes how many bytes we require for the parse
    off_t max_byte_offset = INT_CEIL(mbo[1], sizeof(uint32_t))*sizeof(uint32_t);
    chunk.increase_capacity(max_byte_offset);

    chunk.gram.update_fps(chunk.round);
    finish_byte_parse(chunk, max_byte_offset, phr_with_ovf);
}

template<bool query_sink>
void int_par_l2r(text_chunk& chunk);

template<>
void int_par_l2r<true>(text_chunk& chunk){

    assert(chunk.round>0);
    const uint64_t* fps[2] = {chunk.sink_gram.fps[chunk.round], chunk.gram.fps[chunk.round]};

    uint64_t hash;
    uint32_t next_av_mt_in_sink = chunk.sink_gram.nt_dicts[chunk.round-1].size();
    uint32_t alpha_sink = chunk.sink_gram.alphabet(chunk.round);
    uint64_t sym_offset[2] = {0, alpha_sink};

    uint32_t *text = chunk.parse;
    uint32_t mt_sym, sep_sym=0, txt_size = chunk.parse_size;
    uint32_t left_sym, middle_sym;
    uint64_t left_fp, middle_fp, right_fp;
    off_t i=0, parse_size = 0, phrase_len, lb, rb;
    bool new_str=false, left_is_new, middle_is_new, right_is_new, phrase_is_new;

    lb = 0;
    left_sym = text[i];
    left_is_new = left_sym>alpha_sink;
    phrase_is_new = left_is_new;

    left_fp = fps[left_is_new][left_sym-sym_offset[left_is_new]];

    while(++i<txt_size && text[i]==left_sym);
    assert(i<txt_size);

    middle_sym = text[i];
    middle_is_new = middle_sym>alpha_sink;
    middle_fp = fps[middle_is_new][middle_sym-sym_offset[middle_is_new]];
    rb=i;
    while(++i<txt_size && text[i]==middle_sym);

    while(i<txt_size) {
        right_is_new = text[i]>alpha_sink;
        right_fp = fps[right_is_new][text[i]-sym_offset[right_is_new]];

        if(left_fp>middle_fp && middle_fp<right_fp){//local minimum
            phrase_len = rb-lb;

            //TODO remove this block later
            /*bool tmp=false;
            for(off_t k=lb;k<lb+phrase_len;k++){
                tmp+=text[k]>alpha_sink;
            }
            assert(tmp==phrase_is_new);*/
            //

            hash = XXH3_64bits(&text[lb], phrase_len*sizeof(uint32_t));
            if(phrase_is_new || !chunk.sink_gram.nt_dicts[chunk.round-1].find(&text[lb], phrase_len, mt_sym, hash)){
                mt_sym = next_av_mt_in_sink + chunk.gram.nt_dicts[chunk.round-1].insert(&text[lb], phrase_len, hash);
            }
            text[parse_size]=sep_sym;
            parse_size+=new_str;
            text[parse_size++] = mt_sym+1;
            new_str = text[rb]==sep_sym;
            lb = rb+new_str;
            phrase_is_new=false;
        }

        phrase_is_new |= middle_is_new;
        left_fp = middle_fp;

        middle_fp = right_fp;
        middle_is_new = right_is_new;
        middle_sym = text[i];

        rb = i;
        while(++i<txt_size && text[i]==middle_sym);
    }
    assert(rb==(txt_size-1) && text[rb]==sep_sym);

    phrase_is_new |= middle_is_new;
    phrase_len = rb-lb;

    //TODO remove this block later
    /*bool tmp=false;
    for(off_t k=lb;k<lb+phrase_len;k++){
        tmp+=text[k]>alpha_sink;
    }
    assert(tmp==phrase_is_new);*/
    //

    hash = XXH3_64bits(&text[lb], phrase_len*sizeof(uint32_t));
    if(phrase_is_new || !chunk.sink_gram.nt_dicts[chunk.round-1].find(&text[lb], phrase_len, mt_sym, hash)){
        mt_sym = next_av_mt_in_sink+chunk.gram.nt_dicts[chunk.round-1].insert(&text[lb], phrase_len, hash);
    }

    text[parse_size]=sep_sym;
    parse_size+=new_str;
    text[parse_size++] = mt_sym+1;
    text[parse_size++] = sep_sym;

    chunk.gram.update_fps_with_sink(chunk.round, chunk.sink_gram);
    chunk.parse_size = parse_size;
}

template<>
void int_par_l2r<false>(text_chunk& chunk) {

    uint64_t *fps = chunk.gram.fps[chunk.round];
    uint32_t *text = chunk.parse;

    uint32_t mt_sym, sep_sym=0, txt_size = chunk.parse_size;
    uint32_t left_sym, middle_sym;
    uint64_t left_fp, middle_fp, right_fp;
    off_t i=0, parse_size = 0, phrase_len, lb, rb;
    bool new_str=false;

    lb = 0;
    left_sym = text[i];
    left_fp = fps[left_sym];
    while(++i<txt_size && text[i]==left_sym);
    assert(i<txt_size);

    middle_sym = text[i];
    middle_fp = fps[middle_sym];
    rb=i;
    while(++i<txt_size && text[i]==middle_sym);

    while(i<txt_size) {
        right_fp = fps[text[i]];

        if(left_fp>middle_fp && middle_fp<right_fp){//local minimum
            phrase_len = rb-lb;
            mt_sym = chunk.gram.nt_dicts[chunk.round-1].insert(&text[lb], phrase_len);

            text[parse_size]=sep_sym;
            parse_size+=new_str;
            text[parse_size++] = mt_sym+1;
            new_str = text[rb]==sep_sym;
            lb = rb+new_str;
        }
        left_fp = middle_fp;
        middle_fp = right_fp;
        middle_sym = text[i];
        rb = i;
        while(++i<txt_size && text[i]==middle_sym);
    }
    assert(rb==(txt_size-1) && text[rb]==sep_sym);

    phrase_len = rb-lb;
    mt_sym = chunk.gram.nt_dicts[chunk.round-1].insert(&text[lb], phrase_len);

    text[parse_size]=sep_sym;
    parse_size+=new_str;
    text[parse_size++] = mt_sym+1;
    text[parse_size++] = sep_sym;

    chunk.gram.update_fps(chunk.round);
    chunk.parse_size = parse_size;
}

template<bool query_sink>
void compress_text_chunk(text_chunk& chunk){

    off_t n_strings=0;
    size_t sep_sym = chunk.sep_sym;
    chunk.round = 0;

    if(chunk.format==FASTA){
        chunk.e_bytes = PARSE_FASTA(chunk.text, chunk.e_bytes, chunk.empty_entries);
    }

    byte_par_r2l<query_sink>(chunk, n_strings, sep_sym);

    off_t size_limit = n_strings*2;
    chunk.round++;

    while(chunk.parse_size!=size_limit){
        assert(chunk.parse_size>=size_limit);
        int_par_l2r<query_sink>(chunk);
        chunk.round++;
    }

    //the chunks are small, so they do not can hold many strings
    assert((chunk.parse_size>>1)<=0xFFFFFFF);
    chunk.gram.str_orders.emplace_back(chunk.id, (uint8_t)chunk.round,
                                       chunk.gram.comp_string.size(), (uint32_t)(chunk.parse_size/2));
    size_t pos = chunk.gram.comp_string.size();
    chunk.gram.comp_string.resize(pos+(chunk.parse_size/2));
    for(off_t i=0;i<chunk.parse_size;i+=2){
        assert(i==0 || chunk.parse[i-1]==0);
        chunk.gram.comp_string[pos++] = chunk.parse[i];
    }
}

template<bool query_sink, log_lvl msg_lvl=INFO>
void fill_chunk_grammars(std::vector<text_chunk>& text_chunks, parsing_state& p_state){

    ts_queue<size_t> buffers_to_process;
    ts_queue<size_t> buffers_to_reuse;
    std::atomic<size_t> parser_finished{0};

    auto compressor_worker = [&]() {
        size_t buff_id;
        bool res;

        while (true) {
            res = buffers_to_process.pop(buff_id);
            if (!res){
                parser_finished.fetch_add(1, std::memory_order_acq_rel);
                break;
            }
            assert(text_chunks[buff_id].text_bytes > 0);

            text_chunks[buff_id].t_start = std::chrono::steady_clock::now();
            compress_text_chunk<query_sink>(text_chunks[buff_id]);
            memset(text_chunks[buff_id].text, 0, text_chunks[buff_id].buffer_bytes);

            //note: e_bytes change when the format is fastx
            p_state.f_eff_proc_syms += text_chunks[buff_id].e_bytes;
            p_state.empty_entries += text_chunks[buff_id].empty_entries;

            text_chunks[buff_id].t_end = std::chrono::steady_clock::now();

            buffers_to_reuse.push(buff_id);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(p_state.n_threads);
    for (size_t i = 0; i < p_state.n_threads; i++) {
        threads.emplace_back(compressor_worker);
    }

    auto tmp_ck_size = off_t(INT_CEIL(p_state.chunk_size, sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type));
    size_t buff_id = 0;

    while(buff_id < text_chunks.size() && p_state.f_rem_bytes > 0) {
        tmp_ck_size = std::min(tmp_ck_size, p_state.f_rem_bytes);
        text_chunks[buff_id].text_bytes = tmp_ck_size;
        text_chunks[buff_id].sep_sym = (text_chunk::size_type) p_state.sep_sym;
        text_chunks[buff_id].increase_capacity((tmp_ck_size*115)/100);
        text_chunks[buff_id].id = p_state.chunk_id++;

        read_chunk_from_file(p_state.fd_r, p_state.f_rem_bytes, p_state.f_read_bytes, text_chunks[buff_id]);
        buffers_to_process.push(buff_id);

#ifdef __linux__
        p_state.r_page_cache_bytes+=text_chunks[buff_id].text_bytes;
        if(p_state.r_page_cache_bytes>p_state.page_cache_limit){
            p_state.flush_page_cache();
        }
#endif
        buff_id++;
    }

    std::vector<size_t> byte_counts(text_chunks.size(), 0);
    size_t acc_bytes = p_state.sink_gram_mem_usage;
    float input_frac = float(acc_bytes)/float(p_state.f_size);

    // if the sink grammar already exceeds the input fraction limit,
    // we simulate to add a new text chunk
    if(input_frac>p_state.max_frac){
        float inc_frac = (p_state.max_frac/float(text_chunks.size()))*(float)p_state.factor_inc;
        p_state.max_frac= input_frac+inc_frac;
        p_state.factor_inc++;
    }

    while (p_state.f_rem_bytes > 0 && input_frac<p_state.max_frac) {

        buffers_to_reuse.pop(buff_id);
        p_state.f_proc_syms+=text_chunks[buff_id].text_bytes;

        size_t new_byte_count = text_chunks[buff_id].gram.mem_usage();
        acc_bytes -= byte_counts[buff_id];
        acc_bytes += new_byte_count;
        byte_counts[buff_id] = new_byte_count;
        input_frac = float(acc_bytes)/float(p_state.f_size);

        PARSING_INFO

        text_chunks[buff_id].text_bytes = tmp_ck_size;
        text_chunks[buff_id].id = p_state.chunk_id++;
        read_chunk_from_file(p_state.fd_r, p_state.f_rem_bytes, p_state.f_read_bytes, text_chunks[buff_id]);
        buffers_to_process.push(buff_id);

#ifdef __linux__
        p_state.r_page_cache_bytes+=text_chunks[buff_id].text_bytes;
        if(p_state.r_page_cache_bytes>p_state.page_cache_limit){
            p_state.flush_page_cache();
        }
#endif
    }

    //wait for the queue to be empty and close the input file
    while (!buffers_to_process.empty());
    buffers_to_process.done();

    //wait for all the parsers to finish
    while(parser_finished.load(std::memory_order_acquire)!=threads.size());

    while(!buffers_to_reuse.empty()){
        buffers_to_reuse.pop(buff_id);
        p_state.f_proc_syms+=text_chunks[buff_id].text_bytes;

        size_t new_byte_count = text_chunks[buff_id].gram.mem_usage();
        acc_bytes -= byte_counts[buff_id];
        acc_bytes += new_byte_count;
        byte_counts[buff_id] = new_byte_count;
        input_frac = float(acc_bytes)/float(p_state.f_size);

        PARSING_INFO

#ifdef __linux__
        p_state.r_page_cache_bytes+=text_chunks[buff_id].text_bytes;
        if(p_state.r_page_cache_bytes>p_state.page_cache_limit){
            p_state.flush_page_cache();
        }
#endif
    }
    buffers_to_reuse.done();

    for (auto &thread: threads) {
        thread.join();
    }
}

template<log_lvl msg_lvl>
void process_one_file(parsing_state& par_state, plain_gram& sink_gram){

    std::vector<plain_gram> ck_grams(par_state.n_chunks,
                                     plain_gram(sink_gram.lvl_cap(), sink_gram.sep_sym()));

    std::vector<text_chunk> txt_chunks;
    txt_chunks.reserve(par_state.n_chunks);
    for(size_t i=0;i<par_state.n_chunks;i++){
        txt_chunks.emplace_back(sink_gram, ck_grams[i], par_state.txt_fmt);
    }

    if(sink_gram.empty()){
        fill_chunk_grammars<false, msg_lvl>(txt_chunks, par_state);
    }else{
        fill_chunk_grammars<true, msg_lvl>(txt_chunks, par_state);
    }

    collapse_grams<msg_lvl>(sink_gram, ck_grams);
    sink_gram.update_fps();
    par_state.sink_gram_mem_usage=sink_gram.eff_mem_usage();

    while(par_state.f_rem_bytes>0){
        fill_chunk_grammars<true, msg_lvl>(txt_chunks, par_state);
        collapse_grams<msg_lvl>(sink_gram, ck_grams);
        sink_gram.update_fps();
        par_state.sink_gram_mem_usage=sink_gram.eff_mem_usage();
    }
    logger<msg_lvl, false, true>::info(" ");
    sink_gram.text_size += par_state.f_eff_proc_syms;

    if(par_state.empty_entries){
        logger<msg_lvl>::warning("file contains empty entries");
    }
}

template<log_lvl msg_lvl>
void build_lc_gram(std::vector<std::string>& i_files, text_format& txt_fmt, plain_gram& sink_gram, size_t n_threads, off_t chunk_size, float i_frac) {

    logger<msg_lvl>::info("Building a locally consistent grammar");

    // maximum number of bytes we accept in the page cache when reading the input file.
    // Once we exceed this threshold. we ask the kernel to remove the last file pages
    // from the cache. This option only applies on linux
    off_t page_cache_limit = 1024*1024*1024;
    parsing_state par_state(sink_gram.sep_sym(), n_threads, page_cache_limit);

    for(auto & i_file : i_files){
        auto start = std::chrono::steady_clock::now();
        par_state.new_file(i_file, txt_fmt, chunk_size, i_frac);
        std::string format = (txt_fmt==PLAIN? "Plain" : (txt_fmt==FASTA? "Fasta" : "Fastq"));
        std::string msg ="  File                      : "+i_file+"\n";
        msg+="  Format                    : "+format+"\n";
        msg+="  Parsing threads           : "+std::to_string(par_state.n_threads)+"\n";
        msg+="  Active text chunks in RAM : "+std::to_string(par_state.n_chunks)+"\n";
        msg+="  Size of each chunk        : "+report_space((off_t)par_state.chunk_size)+"\n";
        msg+="  Chunks' approx. mem usage : "+report_space((off_t)par_state.chunks_approx_mem_usage());
        logger<msg_lvl>::info(msg);
        process_one_file<msg_lvl>(par_state, sink_gram);
        auto end = std::chrono::steady_clock::now();
        logger<msg_lvl>::info(report_time(start, end, 2));
    }

    sink_gram.reorder_strings();
    sink_gram.clear_fps();
}
#endif //LCG_LC_PARSING_H
