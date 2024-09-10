//
// Created by Diaz, Diego on 31.3.2023.
//

#ifndef SE_STRAT_TEXT_READER_H
#define SE_STRAT_TEXT_READER_H

#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <atomic>
#include "cds/vbyte_encoding.h"

#ifdef __linux__
#include <malloc.h>
#endif

struct text_stats{
    size_t max_sym{};
    size_t sep_sym{};
    std::vector<off_t> str_ptrs;
    std::vector<uint8_t> alphabet;
    off_t longest_str{};
};

template<class decoder_type,
        bool orig_text=false>
struct text_chunk {

    typedef decoder_type                 decoder_t;
    typedef typename decoder_t::sym_type sym_type;
    static constexpr bool first_round = orig_text;

    size_t id{}; //chunk id
    off_t n_bytes_before{};//number of bytes in the file before this chunk
    sym_type sep_sym{};//symbol in the buffer delimiting consecutive strings

    uint8_t *buffer = nullptr;// chunk's buffer
    off_t bytes_cap{};
    off_t bytes{}; //number of bytes the buffer can hold

    uint8_t *sec_buffer = nullptr; //chunk's parse
    off_t sec_bytes{}; //maximum number of bytes the sec. buffer can hold
    off_t acc_sec_bytes{}; //number of bytes inserted in the sec. buffer so far

    off_t *str_ptr = nullptr; //pointer to the leftmost string overlapping the buffer
    off_t n_str{}; //number of strings overlapping the buffer

    //hashing *perm_func = nullptr;
    uint64_t p_seed=0;
    std::vector<uint32_t>* vbyte_sym_perm=nullptr;
    std::atomic_long scans{0};

    text_chunk() = default;

    ~text_chunk() {
        if (buffer != nullptr) {
            //free(buffer);
            alloc<uint8_t>::deallocate(buffer);
        }

        if (sec_buffer != nullptr) {
            //free(sec_buffer);
            alloc<uint8_t>::deallocate(sec_buffer);
        }
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    [[nodiscard]] inline off_t rm_str() const {
        if (n_str == 0) return -1;
        return str_ptr[n_str - 1] - n_bytes_before;
    }

    inline off_t str_buff_start(off_t str) const {
        if (str < 0) return 0;
        if (str == n_str) return bytes;
        return str_ptr[str] - n_bytes_before;
    }

    inline off_t str_buff_end(off_t str) const {
        str++;
        if (str < 0) return 0;

        off_t pos;
        if (str == n_str){
            pos = bytes;
        }else{
            pos = str_ptr[str] - n_bytes_before;
        }

        if constexpr (orig_text){
            if(pos>0){
                //this code is to remove the separator symbol during the first round of parsing
                const uint8_t * l_boundary = buffer-1;
                uint8_t *ptr = buffer+pos-1;
                sym_type last_sym=0;
                off_t ps = read_backwards(ptr, l_boundary, last_sym)+1;
                if(last_sym==sep_sym) pos-=ps;
            }
        }
        return pos;
    }

    inline off_t eff_buff_bytes() const {
        return bytes;
    }

    inline off_t eff_mt_buff_bytes() const {
        return acc_sec_bytes;
    }

    /*template<class encoder_t>
    inline void append(typename encoder_t::sym_type &sym) {
        off_t b = encoder_t::write_forward(&sec_buffer[acc_sec_bytes], sym);
        acc_sec_bytes+=b;
        assert(acc_sec_bytes<=sec_bytes);
    }*/

    inline off_t read_forward(const uint8_t *ptr, sym_type &sym) const {
        off_t b = decoder_t::read_forward(ptr, sym);
        if constexpr (std::is_same<decoder_t, vbyte_decoder<sym_type>>::value){
            sym = (*vbyte_sym_perm)[sym];
        }
        return b;
    }

    inline off_t read_backwards(uint8_t *&ptr, const uint8_t *& l_boundary, sym_type& sym) const {
        off_t b = decoder_t::read_backwards(ptr, l_boundary, sym);
        if constexpr (std::is_same<decoder_t, vbyte_decoder<sym_type>>::value){
            sym = (*vbyte_sym_perm)[sym];
        }
        return b;
    }

    inline off_t mov_to_prev_diff_sym(uint8_t *&ptr, const uint8_t *l_boundary, sym_type &sym) const {
        off_t b = decoder_t::mov_to_prev_diff_sym(*&ptr, l_boundary, sym);
        if constexpr (std::is_same<decoder_t, vbyte_decoder<sym_type>>::value){
            sym = (*vbyte_sym_perm)[sym];
        }
        return b;
    }

    inline off_t mov_to_next_diff_sym(uint8_t *& ptr, const uint8_t* r_boundary, sym_type& sym) {
        off_t b = decoder_t::mov_to_next_diff_sym(ptr, r_boundary, sym);
        if constexpr (std::is_same<decoder_t, vbyte_decoder<sym_type>>::value){
            sym = (*vbyte_sym_perm)[sym];
        }
        return b;
    }

    void write_chunk_to_file(int fd, off_t& written_bytes){

        off_t block_bytes = 8388608;// 8MiB buffer
        off_t b = acc_sec_bytes;
        uint8_t * data = sec_buffer;

        for(off_t i=0;i<n_str;i++){
            *(str_ptr+i) += written_bytes;
        }

#ifdef __linux__
        off_t bytes_before = written_bytes;
#endif

        off_t n_blocks = b/block_bytes;
        for(off_t i=0;i<n_blocks;i++) {
            written_bytes+=write(fd, data, block_bytes);

#ifdef __linux__
            sync_file_range(fd, bytes_before+(i*block_bytes), block_bytes, SYNC_FILE_RANGE_WRITE);
        if(i){
            sync_file_range(fd, bytes_before+(i-1)*block_bytes, block_bytes, SYNC_FILE_RANGE_WAIT_BEFORE | SYNC_FILE_RANGE_WRITE | SYNC_FILE_RANGE_WAIT_AFTER);
        }
#endif
            data+=block_bytes;
        }

        off_t rem_bytes = b-(block_bytes*n_blocks);
        written_bytes+=write(fd, data, rem_bytes);

#ifdef __linux__
        if(n_blocks){
        sync_file_range(fd, bytes_before+(n_blocks-1)*block_bytes, block_bytes, SYNC_FILE_RANGE_WAIT_BEFORE | SYNC_FILE_RANGE_WRITE | SYNC_FILE_RANGE_WAIT_AFTER);
    }
    sync_file_range(fd, bytes_before + (n_blocks)*block_bytes, rem_bytes, SYNC_FILE_RANGE_WAIT_BEFORE | SYNC_FILE_RANGE_WRITE | SYNC_FILE_RANGE_WAIT_AFTER);
#endif

    }

    template<class parser_t>
    void read_chunk_from_file(int fd, off_t& rem_text_bytes, off_t& read_text_bytes, off_t*& global_str_ptr) {

        off_t chunk_bytes = bytes_cap<rem_text_bytes ? bytes_cap : rem_text_bytes;
        bytes = chunk_bytes;

        off_t acc_bytes = 0;
        off_t read_bytes;
        off_t fd_buff_bytes = 8388608;// 8MB buffer

        auto *data = (uint8_t *)buffer;

        while(chunk_bytes>0) {
            fd_buff_bytes = fd_buff_bytes<chunk_bytes ? fd_buff_bytes : chunk_bytes;
            read_bytes = read(fd, data, fd_buff_bytes);
            assert(read_bytes>0);

            data+=read_bytes;
            chunk_bytes-=read_bytes;
            acc_bytes+=read_bytes;
        }
        assert(bytes==acc_bytes);

        n_bytes_before = read_text_bytes;
        str_ptr = global_str_ptr;
        n_str = 0;

        //find the number of strings this chunk covers
        off_t rb = n_bytes_before+acc_bytes;
        auto * tmp_str_ptr = str_ptr;
        while(*tmp_str_ptr<rb) {
            tmp_str_ptr++;
            n_str++;
        }

        off_t offset, eff_bytes;
        if(*tmp_str_ptr!=rb) {//the buffer ends in the middle of a string, we need to trim

            off_t p_break = parser_t::rm_break(*this); //rightmost break in the buffer (returns -1 if there is none)

            if(p_break<0){ //no breaks in the buffer
                if(n_str==0 || rm_str()==0){
                    //double the buffer size, undo the reading and try again with a bigger buffer size
                    chunk_bytes = 2*bytes<rem_text_bytes ? 2*bytes : rem_text_bytes;
                    assert(chunk_bytes>bytes);
                    bytes = chunk_bytes;
                    //buffer = (uint8_t *)realloc(buffer, bytes);
                    buffer = alloc<uint8_t>::reallocate(buffer, bytes);
                    memset(buffer, 0,  bytes);

                    lseek(fd, acc_bytes*-1, SEEK_CUR);
                    return read_chunk_from_file<parser_t>(fd, rem_text_bytes, read_text_bytes, global_str_ptr);
                }else{
                    eff_bytes = rm_str();
                    bytes = eff_bytes;
                    //go one string back
                    n_str--;
                    tmp_str_ptr--;
                }
            }else{//trim the buffer one byte after the first break of the rightmost string starts
                eff_bytes = p_break+1;
                bytes = eff_bytes;
                //go backward in the buffer to consider the overlap of the parser
                eff_bytes = parser_t::overlap(*this, p_break);
            }
            global_str_ptr = tmp_str_ptr;
            offset = acc_bytes-eff_bytes;
            rem_text_bytes-= eff_bytes;
            read_text_bytes = lseek(fd, offset*-1, SEEK_CUR);
        }else{
            global_str_ptr = tmp_str_ptr;
            rem_text_bytes-=acc_bytes;
            read_text_bytes+=acc_bytes;
        }
    }
};

template<class sym_type>
void compute_text_stats(std::string& input_file, text_stats& txt_stats){

    off_t sym_bytes = sizeof(sym_type);
    off_t buff_size = 8388608; //8 MiB buffer
    off_t n_elms = buff_size/sym_bytes;

    //auto buffer = (sym_type *) calloc(n_elms, sym_bytes);
    auto buffer = alloc<sym_type>::allocate(n_elms);
    memset(buffer, 0, n_elms*sym_bytes);

    int fd = open(input_file.c_str(), O_RDONLY);

    struct stat st{};
    if(stat(input_file.c_str(), &st) != 0)  return;

#ifdef __linux__
    posix_fadvise(fd, 0, st.st_size, POSIX_FADV_SEQUENTIAL);
#endif

    lseek(fd, -1*sym_bytes, SEEK_END);
    read(fd, (char *)&txt_stats.sep_sym, sym_bytes);
    lseek(fd, 0, SEEK_SET);

    size_t read_bytes, read_syms; //str_len, , longest_str=0;
    off_t pos=0, cont=0, str_len;
    std::vector<uint8_t> alph(std::numeric_limits<sym_type>::max()+1, 0);

    while(true){
        read_bytes = read(fd, (char *)buffer, buff_size);
        if(read_bytes>0){
            read_syms = read_bytes/sym_bytes;
            for(size_t i=0;i<read_syms;i++){

                //if(buffer[i]>txt_stats.max_sym) txt_stats.max_sym = buffer[i];
                alph[buffer[i]] = 1;

                cont++;
                if(buffer[i]==txt_stats.sep_sym){
                    txt_stats.str_ptrs.push_back(pos*sym_bytes);
                    str_len = cont - pos;
                    if(str_len>txt_stats.longest_str){
                        txt_stats.longest_str=str_len;
                    }
                    pos = cont;
                }
            }
        }else{
            break;
        }
    }
    txt_stats.str_ptrs.push_back(pos*sym_bytes);

    for(size_t i=0;i<=std::numeric_limits<sym_type>::max();i++){
        if(alph[i]){
            txt_stats.alphabet.push_back(i);
        }
    }
    txt_stats.alphabet.shrink_to_fit();
    txt_stats.max_sym = txt_stats.alphabet.back();
    txt_stats.longest_str /=sizeof(sym_type);
    txt_stats.longest_str--; //we are not considering the separator symbol

#ifdef __linux__
    posix_fadvise(fd, 0, st.st_size, POSIX_FADV_DONTNEED);
#endif
    close(fd);
    alloc<sym_type>::deallocate(buffer);
}
#endif //SE_STRAT_TEXT_READER_H
