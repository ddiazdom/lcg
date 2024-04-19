//
// Created by diego on 17-02-20.
//

#ifndef LMS_COMPRESSOR_FILE_BUFFER_H
#define LMS_COMPRESSOR_FILE_BUFFER_H

#include <fstream>
#include <cassert>
#include <iostream>
#include <cstring>
#include <limits>
#include "bitstream.h"

template<class sym_t>
struct i_file_stream{

    typedef sym_t sym_type;
    std::ifstream text_i;
    size_t tot_cells{};
    size_t block_bg{};
    size_t block_bits{};
    bitstream<sym_t> buffer;

    std::string file;
    static constexpr size_t w_bytes = sizeof(sym_t);
    constexpr static uint8_t w_bits = std::numeric_limits<sym_t>::digits;
    constexpr static uint8_t w_shift = __builtin_ctz(w_bits);

    i_file_stream(i_file_stream<sym_t> &&other) noexcept{
        move(std::forward<i_file_stream<sym_t>>(other));
    }

    i_file_stream(const std::string& i_file, size_t buff_size_){

        file =  i_file;
        text_i = std::ifstream(file, std::ifstream::binary);
        assert(text_i.good());

        text_i.seekg (0, std::ifstream::end);
        size_t file_size =  text_i.tellg();
        text_i.seekg (0, std::ifstream::beg);

        size_t block_bytes = INT_CEIL(buff_size_, w_bytes)*w_bytes;
        buffer.stream_size = block_bytes/w_bytes;
        buffer.stream = new sym_t[buffer.stream_size];
        block_bits = block_bytes*8;

        tot_cells = INT_CEIL(file_size, w_bytes);

        block_bg=0;
        text_i.read((char*) buffer.stream, std::streamsize(std::min<size_t>(file_size, block_bytes)));
    }

    i_file_stream& operator=(i_file_stream<sym_t> &&other) noexcept{
        move(std::forward<i_file_stream<sym_t>>(other));
        return *this;
    }

    void move(i_file_stream<sym_t> &&other){
        std::swap(file, other.file);
        std::swap(text_i, other.text_i);
        std::swap(tot_cells, other.tot_cells);
        std::swap(block_bg, other.block_bg);
        buffer.swap(other.buffer);
    }

    inline size_t size() const{
        return tot_cells;
    }

    inline const std::string& filename() const{
        return file;
    }

    void close(bool rem=false){
        if(buffer.stream!= nullptr){
            delete [] buffer.stream;
            buffer.stream = nullptr;
        }
        if(text_i.is_open()) text_i.close();
        if(rem){
            if(remove(file.c_str())){
                std::cout<<"Error trying to remove file"<<file<<std::endl;
            }
        }
    }

    ~i_file_stream(){
        close();
    }

    inline size_t read(size_t i) {
        assert(i<tot_cells);
        if(i<block_bg || (block_bg+buffer.stream_size)<=i){

            block_bg = (i/buffer.stream_size)*buffer.stream_size;
            text_i.seekg(block_bg*w_bytes);
            text_i.read((char*) buffer.stream, buffer.stream_size*w_bytes);
            if((unsigned int)text_i.gcount()<(buffer.stream_size*w_bytes)){
                text_i.clear();
            }
        }
        return buffer.stream[i-block_bg];
    }

    inline void read_bit_chunk(uint8_t *dst, size_t i, size_t j) {
        assert(i<=j);

        size_t cell_i = i >> w_shift;
        size_t cell_j = j >> w_shift;

        assert(cell_j<tot_cells && (cell_j-cell_i+1)<=buffer.stream_size);

        if(cell_i<block_bg || (block_bg+buffer.stream_size)<=cell_i){
            block_bg = (cell_i/buffer.stream_size)*buffer.stream_size;
            text_i.seekg(block_bg*w_bytes);
            text_i.read((char*) buffer.stream, buffer.stream_size*w_bytes);
            if((unsigned int)text_i.gcount()<(buffer.stream_size*w_bytes)){
                text_i.clear();
            }
        }

        size_t start = i % block_bits;
        if((cell_i/buffer.stream_size) == (cell_j/buffer.stream_size)){ //same block
            return buffer.read_chunk(dst, start, start+(j-i));
        }else{ //different blocks

            buffer.read_chunk(dst, start, block_bits-1);
            size_t bits_read = block_bits-start;
            block_bg+=buffer.stream_size;

            text_i.seekg(block_bg*w_bytes);
            text_i.read((char*) buffer.stream, buffer.stream_size*w_bytes);
            if((unsigned int)text_i.gcount()<(buffer.stream_size*w_bytes)){
                text_i.clear();
            }

            size_t bytes_read = INT_CEIL(bits_read, 8);
            size_t offset = std::min((8*bytes_read), (j-i+1))-bits_read;
            dst+=bytes_read-1;
            if(offset>0){
                uint8_t tail = buffer.read(0, offset-1);
                uint8_t left = bits_read % 8;
                *dst &= ~(bitstream<uint8_t>::masks[offset]<<left);//clean the area
                *dst |= tail << left;//add new data in the area
                bits_read+=offset;
            }

            if(bits_read<(j-i+1)){
                dst++;
                size_t end = j % block_bits;
                buffer.read_chunk(dst, offset, end);
                bits_read+=end-offset+1;
            }
            assert(bits_read==(j-i+1));
        }
    }

    inline size_t read_bits(size_t i, size_t j) {

        assert(i<=j && (j-i+1)<=64);

        size_t cell_i = i >> w_shift;
        size_t cell_j = j >> w_shift;

        assert(cell_j<tot_cells && (cell_j-cell_i+1)<=buffer.stream_size);

        if(cell_i<block_bg || (block_bg+buffer.stream_size)<=cell_i){
            block_bg = (cell_i/buffer.stream_size)*buffer.stream_size;
            text_i.seekg(block_bg*w_bytes);
            text_i.read((char*) buffer.stream, buffer.stream_size*w_bytes);
            if((unsigned int)text_i.gcount()<(buffer.stream_size*w_bytes)){
                text_i.clear();
            }
        }

        size_t start = i % block_bits;
        if((cell_i/buffer.stream_size) == (cell_j/buffer.stream_size)){ //same block
            return buffer.read(start, start+(j-i));
        }else{ //different blocks
            size_t data = buffer.read(start, block_bits-1);
            size_t delta = block_bits-start;
            block_bg+=buffer.stream_size;

            text_i.seekg(block_bg*w_bytes);
            text_i.read((char*) buffer.stream, buffer.stream_size*w_bytes);
            if((unsigned int)text_i.gcount()<(buffer.stream_size*w_bytes)){
                text_i.clear();
            }
            return (buffer.read(0, j % block_bits)<<delta) | data;
        }
    }
};

template<class sym_t>
struct o_file_stream{

    typedef long long int size_type;
    typedef sym_t         sym_type;
    std::ofstream ofs;
    std::ifstream ifs;
    std::string file;

    size_type w_bytes=0;
    size_type block_bg=0;
    size_type rm_pos=-1;
    size_type lm_pos=std::numeric_limits<size_type>::max();
    size_type last_pos=-1;
    bitstream<sym_t> buffer;
    bool modified=false;

    o_file_stream()=default;

    o_file_stream(o_file_stream<sym_t> &&other) noexcept{
        move(std::forward<o_file_stream<sym_t>>(other));
    }

    o_file_stream(const std::string& o_file, size_t buff_size_, std::ios::openmode mode){
        file = o_file;

        ofs.open(file, mode | std::ios::out | std::ios::binary);
        assert(ofs.good());

        ifs.open(file,  std::ios::in | std::ios::binary);
        assert(ifs.good());

        w_bytes = sizeof(sym_t);

        size_t buff_size = buff_size_/w_bytes;
        buffer.stream_size = buff_size;
        buffer.stream = new sym_t[buff_size];

        ifs.seekg (0, std::ifstream::end);
        last_pos = (ifs.tellg()/w_bytes)-1;
        ifs.seekg (0, std::ifstream::beg);

        block_bg=0;
        read_block();
    }

    o_file_stream<sym_t>& operator=(o_file_stream<sym_t> &&other) noexcept{
        if(this!=&other){
            move(std::forward<o_file_stream<sym_t>>(other));
        }
        return *this;
    }

    void move(o_file_stream<sym_t>&& other){
        std::swap(ifs, other.ifs);
        std::swap(ofs, other.ofs);
        std::swap(file, other.file);
        std::swap(w_bytes, other.w_bytes);
        std::swap(block_bg, other.block_bg);
        std::swap(rm_pos, other.rm_pos);
        std::swap(lm_pos, other.lm_pos);
        std::swap(last_pos, other.last_pos);
        std::swap(modified, other.modified);
        buffer.swap(other.buffer);
    }

    void push_back(sym_t val){

        if((block_bg + (size_type)buffer.stream_size) <= (++last_pos)){
            //write data
            write_block();

            //reset buffer
            //memset(buffer.stream, 0, buffer.stream_size*sizeof(sym_t));
            block_bg = (last_pos/buffer.stream_size)*buffer.stream_size;
            rm_pos=0;
            lm_pos=std::numeric_limits<size_type>::max();
        }

        size_type buff_pos = last_pos-block_bg;
        buffer.stream[buff_pos] = val;
        rm_pos = buff_pos;
        if(buff_pos<lm_pos) lm_pos = buff_pos;
        modified=true;
    }

    inline size_t size() const{
        return (last_pos+1);
    }

    const std::string& filename() const{
        return file;
    }

    inline bool copy_to_front(size_type src, size_type len, size_type freq=1){

        assert(src>=0);
        size_type b_size = buffer.stream_size;
        size_type buff_src = src & (b_size-1);
        size_type buff_dst = size() & (b_size-1);

        size_type block_src = src/b_size;
        size_type block_curr = last_pos / b_size;
        size_type block_dst = (last_pos+1) / b_size;

        //corner cases for which we cannot copy
        if((len*freq)>(size_type)b_size) return false; // the number of elements to be copied are more than the buffer size
        if(buff_src<=buff_dst && block_src<block_dst) return false; // the source is not in the buffer anymore
        if(src+len>(size_type)size()) return false; // we are copying a segment that doesn't exist
        if(block_dst-block_src>1) return false; // the source is not in the buffer anymore

        size_type n_syms, rem_syms=len*freq, glob_pos=last_pos+1, left;

        assert(block_src==block_dst || block_src==(block_dst-1));

        if(block_dst != block_curr){
            write_block();
            block_bg = (glob_pos/buffer.stream_size)*buffer.stream_size;
            lm_pos=0;
        }

        bool first_in_front=false;
        size_type i=1, tmp_src = buff_dst;
        while(rem_syms>0){

            n_syms=std::min<size_type>(len*i, rem_syms);

            if((buff_src+n_syms-1)>=b_size){//the source is in the buffer boundary

                assert(buff_dst<buff_src);
                left = b_size-buff_src;
                memcpy(&buffer.stream[buff_dst], &buffer.stream[buff_src], left);

                if((buff_dst+n_syms-1)>=b_size){//the destination is also in the buffer boundary
                    size_type d_left = b_size-(buff_dst+left);
                    memcpy(&buffer.stream[buff_dst+left], &buffer.stream[0], d_left);
                    modified = true;
                    rm_pos = b_size-1;
                    write_block();
                    block_bg = ((glob_pos+left+d_left)/buffer.stream_size)*buffer.stream_size;
                    lm_pos = 0;
                    memcpy(&buffer.stream[0], &buffer.stream[d_left], n_syms - (left+d_left));
                    buff_dst = n_syms - (left+d_left);
                }else{
                    //the source is in the boundary, but the destination is not
                    memcpy(&buffer.stream[buff_dst+left], &buffer.stream[0], n_syms-left);
                    buff_dst +=n_syms;
                }
            } else if((buff_dst+n_syms-1)>=b_size){// the source is not in the boundary, but the destination is

                //complete the buffer and write a copy on disk
                left = b_size-buff_dst;
                memcpy(&buffer.stream[buff_dst], &buffer.stream[buff_src], left);
                modified=true;
                rm_pos = b_size-1;
                write_block();
                block_bg = ((glob_pos+left)/buffer.stream_size)*buffer.stream_size;
                lm_pos = 0;

                //write the rest of the string at the beginning of the buffer
                memcpy(&buffer.stream[0], &buffer.stream[buff_src+left], n_syms-left);
                buff_dst = n_syms-left;
            }else{
                //nor the source or the destination are in the buffer boundary
                memcpy(&buffer.stream[buff_dst], &buffer.stream[buff_src], n_syms);
                buff_dst +=n_syms;
            }

            if(!first_in_front){
                buff_src = tmp_src;
                first_in_front = true;
                i=1;
            }else{
                i*=2;
            }
            glob_pos+=n_syms;
            rem_syms-=n_syms;
        }

        //assert(block_bg==(((glob_pos-1)/buffer.stream_size)*buffer.stream_size));
        assert(lm_pos<=rm_pos);
        last_pos = glob_pos-1;
        rm_pos = last_pos & (b_size-1);
        modified=true;
        return true;
    }

    inline size_type buff_size() const {
        return buffer.stream_size;
    }

    size_t read(size_t i) {
        assert((size_type)i<=last_pos);
        if((size_type)i<block_bg || (block_bg+(size_type)buffer.stream_size)<=(size_type)i){

            if(modified){
                //write data
                write_block();
            }

            block_bg = (i/buffer.stream_size)*buffer.stream_size;
            read_block();
        }
        return buffer.stream[i-block_bg];
    }

    void write_block(){
        assert(lm_pos<=rm_pos);
        ofs.seekp((block_bg+lm_pos)*w_bytes);
        ofs.write((char*) &buffer.stream[lm_pos], (rm_pos-lm_pos+1)*w_bytes);
        //TODO caveat: it might produce desynchronization with the ifs in some platforms (OSx)
        // if I don't flush the ofs buffer. However, flusing the buffer is expensive.
        // Maybe I should define an if macro here.
        assert(ofs.good());
        modified=false;
    }

    void read_block(){
        assert(!modified);
        ifs.seekg(block_bg*w_bytes);
        ifs.read((char*) buffer.stream, buffer.stream_size*w_bytes);
        if(ifs.gcount()<(size_type)buffer.stream_size*w_bytes){
            ifs.clear();
        }
        assert(ifs.good());
        rm_pos=-1;
        lm_pos=std::numeric_limits<size_type>::max();
    }

    void flush(){
        if(modified){
            write_block();
            ofs.flush();
            //init the buffer
            block_bg = 0;
            read_block();
        }else{
            ofs.flush();
        }
    }

    void write(size_type i, sym_t val) {

        if(i<block_bg || (block_bg+(size_type)buffer.stream_size)<=i){
            //write data
            if(modified){
                write_block();
                modified = false;
            }

            //read the new segment
            block_bg = (i/buffer.stream_size)*buffer.stream_size;
            read_block();
        }
        size_type buff_pos = i-block_bg;

        if(i >= last_pos) last_pos = i;
        if(buff_pos>rm_pos) rm_pos = buff_pos;
        if(buff_pos<lm_pos) lm_pos = buff_pos;

        buffer.stream[buff_pos] = val;
        modified = true;
    }

    void close(bool rem=false){
        flush();
        if(buffer.stream!=nullptr){
            delete [] buffer.stream;
            buffer.stream=nullptr;
        }
        if(ofs.is_open()){
            ofs.close();
        }
        if(ifs.is_open()){
            ifs.close();
        }
        if(rem){
            if(remove(file.c_str())){
                std::cout<<"Error trying to remove file"<<file<<std::endl;
            }
        }
    }

    ~o_file_stream(){
        close();
    }
};
#endif //LMS_COMPRESSOR_FILE_BUFFER_H
