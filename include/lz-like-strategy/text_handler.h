//
// Created by Diaz, Diego on 17.9.2023.
//

#ifndef LCG_TEXT_HANDLER_H
#define LCG_TEXT_HANDLER_H

#include "hashing.h"

namespace lzstrat {

   struct text_chunk{

       typedef uint32_t size_type;

       size_t id{}; //chunk id
       off_t n_bytes_before{};//number of bytes in the file before this chunk
       size_type sep_sym{};//symbol in the buffer delimiting consecutive strings
       std::vector<hashing> &hashing_func;

       size_type *buffer = nullptr;// chunk's buffer
       off_t buffer_bytes{};

       uint8_t *data_start = nullptr;
       off_t data_bytes{}; //number of bytes the buffer can hold

       text_chunk(std::vector<hashing>& _hf): hashing_func(_hf){}

       off_t eff_bytes(){
           return data_bytes;
       }

       ~text_chunk(){
           if(buffer!=nullptr){
               free(buffer);
           }
       }
   };

   void read_chunk_from_file(int fd, off_t& rem_text_bytes, off_t& read_text_bytes, text_chunk& chunk){

       off_t chunk_bytes = chunk.data_bytes<rem_text_bytes ? chunk.data_bytes : rem_text_bytes;
       chunk.data_bytes = chunk_bytes;

       off_t acc_bytes = 0;
       off_t read_bytes;
       off_t fd_buff_bytes = 8388608;// 8MB buffer

       chunk.data_start = (uint8_t *) chunk.buffer;
       chunk.data_start = chunk.data_start + chunk.buffer_bytes-chunk.data_bytes;
       uint8_t * data = chunk.data_start;

       while(chunk_bytes>0) {
           fd_buff_bytes = fd_buff_bytes<chunk_bytes ? fd_buff_bytes : chunk_bytes;
           read_bytes = read(fd, data, fd_buff_bytes);
           assert(read_bytes>0);

           data+=read_bytes;
           chunk_bytes-=read_bytes;
           acc_bytes+=read_bytes;
       }
       assert(chunk.data_bytes==acc_bytes);
       chunk.n_bytes_before = read_text_bytes;

       //go to the rightmost separators symbol
       size_t i = chunk.data_bytes-1;
       while(i>0 && chunk.data_start[i]!=chunk.sep_sym) i--;

       off_t eff_bytes = i+1;
       chunk.data_bytes = eff_bytes;

       off_t offset = acc_bytes-eff_bytes;
       rem_text_bytes-= eff_bytes;
       read_text_bytes = lseek(fd, offset*-1, SEEK_CUR);
   }
};
#endif //LCG_TEXT_HANDLER_H
