//
// Created by Diaz, Diego on 17.9.2023.
//

#ifndef LCG_TEXT_HANDLER_H
#define LCG_TEXT_HANDLER_H

#include "partial_gram.h"

namespace lz_like_strat {

   struct text_chunk{

       typedef uint32_t size_type;

       size_t id{}; //chunk id
       off_t n_bytes_before{};//number of bytes in the file before this chunk
       size_type sep_sym{};//symbol in the buffer delimiting consecutive strings

       size_type *buffer = nullptr;// chunk's buffer
       off_t buffer_bytes{};
       off_t e_bytes{};

       uint8_t *text = nullptr;
       off_t text_bytes{}; //number of bytes the buffer can hold

       size_type * parse = nullptr;

       //todo make it compatible with uint16_t
       partial_gram<uint8_t> p_gram;

       [[nodiscard]] off_t eff_buff_bytes() const{
           return e_bytes;
       }

       ~text_chunk(){
           if(buffer!=nullptr){
               free(buffer);
           }
       }
   };

   template<class sym_type>
   void read_chunk_from_file(int fd, off_t& rem_text_bytes, off_t& read_text_bytes, text_chunk& chunk){

       off_t chunk_bytes = chunk.text_bytes<rem_text_bytes ? chunk.text_bytes : rem_text_bytes;
       chunk.text_bytes = chunk_bytes;

       off_t acc_bytes = 0;
       off_t read_bytes;
       off_t fd_buff_bytes = 8388608;// 8MB buffer

       chunk.text = (sym_type *) chunk.buffer;
       sym_type * data = chunk.text;
       chunk.n_bytes_before = read_text_bytes;
       //off_t limit = chunk.text_bytes>>1;
       off_t limit = 0;
       off_t i;

       while(true){

           while(chunk_bytes>0) {
               fd_buff_bytes = fd_buff_bytes<chunk_bytes ? fd_buff_bytes : chunk_bytes;
               read_bytes = read(fd, data, fd_buff_bytes);
               assert(read_bytes>0);

               data+=read_bytes;
               chunk_bytes-=read_bytes;
               acc_bytes+=read_bytes;
           }
           assert(chunk.text_bytes==acc_bytes);

           //go to the rightmost separators symbol
           i = chunk.text_bytes-1;
           while(i>limit && chunk.text[i]!=chunk.sep_sym){
               i--;
           }
           if(i>limit) break;

           off_t tmp_ck_size = INT_CEIL(((chunk.text_bytes*125)/100), sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type);
           tmp_ck_size = std::min(tmp_ck_size, rem_text_bytes);
           chunk_bytes = tmp_ck_size-chunk.text_bytes;
           chunk.text_bytes =  tmp_ck_size;

           //the parse size is (text_len/2)*(sizeof(size_type)/sizeof(sym_type)),
           // where ``text_len'' is the number of input symbols that fits the buffer
           off_t parse_bytes = INT_CEIL((tmp_ck_size/sizeof(sym_type)), 2)*(sizeof(text_chunk::size_type)/sizeof(sym_type));


           chunk.buffer_bytes = off_t(tmp_ck_size + parse_bytes);
           chunk.buffer = (text_chunk::size_type *) realloc(chunk.buffer, chunk.buffer_bytes);
           chunk.text = (sym_type *)chunk.buffer;
           data = &chunk.text[chunk.text_bytes-chunk_bytes];
       }

       off_t eff_bytes = i+1;
       chunk.text_bytes = eff_bytes;
       chunk.e_bytes = eff_bytes;

       off_t offset = acc_bytes-eff_bytes;
       rem_text_bytes-= eff_bytes;

       read_text_bytes = lseek(fd, offset*-1, SEEK_CUR);
   }
}
#endif //LCG_TEXT_HANDLER_H
