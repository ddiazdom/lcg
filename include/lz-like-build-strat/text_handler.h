//
// Created by Diaz, Diego on 17.9.2023.
//

#ifndef LCG_TEXT_HANDLER_H
#define LCG_TEXT_HANDLER_H

#include "partial_gram.h"
#include "cds/mmap_allocator.h"

namespace lz_like_strat {

   struct text_chunk{

       typedef uint32_t size_type;

       size_t id{};//chunk id
       off_t n_bytes_before{};//number of bytes in the file before this chunk
       size_type sep_sym{};//symbol in the buffer delimiting consecutive strings

       off_t buffer_bytes{};
       off_t e_bytes{};

       std::vector<uint64_t *> fps;
       std::vector<size_t> fps_len;

       phrase_set<uint8_t> ter_dict;
       std::vector<phrase_set<uint32_t>> nt_dicts;

       uint8_t *text = nullptr;
       off_t text_bytes{}; //number of bytes the buffer can hold

       size_type * parse = nullptr;
       off_t parse_size=0;
       uint32_t round=0;

       //to measure time
       std::chrono::steady_clock::time_point t_start;
       std::chrono::steady_clock::time_point t_end;

       //area of the buffer free for satellite data
       off_t used_bytes=0;

       explicit text_chunk(uint8_t _sep_sym='\n'): sep_sym(_sep_sym){

           fps.resize(40);
           fps_len.resize(40);
           nt_dicts.resize(39);

           size_t alpha_size = std::numeric_limits<uint8_t>::max()+1;
           fps[0] = mem<uint64_t>::allocate(alpha_size);
           fps_len[0] = alpha_size;
           for(size_t i=0;i<alpha_size;i++){
               fps[0][i] = XXH3_64bits(&i, sizeof(uint8_t));
               assert(fps[0][i]!=0);
           }
           fps[0][sep_sym]=0;//small hack

           //small hack as the metasymbols are one-based
           for(size_t i=1;i<fps.size();i++){
               fps[i] = mem<uint64_t>::allocate(1);
               fps[i][0]=0;
               fps_len[i]=1;
           }
           //chunk.p_gram.max_tsym = std::numeric_limits<sym_type>::max();
           //chunk.p_gram.sep_tsym = chunk.sep_sym;
           //chunk.p_gram.text_size = chunk.text_bytes/sizeof(sym_type);
           //chunk.p_gram.txt_id = chunk.id;
       }

       [[nodiscard]] off_t eff_buff_bytes() const{
           return e_bytes;
       }

       [[nodiscard]] size_t av_bytes() const{
           return buffer_bytes-used_bytes;
       }

       void add_used_bytes(off_t n_bytes){
           used_bytes+=n_bytes;
           assert(used_bytes<=buffer_bytes);
       }

       [[nodiscard]] std::pair<uint8_t*, size_t> get_free_mem_area() const {
           if(used_bytes==buffer_bytes){
               return {nullptr, 0};
           }else{
               return {text+used_bytes, av_bytes()};
           }

           for(auto & fps_set : fps){
               mem<uint64_t>::deallocate(fps_set);
           }
       }

       void update_used_bytes(off_t new_used_bytes){
           //std::cout<<(uintptr_t) text<<" "<<(uintptr_t) next_av_byte<<" "<< (uintptr_t) boundary<<std::endl;
           //assert((uintptr_t) text<= (uintptr_t) new_pos && (uintptr_t) new_pos <= (uintptr_t) boundary);
           used_bytes = new_used_bytes;
       }

       void increase_capacity(off_t new_cap){
           if(new_cap>buffer_bytes){
               buffer_bytes = new_cap;
               if(text== nullptr){
                   text = mem<uint8_t>::allocate(buffer_bytes);
               }else{
                   text = mem<uint8_t>::reallocate(text, buffer_bytes);
               }
               //std::cout<<"now the buffer uses: "<<report_space(buffer_bytes)<<std::endl;
           }
       }

        uintptr_t dist(uint8_t* position) const {
           //std::cout<<(uintptr_t) text<<" "<< (uintptr_t) position <<" "<<(uintptr_t) (text+buffer_bytes)<<std::endl;
           assert((uintptr_t) text<= (uintptr_t) position && (uintptr_t) position <= (uintptr_t) (text+buffer_bytes));
           return (uintptr_t) position - (uintptr_t) text;
       }

       [[nodiscard]] uintptr_t boundary() const {
           return (uintptr_t) (text+buffer_bytes);
       }

       ~text_chunk(){
           if(text!=nullptr){
               mem<uint8_t>::deallocate(text);
           }
           std::cout<<" dict ter "<<ter_dict.size()<<", "<<ter_dict.load_factor()<<std::endl;
           ter_dict.psl_dist();
           size_t i=0;
           for(auto &nt_dict : nt_dicts){
               if(!nt_dict.empty()){
                   std::cout<<i<<" "<<nt_dict.size()<<", "<<nt_dict.load_factor()<<std::endl;
                   nt_dict.psl_dist();
               }
               i++;
           }
       }

       size_t mem_usage(){
           size_t bytes = ter_dict.mem_usage();
           for(auto const& dict: nt_dicts){
               bytes+=dict.mem_usage();
           }

           for(auto const& fps_set : fps_len){
               bytes+=fps_set*sizeof(uint64_t);
           }
           return bytes;
       }

       size_t eff_mem_usage(){
           size_t bytes = ter_dict.eff_mem_usage();
           for(auto const& dict: nt_dicts){
               bytes+=dict.eff_mem_usage();
           }

           for(auto const& fps_set : fps_len){
               bytes+=fps_set*sizeof(uint64_t);
           }
           return bytes;
       }

       size_t dict_buff_av(){
           size_t bytes = ter_dict.buff_bytes_available();
           for(auto const& dict: nt_dicts){
               bytes+=dict.buff_bytes_available();
           }
           return bytes;
       }

   };

   //template<class sym_type>
   void read_chunk_from_file(int fd, off_t& rem_text_bytes, off_t& read_text_bytes, text_chunk& chunk){

       off_t chunk_bytes = chunk.text_bytes<rem_text_bytes ? chunk.text_bytes : rem_text_bytes;
       chunk.text_bytes = chunk_bytes;

       off_t acc_bytes = 0;
       off_t read_bytes;
       off_t fd_buff_bytes = 8388608;// 8MB buffer

       //chunk.text = (sym_type *) chunk.buffer;
       uint8_t * data = chunk.text;
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

           //go to the rightmost separator symbol
           i = chunk.text_bytes-1;
           while(i>limit && chunk.text[i]!=chunk.sep_sym){
               i--;
           }
           if(i>limit) break;

           off_t tmp_ck_size = INT_CEIL(((chunk.text_bytes*125)/100), sizeof(text_chunk::size_type))*sizeof(text_chunk::size_type);
           tmp_ck_size = std::min(tmp_ck_size, rem_text_bytes);
           chunk_bytes = tmp_ck_size-chunk.text_bytes;
           chunk.text_bytes =  tmp_ck_size;

           chunk.increase_capacity(chunk.text_bytes);
           //chunk.buffer_bytes = tmp_ck_size;
           //chunk.text = alloc<uint8_t>::reallocate(chunk.text, chunk.buffer_bytes);

           data = &chunk.text[chunk.text_bytes-chunk_bytes];
       }

       off_t eff_bytes = i+1;
       chunk.text_bytes = eff_bytes;
       chunk.e_bytes = eff_bytes;

       off_t offset = acc_bytes-eff_bytes;
       rem_text_bytes-= eff_bytes;

       read_text_bytes = lseek(fd, offset*-1, SEEK_CUR);
       chunk.update_used_bytes(chunk.text_bytes);
   }
}
#endif //LCG_TEXT_HANDLER_H
