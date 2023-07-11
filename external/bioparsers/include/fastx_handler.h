//
// Created by Diaz, Diego on 28.10.2021.
//

#ifndef RYU_FASTQ_HANDLER_H
#define RYU_FASTQ_HANDLER_H

#include <fstream>
#include <utility>
#include <zlib.h>
#include <thread>
#include "dna_string.h"
#include "file_streams.hpp"
#include "utils.h"
extern"C"{
#include "kseq.h"
}
#include <unistd.h>

KSEQ_INIT(gzFile, gzread)

#define insert_run(head, len)\
if(len>1){\
    len_buffer[u++] = len;\
    if(u==l_buff_size){ \
        len_ofs.write((char *)len_buffer, l_buff_size*sizeof(uint16_t)); \
        u=0; \
    }\
    head += 84;\
    n_hp++;\
    n_s_hp++;\
}\
n_syms+=len;\
head_buffer[k++] = head;\
sym_freqs[head]++;\
if(len>max_hp) max_hp=len;\
if(k==buff_size){\
    head_ofs.write(head_buffer, buff_size);\
    k=0;\
}\

//string collection with homopolymer compressed
struct hp_collection{
    str_collection str_coll;//string collection struct as in grlBWT
    std::string hp_heads;//file with reads with the homopolymer compressed
    std::string hp_lengths; //file storing the homopolymer lengths
    size_t max_hp_len{}; //length of the longest homopolymer
    size_t max_n_s_hp{}; //maximum number of homopolymers in a read
};

class fasta_writer{
    std::ofstream ofs;
    char buffer[1024]={0};
    size_t b_pos=0;
    size_t wrap_width=0;

public:
    void insert_record(const std::string& identifier, const dna_string& string);
    explicit fasta_writer(const std::string& output_file,
                          size_t wrap=0) : ofs(output_file, std::ios::out | std::ios::binary),
                                           wrap_width(wrap){};
    ~fasta_writer(){
        if(b_pos!=0){
            ofs.write((char *)&buffer, (std::streamsize)b_pos);
        }
        ofs.close();
    }
};

template<bool is_gzipped>
struct fasta_parser {
    std::string input_file;
    std::string r_heads_file;
    std::string r_len_file;
    gzFile in_zfile=nullptr;
    std::ifstream ifs;
    std::ofstream head_ofs;
    std::ofstream len_ofs;
public:
    fasta_parser(std::string  ifile, std::string const& ofile): input_file(std::move(ifile)),
                                                                r_heads_file(ofile+".heads"),
                                                                r_len_file(ofile+".len"){};
    hp_collection operator()();
};

template<bool is_gzipped>
class fastq_parser {
    std::string input_file;
    std::string r_heads_file;
    std::string r_len_file;
    gzFile in_zfile=nullptr;
    std::ifstream ifs;
    std::ofstream head_ofs;
    std::ofstream len_ofs;
public:

    fastq_parser(std::string  ifile, std::string const& ofile): input_file(std::move(ifile)),
                                                                r_heads_file(ofile+".heads"),
                                                                r_len_file(ofile+".len"){};
    hp_collection operator()();
};

//transform an input FASTA/Q file in one-string-per-line format (without homopolymer compression)
str_collection fastx2plain_format(const std::string& input_file, std::string const& output_file, bool rc, uint8_t sep);

//transform an input FASTA/Q file in one-string-per-line format (assuming homopolymer compression)
hp_collection fastx2_hp_plain_format(const std::string& input_file, std::string const& output_file);

//get the DNA reverse complements (assuming homopolymer compression)
void get_hp_DNA_rc(std::string const& r_heads_file, std::string const& r_len_file, str_collection& str_coll);



template<bool is_gzipped>
hp_collection fasta_parser<is_gzipped>::operator()() {

    if constexpr(is_gzipped){
        in_zfile = gzopen(input_file.c_str(), "rb");
    }else{
        ifs.open(input_file.c_str(), std::ios::binary | std::ios::in);
        assert(ifs.good());
    }

    std::streamsize k=0, u=0;
    std::streamsize buff_size = 8192, l_buff_size=4096;
    std::streamsize read_bytes;
    char *head_buffer, *in_buffer;
    uint16_t *len_buffer;
    size_t sym_freqs[255]={0};
    hp_collection hp_coll;

    in_buffer = (char *)calloc(buff_size, 1);
    len_buffer = (uint16_t *)calloc(l_buff_size, sizeof(uint16_t));
    head_buffer = (char *)calloc(buff_size, 1);

    head_ofs.open(r_heads_file, std::ios::out | std::ios::binary);
    len_ofs.open(r_len_file, std::ios::out | std::ios::binary);

    bool header=false;
    uint8_t sym, p_sym;

    if constexpr(is_gzipped) {
        read_bytes = (std::streamsize) gzread(in_zfile, in_buffer, buff_size);
    }else{
        ifs.read((char *)in_buffer, buff_size);
        read_bytes = ifs.gcount();
    }

    std::streamsize j=1;
    while(j<read_bytes && in_buffer[j-1]!='\n') j++;
    assert(j+1<read_bytes);

    p_sym = in_buffer[j++];
    size_t r_len=1, n_hp=0, n_syms=0, max_hp=0, n_s_hp=0, max_n_s_hp=0;
    hp_coll.str_coll.n_strings = 1;
    uint8_t r_head;

    while (read_bytes>0) {
        while(j<read_bytes){
            if(in_buffer[j]=='>'){
                insert_run(p_sym, r_len)
                header = true;
                r_len = 1;
                p_sym = '\n';
                hp_coll.str_coll.n_strings++;
                if(n_s_hp>max_n_s_hp) max_n_s_hp=n_s_hp;
                n_s_hp=0;
            }else if (in_buffer[j]=='\n'){
                header = false;
            }else if(!header){
                sym = (uint8_t) toupper(in_buffer[j]);
                if(sym!=p_sym){
                    insert_run(p_sym, r_len)
                    r_len = 0;
                }
                r_len++;
                p_sym = sym;
            }
            j++;
        }

        if constexpr(is_gzipped) {
            read_bytes = (std::streamsize) gzread(in_zfile, in_buffer, buff_size);
        }else{
            ifs.read((char *)in_buffer, buff_size);
            read_bytes = ifs.gcount();
        }
        j=0;
    }

    r_head = p_sym;
    insert_run(r_head, r_len);
    if(n_s_hp>max_n_s_hp) max_n_s_hp=n_s_hp;

    head_buffer[k++] = '\n';
    sym_freqs[10]++;//10 = '\n'

    head_ofs.write(head_buffer, k);
    len_ofs.write((char *)len_buffer, std::streamsize(sizeof(uint16_t)*u));

    if constexpr(is_gzipped){
        gzclose(in_zfile);
    }else{
        ifs.close();
    }
    head_ofs.close();
    len_ofs.close();
    free(in_buffer);
    free(head_buffer);
    free(len_buffer);

    /*for(uint8_t i=0;i<255;i++){
        if(sym_freqs[i]!=0){
            hp_coll.str_coll.alphabet.push_back(i);
            hp_coll.str_coll.sym_freqs.push_back(sym_freqs[i]);
            hp_coll.str_coll.n_syms+=sym_freqs[i];
        }
    }

    hp_coll.str_coll.min_sym = hp_coll.str_coll.alphabet[0];
    hp_coll.str_coll.max_sym = hp_coll.str_coll.alphabet.back();*/
    hp_coll.hp_lengths = r_len_file;
    hp_coll.hp_heads = r_heads_file;
    hp_coll.max_hp_len = max_hp;
    hp_coll.max_n_s_hp = max_n_s_hp;

    get_hp_DNA_rc(r_heads_file, r_len_file, hp_coll.str_coll);

    return hp_coll;
}

template<bool is_gzipped>
hp_collection fastq_parser<is_gzipped>::operator()() {

    if constexpr(is_gzipped){
        in_zfile = gzopen(input_file.c_str(), "rb");
    }else{
        ifs.open(input_file.c_str(), std::ios::binary | std::ios::in);
        assert(ifs.good());
    }

    std::streamsize read_bytes;
    std::streamsize k=0, u=0;
    std::streamsize buff_size = 8192, l_buff_size=4096;
    char *in_buffer, *head_buffer;
    uint16_t* len_buffer;
    in_buffer = (char *)calloc(buff_size, 1);
    head_buffer = (char *)calloc(buff_size, 1);
    len_buffer = (uint16_t *)calloc(l_buff_size, sizeof(uint16_t));
    head_ofs.open(r_heads_file, std::ios::out | std::ios::binary);
    len_ofs.open(r_len_file, std::ios::out | std::ios::binary);
    size_t n_hp=0, n_syms=0, n_s_hp=0, max_n_s_hp=0, max_hp=0;

    hp_collection hp_coll;
    size_t sym_freqs[255]={0};

    if constexpr(is_gzipped) {
        read_bytes = (std::streamsize) gzread(in_zfile, in_buffer, buff_size);
    }else{
        ifs.read((char *)in_buffer, buff_size);
        read_bytes = ifs.gcount();
    }

    std::streamsize i=1;
    while(i<read_bytes && in_buffer[i-1]!='\n') i++;
    assert(i+1<read_bytes);
    size_t r_len=1;
    uint8_t p_sym=in_buffer[i++];
    hp_coll.str_coll.n_strings=1;
    size_t line=1;

    while (read_bytes>0) {
        while(i<read_bytes){
            if(line==1){
                if(p_sym!=in_buffer[i]){
                    //TODO testing
                    //if(r_len>1 && hp_coll.str_coll.n_strings==8000){
                    //    std::cout<<n_hp<<" "<<(char)p_sym<<" "<<r_len<<std::endl;
                    //}
                    //
                    insert_run(p_sym, r_len);
                    r_len = 0;
                    p_sym = in_buffer[i];
                }
                r_len++;
            }

            if(in_buffer[i]=='\n'){
                line = (line+1) & 3UL;
                if(line==1) hp_coll.str_coll.n_strings++;
                if(n_s_hp>max_n_s_hp) max_n_s_hp=n_s_hp;
                n_s_hp=0;
            }
            i++;
        }

        if constexpr(is_gzipped) {
            read_bytes = (std::streamsize) gzread(in_zfile, in_buffer, buff_size);
        }else{
            ifs.read((char *)in_buffer, buff_size);
            read_bytes = ifs.gcount();
        }
        i=0;
    }

    insert_run(p_sym, r_len)

    head_ofs.write( head_buffer, k);
    len_ofs.write((char *)len_buffer, std::streamsize(sizeof(uint16_t)*u));
    head_ofs.close();
    len_ofs.close();

    free(in_buffer);
    free(head_buffer);
    free(len_buffer);

    if constexpr(is_gzipped){
        gzclose(in_zfile);
    }else{
        ifs.close();
    }

    /*for(uint8_t j=0;j<255;j++){
        if(sym_freqs[j]!=0){
            hp_coll.str_coll.alphabet.push_back(j);
            hp_coll.str_coll.sym_freqs.push_back(sym_freqs[j]);
            hp_coll.str_coll.n_syms+=sym_freqs[j];
        }
    }

    hp_coll.str_coll.min_sym = hp_coll.str_coll.alphabet[0];
    hp_coll.str_coll.max_sym = hp_coll.str_coll.alphabet.back();*/
    hp_coll.hp_lengths = r_len_file;
    hp_coll.hp_heads = r_heads_file;
    hp_coll.max_hp_len = max_hp;
    hp_coll.max_n_s_hp = max_n_s_hp;

    std::cout<<"  Input stats:"<<std::endl;
    std::cout<<"    Number of DNA symbols:                            "<<n_syms-hp_coll.str_coll.n_strings<<std::endl;
    std::cout<<"    Number of homopolymers:                           "<<n_hp<<std::endl;
    std::cout<<"    Longest homopolymer:                              "<<hp_coll.max_hp_len<<std::endl;
    std::cout<<"    Maximum number of homopolymers in a read:         "<<hp_coll.max_n_s_hp<<std::endl;
    std::cout<<"    Number of symbols in the input after compression: "<<hp_coll.str_coll.n_syms<<std::endl;
    std::cout<<"    Size reduction:                                   "<<100-(double(hp_coll.str_coll.n_syms)/double(n_syms))*100<<"%"<<std::endl;
    std::cout<<"    Number of strings:                                "<<hp_coll.str_coll.n_strings<<std::endl;

    std::cout<<"  Inserting the reverse complements symbols"<<std::endl;
    get_hp_DNA_rc(r_heads_file, r_len_file, hp_coll.str_coll);

    /*std::cout<<"  Symbol frequencies in the new file: "<<std::endl;
    uint8_t sym;
    for(size_t j=0;j<hp_coll.str_coll.alphabet.size();j++){
        sym = hp_coll.str_coll.alphabet[j];
        if(sym=='\n'){
            std::cout<<"    \\n";
        }else{
            if(sym >84){
                std::cout<<"    *"<<(char)(sym-84);
            }else{
                std::cout<<"    "<<(char)(sym);
            }
        }
        std::cout<<"\t"<<hp_coll.str_coll.sym_freqs[j]<<std::endl;
    }*/
    std::cout<<"  Total symbols: "<<hp_coll.str_coll.n_syms<<std::endl;
    std::cout<<"  Total strings: "<<hp_coll.str_coll.n_strings<<std::endl;
    std::cout<<"  Results were stored in:"<<std::endl;
    std::cout<<"    "<<r_heads_file<<std::endl;
    std::cout<<"    "<<r_len_file<<std::endl;

    return hp_coll;
}

#endif //RYU_FASTQ_HANDLER_H
