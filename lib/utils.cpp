//
// Created by Diaz, Diego on 27.10.2021.
//
#include "cds/utils.h"
#include "cds/macros.h"
#include "cds/memory_handler.hpp"
#include <filesystem>
#include <unistd.h>
#include <random>
#include <sys/resource.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef __linux__
void empty_page_cache(const char *filename) {
    const int fd = open(filename, O_RDWR);
    if (fd == -1) {
        std::perror(filename);
        std::exit(EXIT_FAILURE);
    }
    const off_t length = lseek(fd, 0, SEEK_END);
    lseek(fd, 0L, SEEK_SET);
    posix_fadvise(fd, 0, length, POSIX_FADV_DONTNEED);
    close(fd);
}
#endif


bool file_exists(const std::filesystem::path& p, std::filesystem::file_status const& s){
    if(std::filesystem::status_known(s) ? std::filesystem::exists(s) : std::filesystem::exists(p)){
        return true;
    }else{
        return false;
    }
}

std::string random_string(size_t length){
    const std::string characters = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_int_distribution<size_t> distribution(0, characters.size() - 1);
    std::string random_string;
    random_string.resize(length);
    for(std::size_t i = 0; i < length; i++) {
        random_string[i] = characters[distribution(generator)];
    }
    return random_string;
}

std::vector<std::string> split (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }
    return result;
}

std::string report_space(off_t bytes){
    if(bytes<1000){
        return std::to_string(bytes)+" bytes";
    }else if(bytes<1000000){
        float b = float(bytes)/1000;
        return to_string_with_precision(b, 2)+" KBs";
    } else if(bytes < 1000000000L){
        float b = float(bytes)/1000000;
        return to_string_with_precision(b, 2)+" MBs";
    } else if(bytes < 1000000000000L){
        double b = double(bytes)/1000000000L;
        return to_string_with_precision(b, 2)+" GBs";
    } else {
        double b = double(bytes)/1000000000000L;
        return to_string_with_precision(b, 2)+" TBs";
    }
}

void report_mem_peak(){
    struct rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    std::cout<<"\n peak : "<<report_space(usage.ru_maxrss)<<std::endl;
}

std::ifstream::pos_type file_size(std::string & filename) {
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

str_collection collection_stats(std::string& input_file){

    uint8_t sep_sym;
    {
        std::ifstream ifs(input_file, std::ios::binary);
        ifs.seekg(-1, std::ios::end);
        ifs.read((char *) &sep_sym, 1);
        ifs.seekg(0, std::ios::beg);
    }

    struct stat st{};
    stat(input_file.c_str(), &st);
    off_t fsz = st.st_size;
    int fd_in = open(input_file.c_str(), O_RDONLY);

    off_t buffer_size = 1024*1024*8;
    off_t last_block = INT_CEIL(fsz, buffer_size);
    uint8_t *buffer= nullptr;

#ifdef __linux__
    long page_size = sysconf(_SC_PAGESIZE);
    posix_memalign((void **)&buffer, page_size, buffer_size);
#else
    //buffer = (uint8_t*)malloc(buffer_size);
    buffer = mem<uint8_t>::allocate(buffer_size);
#endif

#ifdef __linux__
    posix_fadvise(fd_in, 0, fsz, POSIX_FADV_SEQUENTIAL);
#endif

    str_collection str_coll;
    size_t str_len;
    uint8_t sym{};
    size_t sym_frq[256] = {0};
    size_t pos = 0, cont=0;

    for (off_t index = 0; index<last_block;index++) {

        ssize_t rd_syms = read(fd_in, buffer, buffer_size);
        assert(rd_syms>0);

        for(off_t i=0;i<rd_syms;i++){
            sym = buffer[i];
            sym_frq[sym]++;
            cont++;
            if(sym==sep_sym){
                str_coll.str_ptrs.push_back((long)pos);
                str_len = cont - pos;
                if(str_len>str_coll.longest_string) str_coll.longest_string=str_len;
                pos = cont;
            }
        }

#ifdef __linux__
        if(index && rd_syms==buffer_size){
            posix_fadvise(fd_in, (index-1)*buffer_size, buffer_size, POSIX_FADV_DONTNEED);
        }
#endif
    }

    //free(buffer);
    mem<uint8_t>::deallocate(buffer);
    close(fd_in);

    /*
    //We assume that the smallest symbol in the collection is the separator symbol.
    // This should be the last character in the file
    str_collection str_coll;

    size_t sym_frq[256] = {0};


    uint8_t sep_sym;
    std::ifstream ifs(input_file, std::ios::binary);
    ifs.seekg(-1, std::ios::end);
    ifs.read((char *)&sep_sym, 1);
    ifs.seekg(0, std::ios::beg);
    size_t pos = 0, cont=0;

    char buffer[BUFSIZ]={0};
    size_t read_bytes, str_len;
    std::streampos buff_size=8192;
    uint8_t sym=255;

    while(true){
        ifs.read((char *)&buffer, buff_size);
        read_bytes = ifs.gcount();
        if(read_bytes>0){
            for(size_t i=0;i<read_bytes;i++){
                sym = (uint8_t)buffer[i];
                sym_frq[sym]++;
                cont++;
                if(sym==sep_sym){
                    str_coll.str_ptrs.push_back((long)pos);
                    str_len = cont - pos;
                    if(str_len>str_coll.longest_string) str_coll.longest_string=str_len;
                    pos = cont;
                }
            }
        }else{
            break;
        }
    }*/

    str_coll.str_ptrs.shrink_to_fit();

    for(size_t i=0;i<256;i++){
        if(sym_frq[i]!=0){
            str_coll.alphabet.push_back(i);
            str_coll.sym_freqs.push_back(sym_frq[i]);
            str_coll.n_char+=sym_frq[i];
        }
    }

    str_coll.min_sym = str_coll.alphabet[0];
    str_coll.max_sym = str_coll.alphabet.back();
    str_coll.n_strings = sym_frq[str_coll.alphabet[0]];
    //ifs.close();

#ifdef __linux__
    empty_page_cache(input_file.c_str());
#endif

    if(sym!=str_coll.min_sym){
        std::cerr<<"Error: the file does not end with the separator symbol"<<std::endl;
        exit(1);
    }
    return str_coll;
}
