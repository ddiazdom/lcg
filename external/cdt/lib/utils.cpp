//
// Created by Diaz, Diego on 27.10.2021.
//
#include "utils.h"
#include <filesystem>
#ifdef __APPLE__
#include <unistd.h>
#endif
#include <random>
#include <zlib.h>
#include <sys/resource.h>

bool is_fastx(const std::string& input_file){

    bool is_gz=check_gzip(input_file);
    char f_sym;

    //read the first symbol in the file to check if it is in FASTA or FASTQ format
    if(is_gz){
        gzFile zfile = gzopen(input_file.c_str(), "rb");
        gzread(zfile, &f_sym, 1);
        gzclose(zfile);
    }else{
        std::ifstream ifs(input_file, std::ios::binary);
        ifs.read(&f_sym, 1);
        ifs.close();
    }
    return f_sym=='>' || f_sym=='@';
}

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

bool check_gzip(const std::string& file){
   std::ifstream ifs(file, std::ios::binary);
   char f_bytes[2]={0,0};
   auto ext = std::filesystem::path(file).extension();
   ifs.read(f_bytes, 2);
   ifs.close();
   return ext==".gz" && f_bytes[0]=='\x1F' && f_bytes[1]=='\x8B';
}

std::string append_to_path(std::string& orig_path, std::string const& suffix){
    return std::filesystem::path(orig_path) / suffix;
}

std::string remove_ext(std::string& input_string, std::vector<std::string>& exts){
    std::string res = input_string;
    for(auto const& ext: exts){
        if(ends_with(res, ext)){
            res.erase(res.size()-ext.size());
            break;
        }
    }
    return res;
}

std::string create_temp_folder(std::string& base_folder, std::string const & prefix){
    std::string tmp_path = base_folder + "/"+prefix+".XXXXXX";
    char temp[200] = {0};
    tmp_path.copy(temp, tmp_path.size() + 1);
    temp[tmp_path.size() + 1] = '\0';
    auto res = mkdtemp(temp);
    if (res == nullptr) {
        std::cout << "Error trying to create a temporal folder" << std::endl;
    }
    return {temp};
}

//a helper function to modify the output
bool ends_with(std::string const & value, std::string const & ending) {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

void report_mem_peak(){
    struct rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    std::cout<<"\n peak : "<<usage.ru_maxrss<<" megabytes "<<std::endl;
}

template<class sym_type>
str_collection collection_stats(std::string& input_file){

    //We assume that the smallest symbol in the collection is the separator symbol.
    // This should be the last character in the file
    str_collection str_coll;

    size_t sym_frq[256] = {0};

    std::ifstream ifs(input_file, std::ios::binary);

    sym_type sep_sym;
    ifs.seekg(-1*sizeof(sym_type), std::ios::end);

    ifs.read((char *)&sep_sym, 1*sizeof(sym_type));
    str_coll.n_syms = ifs.tellg()/sizeof(sym_type);
    assert(str_coll.n_syms>0);
    str_coll.max_sym_freq = str_coll.n_syms;
    ifs.seekg(0, std::ios::beg);
    size_t pos = 0, cont=0;

    std::streampos buff_size = 8388608;
    std::streampos n_elms = buff_size/sizeof(sym_type);
    auto buffer = (sym_type *) calloc(n_elms, sizeof(sym_type));

    size_t read_bytes, str_len, read_syms;
    sym_type sym=0, min_sym=std::numeric_limits<sym_type>::max(), max_sym=0;

    while(true){
        ifs.read((char *)buffer, buff_size);

        read_bytes = ifs.gcount();

        if(read_bytes>0){

            read_syms = read_bytes/sizeof(sym_type);
            for(size_t i=0;i<read_syms;i++){
                sym = buffer[i];

                if constexpr (std::is_same<sym_type, uint8_t>::value){
                    sym_frq[sym]++;
                }else{
                    if(sym>max_sym) max_sym = sym;
                    if(sym<min_sym) min_sym = sym;
                }

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
    }

    str_coll.str_ptrs.shrink_to_fit();

    if constexpr (std::is_same<sym_type, uint8_t>::value) {
        size_t i=0;
        while(sym_frq[i] == 0) i++;
        min_sym = i;
        str_coll.sep_symbol = sep_sym;

        i=255;
        while(sym_frq[i]==0) i--;
        max_sym = i;
        str_coll.alphabet.reserve(max_sym+1);

        size_t max_sym_freq = 0;
        sym_type tmp_sym =0;
        for(unsigned long frq : sym_frq){
            if(frq>0){
                if(frq>max_sym_freq) max_sym_freq = frq;
                str_coll.alphabet.push_back(tmp_sym);
            }
            tmp_sym++;
        }
        str_coll.max_sym_freq = max_sym_freq;
    }

    if(sep_sym!=min_sym || sym!=sep_sym){
        std::cout<<"Error: the file is ill formed"<<std::endl;
        exit(1);
    }

    str_coll.min_sym = min_sym;
    str_coll.max_sym = max_sym;
    str_coll.n_strings = str_coll.str_ptrs.size();

    ifs.close();
    free(buffer);
    return str_coll;
}


template str_collection collection_stats<uint8_t>(std::string& input_file);
template str_collection collection_stats<uint16_t>(std::string& input_file);
template str_collection collection_stats<uint32_t>(std::string& input_file);
template str_collection collection_stats<uint64_t>(std::string& input_file);


