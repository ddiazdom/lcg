//
// Created by Diaz, Diego on 27.10.2021.
//

#ifndef LHTIGS_UTILS_H
#define LHTIGS_UTILS_H
#include <chrono>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <cassert>

#ifdef __APPLE__
#include <unistd.h>
#endif

struct str_collection {
    std::vector<long> str_ptrs;
    std::vector<uint8_t> alphabet;
    uint8_t sep_symbol;
    size_t n_strings=0;
    size_t longest_string=0;
    size_t max_sym=0;
    size_t min_sym=std::numeric_limits<size_t>::max();
    size_t n_syms=0;
    size_t max_sym_freq=0;
};

//check if the file is gzipped
bool check_gzip(const std::string& file);

//check if file exists
bool file_exists(const std::filesystem::path& p, std::filesystem::file_status const& s = std::filesystem::file_status{});

//create temporal folder
std::string create_temp_folder(std::string& base_folder, std::string const & prefix);

//a helper function to modify the output
bool ends_with(std::string const & value, std::string const & ending);

//remove extension
std::string remove_ext(std::string& input_string, std::vector<std::string>& exts);

//produce a random string
std::string random_string(size_t length);

//check if it is in FASTx
bool is_fastx(const std::string& file);


struct tmp_workspace{

    std::string tmp_folder;
    std::string ext;
    bool remove_all;

    explicit tmp_workspace(std::string const& base_folder=std::filesystem::temp_directory_path(),
                           bool rem_all=true,
                           std::string const& prefix="tmp") : remove_all(rem_all) {

        std::string tmp_path = std::filesystem::canonical(std::filesystem::path(base_folder)) / std::string(prefix+".XXXXXX");
        char temp[200] = {0};
        tmp_path.copy(temp, tmp_path.size() + 1);
        temp[tmp_path.size() + 1] = '\0';
        auto res = mkdtemp(temp);
        if (res == nullptr) {
            std::cout << "Error trying to create a temporal folder" << std::endl;
        }
        tmp_folder = std::string(temp);
        ext = random_string(3);
    }

    explicit tmp_workspace(std::string folder,
                           std::string ext_,
                           bool rem_all) : tmp_folder(folder),
                                           ext(ext_),
                                           remove_all(rem_all){
        assert(file_exists(tmp_folder));
    }

    [[nodiscard]] std::string get_file(std::string const& prefix) const {
        return std::filesystem::path(tmp_folder) / std::string(prefix+"_"+ext);
    }

    void remove_file(std::string const& prefix) const {
        std::filesystem::path file =  std::filesystem::path(tmp_folder) / std::string(prefix+"_"+ext);
        bool res = remove(file);
        if(!res){
            std::cout<<"Error trying to remove "<<file<<std::endl;
            exit(1);
        }
    }

    ~tmp_workspace(){
        if(remove_all){
            std::filesystem::remove_all(tmp_folder);
        }
    }

    [[nodiscard]] std::string folder() const{
        return tmp_folder;
    }
};

template<class sym_type>
str_collection collection_stats(std::string& input_file);

template<class time_t>
void report_time(time_t start, time_t end, size_t padding){
    auto dur = end - start;
    auto h = std::chrono::duration_cast<std::chrono::hours>(dur);
    auto m = std::chrono::duration_cast<std::chrono::minutes>(dur -= h);
    auto s = std::chrono::duration_cast<std::chrono::seconds>(dur -= m);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur -= s);

    for(size_t i=0;i<padding;i++) std::cout<<" ";
    if(h.count()>0){
        std::cout<<"Elapsed time (hh:mm:ss.ms): "<<h.count()<<":"<<m.count()<<":"<<s.count()<<"."<<ms.count()<<std::endl;
    }else if(m.count()>0){
        std::cout<<"Elapsed time (mm:ss.ms): "<<size_t(m.count())<<":"<<s.count()<<"."<<ms.count()<<std::endl;
    }else if(s.count()>0){
        std::cout<<"Elapsed time (ss.ms): "<<size_t(s.count())<<"."<<ms.count()<<std::endl;
    }else{
        std::cout<<"Elapsed time (ms): "<<ms.count()<<std::endl;
    }
}

void report_mem_peak();
#endif //LHTIGS_UTILS_H
