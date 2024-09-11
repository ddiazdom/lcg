//
// Created by Diaz, Diego on 27.10.2021.
//

#ifndef CDS_UTILS_H
#define CDS_UTILS_H

#include <chrono>
#include <iostream>
#include <utility>
#include <vector>
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <cassert>
#include <sstream>
#include <iomanip>

#ifdef __APPLE__
#include <unistd.h>
#endif

struct str_collection {
    std::vector<uint8_t> alphabet;
    std::vector<size_t> sym_freqs;
    std::vector<long> str_ptrs;
    size_t n_strings=0;
    size_t longest_string=0;
    uint8_t max_sym=0;
    uint8_t min_sym=255;
    size_t n_char=0;
};

std::ifstream::pos_type file_size(std::string & filename);


//check if file exists
bool file_exists(const std::filesystem::path& p, std::filesystem::file_status const& s = std::filesystem::file_status{});

//produce a random string
std::string random_string(size_t length);

#ifdef __linux__
void empty_page_cache(const char *filename);
#endif

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
                           bool rem_all) : tmp_folder(std::move(folder)),
                                           ext(std::move(ext_)),
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

//https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
std::vector<std::string> split (const std::string &s, char delim);


str_collection collection_stats(std::string& input_file);

//from https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}

template<class time_t>
void report_time(time_t start, time_t end, size_t padding){
    auto dur = end - start;
    auto h = std::chrono::duration_cast<std::chrono::hours>(dur);
    auto m = std::chrono::duration_cast<std::chrono::minutes>(dur -= h);
    auto s = std::chrono::duration_cast<std::chrono::seconds>(dur -= m);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur -= s);

    for(size_t i=0;i<padding;i++) std::cout<<" ";
    if(h.count()>0){
        std::cout<<"Elapsed time (hh:mm:ss.ms): "<<std::setfill('0')<<std::setw(2)<<h.count()<<":"<<std::setfill('0')<<std::setw(2)<<m.count()<<":"<<std::setfill('0')<<std::setw(2)<<s.count()<<"."<<ms.count()<<std::endl;
    }else if(m.count()>0){
        std::cout<<"Elapsed time (mm:ss.ms): "<<std::setfill('0')<<std::setw(2)<<size_t(m.count())<<":"<<std::setfill('0')<<std::setw(2)<<s.count()<<"."<<ms.count()<<std::endl;
    }else if(s.count()>0){
        std::cout<<"Elapsed time (ss.ms): "<<std::setfill('0')<<std::setw(2)<<size_t(s.count())<<"."<<ms.count()<<std::endl;
    }else{
        std::cout<<"Elapsed time (ms): "<<ms.count()<<std::endl;
    }
}

template<class time_t>
std::string report_speed(off_t bytes, time_t start, time_t end){

    auto dur = end - start;
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(dur);
    std::vector<std::string> byte_units = {"B", "KB", "MB", "GB", "TB"};
    std::vector<std::string> time_units = {"us", "ms", "s", "min", "h"};

    float byte_scale = 1024.0;
    std::vector<float> time_scales = {1000, 1000, 60, 60, 24};

    auto b = (float) bytes;
    size_t byte_unit_index =0;
    while(b >= byte_scale && byte_unit_index < (byte_units.size()-1)){
        b /= byte_scale;
        byte_unit_index += 1;
    }

    auto t = (float) us.count();
    size_t time_unit_index=0;
    while(t >= time_scales[time_unit_index] && time_unit_index < (time_scales.size()-1)){
        t /= time_scales[time_unit_index++];
    }

    float speed = b/t;

    std::ostringstream result;
    result <<std::fixed << std::setprecision(2) << speed<< " " << byte_units[byte_unit_index] << "/" << time_units[time_unit_index];
    return result.str();
}

std::string report_space(off_t bytes);

void report_mem_peak();

#endif //CDS_UTILS_H
