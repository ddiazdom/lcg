//
// Created by Diaz, Diego on 3.3.2022.
//

#ifndef CDT_COMMON_H
#define CDT_COMMON_H

#include <iostream>
#include <fstream>

uint8_t sym_width(unsigned long val);

size_t next_power_of_two(unsigned long val);

size_t prev_power_of_two(unsigned long val);

bool is_power_of_two(unsigned long val);

template<class data_type>
void load_from_file(std::string const& file, data_type& dt){
    std::ifstream ifs(file, std::ios::binary);
    dt.load(ifs);
    ifs.close();
}

template<class data_type>
size_t store_to_file(std::string const& file, data_type& dt){
    std::ofstream ofs(file, std::ios::binary);
    size_t written_bytes = dt.serialize(ofs);
    ofs.close();
    return written_bytes;
}

template<class vector_t>
size_t serialize_plain_vector(std::ostream& ofs, vector_t& vector){
    size_t n = vector.size();
    ofs.write((char *)&n, sizeof(n));
    ofs.write((char *)vector.data(), (std::streamsize)(sizeof(typename vector_t::value_type)*n));
    return sizeof(n)+ sizeof(typename vector_t::value_type)*n;
}

template<class vector_t>
size_t basic_store_vector_to_file(std::string const& file, vector_t& vector){
    std::ofstream ofs(file, std::ios::binary);
    ofs.write((char *)vector.data(), (std::streamsize)(sizeof(typename vector_t::value_type)*vector.size()));
    ofs.close();
    return sizeof(typename vector_t::value_type)*vector.size();
}

template<class vector_t>
size_t basic_load_vector_from_file(std::string const& file, vector_t& vector){
    std::ifstream ifs(file, std::ios::binary);
    ifs.seekg(0, std::ios::end);
    std::streamsize n_bytes = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    size_t type_bytes = sizeof(typename vector_t::value_type);
    vector.resize(n_bytes/type_bytes);
    ifs.read((char *)vector.data(), n_bytes);
    ifs.close();
    return n_bytes;
}

template<class size_type>
size_t serialize_raw_vector(std::ostream& ofs, size_type * vector, size_t len){
    ofs.write((char *)&len, sizeof(len));
    ofs.write((char *)vector, (std::streamsize)(sizeof(size_type)*len));
    return sizeof(len)+ sizeof(size_type)*len;
}

template<class val_type>
size_t serialize_elm(std::ostream& ofs, val_type value){
    ofs.write((char *)&value, sizeof(val_type));
    return sizeof(val_type);
}

template<class vector_t>
void load_plain_vector(std::istream& ifs, vector_t& vector){
    size_t n=0;
    ifs.read((char *)&n, sizeof(n));
    vector.resize(n);
    ifs.read((char *)vector.data(), (std::streamsize)(sizeof(typename vector_t::value_type)*n));
}

template<class size_type, class len_type>
void load_raw_vector(std::istream& ifs, size_type*& vector, len_type& len){
    ifs.read((char *)&len, sizeof(len_type));
    if(vector== nullptr){
        vector = (size_type *) malloc(sizeof(size_type)*len);
    }else{
        vector = (size_type *) realloc(vector, sizeof(size_type)*len);
    }
    ifs.read((char *)vector, (std::streamsize)(sizeof(size_type)*len));
}

template<class val_type>
void load_elm(std::istream& ifs, val_type& value){
    ifs.read((char *)&value, sizeof(val_type));
}

template<class vector_type>
void store_pl_vector(std::string const& file, vector_type& vector){
    std::ofstream ofs(file, std::ios::binary);
    serialize_plain_vector<vector_type>(ofs, vector);
    ofs.close();
}

template<class vector_type>
void load_pl_vector(std::string const& file, vector_type& vector){
    std::ifstream ifs(file, std::ios::binary);
    load_plain_vector<vector_type>(ifs, vector);
    ifs.close();
}
#endif //CDT_COMMON_H
