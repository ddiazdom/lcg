//
// Created by Diaz, Diego on 9.12.2024.
//

#ifndef LCG_LOGGER_H
#define LCG_LOGGER_H
#include <iostream>

enum log_lvl{ERROR=0, WARNING=1, INFO=2, DEBUG=3};

template<log_lvl lvl, bool overwrite=false, bool suppress=false>
struct logger{

    inline static void error(const std::string& msg){
        std::cout<<"error: "<<msg<<std::endl;
    }

    inline static void info(const std::string& msg) {
        if constexpr (suppress && lvl!=INFO) return;
        if constexpr (lvl>=INFO){
            if constexpr (overwrite){
                std::cout<< "\r" << msg << std::flush;
            }else{
                std::cout<<msg<<std::endl;
            }
        }
    }

    inline static void debug(const std::string& msg) {
        if constexpr (lvl>=DEBUG){
            std::cout<<msg<<std::endl;
        }
    }

    inline static void warning(const std::string& msg) {
        if constexpr (lvl>=WARNING){
            std::cout<<msg<<std::endl;
        }
    }
};
#endif //LCG_LOGGER_H
