//
// Created by Diaz, Diego on 5.3.2022.
//
#include "cdt_common.hpp"

uint8_t sym_width(unsigned long val){
    if(val==0) return 0;
    return (sizeof(unsigned long)*8) - __builtin_clzl(val);
}

size_t next_power_of_two(unsigned long val){
    uint8_t width = sym_width(val);
    return 1UL<<width;
}

size_t prev_power_of_two(unsigned long val){
    uint8_t width = sym_width(val);
    return 1UL<<(width-1);
}

bool is_power_of_two(unsigned long val){
    return !(val & (val-1));
}
