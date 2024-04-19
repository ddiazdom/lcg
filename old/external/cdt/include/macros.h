//
// Created by diego on 02-09-20.
//

#ifndef LPG_COMPRESSOR_MACROS_H
#define LPG_COMPRESSOR_MACROS_H

#define INT_CEIL(a,b) (a>0? 1+(a-1)/b : 0)

#ifdef __APPLE__
#define GET_MEM_PEAK()\
struct rusage usage{};\
getrusage(RUSAGE_SELF, &usage);\
std::cout<<"Mem. peak:\t"<<(double(usage.ru_maxrss)/1000)/1000<<" MB"<<std::endl;
#else
#define GET_MEM_PEAK()\
struct rusage usage{};\
getrusage(RUSAGE_SELF, &usage);\
std::cout<<"mem peak:\t"<<((usage.ru_maxrss*1000)/(1024.0*1024.0))<<" MiB"<<std::endl;
#endif

#endif //LPG_COMPRESSOR_MACROS_H
