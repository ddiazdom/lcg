//
// Created by Diaz, Diego on 25.3.2024.
//

#ifndef LCG_SPECIAL_ALLOCATOR_H
#define LCG_SPECIAL_ALLOCATOR_H
#include <sys/mman.h>

void * allocate(size_t n_bytes){
    size_t page_size =  sysconf(_SC_PAGE_SIZE);
    n_bytes = INT_CEIL(n_bytes, page_size)*page_size;
    void *ptr = mmap(nullptr, n_bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (ptr == MAP_FAILED){
        perror("Error trying to mmap memory"), exit (EXIT_FAILURE);
        exit(0);
    }
    return ptr;
}

void * deallocate(char * ptr, size_t n_bytes){
    size_t page_size =  sysconf(_SC_PAGE_SIZE);
    n_bytes = INT_CEIL(n_bytes, page_size)*page_size;
    if (munmap(ptr, page_size)){
        perror("Error trying to unmmaping a pointer"), exit(EXIT_FAILURE);
        exit(0);
    }
}
#endif //LCG_SPECIAL_ALLOCATOR_H
