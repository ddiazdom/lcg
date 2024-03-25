//
// Created by Diaz, Diego on 25.3.2024.
//

#ifndef LCG_MMAP_ALLOCATOR_H
#define LCG_MMAP_ALLOCATOR_H
#include <sys/mman.h>
#include <unistd.h>

void * mmap_allocate(size_t n_bytes){
    size_t page_size =  sysconf(_SC_PAGE_SIZE);
    n_bytes = INT_CEIL(n_bytes, page_size)*page_size;
    void *ptr = mmap(nullptr, n_bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (ptr == MAP_FAILED){
        std::cerr << "Memory mapping failed\n";
        return nullptr;
    }
    return ptr;
}

void mmap_deallocate(void * ptr, size_t n_bytes){
    size_t page_size =  sysconf(_SC_PAGE_SIZE);
    n_bytes = INT_CEIL(n_bytes, page_size)*page_size;
    if (munmap(ptr, n_bytes)){
        std::cerr << "Memory mapping failed\n";
    }
}

void* mmap_reallocate(void* ptr, size_t old_size, size_t new_size) {

    size_t page_size =  sysconf(_SC_PAGE_SIZE);
    new_size = INT_CEIL(new_size, page_size)*page_size;
    old_size = INT_CEIL(old_size, page_size)*page_size;

    void* new_ptr = mmap(nullptr, new_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (new_ptr == MAP_FAILED) {
        std::cerr << "Memory mapping failed\n";
        return nullptr;
    }

    // Copy the contents from the old memory region to the new one
    memcpy(new_ptr, ptr, old_size < new_size ? old_size : new_size);

    // Free the old memory region
    mmap_deallocate(ptr, old_size);

    return new_ptr;
}

#endif //LCG_MMAP_ALLOCATOR_H
