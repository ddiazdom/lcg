//
// Created by Diaz, Diego on 10.9.2024.
//

#ifndef LCG_BUFF_INT_VECTOR_H
#define LCG_BUFF_INT_VECTOR_H
#include <cstdlib>
#include <cassert>
#include <algorithm>

template<class type>
struct buff_vector{

    type *data= nullptr;
    size_t len=0;
    size_t cap=0;
    unsigned long offset=0;
    bool mem_alloc=true;

    unsigned long get_aligned_addr(char *& addr, size_t n_bytes){
        auto curr_addr = (uintptr_t)addr;
        uintptr_t alg_addr = ((curr_addr + (sizeof(type) - 1))/sizeof(type))*sizeof(type);
        uintptr_t diff = alg_addr-curr_addr;
        uintptr_t boundary = curr_addr + n_bytes;
        if(alg_addr<boundary){
            addr = (char *)alg_addr;
        }else{
            addr = nullptr;
        }
        return diff;
    }

    buff_vector()= default;

    buff_vector(char* address, size_t n_bytes){
        if(address!= nullptr){
            char * tmp = address;
            offset = get_aligned_addr(tmp, n_bytes);
            if(tmp != nullptr){
                data = (type *) tmp;
                cap = (n_bytes-offset)/sizeof(type);
                mem_alloc = false;
            }
        }
    }

    type& operator[](size_t idx){
        assert(idx<len);
        return data[idx];
    }

    const type& operator[](size_t idx) const {
        assert(idx<len);
        return data[idx];
    }

    [[nodiscard]] inline bool empty() const {
        return len==0;
    }

    inline void clear(){
        len=0;
    }

    [[nodiscard]] inline size_t size() const {
        return len;
    }

    void shrink_to_fit(){
        if(mem_alloc){
            cap = len;
            data = (type *)realloc(data, cap*sizeof(type));
        }
    }

    [[nodiscard]] inline size_t capacity() const{
        return cap;
    }

    void increase_capacity(size_t new_cap){
        if(new_cap>cap){
            cap = new_cap;
            if(!mem_alloc || data==nullptr){
                auto *tmp = (type *) malloc(cap*sizeof(type));
                if(!mem_alloc && len>0){
                    memcpy(tmp, data, len*sizeof(type));
                }
                data = tmp;
                mem_alloc=true;
            }else{
                data = (type *) realloc(data, cap*sizeof(type));
            }
        }
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    inline void push_back(type& new_val){
        if(len==cap){
            increase_capacity(std::max<size_t>(2,cap*2));
        }
        data[len++]=new_val;
    }

    template<typename... Args>
    inline void emplace_back(Args&&... args) {
        if (len == cap) {
            increase_capacity(std::max<size_t>(2,cap*2));
        }
        data[len++] = type(std::forward<Args>(args)...);
    }

    ~buff_vector(){
        if(mem_alloc){
            free(data);
#ifdef __linux__
            malloc_trim(0);
#endif
        }
    }

    [[nodiscard]] inline size_t byte_usage() const {
        return (len*sizeof(type))+offset;
    }

    // Define the iterator class within the container
    class iterator {
    private:
        type * ptr;  // Pointer to the element the iterator points to

    public:
        // Constructor
        explicit iterator(type* p = nullptr) : ptr(p) {}

        // Dereference operator
        type& operator*() const { return *ptr; }

        // Pointer access operator
        type* operator->() const { return ptr; }

        // Pre-increment operator (++it)
        iterator& operator++() {
            ++ptr;
            return *this;
        }

        // Post-increment operator (it++)
        iterator operator++(int) {
            iterator tmp = *this;
            ++ptr;
            return tmp;
        }

        // Equality comparison operator
        inline bool operator==(const iterator& other) const {
            return ptr == other.ptr;
        }

        // Inequality comparison operator
        inline bool operator!=(const iterator& other) const {
            return ptr != other.ptr;
        }
    };

    // Begin iterator
    [[nodiscard]] iterator begin() const {
        return iterator(data);
    }

    // End iterator (points one past the last element)
    iterator end() const {
        return iterator(data + len);
    }
};
#endif //LCG_BUFF_INT_VECTOR_H
