//
// Created by Diaz, Diego on 10.9.2024.
//

#ifndef LCG_BUFF_VECTOR_H
#define LCG_BUFF_VECTOR_H

#include <cstdlib>
#include <cassert>
#include <algorithm>

template<class type, class alloc_t=mem<type>>
struct buff_vector{

    type *data= nullptr;
    size_t len=0;
    size_t cap=0;
    unsigned long offset=0;
    bool mem_alloc=true;

    unsigned long next_aligned_addr(uint8_t *& addr, size_t n_bytes){
        auto curr_addr = (uintptr_t)addr;
        uintptr_t alg_addr = ((curr_addr + (sizeof(type) - 1))/sizeof(type))*sizeof(type);
        uintptr_t diff = alg_addr-curr_addr;
        uintptr_t boundary = curr_addr + n_bytes;
        if(alg_addr<boundary){
            addr = (uint8_t *)alg_addr;
        }else{
            addr = nullptr;
        }
        return diff;
    }

    buff_vector()= default;

    buff_vector(uint8_t * addr, size_t n_bytes){
        if(addr!= nullptr){
            size_t new_offset = next_aligned_addr(addr, n_bytes);
            if(addr != nullptr){
                data = (type *) addr;
                offset = new_offset;
                cap = (n_bytes-offset)/sizeof(type);
                mem_alloc = false;
            }
        }
    }

    void report_buff_area() const{
        if(mem_alloc){
            std::cout<<0<<" -- "<<0<<std::endl;
        }else{
            std::cout<<(uintptr_t)data<<" -- "<<(uintptr_t)(data+len)<<std::endl;
        }
    }

    buff_vector& swap(buff_vector& other){
        std::swap(data, other.data);
        std::swap(len, other.len);
        std::swap(cap, other.cap);
        std::swap(offset, other.offset);
        std::swap(mem_alloc, other.mem_alloc);
        return *this;
    }

    void move_buffer(uint8_t * new_addr, size_t new_bytes){
        bool success_move=false;
        if(new_addr!= nullptr){
            size_t new_offset = next_aligned_addr(new_addr, new_bytes);
            size_t new_cap = (new_bytes-new_offset)/sizeof(type);
            if(new_addr != nullptr && len<=new_cap){
                type * old_data = data;
                data = (type *)new_addr;
                memcpy(data, old_data, len*sizeof(type));
                offset = new_offset;
                cap = new_cap;
                if(mem_alloc){
                    mem_alloc = false;//we are now in static memory
                    alloc_t::deallocate(old_data);
                }
                success_move = true;
            }
        }

        if(!success_move && !mem_alloc){//data was not moved and was originally in static memory
            mem_alloc = true;//move to dynamic memory
            cap = len;
            offset = 0;
            type *tmp = alloc_t::allocate(len);
            memcpy(tmp, data, len*sizeof(type));
            data = tmp;
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
            data = alloc_t::reallocate(data, len);
            cap = len;
        }
    }

    [[nodiscard]] inline size_t capacity() const{
        return cap;
    }

    void increase_capacity(size_t new_cap){
        if(new_cap>cap){
            if(!mem_alloc || data==nullptr){
                auto *tmp = alloc_t::allocate(new_cap);
                if(!mem_alloc && len>0){
                    memcpy(tmp, data, len*sizeof(type));
                }
                data = tmp;
                mem_alloc=true;
                offset=0;
            }else{
                data = alloc_t::reallocate(data, new_cap);
            }
            cap = new_cap;
#ifdef __linux__
            malloc_trim(0);
#endif
        }
    }

    inline type* get_data(){
        return data;
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

    void resize(size_t new_len){
        increase_capacity(new_len);
        len = new_len;
    }

    ~buff_vector(){
        if(mem_alloc){
            alloc_t::deallocate(data);
        }
    }

    [[nodiscard]] inline size_t static_buff_usage() const {
        if(mem_alloc){
            return 0;
        }else{
            return (len*sizeof(type))+offset;
        }
    }

    [[nodiscard]] inline size_t mem_usage() const {
        if(mem_alloc){
            return cap*sizeof(type);
        } else{
            return (len*sizeof(type))+offset;
        }
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

    [[nodiscard]] bool called_malloc() const{
        return mem_alloc;
    }
};
#endif //LCG_BUFF_VECTOR_H
