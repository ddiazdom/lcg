//
// Created by Diego Diaz on 2/10/20.
//

#ifndef LMS_COMPRESSOR_MEMORY_HANDLER_HPP
#define LMS_COMPRESSOR_MEMORY_HANDLER_HPP

#endif //LMS_COMPRESSOR_MEMORY_HANDLER_HPP

class allocator{
public:
    template<class value_t>
    inline static value_t *allocate(size_t n_vals, bool init = true, value_t def_val = 0) {
        auto all_mem = (value_t *)malloc(sizeof(value_t) * n_vals);
        if(init){
            for(size_t i=0; i < n_vals;i++) all_mem[i] = def_val;
        }
        return all_mem;
    }

    template<class value_t>
    static value_t *
    reallocate(value_t *pointer, size_t old_size, size_t new_size, bool init = true, value_t def_val = 0) {
        auto re_all_mem = (value_t*)realloc(pointer, sizeof(value_t)*new_size);
        if(new_size>old_size && init){
            for(size_t i=old_size;i<new_size;i++){
                re_all_mem[i] = def_val;
            }
        }
        return re_all_mem;
    }

    template<class value_t>
    static void deallocate(value_t * pointer){
        if(pointer!= nullptr){
            free(pointer);
        }
    }
};

/* Obtained from
 * https://stackoverflow.com/questions/36517825/is-stephen-lavavejs-mallocator-the-same-in-c11/36521845#36521845
 */
template <class T> struct mallocator {

    typedef T value_type;

    mallocator() noexcept = default; // default ctor not required

    template <class U> explicit mallocator(const mallocator<U>&) noexcept { }

    template <class U>
    bool operator==(const mallocator<U>&) const noexcept {
        return true;
    }

    template <class U> bool operator!=(const mallocator<U>&) const noexcept {
        return false;
    }

    T * allocate(const size_t n) const {
        if (n == 0) {return nullptr;}
        if (n > static_cast<size_t>(-1) / sizeof(T)) {
            throw std::bad_array_new_length();
        }
        void * const pv = malloc(n * sizeof(T));
        if (!pv) { throw std::bad_alloc(); }
        return static_cast<T *>(pv);
    }

    void deallocate(T * const p, size_t) const noexcept {
        free(p);
    }
};
