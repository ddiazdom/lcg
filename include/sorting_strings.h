//
// Created by Diaz, Diego on 16.1.2024.
//

#ifndef LCG_SORT_COLLECTION_H
#define LCG_SORT_COLLECTION_H

struct entry{
    size_t start;
    size_t len;
    size_t symbol;
};

template<class vector_t>
void sort_strings(std::vector<entry>& entries, const vector_t& strings, bool co_lex){

    if(co_lex){
        std::sort(entries.begin(), entries.end(), [&](auto const& a, auto const& b){

            size_t len = std::min(a.len, b.len);
            size_t inv_a = a.start+a.len-1;
            size_t inv_b = b.start+b.len-1;

            for(size_t i=0;i<len;i++){
                if(strings[a.start+i]!=strings[b.start+i]){
                    return strings[inv_a-i]<strings[inv_b-i];
                }
            }
            return a.len<b.len;
        });
    }else{
        std::sort(entries.begin(), entries.end(), [&](auto const& a, auto const& b){
            size_t len = std::min(a.len, b.len);
            for(size_t i=0;i<len;i++){
                if(strings[a.start+i]!=strings[b.start+i]){
                    return strings[a.start+i]<strings[b.start+i];
                }
            }
            return a.len<b.len;
        });
    }
}

#endif //LCG_SORT_COLLECTION_H
