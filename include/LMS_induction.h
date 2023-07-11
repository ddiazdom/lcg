//
// Created by Diaz, Diego on 4.1.2022.
//

#ifndef LPG_COMPRESSOR_LMS_INDUCTION_H
#define LPG_COMPRESSOR_LMS_INDUCTION_H
#include "common.h"
#include "bwt_io.h"

namespace exact_algo {

template <class vector_type, class sa_type>
void induce_L_type(sa_type& sa, const dictionary &dict, vector_type& buckets, bv_t& solved_sym){

    size_t pos, l_sym, bck;
    vector_type ind_bck(dict.alphabet+1, 0);
    bool lcs, first_eq=false;
    size_t lb=1, flag;

    for(size_t i=0;i<sa.size();i++) {

        pos = sa[i];
        lcs = pos & 1UL;
        pos>>=1UL;
        if(pos==0) continue;

        if(!lcs) {
            first_eq = false;
            lb = i+1;
        }
        if(pos==1 || dict.d_lim[pos-2]) continue;

        bck = dict.dict.read(pos-1);
        l_sym = dict.dict.read(pos-2);

        if(!solved_sym[l_sym] && l_sym>=bck){
            if(!first_eq && l_sym==bck){
                first_eq = true;
                lcs = false;
            }

            flag = (ind_bck[l_sym]>=lb && lcs);
            sa[buckets[l_sym]++] = ((pos-1)<<1UL) | flag ;
            ind_bck[l_sym] = i+1;
        }
    }
}

template <class vector_type, class sa_type>
void induce_S_type(sa_type& sa, dictionary& dict, vector_type& buckets, bv_t& solved_sym){

    size_t pos, bck, l_sym, ind_pos;
    vector_type ind_bck(dict.alphabet+1, 0);
    bool lcs, first_eq=true, new_break, p_lcs=false;
    size_t rb=sa.size(), flag;

    for(size_t i=sa.size();i-->0;) {

        pos = sa[i];
        lcs = pos & 1UL;
        pos>>=1UL;
        assert(pos>0);

        if(!p_lcs) {
            first_eq = false;
            rb = i+1;
        }

        if(!(pos==1 || dict.d_lim[pos-2])){
            bck = dict.dict.read(pos-1);
            l_sym = dict.dict.read(pos-2);
            assert(l_sym<dict.end_str_dummy);

            if(!solved_sym[l_sym] && (l_sym < bck || (l_sym==bck && i >= buckets[bck]))){

                ind_pos = buckets[l_sym]--;
                flag = (lcs && ind_pos>0 && sa[ind_pos-1]==0);
                sa[ind_pos] = ((pos-1)<<1UL) | flag;

                new_break = !first_eq && l_sym==bck;
                if(new_break) first_eq = true;

                if(ind_bck[l_sym]==0 || ind_bck[l_sym]>rb || new_break){
                    sa[ind_pos+1] &= ~1UL;
                    if(ind_pos+1==i) lcs=false;
                }
                ind_bck[l_sym] = i+1;
            }
        }
        p_lcs = lcs;
    }
}

template <class vector_type, class sa_type>
void suffix_induction(dictionary &dict, sa_type& sa){

    vector_type buckets(dict.alphabet+1, 0);

    for(size_t i=0;i<dict.dict.size();i++){
        buckets[dict.dict[i]]++;
    }

    size_t acc=0, freq, max_freq=0;
    for(size_t i=0;i<dict.alphabet;i++){
        freq = buckets[i];
        buckets[i] = acc;
        acc +=freq;
        if(freq>max_freq) max_freq = freq;
    }
    buckets[dict.alphabet] = acc;

    vector_t freqs(dict.alphabet+1, 0, sdsl::bits::hi(max_freq)+1);
    bv_t solv_syms(dict.alphabet+1, false);
    size_t sym;

    for(size_t i=dict.dict.size();i-->0;) {
        sym = dict.dict[i];
        freq = buckets[sym+1]-buckets[sym];

        if(freq==1){
            //if it is a suffix of a phrase, we don't need to mark it
            sa[buckets[sym]] = ((i+1)<<1UL);
            solv_syms[sym]= true;
        }else if(dict.d_lim[i]){
            sa[buckets[sym+1]-1-freqs[sym]++] = ((i+1)<<1UL) | 1UL;
        }
    }

    size_t pos;
    for(size_t i=0;i<dict.alphabet;i++){
        if(freqs[i]>0){
            pos = buckets[i+1]-freqs[i];
            sa[pos] &= ~1UL;
        }
    }
    freqs.erase();

    for(size_t i=0;i<solv_syms.size();i++){
        if(solv_syms[i]) buckets[i]++;
    }

    induce_L_type<vector_type>(sa,  dict, buckets, solv_syms);

    //move the pointer of every bucket to the last S symbol
    for(size_t bck=0;bck<dict.alphabet;bck++){
        size_t ptr = buckets[bck];
        size_t limit = buckets[bck+1];
        while(ptr<limit && sa[ptr]==0) ptr++;
        ptr--;
        buckets[bck] = ptr;
    }

    induce_S_type<vector_type>(sa, dict, buckets, solv_syms);

#ifdef __linux__
    malloc_trim(0);
#endif
}

}
#endif //LPG_COMPRESSOR_LMS_INDUCTION_H
