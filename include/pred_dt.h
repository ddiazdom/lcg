//
// Created by Diaz, Diego on 20.9.2024.
//

#ifndef LCG_PRED_DT_H
#define LCG_PRED_DT_H

template<int block=32>
struct pred_dt{

    size_t alph_a;
    size_t alph_b;
    size_t b_start=0;
    size_t nb_start=0;
    size_t len=0;
    size_t prev_elm = std::numeric_limits<size_t>::max();
    size_t str_pos=0;
    size_t str_width;
    size_t smp_pos=0;
    size_t smp_width;

    bitstream<size_t> samples;
    bitstream<size_t> stream;

    pred_dt(size_t _alph_a, size_t _alph_b): alph_a(_alph_a),
                                             alph_b(_alph_b){
        size_t n_samples = INT_CEIL(_alph_a, block);
        str_width = sym_width(block);
        stream.reserve_in_bits(alph_b*str_width);
        smp_width = sym_width(alph_b);
        samples.reserve_in_bits(n_samples*smp_width);
    }

    [[nodiscard]] inline size_t find(size_t val) const {
        size_t block_id = val/block;
        val%=block;

        size_t s = samples.read(block_id*smp_width, (block_id*smp_width)+smp_width-1);
        size_t e = samples.read((block_id+1)*smp_width, ((block_id+1)*smp_width)+smp_width-1);

        for(size_t j=s;j<e;j+=str_width){
            size_t elm = stream.read(j, j+str_width-1);
            if(elm>val){
                break;
            }
        }
        return 0;
    }

    void append(size_t elm){
        assert(elm<=alph_a);
        if(prev_elm!=elm && elm >= nb_start){
            b_start = nb_start;
            nb_start+=block;
            assert(len<alph_b);
            samples.write(smp_pos, smp_pos+smp_width-1, len);
            smp_pos+=smp_width;
        }

        stream.write(str_pos, str_pos+str_width-1, elm-b_start);
        str_pos+=str_width;
        len++;
        prev_elm = elm;
    }

    size_t operator[](size_t val) const {
        assert(val<alph_a);
        return val+find(val);
    }
};
#endif //LCG_PRED_DT_H
