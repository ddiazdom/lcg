//
// Created by Diaz, Diego on 19.7.2023.
//

#ifndef SE_STRAT_LC_PARSING_H
#define SE_STRAT_LC_PARSING_H

#include "gram_buffer.h"
#include "cds/ts_string_map.h"
#include "parse_text.h"
#include "create_dict_full_scan.h"
#include "local_minima_parser.hpp"

#ifdef __linux__
#include <malloc.h>
#endif

template<class parser_t, class map_type, class text_chunk_t>
void par_round(std::string& input_file, std::string& output_file,
               parsing_opts& p_opts, lc_gram_buffer_t& gram_buffer, tmp_workspace& ws){

    p_opts.next_seed();
    map_type map(4, 0.6, p_opts.n_threads);
    std::cout <<"    Creating the dictionary" << std::flush;
    auto start = std::chrono::steady_clock::now();
    p_opts.max_mt_sym_in_buff = create_dict_full_scan<parser_t, map_type, text_chunk_t>()(input_file, map, p_opts);
    map.shrink_to_fit();
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 1);

    std::cout <<"    Assigning metasymbols" << std::flush;
    start = std::chrono::steady_clock::now();
    parser_t::template assign_met<map_type, typename text_chunk_t::decoder_t>(map, p_opts, gram_buffer, ws);
    end = std::chrono::steady_clock::now();
    report_time(start, end, 3);

    uint64_t p_size;
    std::cout <<"    Parsing the text" << std::flush;
    start = std::chrono::steady_clock::now();
    if(p_opts.p_alph_bytes==1){
        if(p_opts.parse_compressible){
            p_size = parse_text<parser_t, map_type, vbyte_decoder<uint8_t>, text_chunk_t>()(input_file, output_file, map, p_opts);
        }else{
            p_size = parse_text<parser_t, map_type, plain_decoder<uint8_t>, text_chunk_t>()(input_file, output_file, map, p_opts);
        }
    }else if(p_opts.p_alph_bytes==2){
        if(p_opts.parse_compressible){
            p_size = parse_text<parser_t, map_type, vbyte_decoder<uint16_t>, text_chunk_t>()(input_file, output_file, map, p_opts);
        }else{
            p_size = parse_text<parser_t, map_type, plain_decoder<uint16_t>, text_chunk_t>()(input_file, output_file, map, p_opts);
        }
    } else if(p_opts.p_alph_bytes<=4){
        if(p_opts.parse_compressible){
            p_size = parse_text<parser_t, map_type, vbyte_decoder<uint32_t>, text_chunk_t>()(input_file, output_file, map, p_opts);
        }else{
            p_size = parse_text<parser_t, map_type, plain_decoder<uint32_t>, text_chunk_t>()(input_file, output_file, map, p_opts);
        }
    } else if(p_opts.p_alph_bytes<=8){
        if(p_opts.parse_compressible){
            p_size = parse_text<parser_t, map_type, vbyte_decoder<uint64_t>, text_chunk_t>()(input_file, output_file, map, p_opts);
        }else{
            p_size = parse_text<parser_t, map_type, plain_decoder<uint64_t>, text_chunk_t>()(input_file, output_file, map, p_opts);
        }
    } else{
        std::cerr<<"Incorrect number of bytes for the parse : "<<p_opts.p_alph_bytes<<std::endl;
        exit(1);
    }
    p_opts.str_ptr->back() = file_size(output_file);

    if(p_opts.parse_compressible){
        p_opts.vbyte_sym_perm.swap(p_opts.new_vbyte_sym_perm);
        std::vector<uint32_t>().swap(p_opts.new_vbyte_sym_perm);
        assert(off_t(p_opts.vbyte_size)==p_opts.str_ptr->back());
    }else{
        std::vector<uint32_t>().swap(p_opts.vbyte_sym_perm);
        std::vector<uint32_t>().swap(p_opts.new_vbyte_sym_perm);
    }
    end = std::chrono::steady_clock::now();
    report_time(start, end, 8);
    assert(p_size == p_opts.tot_phrases);

    std::cout<<"    Stats: "<<std::endl;
    std::cout<<"      Number of phrases:         "<<map.size()<<std::endl;
    std::cout<<"      Parse size:                "<<p_opts.tot_phrases<<std::endl;
    std::cout<<"      Bytes of the parse:        "<<p_opts.str_ptr->back()<<std::endl;
    std::cout<<"      Parse is vbyte-compressed: "<<p_opts.parse_compressible<<std::endl;
    if(p_opts.parse_compressible){
        std::cout<<"        Vbyte comp. ratio:       "<<p_opts.vbyte_comp_ratio<<std::endl;
    }
    p_opts.p_round++;
}

template<class sym_type, class gram_type>
void lc_parsing_algo(std::string& i_file, std::string& gram_file, tmp_workspace& ws, size_t n_threads,
                     size_t active_chunks, size_t chunk_size, uint64_t par_seed) {

    struct stat st{};
    if (stat(i_file.c_str(), &st) != 0) return;

    std::string tmp_i_file = ws.get_file("tmp_input");
    std::string o_file = ws.get_file("tmp_out_file");

    using map_t = par_string_map<size_t>;

    parsing_opts p_opts;
    p_opts.active_chunks = active_chunks!=0? active_chunks : n_threads*2;
    p_opts.chunk_size = chunk_size!=0 ? (off_t)chunk_size : off_t(next_power_of_two(size_t(float(st.st_size) * 0.0025)));
    p_opts.n_threads = n_threads;
    p_opts.vbyte_threshold = 0.4;
    p_opts.vbyte_alphabet_threshold = 16777216;

    //std::random_device rd;  // Will be used to obtain a seed for the random number engine
    //std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::mt19937 gen(par_seed); // Standard mersenne_twister_engine seeded with a fixed value
    std::uniform_int_distribution<uint64_t> distrib(1, std::numeric_limits<uint64_t>::max());
    p_opts.p_seeds.resize(32);
    for(size_t i=0;i<32;i++){
        p_opts.p_seeds[i] = distrib(gen);
    }

    std::cout<<"  Settings"<<std::endl;
    std::cout<<"    Parsing mode              : long strings (> 4 GBs)"<<std::endl;
    std::cout<<"    Parsing threads           : "<<p_opts.n_threads<<std::endl;
    std::cout<<"    Active text chunks in RAM : "<<p_opts.active_chunks<<std::endl;
    std::cout<<"    Chunk size                : "<<report_space(p_opts.chunk_size)<<std::endl;
    std::cout<<"    Chunks' approx. mem usage : "<<report_space(off_t(p_opts.chunk_size*p_opts.active_chunks))<<std::endl;
    std::cout<<"    Vbyte encoding thresholds"<<std::endl;
    std::cout<<"      min. comp. ratio        : "<<p_opts.vbyte_threshold<<std::endl;
    std::cout<<"      max. alphabet size      : "<<p_opts.vbyte_alphabet_threshold<<" symbols "<<std::endl;

    std::cout<<"  Reading the input file before computing the grammar"<<std::endl;
    text_stats txt_stats;
    compute_text_stats<sym_type>(i_file, txt_stats);
    p_opts.str_ptr = &txt_stats.str_ptrs;
    p_opts.sep_sym = txt_stats.sep_sym;

    std::cout<<"    Stats"<<std::endl;
    std::cout<<"      Number of symbols: "<<txt_stats.str_ptrs.back()/sizeof(sym_type)<<std::endl;
    std::cout<<"      Alphabet size:     "<<txt_stats.alphabet.size()<<std::endl;
    std::cout<<"      Smallest symbol:   "<<(int)txt_stats.alphabet[0]<<std::endl;
    std::cout<<"      Greatest symbol:   "<<(int)txt_stats.alphabet.back()<<std::endl;
    std::cout<<"      Sep. symbol:       "<<txt_stats.sep_sym<<std::endl;
    std::cout<<"      Number of strings: "<<txt_stats.str_ptrs.size()-1<<std::endl;
    std::cout<<"      Longest string:    "<<txt_stats.longest_str/sizeof(sym_type)<<std::endl;

    lc_gram_buffer_t gram_buff(gram_file, txt_stats.alphabet, txt_stats.str_ptrs, txt_stats.longest_str, txt_stats.sep_sym);

    size_t iter =0;
    std::cout << "  Parsing round " << ++iter << std::endl;
    auto start = std::chrono::steady_clock::now();
    par_round<lm_parser, map_t, text_chunk<plain_decoder<sym_type>, true>>(i_file, tmp_i_file, p_opts, gram_buff, ws);
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 2);

    size_t n_str = p_opts.str_ptr->size()-1;

    while(p_opts.tot_phrases!=n_str) {

        std::cout << "  Parsing round " << ++iter << std::endl;
        start = std::chrono::steady_clock::now();

        if (p_opts.p_alph_bytes == 1) {
            if(p_opts.parse_compressible){
                par_round<lm_parser, map_t, text_chunk<vbyte_decoder<uint8_t>>>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }else{
                par_round<lm_parser, map_t, text_chunk<plain_decoder<uint8_t>>>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }
        } else if (p_opts.p_alph_bytes == 2) {
            if(p_opts.parse_compressible){
                par_round<lm_parser, map_t, text_chunk<vbyte_decoder<uint16_t>>>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }else{
                par_round<lm_parser, map_t, text_chunk<plain_decoder<uint16_t>>>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }
        } else if (p_opts.p_alph_bytes <= 4){
            if(p_opts.parse_compressible){
                par_round<lm_parser, map_t, text_chunk<vbyte_decoder<uint32_t>>>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }else{
                par_round<lm_parser, map_t, text_chunk<plain_decoder<uint32_t>>>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }
        } else if(p_opts.p_alph_bytes<=8) {
            if(p_opts.parse_compressible){
                par_round<lm_parser, map_t, text_chunk<vbyte_decoder<uint64_t>>>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }else{
                par_round<lm_parser, map_t, text_chunk<plain_decoder<uint64_t>>>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }
        } else{
            std::cerr<<" Incorrect number of bytes of the input text: "<<p_opts.p_alph_bytes<<std::endl;
            exit(1);
        }
        end = std::chrono::steady_clock::now();
        report_time(start, end, 4);

        remove(tmp_i_file.c_str());
        rename(o_file.c_str(), tmp_i_file.c_str());
    }

    if(p_opts.p_alph_bytes<=1){
        gram_buff.insert_comp_string<uint8_t>(tmp_i_file);
    }else if(p_opts.p_alph_bytes<=2){
        gram_buff.insert_comp_string<uint16_t>(tmp_i_file);
    }else if(p_opts.p_alph_bytes<=4){
        gram_buff.insert_comp_string<uint32_t>(tmp_i_file);
    } else if(p_opts.p_alph_bytes<=8){
        gram_buff.insert_comp_string<uint64_t>(tmp_i_file);
    }

    gram_type lc_gram;
    gram_buff.make_gram(lc_gram, par_seed);

    std::cout<<"  Breakdown of the resulting locally-consistent grammar:"<<std::endl;
    lc_gram.breakdown(4);
    store_to_file(gram_file, lc_gram);
}
#endif //SE_STRAT_LC_PARSING_H
