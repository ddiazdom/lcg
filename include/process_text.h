//
// Created by Diaz, Diego on 19.7.2023.
//

#ifndef TEXT_PARSER_PROCESS_TEXT_H
#define TEXT_PARSER_PROCESS_TEXT_H
#include "parsing_strategies.h"
#include "cds/utils.h"
#include "cds/ts_string_map.h"
#include "grammar.h"

#ifdef __linux__
#include <malloc.h>
#endif

template<class parser_t, class map_type>
void par_round(std::string& input_file, std::string& output_file, parsing_opts& p_opts, lc_gram_buffer_t& gram_buffer,
               tmp_workspace& ws){

    map_type map(4, 0.6, p_opts.n_threads);
    std::cout <<"    Creating the dictionary" << std::flush;
    auto start = std::chrono::steady_clock::now();
    p_opts.max_mt_sym_in_buff = create_dict_full_scan<parser_t, map_type>()(input_file, map, p_opts);
    map.shrink_to_fit();
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 1);

    std::cout <<"    Assigning metasymbols" << std::flush;
    start = std::chrono::steady_clock::now();
    parser_t::assign_met(map, p_opts, gram_buffer, ws);
    end = std::chrono::steady_clock::now();
    report_time(start, end, 3);

    uint64_t p_size;
    std::cout <<"    Parsing the text" << std::flush;
    start = std::chrono::steady_clock::now();
    if(p_opts.p_alph_bytes==1){
        if(p_opts.parse_compressible){
            p_size = parse_text<parser_t, map_type, vbyte_decoder<uint8_t>>()(input_file, output_file, map, p_opts);
        }else{
            p_size = parse_text<parser_t, map_type, plain_decoder<uint8_t>>()(input_file, output_file, map, p_opts);
        }
    }else if(p_opts.p_alph_bytes==2){
        if(p_opts.parse_compressible){
            p_size = parse_text<parser_t, map_type, vbyte_decoder<uint16_t>>()(input_file, output_file, map, p_opts);
        }else{
            p_size = parse_text<parser_t, map_type, plain_decoder<uint16_t>>()(input_file, output_file, map, p_opts);
        }
    } else if(p_opts.p_alph_bytes<=4){
        if(p_opts.parse_compressible){
            p_size = parse_text<parser_t, map_type, vbyte_decoder<uint32_t>>()(input_file, output_file, map, p_opts);
        }else{
            p_size = parse_text<parser_t, map_type, plain_decoder<uint32_t>>()(input_file, output_file, map, p_opts);
        }
    } else if(p_opts.p_alph_bytes<=8){
        if(p_opts.parse_compressible){
            p_size = parse_text<parser_t, map_type, vbyte_decoder<uint64_t>>()(input_file, output_file, map, p_opts);
        }else{
            p_size = parse_text<parser_t, map_type, plain_decoder<uint64_t>>()(input_file, output_file, map, p_opts);
        }
    } else{
        std::cerr<<"Incorrect number of bytes for the parse : "<<p_opts.p_alph_bytes<<std::endl;
        exit(1);
    }
    p_opts.str_ptr->back() = file_size(output_file);

    if(p_opts.parse_compressible){
        p_opts.sym_perm.swap(p_opts.new_sym_perm);
        std::vector<uint32_t>().swap(p_opts.new_sym_perm);
        assert(off_t(p_opts.vbyte_size)==p_opts.str_ptr->back());
    }else{
        std::vector<uint32_t>().swap(p_opts.sym_perm);
        std::vector<uint32_t>().swap(p_opts.new_sym_perm);
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


template<class sym_type>
void process_text(std::string& i_file, std::string& o_file, tmp_workspace& ws, size_t n_threads) {

    //std::cout <<"Parallel parsing " << std::endl;

    struct stat st{};
    if (stat(i_file.c_str(), &st) != 0) return;

    std::string tmp_i_file = ws.get_file("tmp_input");
    text_stats txt_stats;

    using map_t = par_string_map<size_t>;

    compute_text_stats<sym_type>(i_file, txt_stats);

    auto chunk_size = off_t(next_power_of_two(size_t(float(st.st_size) * 0.0025)));
    size_t active_chunks = n_threads * 2;
    //auto active_chunks = 1;
    //auto chunk_size = 1024*1024*8;
    //std::cout << "Approximate usage for the buffers " << (active_chunks * chunk_size) / (1024 * 1024) << " MiB "<< std::endl;
    //auto chunk_size = 5912008;
    //active_chunks =1;
    //parsing_opts p_opts = {chunk_size, active_chunks, n_threads, false, 0.7};
    parsing_opts p_opts;
    p_opts.chunk_size = chunk_size;
    p_opts.active_chunks = active_chunks;
    p_opts.n_threads = n_threads;
    p_opts.str_ptr = &txt_stats.str_ptrs;

    p_opts.vbyte_threshold = 0.7;
    p_opts.vbyte_alphabet_threshold = 16777216;

    std::cout<<"Parsing settings"<<std::endl;
    std::cout<<"Parser threads           : "<<p_opts.n_threads<<std::endl;
    std::cout<<"Chunk_size               : "<<p_opts.chunk_size<<" bytes "<<std::endl;
    std::cout<<"Active chunks            : "<<p_opts.active_chunks<<std::endl;
    std::cout<<"Vbyte threshold          : "<<p_opts.vbyte_threshold<<std::endl;
    std::cout<<"Vbyte alphabet threshold : "<<p_opts.vbyte_alphabet_threshold<<" symbols "<<std::endl;


    std::string gram_file = ws.get_file("gram_file");
    lc_gram_buffer_t gram_buff(gram_file, txt_stats.alphabet, txt_stats.str_ptrs, txt_stats.longest_str, txt_stats.sep_sym);

    size_t iter =0;
    std::cout << "  Parsing round " << ++iter << std::endl;
    auto start = std::chrono::steady_clock::now();
    par_round<lms_parsing<plain_decoder<sym_type>>, map_t>(i_file, tmp_i_file, p_opts, gram_buff, ws);
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 2);

    size_t n_str = p_opts.str_ptr->size()-1;
    while(p_opts.tot_phrases!=n_str) {
        std::cout << "  Parsing round " << ++iter << std::endl;
        start = std::chrono::steady_clock::now();
        if (p_opts.p_alph_bytes == 1) {
            if(p_opts.parse_compressible){
                par_round<lms_parsing<vbyte_decoder<uint8_t>>, map_t>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }else{
                par_round<lms_parsing<plain_decoder<uint8_t>>, map_t>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }
        } else if (p_opts.p_alph_bytes == 2) {
            if(p_opts.parse_compressible){
                par_round<lms_parsing<vbyte_decoder<uint16_t>>, map_t>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }else{
                par_round<lms_parsing<plain_decoder<uint16_t>>, map_t>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }
        } else if (p_opts.p_alph_bytes <= 4){
            if(p_opts.parse_compressible){
                par_round<lms_parsing<vbyte_decoder<uint32_t>>, map_t>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }else{
                par_round<lms_parsing<plain_decoder<uint32_t>>, map_t>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }
        } else if(p_opts.p_alph_bytes<=8) {
            if(p_opts.parse_compressible){
                par_round<lms_parsing<vbyte_decoder<uint64_t>>, map_t>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }else{
                par_round<lms_parsing<plain_decoder<uint64_t>>, map_t>(tmp_i_file, o_file, p_opts, gram_buff, ws);
            }
        } else{
            std::cerr<<" Incorrect number of bytes of the input text: "<<p_opts.p_alph_bytes<<std::endl;
            exit(1);
        }
        end = std::chrono::steady_clock::now();
        report_time(start, end, 2);

        remove(tmp_i_file.c_str());
        rename(o_file.c_str(), tmp_i_file.c_str());
    }
    //gram_buff.insert_comp_string<size_t>(tmp_i_file);
}
#endif //TEXT_PARSER_PROCESS_TEXT_H
