#include <iostream>
#include "grammar_algorithms.h"

/*#include "old/cdt/include/file_streams.hpp"
template<class sym_type>
void expand(std::string& input_file){
    i_file_stream<uint8_t> ifs(input_file, 1024*1204*8);
    std::string out_file = input_file+"_"+std::to_string(sizeof(sym_type))+"_bytes";
    o_file_stream<sym_type> ofs(out_file, 1024*1024*8, std::ios::out);
    for(size_t i=0;i<ifs.size();i++){
        ofs.push_back(ifs[i]);
    }
    std::cout<<ofs.size()<<" symbols in "<<out_file<<std::endl;
    ofs.close();
}*/

int main (int argc, char** argv){

    if(argc!=4){
        std::cout<<"usage: ./parallel_parser input_file n_bytes n_threads "<<std::endl;
        std::cout<<"rand_collection: random string collection to test against "<<std::endl;
        std::cout<<"n_threads: number of threads that will work on the hash table "<<std::endl;
        exit(0);
    }

    char * pEnd;
    std::string input_file = std::string(argv[1]);
    size_t n_bytes=strtoull(argv[2], &pEnd, 10);
    size_t n_threads=strtoull(argv[3], &pEnd, 10);

    std::string tmp_dir = "./";
    std::string output_gram = "output_gram.lcg";
    tmp_workspace tmp_ws(tmp_dir, true, "parse");
    std::string pf_file;

    auto start = std::chrono::steady_clock::now();
    switch (n_bytes) {
        case 1:
            gram_algo<uint8_t>(input_file, pf_file, output_gram, tmp_ws, n_threads);
            break;
        case 2:
            gram_algo<uint16_t>(input_file, pf_file, output_gram, tmp_ws, n_threads);
            break;
        case 4:
            gram_algo<uint32_t>(input_file, pf_file, output_gram, tmp_ws, n_threads);
            break;
        case 8:
            gram_algo<uint64_t>(input_file, pf_file, output_gram, tmp_ws, n_threads);
            break;
        default:
            exit(1);
    }
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 0);

    return 0;
}
