/*#include <iostream>
#include "grammar_algorithms.h"

/ *#include "old/cdt/include/file_streams.hpp"
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
}* /

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
}*/

#include <thread>

#include "external/CLI11.hpp"
#include "grammar_algorithms.h"

struct arguments{
    std::string input_file;
    std::string output_file;

    std::string tmp_dir;
    size_t n_threads{};
    size_t n_chunks{};
    size_t chunk_size{};
    size_t n_tries=1;
    bool ver=false;
    bool det=false;
    uint8_t alph_bytes=1;
    std::string p_file;
    std::vector<std::string> grammars_to_merge;
    std::vector<std::string> position_list;
    std::string coord_file;
    std::string rand_acc;
    std::string version= "v1.0.1 alpha";
};

struct CellWidthValidator : public CLI::Validator {
    CellWidthValidator() {
        name_ = "ValidCellWidth";
        func_ = [](const std::string &str) {
            bool valid = (str=="1") || (str=="2") || (str=="4") || (str=="8");
            if(!valid)
                return std::string(str+" is not a valid number of bytes for a native integer type");
            else
                return std::string();
        };
    }
};
const static CellWidthValidator ValidCellWidth;

class MyFormatter : public CLI::Formatter {
public:
    MyFormatter() : Formatter() {}
    std::string make_option_opts(const CLI::Option *) const override { return ""; }
};

static void parse_app(CLI::App& app, struct arguments& args){

    auto fmt = std::make_shared<MyFormatter>();

    fmt->column_width(27);
    app.formatter(fmt);
    app.add_flag("-v,--version", args.ver, "Print the software version and exit");

    CLI::App* comp = app.add_subcommand("comp", "compress text");
    comp->add_option("TEXT", args.input_file, "Input file in one-string-per-line format")->check(CLI::ExistingFile)->required();
    comp->add_option("-o,--output-file", args.output_file, "Output file")->type_name("");
    comp->add_option("-a,--alphabet", args.alph_bytes, "Number of bytes of the input alphabet (def. 1)")->check(CLI::Range(1, 8))->default_val(1)->check(ValidCellWidth);
    comp->add_option("-t,--threads", args.n_threads, "Maximum number of parsing threads")->default_val(1);
    comp->add_option("-c,--text-chunks", args.n_chunks, "Number of text chunks in memory during the parsing (def. n_threads*2)")->default_val(0);
    comp->add_option("-C,--chunk-size", args.chunk_size, "Size in bytes of each text chunk (def. TEXT_SIZE*0.0025)")->default_val(0);
    comp->add_option("-T,--tmp", args.tmp_dir, "Temporary folder (def. /tmp/lcg.xxxx)")-> check(CLI::ExistingDirectory)->default_val("/tmp");
    comp->add_flag("-d,--deterministic", args.det, "The resulting grammar is always the same");
    comp->add_option("-p,--parsing-functions", args.p_file, "File with the hash functions (PF format) to parse the text")->check(CLI::ExistingFile);

    CLI::App* par_func = app.add_subcommand("par", "extract parsing functions from a grammar");
    par_func->add_option("GRAM", args.input_file, "Input grammar in LCG format")->check(CLI::ExistingFile)->required();
    par_func->add_option("-o,--output-file", args.output_file, "Output file (PF format) with the hash functions")->type_name("");

    CLI::App* meta = app.add_subcommand("meta", "get the metadata of a grammar");
    meta->add_option("GRAM", args.input_file, "Input grammar in LCG format")->check(CLI::ExistingFile)->required();

    CLI::App* merge = app.add_subcommand("merge", "merge grammars");
    merge->add_option("GRAM LIST", args.grammars_to_merge, "Grammars to be merged")->check(CLI::ExistingFile)->required();
    merge->add_option("-o,--output-file", args.output_file, "Output grammar");

    CLI::App* access = app.add_subcommand("access", "random access");
    access->add_option("GRAM", args.input_file, "Input grammar in LCG format")->check(CLI::ExistingFile)->required();
    CLI::App* group = access->add_option_group("asdas");
    group->add_option("-p,--position", args.position_list, "Area to access in format str_idx:start-end (0-based)");
    group->add_option("-f,--pos-file", args.coord_file, "File with the list of positions to access");
    access->add_option("-o,--output-file", args.output_file, "Output file");
    group->require_option(1,2);

    app.require_subcommand(1,1);
    app.footer("Report bugs to <diego.diaz@helsinki.fi>");
}

template<class sym_type>
void run_int(std::string input_file, arguments& args) {
    tmp_workspace tmp_ws(args.tmp_dir, true, "lcg");
    std::cout<< "Temporary folder: "<<tmp_ws.folder()<<std::endl;
    gram_algo<uint8_t>(input_file, args.p_file, args.output_file, tmp_ws, args.n_threads, args.n_chunks, args.chunk_size);
}

int main(int argc, char** argv) {

    arguments args;

    CLI::App app("Locally-consistent grammar compression");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);

    if(args.ver){
        std::cout<<args.version<<std::endl;
        exit(0);
    }

    std::cout << "Input file:       " << args.input_file << std::endl;
    if(app.got_subcommand("comp")) {
        if (args.output_file.empty()) args.output_file = std::filesystem::path(args.input_file).filename();
        args.output_file = std::filesystem::path(args.output_file).replace_extension(".lcg");

        std::string input_collection = args.input_file;

        if (args.alph_bytes > 1) {
            std::cout << "Alphabet type:    integer" << std::endl;
        } else {
            std::cout << "Alphabet type:    byte" << std::endl;
        }
        if (args.alph_bytes == 1) {
            run_int<uint8_t>(input_collection, args);
        } else if (args.alph_bytes == 2) {
            run_int<uint16_t>(input_collection, args);
        } else if (args.alph_bytes == 4) {
            run_int<uint32_t>(input_collection, args);
        } else if (args.alph_bytes == 8) {
            run_int<uint64_t>(input_collection, args);
        }
    } else if(app.got_subcommand("par")){
        if (args.output_file.empty()) args.output_file = std::filesystem::path(args.input_file).filename();
        args.output_file = std::filesystem::path(args.output_file).replace_extension(".pf");
        get_par_functions(args.input_file, args.output_file);
    } else if(app.got_subcommand("meta")){
        print_metadata(args.input_file);
    } else if(app.got_subcommand("access")){
        lc_gram_t gram;
        load_from_file(args.input_file, gram);
        estimate_alt_encodings(gram);
    }
    return 0;
}
