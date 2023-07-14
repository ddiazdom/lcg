#include <thread>

#include "CLI11.hpp"
#include "utils.h"
#include "grammar_algorithms.h"

struct arguments{
    std::string input_file;
    std::string output_file;

    std::string tmp_dir;
    size_t n_threads{};
    size_t n_tries=1;
    bool ver=false;
    bool det=false;
    uint8_t alph_bytes=1;
    std::string p_file;
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

    app.add_option("TEXT", args.input_file, "Input file in one-string-per-line format")->check(CLI::ExistingFile)->required();
    app.add_option("-o,--output-file", args.output_file, "Output file")->type_name("");
    app.add_option("-a,--alphabet", args.alph_bytes, "Number of bytes of the input alphabet (def. 1)")->check(CLI::Range(1, 8))->default_val(1)->check(ValidCellWidth);
    auto n_tries = app.add_option("-n,--random-tries", args.n_tries, "Number of random tries to find a small parsing (def. 1)");
    //app.add_option("-t,--threads", args.n_threads, "Maximum number of working threads")->default_val(1);
    app.add_option("-T,--tmp", args.tmp_dir, "Temporary folder (def. /tmp/lcg.xxxx)")-> check(CLI::ExistingDirectory)->default_val("/tmp");
    app.add_flag("-v,--version", args.ver, "Print the software version and exit");
    auto det_flag = app.add_flag("-d,--deterministic", args.det, "The resulting grammar is always the same");
    app.add_flag("-p,--parsing-functions", args.p_file, "List of hash functions to parse the text");
    app.footer("Report bugs to <diego.diaz@helsinki.fi>");
    n_tries->excludes(det_flag);
}

template<class sym_type>
void run_int(std::string input_collection, arguments& args){

    tmp_workspace tmp_ws(args.tmp_dir, true, "lcg");
    std::cout<< "Temporary folder: "<<tmp_ws.folder()<<std::endl;
    gram_algo<sym_type>(input_collection, args.output_file, tmp_ws, args.n_tries, args.n_threads);
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

    std::cout << "Input file:       "<<args.input_file<<std::endl;
    if(args.output_file.empty()) args.output_file = std::filesystem::path(args.input_file).filename();
    args.output_file = std::filesystem::path(args.output_file).replace_extension(".lcg");

    std::string input_collection = args.input_file;

    if(args.alph_bytes>1){
        std::cout<<"Alphabet type:    integer"<<std::endl;
    }else{
        std::cout<<"Alphabet type:    byte"<<std::endl;
    }

    if(args.alph_bytes==1){
        run_int<uint8_t>(input_collection, args);
    }else if(args.alph_bytes==2){
        run_int<uint16_t>(input_collection, args);
    }else if(args.alph_bytes==4){
        run_int<uint32_t>(input_collection, args);
    } else if(args.alph_bytes==8){
        run_int<uint64_t>(input_collection, args);
    }
    return 0;
}
