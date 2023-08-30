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
    CLI::App* comp = app.add_subcommand("comp", "compress text");
    comp->add_option("TEXT", args.input_file, "Input file in one-string-per-line format")->check(CLI::ExistingFile)->required();
    comp->add_option("-o,--output-file", args.output_file, "Output file")->type_name("");
    comp->add_option("-a,--alphabet", args.alph_bytes, "Number of bytes of the input alphabet (def. 1)")->check(CLI::Range(1, 8))->default_val(1)->check(ValidCellWidth);
    auto n_tries = comp->add_option("-n,--random-tries", args.n_tries, "Number of random tries to find a small parsing (def. 1)");
    //app.add_option("-t,--threads", args.n_threads, "Maximum number of working threads")->default_val(1);
    comp->add_option("-T,--tmp", args.tmp_dir, "Temporary folder (def. /tmp/lcg.xxxx)")-> check(CLI::ExistingDirectory)->default_val("/tmp");
    comp->add_flag("-v,--version", args.ver, "Print the software version and exit");
    auto det_flag = comp->add_flag("-d,--deterministic", args.det, "The resulting grammar is always the same");
    auto pf_cmd = comp->add_option("-p,--parsing-functions", args.p_file, "List of hash functions (PF format) to parse the text")->check(CLI::ExistingFile);
    n_tries->excludes(det_flag);
    n_tries->excludes(pf_cmd);

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
void run_int(std::string input_collection, arguments& args){
    tmp_workspace tmp_ws(args.tmp_dir, true, "lcg");
    std::cout<< "Temporary folder: "<<tmp_ws.folder()<<std::endl;
    gram_algo<sym_type>(input_collection, args.p_file, args.output_file, tmp_ws, args.n_tries, args.n_threads);
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

    }
    return 0;
}
