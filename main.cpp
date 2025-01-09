#include "logger.h"
#include "external/CLI11.hpp"
#include "grammar_algorithms.h"

struct arguments{
    std::string file_list;
    std::vector<std::string> input_files;
    std::vector<text_format> input_fmts;
    std::string i_file;
    std::string output_file;

    std::string tmp_dir = std::filesystem::temp_directory_path();
    size_t n_threads=1;
    size_t n_chunks{};
    off_t chunk_size{};
    float i_frac=0.0;
    bool ver=false;
    bool rand_acc=false;
    bool skip_simp=false;
    bool skip_rl=false;
    bool part=false;
    bool check_gram=false;
    bool semi_external=false;
    size_t seed=0;
    log_lvl verbose_lvl=INFO;
    text_format txt_fmt=UNKNOWN;

    std::string p_file;
    std::vector<std::string> grammars_to_merge;
    std::vector<std::string> position_list;

    std::string coord_file;
    std::string version= "v1.0.1 alpha";
    std::vector<std::string> ra_positions;
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

/*
std::function<std::string(std::string)> unit_scale = [](std::string str) -> std::string {

    std::vector<std::string> pats_mb = {"mb", "MB", "mbs", "MBs", "m", "M"};
    std::vector<std::string> pats_kb = {"kb", "KB", "kbs", "KBs", "k", "K"};
    std::vector<std::string> pats_gb = {"gb", "GB", "gbs", "GBs", "g", "G"};
    std::string ext;

    size_t size, scale=0, val;

    if(ends_with(str, pats_mb, ext)){
        scale = 1048576;
    }else if(ends_with(str, pats_gb, ext)){
        scale = 1073741824;
    } else if(ends_with(str, pats_kb, ext)){
        scale = 1024;
    }

    if(!ext.empty()){
        str = str.substr(0, str.length() - ext.length());
        std::stringstream stream(str);
        stream >> val;
        size = size_t(val*scale);
        str = std::to_string(size);
    }
    return str;
};*/

static void parse_app(CLI::App& app, struct arguments& args){

    auto fmt = std::make_shared<MyFormatter>();

    fmt->column_width(27);
    app.formatter(fmt);
    app.add_flag("-v,--version", args.ver, "Print the software version and exit");

    CLI::App* comp = app.add_subcommand("cmp", "Compress one or more text collections");

    //compression
    comp->add_option("TEXTS", args.input_files, "Input files");
    comp->add_option("-o,--output-file", args.output_file, "Output file")->type_name("");

    auto * opt_fmt1 = comp->add_option("-x,--format", args.txt_fmt, "Format of the input files (1=Plain, 2=Fasta, 3=Fastq)")->default_val(UNKNOWN)->check(CLI::Range(1,3));
    auto * opt_fmt2 = comp->add_option("-X,--format-seq", args.input_fmts, "Sequence of input formats (in case they differ)");
    opt_fmt1->excludes(opt_fmt2);

    comp->add_option("-l,--list", args.file_list, "List of input files")->check(CLI::ExistingFile);

    comp->add_flag("-q,--skip-simp", args.skip_simp, "Do not simplify the grammar");
    comp->add_flag("-e,--skip-rl", args.skip_rl, "Do not perform run-length compression");
    comp->add_flag("-r,--random-support", args.rand_acc, "Augment the grammar with random access support");
    comp->add_flag("-g,--check-gram", args.check_gram, "Check that the grammar was compressed correctly");
    auto * se_flag = comp->add_flag("-s,--semi-external", args.semi_external, "Keep some satellite data on disk to reduce RAM usage");
    auto * tmp_fd_opt = comp->add_option("-T,--tmp", args.tmp_dir, "Temporary folder (def. /tmp/lcg.xxxx)")->check(CLI::ExistingDirectory);
    //comp->add_flag("-p,--partial", args.part, "Build a partial grammar representation");
    tmp_fd_opt->needs(se_flag);

    //comp->add_option("-c,--text-chunks", args.n_chunks, "Number of text chunks in memory during the parsing (def. n_threads+1)")->default_val(0);
    comp->add_option("-t,--threads", args.n_threads, "Maximum number of parsing threads")->default_val(1);
    comp->add_option("-f,--fraction", args.i_frac, "The combined threads will try to use TEXT_SIZE*f bytes of RAM to compress TEXT");
    comp->add_option("-c,--chunk-size", args.chunk_size, "Size in bytes of each text chunk (def. min(TEXT_SIZE*0.005, 200MB))")->default_val(0);

    comp->add_option("-v,--verbose-level", args.verbose_lvl, "Verbose level (0=error, 1=warning, 2=info, 3=debug)")->default_val(INFO)->check(CLI::Range(0,3));

    //metadata
    CLI::App* meta = app.add_subcommand("met", "Get the metadata of a grammar");
    meta->add_option("GRAM", args.i_file, "Input grammar in LCG format")->check(CLI::ExistingFile)->required();

    //merge
    //CLI::App* merge = app.add_subcommand("mrg", "Merge grammars");
    //merge->add_option("GRAM LIST", args.grammars_to_merge, "Grammars to be merged")->check(CLI::ExistingFile)->required();
    //merge->add_option("-o,--output-file", args.output_file, "Output grammar")->required(true);
    //merge->add_option("-T,--tmp", args.tmp_dir, "Temporary folder (def. /tmp/lcg.xxxx)")->check(CLI::ExistingDirectory)->default_val("/tmp");

    //access
    CLI::App* access = app.add_subcommand("acc", "Random access");
    access->add_option("GRAM", args.i_file, "Input grammar")->check(CLI::ExistingFile)->required();
    access->add_option("STR:START-END", args.ra_positions, "Coordinates to be accessed");

    //edit
    //CLI::App* edt = app.add_subcommand("edt", "Edit grammar-compressed text");
    //edt->add_option("GRAM", args.input_file, "Grammar to be edited")->check(CLI::ExistingFile)->required();
    //edt->add_option("STR:START-END", args.ra_positions, "Coordinates to be removed");
    //edt->add_option("-o,--output-file", args.output_file, "Output grammar");
    app.require_subcommand(1,1);
}

template<class gram_type>
void comp_int2(arguments& args){
    switch (args.verbose_lvl) {
        case ERROR:
            build_gram<gram_type, ERROR>(args.input_files, args.input_fmts, args.output_file, args.n_threads,
                                         args.chunk_size, args.i_frac, args.skip_simp,
                                         args.part, args.check_gram, args.semi_external, args.tmp_dir);
            break;
        case WARNING:
            build_gram<gram_type, WARNING>(args.input_files, args.input_fmts, args.output_file, args.n_threads,
                                           args.chunk_size, args.i_frac, args.skip_simp,
                                           args.part, args.check_gram, args.semi_external, args.tmp_dir);
            break;
        case INFO:
            build_gram<gram_type, INFO>(args.input_files, args.input_fmts, args.output_file, args.n_threads,
                                        args.chunk_size, args.i_frac, args.skip_simp,
                                        args.part, args.check_gram, args.semi_external, args.tmp_dir);
            break;
        case DEBUG:
            build_gram<gram_type, DEBUG>(args.input_files, args.input_fmts, args.output_file, args.n_threads,
                                         args.chunk_size, args.i_frac, args.skip_simp,
                                         args.part, args.check_gram, args.semi_external, args.tmp_dir);
            break;
        default:
            logger<ERROR>::error("Unknown log level");
            exit(1);
    }
}

void comp_int(arguments& args) {
    if(args.skip_rl){
        if(args.rand_acc){
            using gram_type = lc_gram_t<false, false, true>;
            comp_int2<gram_type>(args);
        }else{
            using gram_type = lc_gram_t<false, false, false>;
            comp_int2<gram_type>(args);
        }
    }else{
        if(args.rand_acc){
            using gram_type = lc_gram_t<false, true, true>;
            comp_int2<gram_type>(args);
        } else {
            using gram_type = lc_gram_t<false, true, false>;
            comp_int2<gram_type>(args);
        }
    }
}

void access_int(std::string& input_file, std::vector<str_coord_type>& query_coords, bool has_rl_rules, bool has_cg_rules, bool has_rand_access) {
    assert(has_rand_access);
    if(has_cg_rules){
        if(has_rl_rules){
            lc_gram_t<true, true, true> gram;
            load_from_file(input_file, gram);
            for(auto const& query : query_coords){
                std::string dc_output;
                gram.im_str_rand_access(query.str, query.start, query.end, dc_output);
                std::cout<<query.start<<":"<<query.start<<"-"<<query.end<<std::endl;
                std::cout<<dc_output<<std::endl;
            }
        }else{
            lc_gram_t<true, false, true> gram;
            load_from_file(input_file, gram);
            for(auto const& query : query_coords){
                std::string dc_output;
                gram.im_str_rand_access(query.str, query.start, query.end, dc_output);
                std::cout<<query.str<<":"<<query.start<<"-"<<query.end<<std::endl;
                std::cout<<dc_output<<std::endl;
            }
        }
    }else{
        if(has_rl_rules){
            lc_gram_t<false, true, true> gram;
            load_from_file(input_file, gram);
            for(auto const& query : query_coords){
                std::string dc_output;
                gram.im_str_rand_access(query.str, query.start, query.end, dc_output);
                std::cout<<query.str<<":"<<query.start<<"-"<<query.end<<std::endl;
                std::cout<<dc_output<<std::endl;
            }
        }else{
            lc_gram_t<false, false, true> gram;
            load_from_file(input_file, gram);
            for(auto const& query : query_coords){
                std::string dc_output;
                gram.im_str_rand_access(query.str, query.start, query.end, dc_output);
                std::cout<<query.str<<":"<<query.start<<"-"<<query.end<<std::endl;
                std::cout<<dc_output<<std::endl;
            }
        }
    }
}

void read_file_list(std::vector<std::string>& input_files, const std::string& list){

    std::ifstream file(list);

    if (!file.is_open()) {
        std::cerr << "Failed to open " << list<<std::endl;
        exit(1);
    }

    std::string line;
    // Iterate over the file line by line
    while (std::getline(file, line)) {
        input_files.push_back(line);
    }

    // Close the file when done
    file.close();
}

std::vector<str_coord_type> parse_query_coords(std::vector<std::string>& str_queries){

    size_t str;
    off_t start, end;
    std::vector<str_coord_type> query_coords;
    query_coords.reserve(str_queries.size());

    for(auto const &coord: str_queries){
        std::vector<std::string> tmp = split(coord, ':');
        if(tmp.size()!=2){
            std::cout<<"Coordinate error: query \""<<coord<<"\" is il-formed"<<std::endl;
            exit(1);
        }
        std::vector<std::string> tmp2 = split(tmp[1], '-');
        if(tmp2.size()!=2){
            std::cout<<"Coordinate error: query \""<<coord<<"\" is il-formed"<<std::endl;
            exit(1);
        }

        try{
            str = stoi(tmp[0]);
        } catch (const std::invalid_argument & e) {
            std::cout << "Coordinate error: \""<<tmp[0] <<"\" in \""<<coord<<"\" is not a valid string ID\n";
            exit(1);
        } catch (const std::out_of_range & e) {
            std::cout << "Coordinate error: \""<<tmp[0] <<"\" in \""<<coord<<"\" is not a valid string ID\n";
            exit(1);
        }

        try{
            start = stoi(tmp2[0]);
        } catch (const std::invalid_argument & e) {
            std::cout << "Coordinate error: \""<<tmp2[0] <<"\" in \""<<coord<<"\" is not a valid start\n";
            exit(1);
        } catch (const std::out_of_range & e) {
            std::cout << "Coordinate error: \""<<tmp2[0] <<"\" in \""<<coord<<"\" is not a valid start\n";
            exit(1);
        }

        try{
            end = stoi(tmp2[1]);
        } catch (const std::invalid_argument & e) {
            std::cout << "Coordinate error: \""<<tmp2[1] <<"\" in \""<<coord<<"\" is not a valid end\n";
            exit(1);
        } catch (const std::out_of_range & e) {
            std::cout << "Coordinate error: \""<<tmp2[1] <<"\" in \""<<coord<<"\" is not a valid end\n";
            exit(1);
        }
        query_coords.emplace_back(str, start, end);
    }
    return query_coords;
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

    if(app.got_subcommand("cmp")) {

        //Note: this is an ugly hack to catch situations where no inputs are given.
        // Oddly, supporting "at least one of these" options in CL11 without affecting
        // the help format is not simple.
        auto * txt_positional = app.get_subcommand(0)->get_option("TEXTS");
        auto * txt_file = app.get_subcommand(0)->get_option("--list");
        if(!txt_positional->count() && !txt_file->count()){
            std::cerr<<"TEXTS/--list: at least one of these parameters is needed"<<std::endl;
            std::cerr<<"Run with --help for more information"<<std::endl;
            exit(1);
        }
        //

        if(!args.file_list.empty()){
            read_file_list(args.input_files, args.file_list);
        }
        assert(!args.input_files.empty());

        if(args.output_file.empty()) args.output_file = std::filesystem::path(args.input_files[0]).filename();
        args.output_file = std::filesystem::path(args.output_file).replace_extension(".lcg");

        if(args.input_fmts.empty()){
            args.input_fmts.resize(args.input_files.size());
            for(auto& fmt : args.input_fmts){
                fmt = args.txt_fmt;
            }
        }

        assert(args.input_files.size()==args.input_fmts.size());
        //TODO only for testing
        /*size_t n_bytes = file_size(args.input_files[0]);
        auto *tmp = allocator::allocate<uint8_t>(n_bytes);
        std::ifstream file(args.input_files[0], std::ios::binary); // `ate` opens at end for size calculation
        file.read((char *)tmp, (std::streamsize)n_bytes);

        auto *tmp2 = allocator::allocate<uint8_t>(n_bytes);
        memcpy(tmp2, tmp, n_bytes);

        auto *tmp3 = allocator::allocate<uint8_t>(n_bytes);
        memcpy(tmp3, tmp, n_bytes);

        auto *tmp4 = allocator::allocate<uint8_t>(n_bytes);
        memcpy(tmp4, tmp, n_bytes);

        size_t n_empty;
        auto start = std::chrono::steady_clock::now();
        size_t parsed_bytes_scalar = fasta_parsing_scalar(tmp, n_bytes, n_empty);
        auto end = std::chrono::steady_clock::now();
        std::cout <<"Time scalar = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - start).count() << "[ns] output_size: "<<parsed_bytes_scalar<<" empty entries "<<n_empty<<std::endl;

#ifdef __ARM_NEON__
        start = std::chrono::steady_clock::now();
        size_t parsed_bytes_neon = fasta_parsing_neon(tmp2, n_bytes, n_empty);
        end = std::chrono::steady_clock::now();
        std::cout <<"Time neon = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - start).count() << "[ns] output_size: "<<parsed_bytes_neon<<" empty entries "<<n_empty<<std::endl;
        assert(parsed_bytes_scalar==parsed_bytes_neon);
        std::cout<<"Equals? "<<(memcmp(tmp,tmp2, parsed_bytes_neon)==0)<<std::endl;
#endif

#ifdef __AVX2__
        start = std::chrono::steady_clock::now();
        size_t parsed_bytes_avx = fasta_parsing_avx2(tmp3, n_bytes, n_empty);
        end = std::chrono::steady_clock::now();
        std::cout <<"Time avx2 = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - start).count() << "[ns] output_size: "<<parsed_bytes_avx<<" empty entries "<<n_empty<<std::endl;
        assert(parsed_bytes_scalar==parsed_bytes_avx);
        std::cout<<"Equals? "<<(memcmp(tmp,tmp3, parsed_bytes_avx)==0)<<std::endl;
#endif

#ifdef __SSE4_2__
        start = std::chrono::steady_clock::now();
        size_t parsed_bytes_sse42 = fasta_parsing_sse42(tmp4, n_bytes, n_empty);
        end = std::chrono::steady_clock::now();
        std::cout <<"Time sse42 = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - start).count() << "[ns] output_size: "<<parsed_bytes_avx<<" empty entries "<<n_empty<<std::endl;
        assert(parsed_bytes_scalar==parsed_bytes_sse42);
        std::cout<<"Equals? "<<(memcmp(tmp,tmp4, parsed_bytes_sse42)==0)<<std::endl;
#endif
        free(tmp);
        free(tmp2);
        free(tmp3);
        free(tmp4);*/
        //
        comp_int(args);
    } else if(app.got_subcommand("met")){
        print_metadata(args.i_file);
    } else if (app.got_subcommand("acc")){
        std::vector<str_coord_type> query_coords = parse_query_coords(args.ra_positions);
        bool has_rl_rules, has_cg_rules, has_rand_access;
        std::tie(has_rl_rules, has_cg_rules, has_rand_access) = read_grammar_flags(args.i_file);
        access_int(args.i_file, query_coords, has_rl_rules, has_cg_rules, has_rand_access);
    } else {
        std::cout<<" Unknown command "<<std::endl;
        return 1;
    }
    return 0;
}
