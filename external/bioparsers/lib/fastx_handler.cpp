//
// Created by Diaz, Diego on 28.10.2021.
//

#include "fastx_handler.h"

str_collection fastx2plain_format(const std::string& input_file, std::string const& output_file, bool rc, uint8_t sep){

    str_collection str_col;
    gzFile fp = gzopen(input_file.c_str(), "r");
    kseq_t *seq = kseq_init(fp);

    o_file_stream<uint8_t> o_buffer(output_file, 1024*1024, std::ios::out);
    uint8_t sym;
    size_t sym_freqs[255]={0};

    str_col.n_strings = 0;
    while(kseq_read(seq) >= 0) {
        for(size_t i=0;i<seq->seq.l;i++){
            sym = seq->seq.s[i];
            sym_freqs[sym]++;
            o_buffer.push_back(sym);
        }
        sym_freqs[sep]++;
        o_buffer.push_back(sep);
        str_col.n_strings++;

        if(rc){
            for(size_t i=seq->seq.l;i-->0;){
                sym = dna_string::comp[(uint8_t)seq->seq.s[i]];
                if(sym==0){
                    std::cerr<<"The input seems not to be DNA (invalid symbol:"<<seq->seq.s[i]<<")"<<std::endl;
                    exit(1);
                }
                sym_freqs[sym]++;
                o_buffer.push_back(sym);
            }
            sym_freqs[sep]++;
            o_buffer.push_back(sep);
            str_col.n_strings++;
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    o_buffer.close();

    /*for(uint8_t i=0;i<255;i++){
        if(sym_freqs[i]!=0){
            str_col.alphabet.push_back(i);
            str_col.sym_freqs.push_back(sym_freqs[i]);
            str_col.n_syms+=sym_freqs[i];
        }
    }
    str_col.min_sym = str_col.alphabet[0];
    str_col.max_sym = str_col.alphabet.back();*/

    return str_col;
}

void fasta_writer::insert_record(const std::string& identifier, const dna_string& seq) {

    assert(!identifier.empty() && identifier[0]!='>' && identifier.back()!='\n');

    buffer[b_pos++] = '>';
    if(b_pos==1024){
        ofs.write((char *)&buffer, 1024);
        b_pos = 0;
    }
    for(auto const& sym : identifier){
        buffer[b_pos++] = sym;
        if(b_pos==1024){
            ofs.write((char *)&buffer, 1024);
            b_pos = 0;
        }
    }
    buffer[b_pos++] = '\n';
    if(b_pos==1024){
        ofs.write((char *)&buffer, 1024);
        b_pos = 0;
    }

    size_t idx=1;
    char in_sym=0;
    auto it = seq.ibegin();
    auto end = seq.iend();
    while(it!=end){
        in_sym = *it;
        buffer[b_pos++] = in_sym;
        if(b_pos==1024){
            ofs.write((char *)&buffer, 1024);
            b_pos = 0;
        }

        if(wrap_width>0 && (idx % wrap_width)==0){ // wrap the sequence
            in_sym = '\n';
            buffer[b_pos++] = in_sym;
            if(b_pos==1024){
                ofs.write((char *)&buffer, 1024);
                b_pos = 0;
            }
        }
        idx++;
        ++it;
    }

    if(in_sym!='\n'){
        buffer[b_pos++] = '\n';
        if(b_pos==1024){
            ofs.write((char *)&buffer, 1024);
            b_pos = 0;
        }
    }
}

hp_collection fastx2_hp_plain_format(const std::string& input_file, std::string const &output_file) {

    bool is_fasta = false, is_fastq = false, is_gz = check_gzip(input_file);
    char f_sym;
    hp_collection hp_coll;

    //read the first symbol in the file to check if it is in FASTA or FASTQ format
    if (is_gz) {
        gzFile zfile = gzopen(input_file.c_str(), "rb");
        gzread(zfile, &f_sym, 1);
        gzclose(zfile);
    } else {
        std::ifstream ifs(input_file, std::ios::binary);
        ifs.read(&f_sym, 1);
        ifs.close();
    }

    if (f_sym == '>') {
        is_fasta = true;
    } else if (f_sym == '@') {
        is_fastq = true;
    }

    if (is_fasta) {
        if (is_gz) {
            fasta_parser<true> fp(input_file, output_file);
            return fp();
        } else {
            fasta_parser<false> fp(input_file, output_file);
            return fp();
        }
    } else if (is_fastq) {
        if (is_gz) {
            fastq_parser<true> fp(input_file, output_file);
            return fp();
        } else {
            fastq_parser<false> fp(input_file, output_file);
            return fp();
        }
    } else {
        std::cout << "Error trying to infer the input" << std::endl;
    }
    return hp_coll;
}

void get_hp_DNA_rc(std::string const& r_heads_file, std::string const& r_len_file, str_collection& str_coll){
    std::ifstream ifs;
    std::ofstream ofs;
    std::streampos buff_size = 8192, read_bytes;

    char *in_h_buffer, *out_h_buffer;
    in_h_buffer = (char *)calloc(buff_size, 1);
    out_h_buffer = (char *)calloc(buff_size, 1);

    //process the run heads
    ifs.open(r_heads_file, std::ios::binary | std::ios::in);
    ofs.open(r_heads_file, std::ios::binary | std::ios::out | std::ios::app);
    ifs.seekg(0, std::ios::end);
    ofs.seekp(0, std::ios::end);
    std::streampos n_bytes = ifs.tellg() - (std::streampos)1;
    uint8_t sym;
    while(true){
        read_bytes = std::min<std::streampos>(n_bytes, buff_size);
        ifs.seekg(n_bytes - (std::streampos)read_bytes);
        ifs.read((char *)in_h_buffer, read_bytes);

        for(size_t j=read_bytes, i=0;j-->0;i++){
            sym = in_h_buffer[j];
            if(sym<85){
                out_h_buffer[i] = (char)dna_string::comp[sym];
            }else{
                out_h_buffer[i] = (char) (dna_string::comp[sym-84]+84);
            }
        }
        ofs.write(out_h_buffer, read_bytes);
        n_bytes -= read_bytes;
        if(n_bytes==0) break;
    }
    out_h_buffer[0] = '\n';
    ofs.write(out_h_buffer, 1);
    ofs.close();
    ifs.close();
    free(in_h_buffer);
    free(out_h_buffer);

    //process the run lengths
    uint16_t *in_l_buffer, *out_l_buffer;
    in_l_buffer = (uint16_t *)calloc(buff_size/sizeof(uint16_t), sizeof(uint16_t));
    out_l_buffer = (uint16_t *)calloc(buff_size/sizeof(uint16_t), sizeof(uint16_t));
    ifs.open(r_len_file, std::ios::binary | std::ios::in);
    ofs.open(r_len_file, std::ios::binary | std::ios::out | std::ios::app);
    ifs.seekg(0, std::ios::end);
    ofs.seekp(0, std::ios::end);
    n_bytes = ifs.tellg();
    size_t n_elms;
    while(true){
        read_bytes = std::min<std::streampos>(n_bytes, buff_size);
        ifs.seekg(n_bytes - (std::streampos)read_bytes);
        ifs.read((char *)in_l_buffer, read_bytes);
        n_elms = read_bytes/sizeof(uint16_t);

        for(size_t j=n_elms, i=0;j-->0;i++){
            out_l_buffer[i] = in_l_buffer[j];
        }
        ofs.write((char *)out_l_buffer, read_bytes);
        n_bytes -= read_bytes;
        if(n_bytes==0) break;
    }

    ofs.close();
    ifs.close();
    free(in_l_buffer);
    free(out_l_buffer);

    str_coll.n_strings*=2;
    str_coll.n_syms*=2;
    //size_t inv_alph[255]={0};
    /*for(size_t i=0;i<str_coll.alphabet.size();i++){
        inv_alph[str_coll.alphabet[i]] = i;
    }

    std::vector<size_t> new_freqs(str_coll.alphabet.size(), 0);
    uint8_t rc_sym;
    for(size_t i=0;i<str_coll.alphabet.size();i++){
        sym = str_coll.alphabet[i];
        rc_sym = dna_string::comp[sym];
        if(sym>=85) rc_sym =dna_string::comp[sym-84]+84;
        new_freqs[i] = str_coll.sym_freqs[i] + str_coll.sym_freqs[inv_alph[rc_sym]];
    }
    std::swap(new_freqs, str_coll.sym_freqs);*/
}
