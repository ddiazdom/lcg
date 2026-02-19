//
// FASTA file list reader for pangenome-scale workflows
//

#ifndef LCG_FASTA_READER_H
#define LCG_FASTA_READER_H

#include "input_reader.h"
#include <vector>
#include <string>
#include <memory>

// Parse a file list: one FASTA path per line, # comments, blank lines skipped.
std::vector<std::string> parse_fasta_file_list(const std::string& list_path);

// fasta_reader streams unwrapped sequence data from a list of FASTA files,
// emitting sequences separated by '\n'. Implements the input_reader interface
// so the downstream chunked parser sees the same byte format as plain text.
//
// Pass 1 (constructor): scans all files to collect sequence names and compute
// total uncompressed data size.
// Pass 2 (read_data): re-opens files sequentially, skips '>' header lines,
// strips internal line wrapping, emits '\n' after each complete sequence.
class fasta_reader : public input_reader {
public:
    explicit fasta_reader(std::vector<std::string> file_paths);
    ~fasta_reader() override = default;

    ssize_t read_data(void* buf, size_t len) override;
    size_t uncompressed_size() const override;

    void advise_sequential() override;
    void advise_dontneed(off_t offset, off_t len) override;
    void advise_dontneed_all() override;

    const std::vector<std::string>& sequence_names() const { return seq_names_; }

private:
    std::vector<std::string> file_paths_;
    std::vector<std::string> seq_names_;
    size_t total_data_bytes_ = 0;

    // read_data state
    size_t current_file_idx_ = 0;
    std::unique_ptr<input_reader> current_reader_;
    bool in_header_ = false;
    bool need_separator_ = false; // emit '\n' before next sequence data
    bool first_seq_in_stream_ = true; // first sequence overall, no leading '\n'
    bool at_eof_ = false;

    // internal raw buffer for reading from current_reader_
    static constexpr size_t RAW_BUF_SIZE = 8 * 1024 * 1024; // 8MB
    std::vector<uint8_t> raw_buf_;
    size_t raw_pos_ = 0;
    size_t raw_end_ = 0;

    bool open_next_file();
    bool fill_raw_buf();
};

#endif //LCG_FASTA_READER_H
