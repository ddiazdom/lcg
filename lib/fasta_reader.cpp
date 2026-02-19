//
// FASTA file list reader implementation
//

#include "build_gram/fasta_reader.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_set>
#include <algorithm>
#include <cstring>

// --- parse_fasta_file_list ---

std::vector<std::string> parse_fasta_file_list(const std::string& list_path) {
    std::ifstream ifs(list_path);
    if (!ifs.is_open()) {
        throw std::runtime_error("Cannot open file list: " + list_path);
    }

    std::vector<std::string> paths;
    std::string line;
    while (std::getline(ifs, line)) {
        // strip trailing \r
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        // strip leading/trailing whitespace
        size_t start = line.find_first_not_of(" \t");
        if (start == std::string::npos) continue; // blank line
        size_t end = line.find_last_not_of(" \t");
        line = line.substr(start, end - start + 1);
        // skip comments
        if (line.empty() || line[0] == '#') continue;
        paths.push_back(line);
    }
    return paths;
}

// --- fasta_reader Pass 1: scan ---

// Helper: scan a single FASTA file to collect sequence names and data bytes.
// Uses the input_reader interface so it handles both plain and BGZF files.
static void scan_fasta_file(const std::string& path,
                            std::vector<std::string>& seq_names,
                            std::unordered_set<std::string>& seen_names,
                            size_t& total_data_bytes) {
    auto reader = make_input_reader(path);

    static constexpr size_t SCAN_BUF_SIZE = 4 * 1024 * 1024;
    std::vector<uint8_t> buf(SCAN_BUF_SIZE);
    size_t buf_pos = 0, buf_end = 0;

    bool in_header = false;
    bool first_line = true;
    bool seen_any_header = false;
    std::string current_header;
    size_t current_seq_bytes = 0;

    auto fill = [&]() -> bool {
        ssize_t n = reader->read_data(buf.data(), SCAN_BUF_SIZE);
        if (n <= 0) return false;
        buf_pos = 0;
        buf_end = (size_t)n;
        return true;
    };

    // Process byte by byte through the buffer (fast inner loop)
    while (true) {
        if (buf_pos >= buf_end) {
            if (!fill()) break;
        }

        // Find end of current line
        size_t line_start = buf_pos;
        while (buf_pos < buf_end && buf[buf_pos] != '\n') {
            buf_pos++;
        }

        if (buf_pos < buf_end) {
            // Found '\n' — we have a complete line segment
            size_t line_end = buf_pos;
            buf_pos++; // skip '\n'

            // strip \r
            if (line_end > line_start && buf[line_end - 1] == '\r') {
                line_end--;
            }

            if (in_header) {
                // Continuation of a header line that was split across buffers
                // (this shouldn't happen often since headers are short)
                for (size_t i = line_start; i < line_end; i++) {
                    current_header.push_back((char)buf[i]);
                }
                in_header = false;

                // Parse the name from the header
                std::string name;
                size_t name_start = 0;
                // skip leading whitespace after '>'
                while (name_start < current_header.size() &&
                       (current_header[name_start] == ' ' || current_header[name_start] == '\t')) {
                    name_start++;
                }
                size_t name_end = name_start;
                while (name_end < current_header.size() &&
                       current_header[name_end] != ' ' && current_header[name_end] != '\t') {
                    name_end++;
                }
                name = current_header.substr(name_start, name_end - name_start);

                if (name.empty()) {
                    throw std::runtime_error("Empty sequence name in FASTA header in " + path);
                }
                if (seen_names.count(name)) {
                    throw std::runtime_error("Duplicate sequence name '" + name + "' in " + path);
                }
                seen_names.insert(name);
                seq_names.push_back(name);
                seen_any_header = true;
                current_header.clear();
                continue;
            }

            // Start of a new line
            if (line_start < line_end && buf[line_start] == '>') {
                // This is a header line
                if (!first_line && seen_any_header) {
                    // End of previous sequence
                    if (current_seq_bytes == 0) {
                        std::cerr << "Warning: empty sequence '" << seq_names.back() << "' in " << path << std::endl;
                    }
                    total_data_bytes += current_seq_bytes + 1; // +1 for '\n' separator
                    current_seq_bytes = 0;
                }
                first_line = false;

                // Parse header: extract first whitespace-delimited token after '>'
                current_header.clear();
                for (size_t i = line_start + 1; i < line_end; i++) {
                    current_header.push_back((char)buf[i]);
                }

                std::string name;
                size_t name_start = 0;
                while (name_start < current_header.size() &&
                       (current_header[name_start] == ' ' || current_header[name_start] == '\t')) {
                    name_start++;
                }
                size_t name_end = name_start;
                while (name_end < current_header.size() &&
                       current_header[name_end] != ' ' && current_header[name_end] != '\t') {
                    name_end++;
                }
                name = current_header.substr(name_start, name_end - name_start);

                if (name.empty()) {
                    throw std::runtime_error("Empty sequence name in FASTA header in " + path);
                }
                if (seen_names.count(name)) {
                    throw std::runtime_error("Duplicate sequence name '" + name + "' in " + path);
                }
                seen_names.insert(name);
                seq_names.push_back(name);
                seen_any_header = true;
                current_header.clear();
            } else {
                // Sequence data line
                if (!seen_any_header && line_start < line_end) {
                    throw std::runtime_error("File does not appear to be FASTA (data before first header): " + path);
                }
                current_seq_bytes += (line_end - line_start);
            }
            first_line = false;
        } else {
            // Reached end of buffer without finding '\n' — partial line
            if (line_start < buf_end) {
                if (!in_header && buf[line_start] == '>') {
                    // Start of a header that's split across buffers
                    if (!first_line && seen_any_header) {
                        if (current_seq_bytes == 0) {
                            std::cerr << "Warning: empty sequence '" << seq_names.back() << "' in " << path << std::endl;
                        }
                        total_data_bytes += current_seq_bytes + 1;
                        current_seq_bytes = 0;
                    }
                    first_line = false;
                    in_header = true;
                    current_header.clear();
                    for (size_t i = line_start + 1; i < buf_end; i++) {
                        current_header.push_back((char)buf[i]);
                    }
                } else if (in_header) {
                    for (size_t i = line_start; i < buf_end; i++) {
                        current_header.push_back((char)buf[i]);
                    }
                } else {
                    // Partial sequence data line
                    if (!seen_any_header) {
                        throw std::runtime_error("File does not appear to be FASTA (data before first header): " + path);
                    }
                    current_seq_bytes += (buf_end - line_start);
                }
            }
            // Need more data
        }
    }

    // Handle last sequence
    if (seen_any_header) {
        if (current_seq_bytes == 0) {
            std::cerr << "Warning: empty sequence '" << seq_names.back() << "' in " << path << std::endl;
        }
        total_data_bytes += current_seq_bytes + 1; // +1 for '\n' separator
    }

    if (!seen_any_header) {
        throw std::runtime_error("File does not appear to be FASTA (no headers found): " + path);
    }
}

fasta_reader::fasta_reader(std::vector<std::string> file_paths)
    : file_paths_(std::move(file_paths)) {

    if (file_paths_.empty()) {
        throw std::runtime_error("fasta_reader: empty file list");
    }

    // Pass 1: scan all files
    std::unordered_set<std::string> seen_names;
    for (const auto& path : file_paths_) {
        scan_fasta_file(path, seq_names_, seen_names, total_data_bytes_);
    }

    // Allocate raw buffer for Pass 2
    raw_buf_.resize(RAW_BUF_SIZE);
}

// --- fasta_reader Pass 2: read_data ---

size_t fasta_reader::uncompressed_size() const {
    return total_data_bytes_;
}

void fasta_reader::advise_sequential() {}
void fasta_reader::advise_dontneed(off_t, off_t) {}
void fasta_reader::advise_dontneed_all() {}

bool fasta_reader::open_next_file() {
    if (current_file_idx_ >= file_paths_.size()) {
        at_eof_ = true;
        return false;
    }
    current_reader_ = make_input_reader(file_paths_[current_file_idx_++]);
    raw_pos_ = 0;
    raw_end_ = 0;
    in_header_ = false;
    return true;
}

bool fasta_reader::fill_raw_buf() {
    if (!current_reader_) return false;
    ssize_t n = current_reader_->read_data(raw_buf_.data(), RAW_BUF_SIZE);
    if (n <= 0) return false;
    raw_pos_ = 0;
    raw_end_ = (size_t)n;
    return true;
}

ssize_t fasta_reader::read_data(void* buf, size_t len) {
    if (at_eof_ || len == 0) return 0;

    uint8_t* out = static_cast<uint8_t*>(buf);
    size_t written = 0;

    while (written < len) {
        // Ensure we have a reader open
        if (!current_reader_) {
            if (!open_next_file()) break;
        }

        // Ensure we have data in raw buffer
        if (raw_pos_ >= raw_end_) {
            if (!fill_raw_buf()) {
                // Current file exhausted. Emit final '\n' for last sequence if needed.
                // (We'll handle this via need_separator_ on next file open)
                current_reader_.reset();
                continue;
            }
        }

        // Process raw buffer bytes → output buffer
        while (written < len && raw_pos_ < raw_end_) {
            uint8_t ch = raw_buf_[raw_pos_];

            if (in_header_) {
                // Skip until end of header line
                if (ch == '\n') {
                    in_header_ = false;
                }
                raw_pos_++;
                continue;
            }

            if (ch == '>') {
                // Start of a new header
                // If we have a previous sequence, emit separator
                if (need_separator_) {
                    out[written++] = '\n';
                    need_separator_ = false;
                }
                in_header_ = true;
                // After this header, the next sequence data should set need_separator_
                raw_pos_++;
                continue;
            }

            if (ch == '\n' || ch == '\r') {
                // Skip line wrapping within sequence data
                raw_pos_++;
                continue;
            }

            // Sequence data character
            // Before first char of a sequence, mark that we'll need a separator later
            if (!need_separator_ && !first_seq_in_stream_) {
                // We just finished a header, about to start data
            }
            first_seq_in_stream_ = false;
            need_separator_ = true;

            out[written++] = ch;
            raw_pos_++;
        }
    }

    // If we've exhausted all files, emit final '\n'
    if (written < len && !current_reader_ && current_file_idx_ >= file_paths_.size() && need_separator_) {
        out[written++] = '\n';
        need_separator_ = false;
        at_eof_ = true;
    }

    return (ssize_t)written;
}
