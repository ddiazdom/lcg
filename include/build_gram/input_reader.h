//
// BGZF-aware input reader abstraction for LCG
//

#ifndef LCG_INPUT_READER_H
#define LCG_INPUT_READER_H

#include <cstdint>
#include <cstddef>
#include <string>
#include <memory>
#include <sys/types.h>

class input_reader {
public:
    virtual ~input_reader() = default;

    // Read up to len bytes into buf. Returns number of bytes read (0 at EOF).
    virtual ssize_t read_data(void* buf, size_t len) = 0;

    // Return the uncompressed file size.
    virtual size_t uncompressed_size() const = 0;

    // Hint to OS for sequential access (no-op for BGZF).
    virtual void advise_sequential() = 0;

    // Hint to OS to drop pages from cache (no-op for BGZF).
    virtual void advise_dontneed(off_t offset, off_t len) = 0;

    // Drop all pages from cache (no-op for BGZF).
    virtual void advise_dontneed_all() = 0;
};

class plain_reader : public input_reader {
    int fd_;
    size_t file_size_;
public:
    explicit plain_reader(const std::string& path);
    ~plain_reader() override;

    ssize_t read_data(void* buf, size_t len) override;
    size_t uncompressed_size() const override;
    void advise_sequential() override;
    void advise_dontneed(off_t offset, off_t len) override;
    void advise_dontneed_all() override;
};

class bgzf_reader : public input_reader {
    void* bgzf_handle_; // BGZF* -- opaque to avoid exposing htslib headers
    size_t uncomp_size_;
public:
    explicit bgzf_reader(const std::string& path);
    ~bgzf_reader() override;

    ssize_t read_data(void* buf, size_t len) override;
    size_t uncompressed_size() const override;
    void advise_sequential() override;
    void advise_dontneed(off_t offset, off_t len) override;
    void advise_dontneed_all() override;
};

// Detect whether a file is BGZF-compressed by checking magic bytes.
bool is_bgzf_file(const std::string& path);

// Compute uncompressed size of a BGZF file by scanning block headers.
size_t compute_bgzf_uncompressed_size(const std::string& path);

// Return uncompressed file size regardless of format (plain or BGZF).
size_t input_file_size(const std::string& path);

// Factory: auto-detect format and return the appropriate reader.
std::unique_ptr<input_reader> make_input_reader(const std::string& path);

#endif //LCG_INPUT_READER_H
