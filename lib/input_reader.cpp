//
// BGZF-aware input reader implementation for LCG
//

#include "build_gram/input_reader.h"

#include <fcntl.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>

#ifdef __linux__
#include <sys/types.h>
#endif

#include <htslib/bgzf.h>

// --- plain_reader ---

plain_reader::plain_reader(const std::string& path) {
    fd_ = open(path.c_str(), O_RDONLY);
    if (fd_ < 0) {
        throw std::runtime_error("plain_reader: cannot open " + path);
    }
    std::ifstream in(path, std::ifstream::ate | std::ifstream::binary);
    file_size_ = in.tellg();
}

plain_reader::~plain_reader() {
    if (fd_ >= 0) close(fd_);
}

ssize_t plain_reader::read_data(void* buf, size_t len) {
    return ::read(fd_, buf, len);
}

size_t plain_reader::uncompressed_size() const {
    return file_size_;
}

void plain_reader::advise_sequential() {
#ifdef __linux__
    posix_fadvise(fd_, 0, file_size_, POSIX_FADV_SEQUENTIAL);
#endif
}

void plain_reader::advise_dontneed(off_t offset, off_t len) {
#ifdef __linux__
    posix_fadvise(fd_, offset, len, POSIX_FADV_DONTNEED);
#else
    (void)offset; (void)len;
#endif
}

void plain_reader::advise_dontneed_all() {
#ifdef __linux__
    posix_fadvise(fd_, 0, file_size_, POSIX_FADV_DONTNEED);
#endif
}

// --- bgzf_reader ---

bgzf_reader::bgzf_reader(const std::string& path) {
    uncomp_size_ = compute_bgzf_uncompressed_size(path);
    BGZF* fp = bgzf_open(path.c_str(), "r");
    if (!fp) {
        throw std::runtime_error("bgzf_reader: cannot open " + path);
    }
    bgzf_handle_ = fp;
}

bgzf_reader::~bgzf_reader() {
    if (bgzf_handle_) {
        bgzf_close(static_cast<BGZF*>(bgzf_handle_));
    }
}

ssize_t bgzf_reader::read_data(void* buf, size_t len) {
    ssize_t ret = bgzf_read(static_cast<BGZF*>(bgzf_handle_), buf, len);
    return ret;
}

size_t bgzf_reader::uncompressed_size() const {
    return uncomp_size_;
}

void bgzf_reader::advise_sequential() {
    // No meaningful fadvise for BGZF
}

void bgzf_reader::advise_dontneed(off_t, off_t) {
    // No meaningful fadvise for BGZF
}

void bgzf_reader::advise_dontneed_all() {
    // No meaningful fadvise for BGZF
}

// --- Detection and sizing utilities ---

bool is_bgzf_file(const std::string& path) {
    uint8_t buf[18];
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) return false;

    ssize_t n = ::read(fd, buf, 18);
    close(fd);

    if (n < 18) return false;

    // Check gzip magic
    if (buf[0] != 0x1f || buf[1] != 0x8b) return false;
    // Check deflate method
    if (buf[2] != 0x08) return false;
    // FLG.FEXTRA must be set (bit 2)
    if (!(buf[3] & 0x04)) return false;
    // XLEN at offset 10 (little-endian) must be >= 6
    uint16_t xlen = buf[10] | (buf[11] << 8);
    if (xlen < 6) return false;
    // Check BGZF extra field: SI1='B', SI2='C' at offset 12-13
    if (buf[12] != 'B' || buf[13] != 'C') return false;

    return true;
}

size_t compute_bgzf_uncompressed_size(const std::string& path) {
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        throw std::runtime_error("compute_bgzf_uncompressed_size: cannot open " + path);
    }

#ifdef __linux__
    // Get file size for fadvise
    off_t fsize = lseek(fd, 0, SEEK_END);
    lseek(fd, 0, SEEK_SET);
    posix_fadvise(fd, 0, fsize, POSIX_FADV_SEQUENTIAL);
#endif

    // Read file in large chunks, parse BGZF block headers from the buffer
    static constexpr size_t READ_BUF_SIZE = 1024 * 1024; // 1MB
    std::vector<uint8_t> read_buf(READ_BUF_SIZE);
    size_t total_uncompressed = 0;

    // Buffered reading state
    size_t buf_pos = 0;   // current position in read_buf
    size_t buf_end = 0;   // valid bytes in read_buf

    auto fill_buffer = [&]() -> bool {
        // Move unconsumed data to front
        if (buf_pos > 0 && buf_pos < buf_end) {
            memmove(read_buf.data(), read_buf.data() + buf_pos, buf_end - buf_pos);
            buf_end -= buf_pos;
            buf_pos = 0;
        } else if (buf_pos >= buf_end) {
            buf_pos = 0;
            buf_end = 0;
        }
        // Read more data
        ssize_t n = ::read(fd, read_buf.data() + buf_end, READ_BUF_SIZE - buf_end);
        if (n <= 0) return false;
        buf_end += n;
        return true;
    };

    // Ensure at least `needed` bytes are available starting at buf_pos
    auto ensure_bytes = [&](size_t needed) -> bool {
        while ((buf_end - buf_pos) < needed) {
            if (!fill_buffer()) return false;
        }
        return true;
    };

    while (true) {
        // Need at least 18 bytes for the BGZF block header
        if (!ensure_bytes(18)) break;

        uint8_t* hdr = read_buf.data() + buf_pos;

        // Validate gzip magic + BGZF signature
        if (hdr[0] != 0x1f || hdr[1] != 0x8b || hdr[2] != 0x08) break;
        if (!(hdr[3] & 0x04)) break;

        // BSIZE is at offset 16-17 in the header (within the BC extra field)
        // hdr[10..11] = XLEN, hdr[12]='B', hdr[13]='C', hdr[14..15]=SLEN(2), hdr[16..17]=BSIZE
        uint16_t bsize = hdr[16] | (hdr[17] << 8);
        size_t block_size = (size_t)bsize + 1; // total block size

        // The ISIZE (uncompressed size of this block) is the last 4 bytes of the block
        // We need to read up to offset (block_size - 4) from block start to get ISIZE
        if (!ensure_bytes(block_size)) break;

        uint8_t* block = read_buf.data() + buf_pos;
        uint32_t isize = (uint32_t)block[block_size - 4]
                       | ((uint32_t)block[block_size - 3] << 8)
                       | ((uint32_t)block[block_size - 2] << 16)
                       | ((uint32_t)block[block_size - 1] << 24);

        total_uncompressed += isize;
        buf_pos += block_size;
    }

    close(fd);
    return total_uncompressed;
}

size_t input_file_size(const std::string& path) {
    if (is_bgzf_file(path)) {
        return compute_bgzf_uncompressed_size(path);
    }
    std::ifstream in(path, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

std::unique_ptr<input_reader> make_input_reader(const std::string& path) {
    if (is_bgzf_file(path)) {
        return std::make_unique<bgzf_reader>(path);
    }
    return std::make_unique<plain_reader>(path);
}
