//
// Created by Diaz, Diego on 13.2.2023.
//

#ifndef GRLBWT_MXM3_H
#define GRLBWT_MXM3_H
#pragma once

// author: Jon Maiga, 2020-08-03, jonkagstrom.com, @jonkagstrom
// license: CC0 license
//
// v1.0.0, 2020-08-13: https://github.com/jonmaiga/mx3/releases/tag/v1.0.0
// v2.0.0, 2020-10-22: https://github.com/jonmaiga/mx3/releases/tag/v2.0.0
// v3.0.0, 2022-04-19: https://github.com/jonmaiga/mx3/releases/tag/v3.0.0

#include <stddef.h>
#include <stdint.h>

namespace mx3 {

    static const uint64_t C = 0xbea225f9eb34556d;

    inline uint64_t mix(uint64_t x) {
        x ^= x >> 32;
        x *= C;
        x ^= x >> 29;
        x *= C;
        x ^= x >> 32;
        x *= C;
        x ^= x >> 29;
        return x;
    }

    class random {
    public:
        explicit random(uint64_t seed) : _counter(mix(seed + C)) {}
        uint64_t operator()() { return mix(_counter++); }
    private:
        uint64_t _counter;
    };

    inline uint64_t mix_stream(uint64_t h, uint64_t x) {
        x *= C;
        x ^= x >> 39;
        h += x * C;
        h *= C;
        return h;
    }

    inline uint64_t mix_stream(uint64_t h, uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
        a *= C;
        b *= C;
        c *= C;
        d *= C;
        a ^= a >> 39;
        b ^= b >> 39;
        c ^= c >> 39;
        d ^= d >> 39;
        h += a * C;
        h *= C;
        h += b * C;
        h *= C;
        h += c * C;
        h *= C;
        h += d * C;
        h *= C;
        return h;
    }

    inline uint64_t hash(const uint8_t* buf, size_t len, uint64_t seed) {
        const auto* buf64 = reinterpret_cast<const uint64_t*>(buf);
        uint64_t h = mix_stream(seed, len + 1);
        while (len >= 64) {
            len -= 64;
            h = mix_stream(h, buf64[0], buf64[1], buf64[2], buf64[3]);
            h = mix_stream(h, buf64[4], buf64[5], buf64[6], buf64[7]);
            buf64 += 8;
        }

        while (len >= 8) {
            len -= 8;
            h = mix_stream(h, buf64[0]);
            buf64 += 1;
        }

        const auto* const tail8 = reinterpret_cast<const uint8_t*>(buf64);
        const auto* const tail16 = reinterpret_cast<const uint16_t*>(buf64);
        const auto* const tail32 = reinterpret_cast<const uint32_t*>(buf64);
        switch (len) {
            case 0: return mix(h);
            case 1: return mix(mix_stream(h, tail8[0]));
            case 2: return mix(mix_stream(h, tail16[0]));
            case 3: return mix(mix_stream(h, tail16[0] | static_cast<uint64_t>(tail8[2]) << 16));
            case 4: return mix(mix_stream(h, tail32[0]));
            case 5: return mix(mix_stream(h, tail32[0] | static_cast<uint64_t>(tail8[4]) << 32));
            case 6: return mix(mix_stream(h, tail32[0] | static_cast<uint64_t>(tail16[2]) << 32));
            case 7: return mix(mix_stream(h, tail32[0] | static_cast<uint64_t>(tail16[2]) << 32 | static_cast<uint64_t>(tail8[6]) << 48));
            default: ;
        }
        return mix(h);
    }

    inline uint64_t hash_uint64(uint64_t* buf, size_t len, uint64_t seed) {

        uint64_t h = mix_stream(seed, len + 1);

        while (len >= 8) {
            len -= 8;
            h = mix_stream(h, buf[0], buf[1], buf[2], buf[3]);
            h = mix_stream(h, buf[4], buf[5], buf[6], buf[7]);
            buf += 8;
        }

        while (len--> 0) {
            h = mix_stream(h, buf[0]);
            buf += 1;
        }
        return mix(h);
    }

}

#endif //GRLBWT_MXM3_H
