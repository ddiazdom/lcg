# LCG : Locally-consistent grammar compression 

This repository contains the implementation of **LCG**, a scalable text compressor that relies on the concept of
locally-consistent grammars. The main features of **LCG** are:

1. Compression of TBs of text efficiently 
2. Random access to the compressed strings  
3. Insertion or deletion of symbols in the compressed strings 
4. Merge multiple compressed files into one single compressed file 

# Third-party libraries

1. [xxHash](https://github.com/Cyan4973/xxHash)
2. [CLI11](https://github.com/CLIUtils/CLI11)

# Prerequisites

1. C++ >= 17
2. CMake >= 3.7

The xxHash and CLI11 libraries are already included in the source files of this repository.

# Installation

Clone repository, enter the project folder and execute the following commands:

```
mkdir build
cd build
cmake ..
make
```

# Compressing text 

```
./lcg comp sample_file.txt
```

## Input 

Our tool currently assumes the input is a concatenated collection of one or more strings, where every string ends with
the same separator symbol. The tools assumes the last symbol in the file is the separator.

For collections of ASCII characters (i.e, regular text, DNA, protein, etc), inputs in one-string-per-line format should
work just fine. 

### Integer alphabets

LCG assumes by default that the input has a byte alphabet (<=256 symbols). However, it is also possible to compress
collections with an arbitrarily large integer alphabet (up to alphabets of 4 bytes). You just need to pass the
``-a/--alphabet`` parameter. 

In this case, the program expects the collection to be encoded in a serialized integer vector (the restrictions
on the separator symbol remain the same). The value for ``-a`` tells the program the native integer type
the input uses. These are all the possible values:

* 1 : 1-byte cells (uint8_t)
* 2 : 2-byte cells (uint16_t)
* 4 : 4-byte cells (uint32_t)

For instance, if your collection has an alphabet of *at most* $2^{16} = 65536$ symbols, then store it using 2-byte cells
(i.e., short integers) and run LCG as:

```
./lcg comp sample_file.txt -a 2
```

**Disclaimer**: the performance of **LCG** has not been extensively tested in inputs with integer alphabets. From a practical point of view,
the algorithms makes non difference between byte and integer alphabets, but bugs can always be present.

## Long and short strings 

The compression algorithm of **LCG** is optimized to work with collections of strings that do not exceed the 4 GBs in length.
This cap in enough for most practical applications. However, if your collection contains strings longer than 
that value, you can pass the flag ``-b/--long-strings``.

## Determinism

The compression algorithm is also randomized. Specifically, different runs with the same input yield different grammars.
However, you can tell **LCG** to build always the same output file for a given input using the flag ``--d/--deterministic``

```
./lcg comp sample_file.txt -d  
```

It is also possible to extract the hash functions that were used to construct a specific compressed file.

## Output

Our tool outputs the BCR BWT as a run-length compressed array. Concretely, it produces a sequence of equal-symbol runs
encoded as pairs $(a,\ell)$, where $a$ is the run symbol and $\ell$ is its length. We represent the run symbols and the run
lengths using cells of different widths to reduce space usage. The width for the symbols is the smallest number of bytes
that fits the alphabet. On the other hand, the width for the lengths is the smallest number of bytes that fits the
length of the longest equal-symbol run in the BCR BWT.

# Decompression and random access

# Performing text updates 

# Merging compressed files 

# Bugs

This tool still contains experimental code, so you will probably find bugs. Please report them in this repository.

# Author

This implementation was written by [Diego Díaz](https://github.com/ddiazdom) .
