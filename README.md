# LCG : Locally-consistent grammar compression 

This repository contains the implementation of **LCG**, a scalable text compressor that relies on the concept of
locally-consistent grammars. The main features of **LCG** are:

1. Compression of TBs of text efficiently 

# Third-party libraries

1. [xxHash](https://github.com/Cyan4973/xxHash)
2. [CLI11](https://github.com/CLIUtils/CLI11)

# Prerequisites

1. C++ >= 17
2. CMake >= 3.7

The xxHash and CLI11 libraries are already included in the source files of this repository.

## DISCLAIMER

This implementation is still under development. Not all the features described in the help have been tested, and
some of them are partially implemented.

For the moment, use LCG only to measure compression ratios.

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
the same separator symbol. The tool assumes the last symbol in the file is the separator.

For collections of ASCII characters (i.e, regular text, DNA, protein, etc), inputs in one-string-per-line format should
work just fine. 

## Long and short strings 

The compression algorithm of **LCG** is optimized to work with collections of strings that do not exceed the 4 GBs in length.
This cap in enough for most practical applications. However, if your collection contains strings longer than 
that value, you can pass the flag ``-l/--long-strings``.

**NOTE**: The current implementation only works in "long-string mode". Please, always use the ``-l`` flag regardless
if the strings are short. The output should be the same.


# Bugs

This tool still contains experimental code, so you will probably find bugs. Please report them in this repository.

# Author

This implementation was written by [Diego Díaz](https://github.com/ddiazdom) .
