## grlBWT : Efficient construction of the BWT using string compression

This repository contains the implementation of **grlBWT**, a linear-time semi-external algorithm for building the BCR BWT
of a string collection. The novelty of grlBWT is that it resorts to compression to maintain the intermediate data of
the BWT construction in compact form. The compact format used by grlBWT reduces working memory and avoids redundant
computations, thus decreasing computing time too.

## Third-party libraries

1. [SDSL-lite](https://github.com/simongog/sdsl-lite)
2. [xxHash](https://github.com/Cyan4973/xxHash)
3. [CLI11](https://github.com/CLIUtils/CLI11)

## Prerequisites

1. C++ >= 17
2. CMake >= 3.7
3. SDSL-lite

The xxHash and CLI11 libraries are already included in the source files of this repository. The SDSL-lite has to be installed beforehand.
However, we include a CMake module that will search for its local
installation. No need to indicate the path during the compilation.

## Installation

Clone repository, enter the project folder and execute the following commands:

```
mkdir build
cd build
cmake ..
make
```

## Example execution

```
./grlbwt sample_file.txt
```

## Input 

Our tool currently assumes the input is a concatenated collection of one or more strings, where every string ends with
the same separator symbol. The only restriction on the separator is that it has be the smallest symbol in the input's
alphabet.

For collections of ASCII characters (i.e, regular text, DNA, protein, etc), inputs in one-string-per-line format should
work just fine. 

Notice that, when the collection has only one string, grlBWT produces the standard BWT. Still, one can produce the standard
BWT of a multi-string collection by tweaking the input as follows: concatenate all the strings separated 
by a special symbol (any). Then append another special symbol at the end of the resulting string, which in this case, has to be
the smallest in the alphabet.

### Inputs with integer alphabets

grlBWT assumes by default that the input has a byte alphabet (<256 symbols). However, it is also possible to build the
BCR BWT of a collection that has an arbitrarily large integer alphabet by using the ``-a/--alphabet`` parameter. 

In this case, the program expects the collection to be encoded in a serialized integer vector (the restrictions
on the separator symbol remain the same). The value for ``-a`` tells the program the native integer type
the input uses. These are all the possible values:

* 1 : 1-byte cells (uint8_t)
* 2 : 2-byte cells (uint16_t)
* 4 : 4-byte cells (uint32_t)
* 8 : 8-byte cells (uint64_t)

For instance, if your collection has an alphabet of *at most* $2^{16} = 65536$ symbols, then store it using 2-byte cells
(i.e., short integers) and run grlBWT as:

```
./grlbwt sample_file.txt -a 2
```

**Disclaimer**: the performance of grlBWT has not been extensively tested in inputs with integer alphabets. From a practical point of view,
the algorithms makes non difference between byte and integer alphabets, but bugs can always be present.

## Output

Our tool outputs the BCR BWT as a run-length compressed array. Concretely, it produces a sequence of equal-symbol runs
encoded as pairs $(a,\ell)$, where $a$ is the run symbol and $\ell$ is its length. We represent the run symbols and the run
lengths using cells of different widths to reduce space usage. The width for the symbols is the smallest number of bytes
that fits the alphabet. On the other hand, the width for the lengths is the smallest number of bytes that fits the
length of the longest equal-symbol run in the BCR BWT.

## Reversing the BCR BWT to the original text (only byte alphabet)

If you want to print the original strings, please use the ``reverse_bwt`` binary:

```
./reverse_bwt sample_file.txt.rl_bwt reversed_strings.txt 10
```

This command will store the first 10 strings of the original input in the file ``reversed_strings``. If you don't specify
a number of strings, the script will decompress all the strings.

## Store the output BCR BWT in other formats (only byte alphabet)  

The output of our tool is in a custom run-length-encoded format. However, you can easily transform it to
plain format.  

```
./grl2plain sample_file.txt.rl_bwt sample_file_bwt_plain.txt 
```

It is also possible to have the BWT in standard run-length-encoding:

```
./grl2rle sample_file.txt.rl_bwt sample_file_bwt_plain.rle 
```

The file `sample_file_bwt_plain.rle` is a sequence of pairs $(a_1, \ell_1), \ldots , (a_r, \ell_r)$ that encodes the 
$r$ runs of the BCR BWT. Each pair $(a_i, \ell_i)$ uses 16 bytes, 8 for $a_i$ and 8 for $\ell_i$. 

## Multithreading

Our algorithm has a hashing step that supports multithreading to some extent. If you want to use it, pass the parameter
`-t` to the execution of grlBWT.

## Controlling the space usage of the runs

Our algorithm produces partial BCR BWTs so that each new BWT induces the next one. A relevant issue in this scheme is how to
encode the runs of those BWTs. Our default strategy is as follows: if a BWT has an alphabet of length sigma and the longest
run has l symbols, then we use ceil(log(sigma))+ ceil(log(2)) bits for every run. However, in most of the cases, this
value is an overestimation. In non-repetitive instances,
over ~97% of the runs are shorter than 256 (1 byte). This percentage might decrease to 80-60% in highly repetitive datasets,
but ~95% of the run lengths still fit two bytes. In contrast, the number of bytes to encode l can easily exceed the 4 bytes
in any kind of dataset.

We now have an alternative strategy to control the space usage of the runs. You can pass a fix value x (1-4) by the command
line (`-b` parameter) to tell our method that it can't use more than x bytes to encode run lengths. The algorithm will keep the runs with overflow
in a separate data structure. In most of the cases, choosing 1 or 2 should greatly reduce the space usage
during the induction of the BWTs, with almost not impact in the performance.


The general rule is simple: 

- non-repetitive input = `-b 1`
- highly repetitive input = `-b 2`
- insanely repetitive and massive input = `-b 3` (we suppose, but we haven't found an input for this)
- you don't know if the input is repetitive = don't use `-b`. The algorithm will use the default strategy

Example:

```
./grlbwt ../tests/sample_file.txt -b 1
```

If you chose a bad number, don't worry, it won't crash. It might run a bit slower or use more RAM than
it actually needed.

## How to cite

Please use the following bibtex entry to cite this repository:

```tex
@inproceedings{ddgncpm22,
  author ={D{\'\i}az-Dom{\'\i}nguez, Diego and Navarro, Gonzalo},
  title ={Efficient Construction of the {BWT} for Repetitive Text Using String Compression},
  booktitle ={33rd Annual Symposium on Combinatorial Pattern Matching (CPM 2022)},
  pages ={article 29},
  year ={2022},
  volume ={223}
}
```

We hope to update this cite to a more recent version of our work soon.

## Bugs

This tool still contains experimental code, so you will probably find bugs. Please report them in this repository.

## Author

This implementation was written by [Diego Díaz](https://github.com/ddiazdom) .
