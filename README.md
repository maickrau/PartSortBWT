## PartSortBWT

Low memory BWT computation. Usage:

```
std::string input = ...
std::string output;
PartSortBWT(input, output); // computes BWT, stores in output
PartSortBWT(input, input); // computes BWT, stores in the same string and overwrites
```

The input and output can be the same string. For a real example see [src/main.cpp](src/main.cpp). The input string must end with 0, and all other characters must be between 1-5 representing the alphabet `$NACGT`. Not recommended for strings which have long runs of N's, eg. GRCh38.

### Benchmark

Measured with `/usr/bin/time -v bin/main input_file > /dev/null` using commit [19bb9f4](https://github.com/maickrau/PartSortBWT/commit/19bb9f42b6892d62f7952889ad96fb475c5f35a5) on a laptop.

| Dataset | Size | Time (h:mm:ss) | Memory |
| --- | --- | --- | --- |
| CHM13 chr1 hpc | 175 Mbp | 0:00:44 | 291 Mb |
| CHM13 hpc | 2.1 Gbp | 0:11:35 | 3.6 Gb |
| CHM13 | 3.1 Gbp | 0:30:03 | 4.9 Gb |
| Assembly graph hpc | 3.8 Gbp | 0:37:26 | 6.5 Gb |

### Method

Classify suffixes according to their first 4 characters. First collect and sort all suffixes starting with AAAA to get the first x characters of the BWT, then suffixes with AAAC to get the next y characters, etc. Repeat for all 4-character prefixes.

Expected memory use for a random DNA string is `input + output + ~0.875n` bytes. If input == output then total memory use is ~1.875n bytes. Less random strings have higher memory use and worst case is one character repeating (eg `NNNNNNNN...`) using ~128n bytes. Expected runtime for a random string is `O(n log^2 n)` and worst case is `O(n^2 log n)`.
