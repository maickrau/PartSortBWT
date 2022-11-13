## PartSortBWT

Low memory BWT computation. Usage:

```
std::string input;
input.push_back(1); // N
input.push_back(2); // A
input.push_back(3); // C
input.push_back(4); // G
input.push_back(5); // T
input.push_back(0); // $
std::string output;
PartSortBWT(input, output); // computes BWT, stores in output
PartSortBWT(input, input); // computes BWT, stores in the same string and overwrites
```

The input and output can be the same string. For a real example see [src/main.cpp](src/main.cpp). The input string must end with 0, and all other characters must be between 1-5 representing the alphabet `$NACGT`.

### Benchmark

Measured with `/usr/bin/time -v bin/main input_file > /dev/null` using commit [489c112](https://github.com/maickrau/PartSortBWT/commit/489c112769bbd7f4f40aa70021ca52794652c1cb) on a laptop.

| Dataset | Size | Time (h:mm:ss) | Memory |
| --- | --- | --- | --- |
| Human assembly | 6.0 Gbp | 1:15:51 | 9.0 Gb |
| GRCh38+alts | 3.3 Gbp | 0:40:01 | 4.9 Gb |
| Assembly graph hpc | 3.8 Gbp | 0:31:11 | 6.5 Gb |
| CHM13 | 3.1 Gbp | 0:23:45 | 4.6 Gb |
| CHM13 hpc | 2.1 Gbp | 0:12:02 | 3.6 Gb |
| CHM13 chr1 | 248 Mbp | 0:01:40 | 375 Mb |
| CHM13 chr1 hpc | 175 Mbp | 0:00:49 | 291 Mb |
| ACACACAC... | 1 Mbp | 0:05:25 | 11 Mbp |
| ACACACAC... | 100 kbp | 0:00:02 | 3 Mbp |

### Method

Classify suffixes according to their first 4 characters. First collect and sort all suffixes starting with AAAA to get the first x characters of the BWT, then suffixes with AAAC to get the next y characters, etc. Repeat for all 4-character prefixes. Special case for homopolymer runs of >=4 chars: sort the suffixes starting at the ends of the homopolymer runs, then do induced sorting to sort the rest of the homopolymer runs.

Expected memory use for a random DNA string is `input + output + ~0.875n` bytes. If input == output then total memory use is ~1.875n bytes. Less random strings have higher memory use and worst case is dinucleotide repeating (eg `ACACACAC...`) using ~8n bytes. Expected runtime for a random string is `O(n (log n)^2)` and worst case is `O(n^2 log n)`.
