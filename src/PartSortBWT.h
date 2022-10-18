#ifndef PartSortBWT_h
#define PartSortBWT_h

#include <string>

// compute BWT, store in output. Input and output strings can be the same.
// Uses expected input + output + ~0.875n bytes for random strings.
// If input == output then ~1.875n bytes total.
void partSortBWT(const std::string& input, std::string& output);

#endif
