#ifndef PartSortInternals_h
#define PartSortInternals_h

#include <cstdint>
#include <cassert>
#include <vector>
#include <cmath>
#include <algorithm>

const size_t PREFIX_LENGTH = 4;
const size_t ALPHABET_SIZE = 8;
constexpr size_t MAX_PREFIX = pow(ALPHABET_SIZE, PREFIX_LENGTH);
constexpr uint64_t MASK = MAX_PREFIX-1;

template <int length, typename F>
void iteratePrefixSuffixes(const std::vector<uint64_t>& bitseq, const size_t realSize, F callback)
{
	for (size_t i = 0; i < realSize/21; i++)
	{
		uint64_t chunk = bitseq[i];
		for (size_t j = 0; j < 21-PREFIX_LENGTH; j++)
		{
			callback((chunk >> ((21-PREFIX_LENGTH-j)*3)) & MASK, i*21+j);
		}
		chunk <<= PREFIX_LENGTH*3;
		chunk += bitseq[i+1] >> (21-PREFIX_LENGTH)*3;
		for (size_t j = 21-PREFIX_LENGTH; j < 21; j++)
		{
			callback((chunk >> ((21-j)*3)) & MASK, i*21+j);
		}
	}
	if (realSize % 21 != 0)
	{
		size_t i = realSize/21;
		uint64_t chunk = bitseq[i];
		for (size_t j = 0; j < 21-PREFIX_LENGTH && j < realSize % 21; j++)
		{
			callback((chunk >> ((21-PREFIX_LENGTH-j)*3)) & MASK, i*21+j);
		}
		chunk <<= PREFIX_LENGTH*3;
		chunk += bitseq[i+1] >> (21-PREFIX_LENGTH)*3;
		for (size_t j = 21-PREFIX_LENGTH; j < 21 && j < realSize % 21; j++)
		{
			callback((chunk >> ((21-j)*3)) & MASK, i*21+j);
		}
	}
}

uint8_t getChar(const std::vector<uint64_t>& bitseq, const size_t realSize, size_t i);
uint64_t getChunk(const std::vector<uint64_t>& bitseq, const size_t realSize, size_t i);
void chunkRadixSortSuffixesInPlace(const std::vector<uint64_t>& bitseq, const size_t realSize, std::vector<std::pair<size_t, uint64_t>>& positions, size_t offset, const size_t vecStart, const size_t vecEnd);
void chunkRadixSortSuffixesInPlaceNoEscape(const std::vector<uint64_t>& bitseq, const size_t realSize, std::vector<std::pair<size_t, uint64_t>>& positions, size_t offset, const size_t vecStart, const size_t vecEnd);

#endif
