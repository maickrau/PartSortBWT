#include <cassert>
#include <vector>
#include <cmath>
#include <algorithm>
#include "PartSortInternals.h"

uint8_t getChar(const std::vector<uint64_t>& bitseq, const size_t realSize, size_t i)
{
	i = i % realSize;
	size_t chunk = i / 21;
	size_t offset = 20 - (i%21);
	uint8_t result = (bitseq[chunk] >> (offset * 3)) & 0x7;
	assert(result <= 5);
	return result;
}

uint64_t getChunk(const std::vector<uint64_t>& bitseq, const size_t realSize, size_t i)
{
	i = i % realSize;
	size_t chunk = i/21;
	size_t offset = i%21;
	if (offset == 0)
	{
		return bitseq[chunk];
	}
	uint64_t high = (bitseq[chunk] << (offset*3 + 1)) >> 1; // +1 and -1 to handle the extra 64th bit
	uint64_t low = bitseq[chunk+1] >> ((21 - offset)*3);
	return low + high;
}

void chunkRadixSortSuffixesInPlace(const std::vector<uint64_t>& bitseq, const size_t realSize, std::vector<std::pair<size_t, uint64_t>>& positions, size_t offset, const size_t vecStart, const size_t vecEnd)
{
	if (vecEnd - vecStart < 2) return;
	const size_t posMask = ((size_t)1 << (size_t)61) - (size_t)1;
	for (size_t i = vecStart; i < vecEnd; i++)
	{
		positions[i].second = getChunk(bitseq, realSize, (positions[i].first & posMask) + offset);
	}
	std::sort(positions.begin()+vecStart, positions.begin()+vecEnd, [](std::pair<size_t, uint64_t> left, std::pair<size_t, uint64_t> right) { return left.second < right.second; });
	size_t start = vecStart;
	uint8_t uniquePrefixChar = positions[vecStart].first >> 61;
	for (size_t i = vecStart+1; i < vecEnd; i++)
	{
		if (positions[i].second != positions[i-1].second)
		{
			if (i > start+1 && uniquePrefixChar == 7) chunkRadixSortSuffixesInPlace(bitseq, realSize, positions, offset+21, start, i);
			start = i;
			uniquePrefixChar = positions[i].first >> 61;
		}
		else if (positions[i].first >> 61 != uniquePrefixChar)
		{
			uniquePrefixChar = 7;
		}
	}
	if (start != vecEnd-1 && uniquePrefixChar == 7)
	{
		chunkRadixSortSuffixesInPlace(bitseq, realSize, positions, offset+21, start, vecEnd);
	}
}

void chunkRadixSortSuffixesInPlaceNoEscape(const std::vector<uint64_t>& bitseq, const size_t realSize, std::vector<std::pair<size_t, uint64_t>>& positions, size_t offset, const size_t vecStart, const size_t vecEnd)
{
	if (vecEnd - vecStart < 2) return;
	const size_t posMask = ((size_t)1 << (size_t)61) - (size_t)1;
	for (size_t i = vecStart; i < vecEnd; i++)
	{
		positions[i].second = getChunk(bitseq, realSize, (positions[i].first & posMask) + offset);
	}
	std::sort(positions.begin()+vecStart, positions.begin()+vecEnd, [](std::pair<size_t, uint64_t> left, std::pair<size_t, uint64_t> right) { return left.second < right.second; });
	size_t start = vecStart;
	for (size_t i = vecStart+1; i < vecEnd; i++)
	{
		if (positions[i].second != positions[i-1].second)
		{
			if (i > start+1) chunkRadixSortSuffixesInPlaceNoEscape(bitseq, realSize, positions, offset+21, start, i);
			start = i;
		}
	}
	if (start != vecEnd-1)
	{
		chunkRadixSortSuffixesInPlaceNoEscape(bitseq, realSize, positions, offset+21, start, vecEnd);
	}
}
