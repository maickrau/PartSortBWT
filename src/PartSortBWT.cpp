#include <cassert>
#include <vector>
#include <cmath>
#include <algorithm>
#include "PartSortBWT.h"

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
	for (size_t i = vecStart; i < vecEnd; i++)
	{
		positions[i].second = getChunk(bitseq, realSize, positions[i].first+offset);
	}
	std::sort(positions.begin()+vecStart, positions.begin()+vecEnd, [](std::pair<size_t, uint64_t> left, std::pair<size_t, uint64_t> right) { return left.second < right.second; });
	size_t start = vecStart;
	for (size_t i = vecStart+1; i < vecEnd; i++)
	{
		if (positions[i].second != positions[i-1].second)
		{
			if (i > start+1) chunkRadixSortSuffixesInPlace(bitseq, realSize, positions, offset+21, start, i);
			start = i;
		}
	}
	if (start != vecEnd-1)
	{
		chunkRadixSortSuffixesInPlace(bitseq, realSize, positions, offset+21, start, vecEnd);
	}
}

size_t sortSuffixesWithPrefix(const std::vector<uint64_t>& bitseq, const size_t realSize, std::string& bwt, const uint64_t refPrefix, const size_t suffixCount, const size_t doneAlready, std::vector<std::pair<size_t, uint64_t>>& positions)
{
	assert(positions.capacity() >= suffixCount);
	iteratePrefixSuffixes<PREFIX_LENGTH>(bitseq, realSize, [refPrefix, &positions, &bitseq, realSize](uint64_t prefix, size_t pos)
	{
		if (prefix != refPrefix) return;
		uint64_t suffixPrefix = getChunk(bitseq, realSize, pos+PREFIX_LENGTH);
		positions.emplace_back(pos, suffixPrefix);
	});
	std::sort(positions.begin(), positions.end(), [](std::pair<size_t, uint64_t> left, std::pair<size_t, uint64_t> right) { return left.second < right.second; });
	size_t start = 0;
	for (size_t i = 1; i < positions.size(); i++)
	{
		if (positions[i].second != positions[i-1].second)
		{
			if (i > start+1)
			{
				chunkRadixSortSuffixesInPlace(bitseq, realSize, positions, PREFIX_LENGTH+21, start, i);
			}
			start = i;
		}
	}
	if (start < positions.size()-1)
	{
		chunkRadixSortSuffixesInPlace(bitseq, realSize, positions, PREFIX_LENGTH+21, start, positions.size());
	}
	assert(doneAlready + positions.size() <= bwt.size());
	for (size_t i = 0; i < positions.size(); i++)
	{
		assert(positions[i].first < realSize);
		size_t pos = 0;
		if (positions[i].first == 0)
		{
			pos = realSize-1;
		}
		else
		{
			pos = positions[i].first-1;
		}
		assert(bwt[doneAlready+i] == 7);
		bwt[doneAlready+i] = getChar(bitseq, realSize, pos);
	}
	return positions.size();
}

void partSortBWT(std::string& result, const std::vector<uint64_t>& bitseq, const size_t realSize)
{
	result.assign(realSize, 7);
	size_t doneCount = 0;
	std::vector<size_t> prefixCount;
	prefixCount.resize(MAX_PREFIX, 0);
	size_t counted = 0;
	iteratePrefixSuffixes<PREFIX_LENGTH>(bitseq, realSize, [&counted, &prefixCount](uint64_t prefix, size_t pos)
	{
		prefixCount[prefix] += 1;
		counted += 1;
	});
	assert(counted == realSize);
	std::vector<std::pair<size_t, uint64_t>> tmpData;
	size_t maxCount = 0;
	for (size_t i = 0; i < MAX_PREFIX; i++)
	{
		maxCount = std::max(maxCount, prefixCount[i]);
	}
	tmpData.reserve(maxCount);
	for (uint64_t i = 0; i < MAX_PREFIX; i++)
	{
		if (prefixCount[i] == 0) continue;
		uint64_t prefix = i;
		size_t sorted = sortSuffixesWithPrefix(bitseq, realSize, result, prefix, prefixCount[i], doneCount, tmpData);
		tmpData.resize(0);
		doneCount += sorted;
		assert(sorted == prefixCount[i]);
	}
	assert(doneCount == realSize);
}

void partSortBWT(const std::string& input, std::string& output)
{
	for (size_t i = 0; i < input.size(); i++)
	{
		assert(input[i] <= 5);
		assert(input[i] >= 1 || (i == input.size()-1 && input[i] == 0));
	}
	std::vector<uint64_t> bitseq;
	bitseq.resize(input.size() / 21 + 2);
	for (size_t i = 0; i < bitseq.size(); i++)
	{
		for (size_t j = 0; j < 21; j++)
		{
			bitseq[i] <<= 3;
			bitseq[i] += input[(i*21+j) % input.size()];
		}
	}
	size_t realSize = input.size();
	partSortBWT(output, bitseq, realSize);
}
