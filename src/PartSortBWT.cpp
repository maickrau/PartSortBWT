#include <cassert>
#include <vector>
#include <cmath>
#include <algorithm>
#include "PartSortBWT.h"
#include "PartSortInternals.h"

size_t sortSuffixesWithPrefix(const std::vector<uint64_t>& bitseq, const size_t realSize, std::string& bwt, const uint64_t refPrefix, const size_t suffixCount, const size_t doneAlready, std::vector<std::pair<size_t, uint64_t>>& positions)
{
	assert(positions.capacity() >= suffixCount);
	const size_t posMask = ((size_t)1 << (size_t)61) - 1;
	iteratePrefixSuffixes<PREFIX_LENGTH>(bitseq, realSize, [refPrefix, &positions, &bitseq, realSize](uint64_t prefix, size_t pos)
	{
		if (prefix != refPrefix) return;
		uint8_t charBefore = 0;
		if (pos > 0) charBefore = getChar(bitseq, realSize, pos-1);
		assert(charBefore < 8);
		assert(pos < (size_t)1 << (size_t)61);
		uint64_t suffixPrefix = getChunk(bitseq, realSize, pos+PREFIX_LENGTH);
		positions.emplace_back(pos + ((size_t)charBefore << (size_t)61), suffixPrefix);
	});
	std::sort(positions.begin(), positions.end(), [](std::pair<size_t, uint64_t> left, std::pair<size_t, uint64_t> right) { return left.second < right.second; });
	size_t start = 0;
	uint8_t uniquePrefixChar = positions[0].first >> 61;
	for (size_t i = 1; i < positions.size(); i++)
	{
		if (positions[i].second != positions[i-1].second)
		{
			if (i > start+1 && uniquePrefixChar == 7)
			{
				chunkRadixSortSuffixesInPlace(bitseq, realSize, positions, PREFIX_LENGTH+21, start, i);
			}
			start = i;
			uniquePrefixChar = positions[i].first >> 61;
		}
		else if ((positions[i].first >> 61) != uniquePrefixChar)
		{
			uniquePrefixChar = 7;
		}
	}
	if (start < positions.size()-1 && uniquePrefixChar == 7)
	{
		chunkRadixSortSuffixesInPlace(bitseq, realSize, positions, PREFIX_LENGTH+21, start, positions.size());
	}
	assert(doneAlready + positions.size() <= bwt.size());
	for (size_t i = 0; i < positions.size(); i++)
	{
		assert(bwt[doneAlready+i] == 7);
		bwt[doneAlready+i] = positions[i].first >> 61;
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
