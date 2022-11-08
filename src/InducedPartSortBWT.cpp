#include <cassert>
#include <vector>
#include <cmath>
#include <algorithm>
#include "InducedPartSortBWT.h"
#include "PartSortInternals.h"

size_t sortSuffixesWithPrefix(const std::vector<uint64_t>& bitseq, const size_t realSize, std::string& bwt, const uint64_t refPrefix, const size_t suffixCount, const size_t doneAlready, std::vector<std::pair<size_t, uint64_t>>& positions, std::vector<uint8_t>& inducedHpcQueue, bool reverseOrder)
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
	assert(positions.size() == suffixCount);
	std::sort(positions.begin(), positions.end(), [](std::pair<size_t, uint64_t> left, std::pair<size_t, uint64_t> right) { return left.second < right.second; });
	size_t start = 0;
	for (size_t i = 1; i < positions.size(); i++)
	{
		if (positions[i].second != positions[i-1].second)
		{
			if (i > start+1)
			{
				chunkRadixSortSuffixesInPlaceNoEscape(bitseq, realSize, positions, PREFIX_LENGTH+21, start, i);
			}
			start = i;
		}
	}
	if (start < positions.size()-1)
	{
		chunkRadixSortSuffixesInPlaceNoEscape(bitseq, realSize, positions, PREFIX_LENGTH+21, start, positions.size());
	}
	assert(doneAlready + positions.size() <= bwt.size());
	for (size_t i = 0; i < positions.size(); i++)
	{
		assert(bwt[doneAlready+i] == 7);
		bwt[doneAlready+i] = positions[i].first >> 61;
		uint8_t hpcchar = refPrefix >> ((PREFIX_LENGTH-1)*3);
		size_t hpcBefore = 0;
		size_t startPos = positions[i].first & posMask;
		while (hpcBefore+1 <= startPos && getChar(bitseq, realSize, startPos-1-hpcBefore) == hpcchar) hpcBefore += 1;
		if (hpcBefore == 0)
		{
			// no homopolymer to induce
			continue;
		}
		uint8_t charBefore = 0;
		if (hpcBefore+1 <= startPos) charBefore = getChar(bitseq, realSize, startPos-1-hpcBefore);
		assert(charBefore != hpcchar);
		hpcBefore -= 1; // the second-to-last char of the homopolymer run was already inserted to result
		if (hpcBefore < 16)
		{
			inducedHpcQueue.push_back((charBefore << 4) + hpcBefore);
		}
		else if (hpcBefore < 2048)
		{
			inducedHpcQueue.push_back(128 + (charBefore << 4) + (hpcBefore & 15));
			inducedHpcQueue.push_back(hpcBefore >> 4);
			if (reverseOrder) std::reverse(inducedHpcQueue.end()-2, inducedHpcQueue.end());
		}
		else
		{
			assert(hpcBefore < 524288); // 2^19, biggest num which can be represented by inducedHpcQueue
			inducedHpcQueue.push_back(128 + (charBefore << 4) + (hpcBefore & 15));
			inducedHpcQueue.push_back(128 + ((hpcBefore >> 4) & 127));
			inducedHpcQueue.push_back(hpcBefore >> 11);
			if (reverseOrder) std::reverse(inducedHpcQueue.end()-3, inducedHpcQueue.end());
		}
	}
	return positions.size();
}

size_t getHomopolymerCounts(const std::vector<uint64_t>& bitseq, const size_t realSize, std::vector<size_t>& homopolymerLCounts, std::vector<size_t>& homopolymerSCounts)
{
	bool nowS = true;
	uint8_t prevChar = 0;
	std::vector<size_t> LqueueSize;
	std::vector<size_t> SqueueSize;
	LqueueSize.resize(8, 0);
	SqueueSize.resize(8, 0);
	size_t currentLength = 0;
	bool prevS = true;
	for (size_t pos = realSize-2; pos < realSize; pos--)
	{
		uint8_t charHere = getChar(bitseq, realSize, pos);
		if (charHere < prevChar) nowS = true;
		if (charHere > prevChar) nowS = false;
		if (charHere == prevChar)
		{
			currentLength += 1;
			if (nowS)
			{
				homopolymerSCounts[charHere] += 1;
			}
			else
			{
				homopolymerLCounts[charHere] += 1;
			}
		}
		else
		{
			size_t addHere = 0;
			if (currentLength > 0 && currentLength < 16)
			{
				addHere = 1;
			}
			else if (currentLength > 0 && currentLength < 2048)
			{
				addHere = 2;
			}
			else if (currentLength > 0)
			{
				addHere = 3;
			}
			if (prevS)
			{
				SqueueSize[prevChar] += addHere;
			}
			else
			{
				LqueueSize[prevChar] += addHere;
			}
			currentLength = 0;
		}
		prevS = nowS;
		prevChar = charHere;
	}
	size_t addHere = 0;
	if (currentLength > 0 && currentLength < 16)
	{
		addHere = 1;
	}
	else if (currentLength > 0 && currentLength < 2048)
	{
		addHere = 2;
	}
	else if (currentLength > 0)
	{
		addHere = 3;
	}
	if (prevS)
	{
		SqueueSize[prevChar] += addHere;
	}
	else
	{
		LqueueSize[prevChar] += addHere;
	}
	size_t result = 0;
	for (size_t i = 0; i < 8; i++)
	{
		result = std::max(result, LqueueSize[i]);
		result = std::max(result, SqueueSize[i]);
	}
	return result;
}

template <int delta>
size_t induceHomopolymers(std::string& result, std::vector<uint8_t>& inducedHpcQueue, uint8_t homopolymerChar, size_t startPos)
{
	size_t inducedCount = 0;
	size_t round = 0;
	while (inducedHpcQueue.size() > 0)
	{
		size_t queueInsertPos = 0;
		size_t queueReadPos = 0;
		size_t induced = 0;
		while (queueReadPos < inducedHpcQueue.size())
		{
			uint8_t charBefore = (inducedHpcQueue[queueReadPos] >> 4) & 7;
			size_t homopolymerLength = inducedHpcQueue[queueReadPos] & 15;
			bool more = inducedHpcQueue[queueReadPos] >= 128;
			queueReadPos += 1;
			if (more)
			{
				assert(queueReadPos < inducedHpcQueue.size());
				homopolymerLength += 16 * (size_t)(inducedHpcQueue[queueReadPos] & 127);
				bool evenMore = inducedHpcQueue[queueReadPos] >= 128;
				queueReadPos += 1;
				if (evenMore)
				{
					assert(queueReadPos < inducedHpcQueue.size());
					homopolymerLength += 2048 * (size_t)inducedHpcQueue[queueReadPos];
					queueReadPos += 1;
				}
			}
			if (homopolymerLength > 0)
			{
				assert(result[startPos] == 7);
				result[startPos] = homopolymerChar;
				homopolymerLength -= 1;
				if (homopolymerLength < 16)
				{
					assert(queueInsertPos < queueReadPos);
					inducedHpcQueue[queueInsertPos] = (charBefore << 4) + homopolymerLength;
					queueInsertPos += 1;
				}
				else if (homopolymerLength < 2048)
				{
					assert(queueInsertPos+1 < queueReadPos);
					inducedHpcQueue[queueInsertPos] = 128 + (charBefore << 4) + (homopolymerLength & 15);
					inducedHpcQueue[queueInsertPos+1] = homopolymerLength >> 4;
					queueInsertPos += 2;
				}
				else
				{
					assert(queueInsertPos+2 < queueReadPos);
					assert(homopolymerLength < 524288);
					inducedHpcQueue[queueInsertPos] = 128 + (charBefore << 4) + (homopolymerLength & 15);
					inducedHpcQueue[queueInsertPos+1] = 128 + ((homopolymerLength >> 4) & 127);
					inducedHpcQueue[queueInsertPos+2] = homopolymerLength >> 11;
					queueInsertPos += 3;
				}
			}
			else
			{
				assert(result[startPos] == 7);
				result[startPos] = charBefore;
			}
			assert(queueInsertPos <= queueReadPos);
			startPos += delta;
			inducedCount += 1;
		}
		if (queueInsertPos == 0) break;
		assert(queueInsertPos <= inducedHpcQueue.size());
		inducedHpcQueue.resize(queueInsertPos);
	}
	return inducedCount;
}

void inducedPartSortBWT(std::string& result, const std::vector<uint64_t>& bitseq, const size_t realSize)
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
	size_t hasCount = 0;
	for (size_t i = 0; i < MAX_PREFIX; i++)
	{
		maxCount = std::max(maxCount, prefixCount[i]);
		if (prefixCount[i] > 0) hasCount += 1;
	}
	std::vector<size_t> homopolymerLCounts;
	std::vector<size_t> homopolymerSCounts;
	homopolymerLCounts.resize(8, 0);
	homopolymerSCounts.resize(8, 0);
	size_t biggestHpcQueue = getHomopolymerCounts(bitseq, realSize, homopolymerLCounts, homopolymerSCounts);
	tmpData.reserve(maxCount);
	std::vector<uint8_t> inducedHpcQueue;
	inducedHpcQueue.reserve(biggestHpcQueue);
	for (uint8_t firstChar = 0; firstChar < 8; firstChar++)
	{
		inducedHpcQueue.clear();
		for (uint8_t secondChar = 0; secondChar < firstChar; secondChar++)
		{
			for (uint64_t i = 0; i < MAX_PREFIX / 8 / 8; i++)
			{
				uint64_t prefix = firstChar << ((PREFIX_LENGTH-1) * 3);
				prefix += secondChar << ((PREFIX_LENGTH-2) * 3);
				prefix += i;
				if (prefixCount[prefix] == 0) continue;
				size_t sorted = sortSuffixesWithPrefix(bitseq, realSize, result, prefix, prefixCount[prefix], doneCount, tmpData, inducedHpcQueue, false);
				tmpData.resize(0);
				doneCount += sorted;
				assert(sorted == prefixCount[prefix]);
			}
		}
		assert(inducedHpcQueue.size() <= biggestHpcQueue);
		size_t induced = induceHomopolymers<1>(result, inducedHpcQueue, firstChar, doneCount);
		assert(induced == homopolymerLCounts[firstChar]);
		doneCount += induced;
		size_t StypeStart = doneCount;
		inducedHpcQueue.clear();
		for (uint8_t secondChar = firstChar+1; secondChar < 8; secondChar++)
		{
			for (uint64_t i = 0; i < MAX_PREFIX / 8 / 8; i++)
			{
				uint64_t prefix = firstChar << ((PREFIX_LENGTH-1) * 3);
				prefix += secondChar << ((PREFIX_LENGTH-2) * 3);
				prefix += i;
				if (prefixCount[prefix] == 0) continue;
				size_t sorted = sortSuffixesWithPrefix(bitseq, realSize, result, prefix, prefixCount[prefix], doneCount+homopolymerSCounts[firstChar], tmpData, inducedHpcQueue, true);
				tmpData.resize(0);
				doneCount += sorted;
				assert(sorted == prefixCount[prefix]);
			}
		}
		assert(inducedHpcQueue.size() <= biggestHpcQueue);
		std::reverse(inducedHpcQueue.begin(), inducedHpcQueue.end());
		induced = induceHomopolymers<-1>(result, inducedHpcQueue, firstChar, StypeStart+homopolymerSCounts[firstChar]-1);
		assert(induced == homopolymerSCounts[firstChar]);
		doneCount += induced;
	}
}

void inducedPartSortBWT(const std::string& input, std::string& output)
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
	inducedPartSortBWT(output, bitseq, realSize);
}
