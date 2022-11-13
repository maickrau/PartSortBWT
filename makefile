GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g

all: bin/main

clean:
	rm -rf bin

bin/main: src/main.cpp src/PartSortBWT.cpp src/PartSortBWT.h
	$(shell mkdir -p bin)
	$(GPP) -o $@ $^ $(CPPFLAGS) $(LINKFLAGS)
