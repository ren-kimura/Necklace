CXX      := g++
CXXFLAGS := -std=c++23 -g -Wall -Wextra -fsanitize=undefined -Wno-deprecated -O3

all: rep birep

rep: rep.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

birep: birep.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

memcheck-rep: rep
	valgrind --leak-check=full --show-leak-kinds=all ./rep

memcheck-birep: birep
	valgrind --leak-check=full --show-leak-kinds=all ./birep

.PHONY: all clean memcheck-rep memcheck-birep
clean:
	rm -f rep birep
