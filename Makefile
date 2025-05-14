SDSL_INCLUDE = $(HOME)/local/include
SDSL_LIB = $(HOME)/local/lib

all: rep

rep: rep.cpp
	g++ -std=c++23 -o rep rep.cpp \
	    -I$(SDSL_INCLUDE) -L$(SDSL_LIB) \
	    -lsdsl -ldivsufsort -ldivsufsort64 \
	    -g -Wall -Wextra -fsanitize=undefined -Wno-deprecated -O3

memcheck: rep
	valgrind --leak-check=full --show-leak-kinds=all ./rep
