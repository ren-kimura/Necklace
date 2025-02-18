all: max.out

max.out: max.cpp
	g++ -std=c++23 -o max.out max.cpp -g -Wall -Wextra -fsanitize=undefined -Wno-deprecated -O3