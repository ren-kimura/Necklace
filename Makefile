all: rep.out

rep.out: rep.cpp
	g++ -std=c++23 -o rep.out rep.cpp -g -Wall -Wextra -fsanitize=undefined -Wno-deprecated -O3

memcheck: rep.out
	valgrind --leak-check=full --show-leak-kinds=all ./rep.out