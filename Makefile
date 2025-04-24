all: rep

rep: rep.cpp
	g++ -std=c++23 -o rep rep.cpp -g -Wall -Wextra -fsanitize=undefined -Wno-deprecated -O3

memcheck: rep
	valgrind --leak-check=full --show-leak-kinds=all ./rep