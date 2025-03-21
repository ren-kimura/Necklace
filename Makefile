all: max.out

max.out: max.cpp
	g++ -std=c++23 -o max.out max.cpp -g -Wall -Wextra -fsanitize=undefined -Wno-deprecated -O3

run: max.out
	./max.out

memcheck: max.out
	valgrind --leak-check=full --show-leak-kinds=all ./max.out

massif: max.out
	valgrind --tool=massif --massif-out-file=massif.out ./max.out dat/c/c.fa 17 2
	ms_print massif.out
