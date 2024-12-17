all: max.out, gen.out

max.out: max.cpp
	g++ -std=c++23 -o max.out max.cpp -g -Wall -Wextra -fsanitize=undefined -Wno-deprecated -O3

gen.out: sequence_gen.cpp
	g++ -std=c++23 -o gen.out sequence_gen.cpp -g -Wall -Wextra -fsanitize=undefined -Wno-deprecated -O3
