all: max.out, gen.out

max.out: max.cpp
	g++ -o max.out max.cpp -g -O3

gen.out: sequence_gen.cpp
	g++ -o gen.out sequence_gen.cpp -g -O3