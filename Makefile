all: max.out

max.out: max.cpp
	g++ -o max.out max.cpp -g -O3
