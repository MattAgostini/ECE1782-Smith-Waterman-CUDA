#CC=g++
#CFLAGS=
#
#make:
#	$(CC) -o main main.cpp $(CFLAGS)
#
#
CC=nvcc
CFLAGS=-ccbin clang++-3.8 -arch=sm_52

make:
	$(CC) -o main main.cu $(CFLAGS)