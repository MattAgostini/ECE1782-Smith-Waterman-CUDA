#CC=g++
#CFLAGS=
#
#make:
#	$(CC) -o main main.cpp $(CFLAGS)
#
#
CC=nvcc
CFLAGS=-ccbin clang++-3.8 -arch=sm_52

.PHONY: all clean
all: main testchar

main:
	$(CC) -o main main.cu $(CFLAGS)

testchar:
	$(CC) -o testchar testchar.cu $(CFLAGS)

clean:
	rm main testchar
