#CC=g++
#CFLAGS=
#
#make:
#	$(CC) -o main main.cpp $(CFLAGS)
#
#
CC=nvcc
CFLAGS=-ccbin clang++-3.8 -arch=sm_52 -lboost_program_options
TARGETPATH = ./bin
SOURCEPATH = ./src
MKDIR_P = mkdir -p

.PHONY: all directories clean
all: directories $(TARGETPATH)/main $(TARGETPATH)/testchar

directories:
	$(MKDIR_P) $(TARGETPATH)

$(TARGETPATH)/main: $(SOURCEPATH)/main.cu
	$(CC) -o $(TARGETPATH)/main $(SOURCEPATH)/main.cu $(CFLAGS)

$(TARGETPATH)/testchar: $(SOURCEPATH)/testchar.cu
	$(CC) -o $(TARGETPATH)/testchar $(SOURCEPATH)/testchar.cu $(CFLAGS)

clean:
	rm -rf ./bin
