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
TESTPATH = ./test
MKDIR_P = mkdir -p

.PHONY: all directories clean
all: directories $(TARGETPATH)/SWSolver.o $(TARGETPATH)/main $(TARGETPATH)/swissprot_tests

directories:
	$(MKDIR_P) $(TARGETPATH)

$(TARGETPATH)/SWSolver.o: $(SOURCEPATH)/SWSolver.cu
	$(CC) -O2 -c $(SOURCEPATH)/SWSolver.cu -o $(TARGETPATH)/SWSolver.o $(CFLAGS)

$(TARGETPATH)/main: $(SOURCEPATH)/main.cpp ${SOURCEPATH}/SWSolver.cu
	$(CC) -O2 -o $(TARGETPATH)/main $(SOURCEPATH)/main.cpp $(TARGETPATH)/SWSolver.o $(CFLAGS)

$(TARGETPATH)/swissprot_tests: $(SOURCEPATH)/SWSolver.cu $(TESTPATH)/swissprot_tests.cpp
	$(CC) -O2 -c $(TESTPATH)/swissprot_tests.cpp -o $(TARGETPATH)/swissprot_tests.o $(CFLAGS)
	$(CC) -O2 -o $(TARGETPATH)/swissprot_tests $(TARGETPATH)/SWSolver.o $(TARGETPATH)/swissprot_tests.o $(CFLAGS) -lboost_unit_test_framework

clean:
	rm -rf ./bin
