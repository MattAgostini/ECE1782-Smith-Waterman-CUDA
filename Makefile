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
all: directories $(TARGETPATH)/SWSolver.o $(TARGETPATH)/SWSolver_char.o $(TARGETPATH)/main $(TARGETPATH)/testchar $(TARGETPATH)/swissprot_tests

directories:
	$(MKDIR_P) $(TARGETPATH)

$(TARGETPATH)/SWSolver.o: $(SOURCEPATH)/SWSolver.cu
	$(CC) -O2 -c $(SOURCEPATH)/SWSolver.cu -o $(TARGETPATH)/SWSolver.o $(CFLAGS)

$(TARGETPATH)/SWSolver_char.o: $(SOURCEPATH)/SWSolver_char.cu
	$(CC) -O2 -c $(SOURCEPATH)/SWSolver_char.cu -o $(TARGETPATH)/SWSolver_char.o $(CFLAGS)

$(TARGETPATH)/main: $(SOURCEPATH)/main.cpp ${SOURCEPATH}/SWSolver.cu
	$(CC) -O2 -o $(TARGETPATH)/main $(SOURCEPATH)/main.cpp $(TARGETPATH)/SWSolver.o $(CFLAGS)

$(TARGETPATH)/testchar: $(SOURCEPATH)/testchar.cpp ${SOURCEPATH}/SWSolver_char.cu
	$(CC) -O2 -o $(TARGETPATH)/testchar $(SOURCEPATH)/testchar.cpp $(TARGETPATH)/SWSolver_char.o $(CFLAGS)

$(TARGETPATH)/swissprot_tests: $(SOURCEPATH)/SWSolver.cu $(SOURCEPATH)/SWSolver_char.cu $(TESTPATH)/swissprot_tests.cpp
	$(CC) -O2 -c $(TESTPATH)/swissprot_tests.cpp -o $(TARGETPATH)/swissprot_tests.o $(CFLAGS)
	$(CC) -O2 -o $(TARGETPATH)/swissprot_tests $(TARGETPATH)/SWSolver.o $(TARGETPATH)/SWSolver_char.o $(TARGETPATH)/swissprot_tests.o $(CFLAGS) -lboost_unit_test_framework

clean:
	rm -rf ./bin
