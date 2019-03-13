#include <sys/time.h>
#include <stdio.h>
#include <iostream>

#define SEQ_EQUAL 3
#define SEQ_DIFF -3
#define GAP_PENALTY 2

#define FROM_LEFT 1
#define FROM_TOP 2
#define FROM_TOP_LEFT 3

#define BLOCK_X_DIM 32
#define BLOCK_Y_DIM 32
#define MAX_BLOCK_SIZE 1024
#define MAX_GRID_DIM 65535

using namespace std;

// Time stamp function
double getTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_usec/1000000 + tv.tv_sec;
}

__global__ void f_scoreSequence(string seqA, string seqB, float* scoringMatrix, float* tracebackMatrix, int width, int height) {
    // Do the scoring
    int subMatrix[2] = {SEQ_EQUAL, SEQ_DIFF};
    
    int maxScore = 0;
    int maxIndexI = 0;
    int maxIndexJ = 0;
    for (int i = 1; i < (seqA.length() + 1); i++) {
        for (int j = 1; j < (seqB.length() + 1); j++) {
            int score = 0;
            
            if (scoringMatrix[i][j - 1] - GAP_PENALTY > score) {
                score = scoringMatrix[i][j - 1] - GAP_PENALTY;
                tracebackMatrix[i][j] = FROM_LEFT;
            }
            
            if (scoringMatrix[i - 1][j] - GAP_PENALTY > score) {
                score = scoringMatrix[i - 1][j] - GAP_PENALTY;
                tracebackMatrix[i][j] = FROM_TOP;
            }
            
            int similarityScore = 0;
            if (seqA[i - 1] == seqB[j - 1]) similarityScore = subMatrix[0];
            else similarityScore = subMatrix[1];
            
            if (scoringMatrix[i - 1][j - 1] + similarityScore > score) {
                score = scoringMatrix[i - 1][j - 1] + similarityScore;
                tracebackMatrix[i][j] = FROM_TOP_LEFT;
            }
            
            if (score > maxScore) {
                maxScore = score;
                maxIndexI = i;
                maxIndexJ = j;
            }
                    
            scoringMatrix[i][j] = score;
        }
    }
}

int main( int argc, char *argv[] ) {
    // get program arguments
    if (argc != 3) {
        printf("Error: wrong number of args\n");
        exit(1);
    }
    
    // Need to take a database as input, read in a number of sequences, then parallelize the queries per thread/block

    string seqA = argv[1];
    string seqB = argv[2];
    
    int scoringMatrix[(seqA.length() + 1)][(seqB.length() + 1)];
    int tracebackMatrix[(seqA.length() + 1)][(seqB.length() + 1)];
    
    // Initialize Matrix
    for (int i = 0; i < (seqA.length() + 1); i++) {
        for (int j = 0; j < (seqB.length() + 1); j++) {
            scoringMatrix[i][j] = 0;
            tracebackMatrix[i][j] = 0;
        }
    }
    
    int numElems = (seqA.length() + 1) * (seqB.length() + 1)
    int bytes = numElems * sizeof(float);
    
    cudaError_t status;
    
    // alloc memory on GPU
    float* d_output_scoring;
    cudaMallocManaged((void**) &d_output_scoring, bytes);
    
        // alloc memory on GPU
    float* d_output_traceback;
    cudaMallocManaged((void**) &d_output_traceback, bytes);
    
    // Call GPU
    dim3 block(1, 1);
    dim3 grid(1, 1);
 
    f_scoreSequence<<<grid, block>>>(seqA, seqB, d_output_scoring, d_output_traceback, 1, 1);
    
    cudaDeviceSynchronize();
    
    string resultA;
    string resultB;
    
    // Do the traceback to get the sequence
    int iTrace = maxIndexI;
    int jTrace = maxIndexJ;
    int value = scoringMatrix[iTrace][jTrace];
    while (value != 0) {
        if (tracebackMatrix[iTrace][jTrace] == FROM_LEFT) {
            jTrace--;
            resultA.append("-");
            resultB.append(seqB, jTrace, 1);
            value = scoringMatrix[iTrace][jTrace];
        } 
        else if (tracebackMatrix[iTrace][jTrace] == FROM_TOP) {
            iTrace--;
            resultA.append(seqA, iTrace, 1);
            resultB.append("-");
            value = scoringMatrix[iTrace][jTrace];
        } 
        else if (tracebackMatrix[iTrace][jTrace] == FROM_TOP_LEFT) {
            jTrace--;
            iTrace--;
            resultA.append(seqA, iTrace, 1);
            resultB.append(seqB, jTrace, 1);
            value = scoringMatrix[iTrace][jTrace];
        }
    }
    
    reverse(resultA.begin(), resultA.end());
    reverse(resultB.begin(), resultB.end());
    cout << resultA << endl;
    cout << resultB << endl;
    
    // Free device memory
    cudaFree(d_output_scoring);
    cudaFree(d_output_traceback);
    cudaDeviceReset();
}
