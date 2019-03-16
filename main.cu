#include <sys/time.h>
#include <stdio.h>
#include <iostream>

#define SEQ_EQUAL 3
#define SEQ_DIFF -3
#define GAP_PENALTY 2
// define affine penalty ?

#define FROM_LEFT 1
#define FROM_TOP 2
#define FROM_TOP_LEFT 3

#define MAX_BLOCK_SIZE 1024
#define MAX_GRID_DIM 65535

#define A 1
#define G 2
#define C 3
#define T 4

using namespace std;

// Time stamp function
double getTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_usec/1000000 + tv.tv_sec;
}

__global__ void f_scoreSequence(float* seqA, float* seqB, float* scoringMatrix, float* tracebackMatrix, int width, int height) {
    // Do the scoring
    int substitutionMatrix[2] = {SEQ_EQUAL, SEQ_DIFF};
    
    int maxScore = 0;
    for (int i = 1; i < (height + 1); i++) {
        for (int j = 1; j < (width + 1); j++) {
            int score = 0;
            
            if (scoringMatrix[(i * (width + 1)) + j - 1] - GAP_PENALTY > score) {
                score = scoringMatrix[(i * (width + 1)) + j - 1] - GAP_PENALTY;
                tracebackMatrix[(i * (width + 1)) + j] = FROM_LEFT;
            }
            
            if (scoringMatrix[((i - 1) * (width + 1)) + j] - GAP_PENALTY > score) {
                score = scoringMatrix[((i - 1) * (width + 1)) + j] - GAP_PENALTY;
                tracebackMatrix[(i * (width + 1)) + j] = FROM_TOP;
            }
            
            int similarityScore = 0;
            if (seqA[i - 1] == seqB[j - 1]) similarityScore = substitutionMatrix[0];
            else similarityScore = substitutionMatrix[1];
            
            if (scoringMatrix[((i - 1) * (width + 1)) + j - 1] + similarityScore > score) {
                score = scoringMatrix[((i - 1) * (width + 1)) + j - 1] + similarityScore;
                tracebackMatrix[(i * (width + 1)) + j] = FROM_TOP_LEFT;
            }
            
            if (score > maxScore) {
                maxScore = score;
            }
                    
            scoringMatrix[(i * (width + 1)) + j] = score;
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
    
    int numElems = (seqA.length() + 1) * (seqB.length() + 1);
    int bytes = numElems * sizeof(float);
    
    // alloc memory on GPU
    float* d_input_query;
    cudaMallocManaged((void**) &d_input_query, bytes);
    
    float* d_input_subject;
    cudaMallocManaged((void**) &d_input_subject, bytes);
    
    float* d_output_scoring;
    cudaMallocManaged((void**) &d_output_scoring, bytes);
    
    float* d_output_traceback;
    cudaMallocManaged((void**) &d_output_traceback, bytes);
    
    // Convert string to float representation (can't really use strings on the GPU)
    for (int i = 0; i < seqA.length();i++) {
        switch(seqA[i])
        {
            case 'A': { d_input_query[i] = A;
                        break;
                    }
            case 'G': { d_input_query[i] = G;
                        break;
                    }
            case 'C': { d_input_query[i] = C;
                        break;
                    }
            case 'T': { d_input_query[i] = T;
                        break;
                    }
        }
    }
    
    for (int i = 0; i < seqB.length(); i++) {
        switch(seqB[i])
        {
            case 'A': { d_input_subject[i] = A;
                        break;
                    }
            case 'G': { d_input_subject[i] = G;
                        break;
                    }
            case 'C': { d_input_subject[i] = C;
                        break;
                    }
            case 'T': { d_input_subject[i] = T;
                        break;
                    }
        }
    }
    
    // Call GPU
    dim3 block(1, 1);
    dim3 grid(1, 1);
 
    f_scoreSequence<<<grid, block>>>(d_input_query, d_input_subject, d_output_scoring, d_output_traceback, seqB.length(), seqA.length());
    
    cudaDeviceSynchronize();
    
    /*
    
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
    
    */
    
    // Print results
    cout << "    ";
    for (int j = 0; j < (seqB.length() + 1); j++) {
        cout << seqB[j] << " "; 
    }
    cout << endl;
    
    for (int i = 0; i < (seqA.length() + 1); i++) {
        if (i != 0) cout << seqA[i - 1] << " ";
        else cout << "  ";
        for (int j = 0; j < (seqB.length() + 1); j++) {
            cout << d_output_scoring[(i * (seqB.length() + 1)) + j] << " "; 
        }
        cout << endl;
    }
    
    // Free device memory
    cudaFree(d_input_query);
    cudaFree(d_input_subject);
    cudaFree(d_output_scoring);
    cudaFree(d_output_traceback);
    cudaDeviceReset();
}
