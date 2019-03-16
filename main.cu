#include <sys/time.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

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

__global__ void f_scoreSequence(float* seqA, float* seqB, float* scoringMatrix, int width, int height) {
    // Do the scoring
    int substitutionMatrix[2] = {SEQ_EQUAL, SEQ_DIFF};
    
    register int xIndex = threadIdx.x + blockIdx.x * blockDim.x;
    //register int yIndex = threadIdx.y + blockIdx.y * blockDim.y;
    
    int maxScore = 0;
    for (int i = 1; i < (height + 1); i++) {
        for (int j = 1; j < (width + 1); j++) {
            int score = 0;
            
            if (scoringMatrix[(i * (width + 1)) + j - 1] - GAP_PENALTY > score) {
                score = scoringMatrix[(i * (width + 1)) + j - 1] - GAP_PENALTY;
            }
            
            if (scoringMatrix[((i - 1) * (width + 1)) + j] - GAP_PENALTY > score) {
                score = scoringMatrix[((i - 1) * (width + 1)) + j] - GAP_PENALTY;
            }
            
            int similarityScore = 0;
            if (seqA[i - 1] == seqB[j - 1]) similarityScore = substitutionMatrix[0];
            else similarityScore = substitutionMatrix[1];
            
            if (scoringMatrix[((i - 1) * (width + 1)) + j - 1] + similarityScore > score) {
                score = scoringMatrix[((i - 1) * (width + 1)) + j - 1] + similarityScore;
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

    string querySequence = argv[1];
    
    // Parse query file
    ifstream datafile;
    datafile.open(argv[2]);
    
    string temp;
    vector<string> subjectSequences;
    while (datafile >> temp) {
       subjectSequences.push_back(temp);
    }
    
    // Just do the first 32 elements for a test
    int largestSubjectLength = subjectSequences[31].length();
    
    datafile.close();
    
    // alloc memory on GPU
    float* d_input_query;
    cudaMallocManaged((void**) &d_input_query, querySequence.length() * sizeof(float));
    
    float* d_input_subject;
    cudaMallocManaged((void**) &d_input_subject, (largestSubjectLength * 32) * sizeof(float));
    
    float* d_output_scoring;
    cudaMallocManaged((void**) &d_output_scoring, ((querySequence.length() + 1) * (largestSubjectLength + 1) * 32) * sizeof(float));
    
    // Convert string to float representation (can't really use strings on the GPU)
    for (int i = 0; i < querySequence.length();i++) {
        switch(querySequence[i])
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
    
    for (int j = 0; j < 32; j++) {
        for (int i = 0; i < largestSubjectLength; i++) {
            switch(subjectSequences[j][i])
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
    }
    
    // Call GPU
    dim3 block(32, 1);
    dim3 grid(1, 1);
 
    f_scoreSequence<<<grid, block>>>(d_input_query, d_input_subject, d_output_scoring, largestSubjectLength, querySequence.length());
    
    cudaDeviceSynchronize();
    
    // Print results
    string seqA = querySequence;
    string seqB = subjectSequences[0];
    
    cout << "    ";
    for (int j = 0; j < (seqB.length() + 1); j++) {
        cout << seqB[j] << " "; 
    }
    cout << endl;
    
    for (int i = 0; i < (seqA.length() + 1); i++) {
        if (i != 0) cout << seqA[i - 1] << " ";
        else cout << "  ";
        for (int j = 0; j < (seqB.length() + 1); j++) {
            cout << d_output_scoring[i * (seqB.length() + 1) + j] << " "; 
        }
        cout << endl;
    }
    
    // Free device memory
    cudaFree(d_input_query);
    cudaFree(d_input_subject);
    cudaFree(d_output_scoring);
    cudaDeviceReset();
}
