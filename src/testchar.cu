#include <sys/time.h>
#include <stdio.h>
#include <iostream>
#include <numeric>
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

__constant__ char constQuery[1024];

// Time stamp function
double getTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_usec/1000000 + tv.tv_sec;
}

__global__ void f_scoreSequence(char* subject, float* scoringMatrix, int width, int height) {
    // Do the scoring
    int substitutionMatrix[2] = {SEQ_EQUAL, SEQ_DIFF};

    register int xIndex = threadIdx.x + blockIdx.x * blockDim.x;
    register int yIndex = threadIdx.y + blockIdx.y * blockDim.y;

    float maxScore = 0;
    for (int i = 1; i < (height + 1); i++) {
        for (int j = 1; j < (width + 1); j++) {
            float score = 0;

            score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + (i * (width + 1)) + j - 1] - GAP_PENALTY);
            score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + ((i - 1) * (width + 1)) + j] - GAP_PENALTY);

            int similarityScore = 0;

            // Just index scoring matrix from shared/constant memory in the future
            if (constQuery[i - 1] == subject[width*yIndex + j - 1]) similarityScore = substitutionMatrix[0];
            else similarityScore = substitutionMatrix[1];

            score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + ((i - 1) * (width + 1)) + j - 1] + similarityScore);
            
            maxScore = max(maxScore, score);

            scoringMatrix[(width + 1)*(height + 1)*yIndex + (i * (width + 1)) + j] = score;
        }
    }
}

int main( int argc, char *argv[] ) {
    double time_start = getTimeStamp();

    // get program arguments
    if (argc != 3) {
        printf("Error: wrong number of args\n");
        exit(1);
    }

    char* querySequence = argv[1];

    // Parse query file
    ifstream datafile;
    datafile.open(argv[2]);

    int subjectLengthSum = 0;



    char* temp;
    vector<char*> subjectSequences;
    while (datafile >> temp) {
       subjectSequences.push_back(temp);
       subjectLengthSum += strlen(temp);
    }

    // Just do the first 32 elements for a test
    int largestSubjectLength = strlen(subjectSequences[31]);

    datafile.close();

    // alloc memory on GPU
    float* d_input_query = new float[strlen(querySequence)];
	memset(d_input_query, 0, sizeof(float) * strlen(querySequence));
 
    char* d_input_subject;
    cudaMallocManaged((void**) &d_input_subject, (largestSubjectLength * 32) * sizeof(char));
	memcpy(d_input_subject, subjectSequences[0], ((largestSubjectLength * 32) + 1) * sizeof(char));
	
	for (int i = 1; i < subjectSequences.size(); i++)
		strcat(d_input_subject, subjectSequences[i]);

    float* d_output_scoring;
    cudaMallocManaged((void**) &d_output_scoring, ((strlen(querySequence) + 1) * (largestSubjectLength + 1) * 32) * sizeof(float));

    // Convert string to float representation (can't really use strings on the GPU)
    for (int i = 0; i < strlen(querySequence);i++) {
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
	/*
    for (int j = 0; j < 32; j++) {
        for (int i = 0; i < largestSubjectLength; i++) {
            switch(subjectSequences[j][i])
            {
                case 'A': { d_input_subject[i] = 'A';
                            break;
                        }
                case 'G': { d_input_subject[i] = 'G';
                            break;
                        }
                case 'C': { d_input_subject[i] = 'C';
                            break;
                        }
                case 'T': { d_input_subject[i] = 'T';
                            break;
                        }
            }
        }
    }
	*/
	cudaMemcpyToSymbol(constQuery, querySequence, sizeof(char)*strlen(querySequence));

    // Call GPU
    dim3 block(1, 32);
    dim3 grid(1, 1);
 
    f_scoreSequence<<<grid, block>>>(d_input_subject, d_output_scoring, largestSubjectLength, strlen(querySequence));

    cudaDeviceSynchronize();

    // Print results for 1 subject query
    for (int subject = 0; subject < 32; subject++) {
        char* seqA = querySequence;
        char* seqB = subjectSequences[subject];

        cout << "    ";
        for (int j = 0; j < (strlen(seqB) + 1); j++) {
            cout << seqB[j] << " ";
        }
        cout << endl;

        for (int i = 0; i < (strlen(seqA) + 1); i++) {
            if (i != 0) cout << seqA[i - 1] << " ";
            else cout << "  ";
            for (int j = 0; j < (strlen(seqB) + 1); j++) {
                cout << d_output_scoring[((largestSubjectLength + 1) * (strlen(querySequence) + 1) * subject) + (i * (strlen(seqB) + 1)) + j] << " ";
            }
            cout << endl;
        }
    }

    double time_end = getTimeStamp();
    double seconds_elapsed = time_end - time_start;

    std::cout << std::string(80, '=') << std::endl;
    cout << "METRICS:" << endl;
    cout << "Query length: " << strlen(querySequence) << " chars." << endl;
    cout << "Sum of DB length: " << subjectLengthSum << " chars." << endl;
    cout << "Time elapsed: " << seconds_elapsed << " seconds." << endl;
    cout << "Performance: " << 1E-9 * (strlen(querySequence) * subjectLengthSum)
            / seconds_elapsed << " GCUPS." << endl;

    // Free device memory
    cudaFree(d_input_query);
    cudaFree(d_input_subject);
    cudaFree(d_output_scoring);
    cudaDeviceReset();
}
