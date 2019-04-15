#include <iostream>
#include <algorithm>

#include "FASTAParsers.h"
#include "SWSolver.h"

#define GAP_PENALTY 2
// define affine penalty ?

#define FROM_LEFT 1
#define FROM_TOP 2
#define FROM_TOP_LEFT 3

#define MAX_BLOCK_SIZE 1024
#define MAX_GRID_DIM 65535

#define A 0
#define R 1
#define N 2
#define D 3
#define C 4
#define Q 5
#define E 6
#define G 7
#define H 8
#define I 9
#define L 10
#define K 11
#define M 12
#define F 13
#define P 14
#define S 15
#define T 16
#define W 17
#define Y 18
#define V 19
#define B 20
#define J 21
#define Z 22
#define X 23
#define STAR 24

#define BLOCK_Y_DIM 32.0

#define MAX_LENGTH 105
#define CONSTANT_SIZES 4096
#define TILE_SIZE 8

#define GPU_MEM_THRESH 3200000000
#define CPU_MEM_THRESH  150000000

// first is sequence ID, second is max score

short blosum50[25][25] = {
    //        A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
    /* A */ { 5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-2,-1,-1,-3,-1, 1, 0,-3,-2, 0,-2,-2,-1,-1, 0 },
    /* R */ {-2, 7,-1,-2,-4, 1, 0,-3, 0,-4,-3, 3,-2,-3,-3,-1,-1,-3,-1,-3,-1,-3, 0,-1, 0 },
    /* N */ {-1,-1, 7, 2,-2, 0, 0, 0, 1,-3,-4, 0,-2,-4,-2, 1, 0,-4,-2,-3, 5,-4, 0,-1, 0 },
    /* D */ {-2,-2, 2, 8,-4, 0, 2,-1,-1,-4,-4,-1,-4,-5,-1, 0,-1,-5,-3,-4, 6,-4, 1,-1, 0 },
    /* C */ {-1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1,-3,-2,-3,-1, 0 }, 
    /* Q */ {-1, 1, 0, 0,-3, 7, 2,-2, 1,-3,-2, 2, 0,-4,-1, 0,-1,-1,-1,-3, 0,-3, 4,-1, 0 },  
    /* E */ {-1, 0, 0, 2,-3, 2, 6,-3, 0,-4,-3, 1,-2,-3,-1,-1,-1,-3,-2,-3, 1,-3, 5,-1, 0 }, 
    /* G */ { 0,-3, 0,-1,-3,-2,-3, 8,-2,-4,-4,-2,-3,-4,-2, 0,-2,-3,-3,-4,-1,-4,-2,-1, 0 },
    /* H */ {-2, 0, 1,-1,-3, 1, 0,-2,10,-4,-3, 0,-1,-1,-2,-1,-2,-3, 2,-4, 0,-3, 0,-1, 0 },
    /* I */ {-1,-4,-3,-4,-2,-3,-4,-4,-4, 5, 2,-3, 2, 0,-3,-3,-1,-3,-1, 4,-4, 4,-3,-1, 0 },
    /* L */ {-2,-3,-4,-4,-2,-2,-3,-4,-3, 2, 5,-3, 3, 1,-4,-3,-1,-2,-1, 1,-4, 4,-3,-1, 0 },
    /* K */ {-1, 3, 0,-1,-3, 2, 1,-2, 0,-3,-3, 6,-2,-4,-1, 0,-1,-3,-2,-3, 0,-3, 1,-1, 0 },
    /* M */ {-1,-2,-2,-4,-2, 0,-2,-3,-1, 2, 3,-2, 7, 0,-3,-2,-1,-1, 0, 1,-3, 2,-1,-1, 0 },
    /* F */ {-3,-3,-4,-5,-2,-4,-3,-4,-1, 0, 1,-4, 0, 8,-4,-3,-2, 1, 4,-1,-4, 1,-4,-1, 0 },
    /* P */ {-1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3,-2,-3,-1,-1, 0 },
    /* S */ { 1,-1, 1, 0,-1, 0,-1, 0,-1,-3,-3, 0,-2,-3,-1, 5, 2,-4,-2,-2, 0,-3, 0,-1, 0 },
    /* T */ { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 2, 5,-3,-2, 0, 0,-1,-1,-1, 0 },
    /* W */ {-3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1, 1,-4,-4,-3,15, 2,-3,-5,-2,-2,-1, 0 },
    /* Y */ {-2,-1,-2,-3,-3,-1,-2,-3, 2,-1,-1,-2, 0, 4,-3,-2,-2, 2, 8,-1,-3,-1,-2,-1, 0 },
    /* V */ { 0,-3,-3,-4,-1,-3,-3,-4,-4, 4, 1,-3, 1,-1,-3,-2, 0,-3,-1, 5,-3, 2,-3,-1, 0 },
    /* B */ {-2,-1, 5, 6,-3, 0, 1,-1, 0,-4,-4, 0,-3,-4,-2, 0, 0,-5,-3,-3, 6,-4, 1,-1, 0 },
    /* J */ {-2,-3,-4,-4,-2,-3,-3,-4,-3, 4, 4,-3, 2, 1,-3,-3,-1,-2,-1, 2,-4, 4,-3,-1, 0 },
    /* Z */ {-1, 0, 0, 1,-3, 4, 5,-2, 0,-3,-3, 1,-1,-4,-1, 0,-1,-2,-2,-3, 1,-3, 5,-1, 0 }, 
    /* X */ {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0 },
    /* * */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
};

using namespace std;

__constant__ short constQuery[1024];
__constant__ short constSubstitutionMatrix[625];
__constant__ unsigned int constSubjectLengths[CONSTANT_SIZES];
__constant__ unsigned int constSubjectOffsets[CONSTANT_SIZES];
__constant__ unsigned int constScoringOffsets[CONSTANT_SIZES];

inline float convertStringToFloat(char character) {
    switch(character)
    {
        case 'A': { return A; }
        case 'R': { return R; }
        case 'N': { return N; }
        case 'D': { return D; }
        case 'C': { return C; }
        case 'Q': { return Q; }
        case 'E': { return E; }
        case 'G': { return G; }
        case 'H': { return H; }
        case 'I': { return I; }
        case 'L': { return L; }
        case 'K': { return K; }
        case 'M': { return M; }
        case 'F': { return F; }
        case 'P': { return P; }
        case 'S': { return S; }
        case 'T': { return T; }
        case 'W': { return W; }
        case 'Y': { return Y; }
        case 'V': { return V; }
        case 'B': { return B; }
        case 'J': { return J; }
        case 'Z': { return Z; }
        case 'X': { return X; }
    }
    return STAR;
}

// Kernel function for computing the scoring matrix of a sequence
__global__ void f_scoreSequence(float* subject, float* scoringMatrix, float* maxScoreList, 
        int width /*largestSubjectLength*/, int height /*querySequence.length()*/, int numSubjects) {


    //register int xIndex = threadIdx.x + blockIdx.x * blockDim.x;
    register int yIndex = threadIdx.y + blockIdx.y * blockDim.y;

    if (yIndex == 0) {
        //printf("GPU: %d %f\n", constSubstitutionMatrix[24], constQuery[0]);
    }

    float maxScore = 0;
    if (yIndex < numSubjects) {
        for (int i = 1; i < (height + 1); i += 8) {
            for (int j = 1; j < (width + 1); j += 8) {
                for (int k = 0; k < 8; k++) {
                    float score = 0;

                    score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + ((i + k) * (width + 1)) + (j + k) - 1] - GAP_PENALTY);
                    score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + (((i + k) - 1) * (width + 1)) + (j + k)] - GAP_PENALTY);

                    int similarityScore;
                    // check for padded /'s, (dummy symbol) -> score of zero with a real symbol
                    if ((int)constQuery[(i + k) - 1] == 25 || (int)subject[width*yIndex + (j + k) - 1] == 25)
                        similarityScore = 0;
                    else
                        similarityScore = constSubstitutionMatrix[((int)constQuery[(i + k) - 1] * 25) + (int)subject[width*yIndex + (j + k) - 1]];

                    score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + (((i + k) - 1) * (width + 1)) + (j + k) - 1] + similarityScore);

                    maxScore = max(maxScore, score);

                    scoringMatrix[(width + 1)*(height + 1)*yIndex + ((i + k) * (width + 1)) + (j + k)] = score;
                }
            }
        }
        maxScoreList[yIndex] = maxScore;
    }
}

// Kernel function for computing the scoring matrix of a sequence
__global__ void f_scoreSequenceCoalesced(short* subject, short* scoringMatrix, short* maxScoreList, int height /*querySequence.length()*/, int scoreOffset) {

    //register int xIndex = threadIdx.x + blockIdx.x * blockDim.x;
    register int yIndex = threadIdx.y + blockIdx.y * blockDim.y;

    // Use map for different offsets (Change the width)
    unsigned int width = constSubjectLengths[blockIdx.y];
    unsigned int subjectOffset = constSubjectOffsets[blockIdx.y];
    unsigned int blockOffset = constScoringOffsets[blockIdx.y];

    for (int i = 0; i < (height + 1); ++i) {
        scoringMatrix[blockOffset + (threadIdx.y + (blockDim.y * (i)))] = 0;
    }

    for (int j = 0; j < (width + 1); ++j) {
        scoringMatrix[blockOffset + (threadIdx.y + ((j) * blockDim.y * (height + 1)))] = 0;
    }

    int maxScore = 0;
    for (int i = 1; i < (height + 1); i += TILE_SIZE) {
        for (int j = 1; j < (width + 1); j += TILE_SIZE) {
            for (int k = 0; k < TILE_SIZE; k++) {
                for (int m = 0; m < TILE_SIZE; m++) {
                    int score = 0;

                    score = max(score, (int)scoringMatrix[blockOffset + (threadIdx.y + (((j + k) - 1) * blockDim.y * (height + 1))) + (blockDim.y * (i + m))] - GAP_PENALTY); // F[i, j]
                    score = max(score, (int)scoringMatrix[blockOffset + (threadIdx.y + ((j + k) * blockDim.y * (height + 1))) + (blockDim.y * ((i + m) - 1))] - GAP_PENALTY); // E[i, j]

                    int similarityScore = constSubstitutionMatrix[((int)constQuery[(i + m) - 1] * 25) + (int)subject[subjectOffset + threadIdx.y + (((j + k) - 1) * blockDim.y)]];
                    score = max(score, (int)scoringMatrix[blockOffset + (threadIdx.y + (((j + k) - 1) * blockDim.y * (height + 1))) + (blockDim.y * ((i + m) - 1))] + similarityScore); // H(i-1, j-1) + sbt(Sa[i], Sb[j])

                    maxScore = max(maxScore, score); // H[i, j]

               	    scoringMatrix[blockOffset + (threadIdx.y + ((j + k) * blockDim.y * (height + 1))) + (blockDim.y * (i + m))] = score;
                }
            }
        }
    }

    maxScoreList[scoreOffset + yIndex] = maxScore;
}

__global__ void f_scoreSequenceTiledCoalesced(short* subject, short* scoringMatrix, short* maxScoreList, int height /*querySequence.length()*/, int scoreOffset) {

    __shared__ int left_tile[TILE_SIZE];

    //register int xIndex = threadIdx.x + blockIdx.x * blockDim.x;
    register int yIndex = threadIdx.y + blockIdx.y * blockDim.y;

    // Use map for different offsets (Change the width)
    unsigned int width = constSubjectLengths[blockIdx.y];
    unsigned int subjectOffset = constSubjectOffsets[blockIdx.y];
    unsigned int blockOffset = constScoringOffsets[blockIdx.y];

    for (int i = 0; i < (height + 1); ++i) {
        scoringMatrix[blockOffset + (threadIdx.y + (blockDim.y * (i)))] = 0;
    }

    for (int j = 0; j < (width + 1); ++j) {
        scoringMatrix[blockOffset + (threadIdx.y + ((j) * blockDim.y * (height + 1)))] = 0;
    }

    int maxScore = 0;
    int up_data, diagonal_nosim_data, diagonal_sim_data = 0;
    int similarityScore, score;
    for (int i = 1; i < (height + 1); i += TILE_SIZE) {
        // set all values in left_tile to 0
        for (int p = 0; p < TILE_SIZE; p++)
            left_tile[p] = 0;

        for (int j = 1; j < (width + 1); j += TILE_SIZE) {
            for (int k = 0; k < TILE_SIZE; k++) {
		// load up and diagonal for the first access
		up_data = (int)scoringMatrix[blockOffset + (threadIdx.y + ((j + k) * blockDim.y * (height + 1))) + (blockDim.y * (i - 1))]; // E[]

		diagonal_nosim_data = (int)scoringMatrix[blockOffset + (threadIdx.y + (((j + k) - 1) * blockDim.y * (height + 1))) + (blockDim.y * (i - 1))];

                for (int m = 0; m < TILE_SIZE; m++) {
                    int left_data;
                    score = 0;

		    // calculate diagonal_sim_data
		    similarityScore = constSubstitutionMatrix[((int)constQuery[(i + m) - 1] * 25) + (int)subject[subjectOffset + threadIdx.y + (((j + k) - 1) * blockDim.y)]];
		    diagonal_sim_data = diagonal_nosim_data + similarityScore;

		    left_data = left_tile[m];

		    // calculate the new score for cell H[i+m, j+k]
		    score = max(max(max(score, left_data - GAP_PENALTY), up_data - GAP_PENALTY), diagonal_sim_data);
		    
                    maxScore = max(maxScore, score); // H[i, j]

                    left_tile[m] = score;
               	    //scoringMatrix[blockOffset + (threadIdx.y + ((j + k) * blockDim.y * (height + 1))) + (blockDim.y * (i + m))] = score;

		    // set next up_data to H value and next diagonal_data to left_data
		    up_data = score;
		    diagonal_nosim_data = left_data;

                }
                scoringMatrix[blockOffset + (threadIdx.y + ((j + k) * blockDim.y * (height + 1))) + (blockDim.y * (i + 7))] = score;
            }
        }
    }

    maxScoreList[scoreOffset + yIndex] = maxScore;
}

void smith_waterman_cuda(FASTAQuery &query, FASTADatabase &db, vector<seqid_score> &scores) {
    string querySequence = query.get_buffer();
    while (querySequence.length() % TILE_SIZE != 0) // pad to nearest 8
        querySequence = querySequence + "/";

    short* d_input_query = new short[querySequence.length()];
    unsigned int* subject_lengths = new unsigned int[CONSTANT_SIZES];
    unsigned int* subject_offsets = new unsigned int[CONSTANT_SIZES];
    unsigned int* scoring_offsets = new unsigned int[CONSTANT_SIZES];

    // TODO: Should probably put error checking and cleanup here.

    int paddedSubjects = ceil(db.numSubjects / BLOCK_Y_DIM) * BLOCK_Y_DIM;

    short* d_input_subject;
    //cudaMallocManaged((void**) &d_input_subject, (db.largestSubjectLength * paddedSubjects) * sizeof(short));
    cudaMallocManaged((void**) &d_input_subject, 200000000 * sizeof(short));

    short* d_output_max_score;
    cudaMallocManaged((void**) &d_output_max_score, paddedSubjects * sizeof(short));

    short* d_output_scoring;
    cudaMalloc((void**) &d_output_scoring, 3720000000);

    // Convert string to float representation (can't really use strings on the GPU)
    for (int i = 0; i < querySequence.length();i++) { // Pad to nearest 8 eventually here
        d_input_query[i] = convertStringToFloat(querySequence[i]);
    }

    //cout << cudaGetErrorString(cudaGetLastError()) << endl; 

    // Load in the constant memory
    cudaMemcpyToSymbol(constQuery, d_input_query, sizeof(short)*querySequence.length());
    cudaMemcpyToSymbol(constSubstitutionMatrix, blosum50, sizeof(short)*625);

    int blockPop = 0;
    int blockNum = 1;
    subject_offsets[0] = 0;
    scoring_offsets[0] = 0;

    int resultOffset = 0;

    int blockWidth = 0;
    for (map<int, vector<subject_sequence> >::reverse_iterator it = db.parsedDB.rbegin(); it != db.parsedDB.rend(); ++it) {

        blockWidth = max(blockWidth, it->first);

        for (int i = 0; i < it->second.size(); ++i) {
            for (int j = 0; j < blockWidth; j++) {
                if (j < it->second[i].sequence.length()) {
                    d_input_subject[subject_offsets[blockNum - 1] + (j * (int)BLOCK_Y_DIM) + blockPop] = convertStringToFloat(it->second[i].sequence[j]);
                }
                else d_input_subject[subject_offsets[blockNum - 1] + (j * (int)BLOCK_Y_DIM) + blockPop] = STAR;
            }

            blockPop++;

            if (blockPop >= BLOCK_Y_DIM) {
                subject_lengths[blockNum - 1] = blockWidth;
                subject_offsets[blockNum] = subject_offsets[blockNum - 1] + ((BLOCK_Y_DIM)*(blockWidth));
                scoring_offsets[blockNum] = scoring_offsets[blockNum - 1] + ((BLOCK_Y_DIM)*(blockWidth + 1)*(querySequence.length() + 1));

                blockPop = 0;
                blockWidth = it->first;

                // If are going to exceed our resources we need to run a kernel and clean up
                if ((subject_offsets[blockNum] * sizeof(short)) > CPU_MEM_THRESH || 
                        ((scoring_offsets[blockNum] + scoring_offsets[blockNum] - scoring_offsets[blockNum- 1])  * sizeof(short)) > GPU_MEM_THRESH || (blockNum + 1) >= 4000) {

                    // Load in the constant memory
                    cudaMemcpyToSymbol(constSubjectLengths, subject_lengths, sizeof(unsigned int)*CONSTANT_SIZES);
                    cudaMemcpyToSymbol(constSubjectOffsets, subject_offsets, sizeof(unsigned int)*CONSTANT_SIZES);
                    cudaMemcpyToSymbol(constScoringOffsets, scoring_offsets, sizeof(unsigned int)*CONSTANT_SIZES);

                    // Call GPU
                    int grid_y_dim = blockNum;

                    dim3 block(1, BLOCK_Y_DIM);
                    dim3 grid(1, grid_y_dim);

                    f_scoreSequenceTiledCoalesced<<<grid, block>>>(d_input_subject, d_output_scoring, d_output_max_score, querySequence.length(), resultOffset);
                    resultOffset = resultOffset + (blockNum) * BLOCK_Y_DIM;

                    cudaDeviceSynchronize();

                    subject_offsets[0] = 0;
                    scoring_offsets[0] = 0;
                    blockNum = 0;
                }

                blockNum++;
            }
        }
    }

    // Need to fill in last block
    if (blockPop != 0) {
        subject_lengths[blockNum - 1] = blockWidth;
        subject_offsets[blockNum] = subject_offsets[blockNum - 1] + ((BLOCK_Y_DIM)*(blockWidth));
        scoring_offsets[blockNum] = scoring_offsets[blockNum - 1] + ((BLOCK_Y_DIM)*(blockWidth + 1)*(querySequence.length() + 1));
    }

    // Load in the constant memory
    cudaMemcpyToSymbol(constSubjectLengths, subject_lengths, sizeof(unsigned int)*CONSTANT_SIZES);
    cudaMemcpyToSymbol(constSubjectOffsets, subject_offsets, sizeof(unsigned int)*CONSTANT_SIZES);
    cudaMemcpyToSymbol(constScoringOffsets, scoring_offsets, sizeof(unsigned int)*CONSTANT_SIZES);

    // Call GPU
    int grid_y_dim = blockNum;

    dim3 block(1, BLOCK_Y_DIM);
    dim3 grid(1, grid_y_dim);

    f_scoreSequenceTiledCoalesced<<<grid, block>>>(d_input_subject, d_output_scoring, d_output_max_score, querySequence.length(), resultOffset);

    cudaDeviceSynchronize();

    int subject = 0;
    map<int, vector<subject_sequence> >::reverse_iterator it = db.parsedDB.rbegin();
    for (; it != db.parsedDB.rend(); ++it) {
        for (int i = 0; i < it->second.size(); ++i) { 
            scores.push_back(make_pair(it->second[i].id, d_output_max_score[subject])); // change this
            subject++;
        }
    }

    delete[] d_input_query;
    delete[] subject_lengths;
    delete[] subject_offsets;
    delete[] scoring_offsets;

    // Free device memory
    cudaFree(d_input_subject);
    cudaFree(d_output_scoring);
    cudaFree(d_output_max_score);

    cudaDeviceReset(); // Comment for performance later
    return;
}

