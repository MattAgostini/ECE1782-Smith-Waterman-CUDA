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

// first is sequence ID, second is max score

int blosum50[25][25] = {
//        A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
/* A */ { 5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-2,-1,-1,-3,-1, 1, 0,-3,-2, 0,-2,-2,-1,-1,-5 },
/* R */ {-2, 7,-1,-2,-4, 1, 0,-3, 0,-4,-3, 3,-2,-3,-3,-1,-1,-3,-1,-3,-1,-3, 0,-1,-5 },
/* N */ {-1,-1, 7, 2,-2, 0, 0, 0, 1,-3,-4, 0,-2,-4,-2, 1, 0,-4,-2,-3, 5,-4, 0,-1,-5 },
/* D */ {-2,-2, 2, 8,-4, 0, 2,-1,-1,-4,-4,-1,-4,-5,-1, 0,-1,-5,-3,-4, 6,-4, 1,-1,-5 },
/* C */ {-1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1,-3,-2,-3,-1,-5 }, 
/* Q */ {-1, 1, 0, 0,-3, 7, 2,-2, 1,-3,-2, 2, 0,-4,-1, 0,-1,-1,-1,-3, 0,-3, 4,-1,-5 },  
/* E */ {-1, 0, 0, 2,-3, 2, 6,-3, 0,-4,-3, 1,-2,-3,-1,-1,-1,-3,-2,-3, 1,-3, 5,-1,-5 }, 
/* G */ { 0,-3, 0,-1,-3,-2,-3, 8,-2,-4,-4,-2,-3,-4,-2, 0,-2,-3,-3,-4,-1,-4,-2,-1,-5 },
/* H */ {-2, 0, 1,-1,-3, 1, 0,-2,10,-4,-3, 0,-1,-1,-2,-1,-2,-3, 2,-4, 0,-3, 0,-1,-5 },
/* I */ {-1,-4,-3,-4,-2,-3,-4,-4,-4, 5, 2,-3, 2, 0,-3,-3,-1,-3,-1, 4,-4, 4,-3,-1,-5 },
/* L */ {-2,-3,-4,-4,-2,-2,-3,-4,-3, 2, 5,-3, 3, 1,-4,-3,-1,-2,-1, 1,-4, 4,-3,-1,-5 },
/* K */ {-1, 3, 0,-1,-3, 2, 1,-2, 0,-3,-3, 6,-2,-4,-1, 0,-1,-3,-2,-3, 0,-3, 1,-1,-5 },
/* M */ {-1,-2,-2,-4,-2, 0,-2,-3,-1, 2, 3,-2, 7, 0,-3,-2,-1,-1, 0, 1,-3, 2,-1,-1,-5 },
/* F */ {-3,-3,-4,-5,-2,-4,-3,-4,-1, 0, 1,-4, 0, 8,-4,-3,-2, 1, 4,-1,-4, 1,-4,-1,-5 },
/* P */ {-1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3,-2,-3,-1,-1,-5 },
/* S */ { 1,-1, 1, 0,-1, 0,-1, 0,-1,-3,-3, 0,-2,-3,-1, 5, 2,-4,-2,-2, 0,-3, 0,-1,-5 },
/* T */ { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 2, 5,-3,-2, 0, 0,-1,-1,-1,-5 },
/* W */ {-3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1, 1,-4,-4,-3,15, 2,-3,-5,-2,-2,-1,-5 },
/* Y */ {-2,-1,-2,-3,-3,-1,-2,-3, 2,-1,-1,-2, 0, 4,-3,-2,-2, 2, 8,-1,-3,-1,-2,-1,-5 },
/* V */ { 0,-3,-3,-4,-1,-3,-3,-4,-4, 4, 1,-3, 1,-1,-3,-2, 0,-3,-1, 5,-3, 2,-3,-1,-5 },
/* B */ {-2,-1, 5, 6,-3, 0, 1,-1, 0,-4,-4, 0,-3,-4,-2, 0, 0,-5,-3,-3, 6,-4, 1,-1,-5 },
/* J */ {-2,-3,-4,-4,-2,-3,-3,-4,-3, 4, 4,-3, 2, 1,-3,-3,-1,-2,-1, 2,-4, 4,-3,-1,-5 },
/* Z */ {-1, 0, 0, 1,-3, 4, 5,-2, 0,-3,-3, 1,-1,-4,-1, 0,-1,-2,-2,-3, 1,-3, 5,-1,-5 }, 
/* X */ {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-5 },
/* * */ {-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5, 1 }
};

using namespace std;

__constant__ float constQuery[1024];
__constant__ int constSubstitutionMatrix[625];
__constant__ int constMemoryOffsets[2048];

float convertStringToFloat(char character) {
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
        for (int i = 1; i < (height + 1); i++) {
            for (int j = 1; j < (width + 1); j++) {
                float score = 0;

                score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + (i * (width + 1)) + j - 1] - GAP_PENALTY);
                score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + ((i - 1) * (width + 1)) + j] - GAP_PENALTY);

                int similarityScore = constSubstitutionMatrix[((int)constQuery[i - 1] * 25) + (int)subject[width*yIndex + j - 1]];
                score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + ((i - 1) * (width + 1)) + j - 1] + similarityScore);

                maxScore = max(maxScore, score);

                scoringMatrix[(width + 1)*(height + 1)*yIndex + (i * (width + 1)) + j] = score;
            }
        }
        maxScoreList[yIndex] = maxScore;
    }
}

// Kernel function for computing the scoring matrix of a sequence
__global__ void f_scoreSequenceCoalesced(float* subject, float* scoringMatrix, float* maxScoreList, 
                int width /*largestSubjectLength*/, int height /*querySequence.length()*/, int numSubjects) {

    //register int xIndex = threadIdx.x + blockIdx.x * blockDim.x;
    register int yIndex = threadIdx.y + blockIdx.y * blockDim.y;
    
    // Use map for different offsets (Change the width)
    int subjectOffset = blockIdx.y * ((blockDim.y)*(width));
    int blockOffset = blockIdx.y * ((blockDim.y)*(width + 1)*(height + 1));

    float maxScore = 0;
    for (int i = 1; i < (height + 1); i++) {
        for (int j = 1; j < (width + 1); j++) {
            float score = 0;

            score = max(score, scoringMatrix[blockOffset + (threadIdx.y + ((j - 1) * blockDim.y * (height + 1))) + (blockDim.y * i)] - GAP_PENALTY);
            score = max(score, scoringMatrix[blockOffset + (threadIdx.y + (j * blockDim.y * (height + 1))) + (blockDim.y * (i - 1))] - GAP_PENALTY);

            int similarityScore = constSubstitutionMatrix[((int)constQuery[i - 1] * 25) + (int)subject[subjectOffset + threadIdx.y + ((j - 1) * blockDim.y)]];
            score = max(score, scoringMatrix[blockOffset + (threadIdx.y + ((j - 1) * blockDim.y * (height + 1))) + (blockDim.y * (i - 1))] + similarityScore);

            maxScore = max(maxScore, score);

            scoringMatrix[blockOffset + (threadIdx.y + (j * blockDim.y * (height + 1))) + (blockDim.y * i)] = score;
        }
    }
    
    maxScoreList[yIndex] = maxScore;
}

vector<seqid_score> smith_waterman_cuda(FASTAQuery &query, FASTADatabase &db) {
    string querySequence = query.get_buffer();
    vector<seqid_score> scores;
    
    // alloc memory on GPU
    float* d_input_query = new float[querySequence.length()];
    int* memory_offsets = new int[2048];
    
    // TODO: Should probably put error checking and cleanup here.
    
    int paddedSubjects = ceil(db.numSubjects / BLOCK_Y_DIM) * BLOCK_Y_DIM;
    
    float* d_input_subject;
    cudaMallocManaged((void**) &d_input_subject, (db.largestSubjectLength * paddedSubjects) * sizeof(float));
    
    // Set up offsets 
    int grid_y_dim = ceil(db.numSubjects / BLOCK_Y_DIM);
    
    float* d_output_scoring;
    cudaMalloc((void**) &d_output_scoring, ((querySequence.length() + 1) *
                (db.largestSubjectLength + 1) * paddedSubjects) * sizeof(float));
    
    float* d_output_max_score;
    cudaMallocManaged((void**) &d_output_max_score, paddedSubjects * sizeof(float));
    
    // Convert string to float representation (can't really use strings on the GPU)
    for (int i = 0; i < querySequence.length();i++) { // Pad to nearest 8 eventually here
        d_input_query[i] = convertStringToFloat(querySequence[i]);
    }
    
    /*
    int blockPop = 0;
    int blockNum = 1;
    int blockWidth = 0;
    for (map<int, vector<subject_sequence> >::reverse_iterator it = parsedDB.rbegin(); it != parsedDB.rend(); ++it) {
        blockWidth = max(blockWidth, it->first);
        for (int i = 0; i < it->second.size(); ++i) {
            if (blockPop >= BLOCK_Y_DIM) {
                blockPop = 0;
                d_input_offsets[blockNum] = d_input_offsets[blockNum - 1] + (BLOCK_Y_DIM * blockWidth); // Need to include the query length for scoring matrix
                blockNum++;
                blockWidth = it->first;
            }
            
            
            blockPop++;
        }
    }
    */
    
    int blockOffset = BLOCK_Y_DIM * db.largestSubjectLength;
    
    for (int block = 0; block < ceil(db.numSubjects / BLOCK_Y_DIM); block++) {
        for (int i = 0; i < BLOCK_Y_DIM; i++) {
            for (int j = 0; j < db.largestSubjectLength; j++) { // Will need to pad here
                if (j < db.subjectSequences[block*BLOCK_Y_DIM + i].sequence.length()) {
                    //d_input_subject[i*db.largestSubjectLength + j] = convertStringToFloat(db.subjectSequences[i].sequence[j]);
                    d_input_subject[blockOffset*block + (j * (int)BLOCK_Y_DIM) + i] = convertStringToFloat(db.subjectSequences[block*BLOCK_Y_DIM + i].sequence[j]);
                }
                else d_input_subject[blockOffset*block + (j * (int)BLOCK_Y_DIM) + i] = STAR;
            }
        }
    }
    
    // Load in the constant memory
    cudaMemcpyToSymbol(constQuery, d_input_query, sizeof(float)*querySequence.length());
    cudaMemcpyToSymbol(constSubstitutionMatrix, blosum50, sizeof(int)*625);
    cudaMemcpyToSymbol(constMemoryOffsets, memory_offsets, sizeof(int)*2048);
    
    // Call GPU
    dim3 block(1, BLOCK_Y_DIM);
    dim3 grid(1, grid_y_dim);
    
    f_scoreSequenceCoalesced<<<grid, block>>>(d_input_subject, d_output_scoring, d_output_max_score, db.largestSubjectLength, querySequence.length(), db.numSubjects);
    
    cudaDeviceSynchronize();
    
    for (int subject = 0; subject < db.numSubjects; subject++) {
        scores.push_back(make_pair(db.subjectSequences[subject].id, d_output_max_score[subject])); // change this
    }
    
    delete[] d_input_query;
    delete[] memory_offsets;
    
    // Free device memory
    cudaFree(d_input_subject);
    cudaFree(d_output_scoring);
    cudaFree(d_output_max_score);
    cudaDeviceReset();
    
    return scores;
}

