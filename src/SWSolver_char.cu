#include <iostream>
#include <algorithm>

#include "FASTAParsers.h"
#include "SWSolver_char.h"

#define GAP_PENALTY 2
// define affine penalty ?

#define FROM_LEFT 1
#define FROM_TOP 2
#define FROM_TOP_LEFT 3

#define MAX_BLOCK_SIZE 1024
#define MAX_GRID_DIM 65535

#define BLOCK_Y_DIM 32.0

// first is sequence ID, second is max score


int blosum50_alpha[25][25] = {
//        A  B  C  D  E  F  G  H  I  J  K  L  M  N  P  Q  R  S  T  V  W  X  Y  Z  *
/* A */ { 5,-2,-1,-2,-1,-3, 0,-2,-1,-2,-1,-2,-1,-1,-1,-1,-2, 1, 0, 0,-3,-1,-2,-1,-5 },
/* B */ {-2, 6,-3, 6, 1,-4,-1, 0,-4,-4, 0,-4,-3, 5,-2, 0,-1, 0, 0,-3,-5,-1,-3, 1,-5 },
/* C */ {-1,-3,13,-4,-3,-2,-3,-3,-2,-2,-3,-2,-2,-2,-4,-3,-4,-1,-1,-1,-5,-1,-3,-3,-5 },
/* D */ {-2, 6,-4, 8, 2,-5,-1,-1,-4,-4,-1,-4,-4, 2,-1, 0,-2, 0,-1,-4,-5,-1,-3, 1,-5 },
/* E */ {-1, 1,-3, 2, 6,-3,-3, 0,-4,-3, 1,-3,-2, 0,-1, 2, 0,-1,-1,-3,-3,-1,-2, 5,-5 }, 
/* F */ {-3,-4,-2,-5,-3, 8,-4,-1, 0, 1,-4, 1, 0,-4,-4,-4,-3,-3,-2,-1, 1,-1, 4,-4,-5 },  
/* G */ { 0,-1,-3,-1,-3,-4, 8,-2,-4,-4,-2,-4,-3, 0,-2,-2,-3, 0,-2,-4,-3,-1,-3,-2,-5 }, 
/* H */ {-2, 0,-3,-1, 0,-1,-2,10,-4,-3, 0,-3,-1, 1,-2, 1, 0,-1,-2,-4,-3,-1, 2, 0,-5 },
/* I */ {-1,-4,-2,-4,-4, 0,-4,-4, 5, 4,-3, 2, 2,-3,-3,-3,-4,-3,-1, 4,-3,-1,-1,-3,-5 },
/* J */ {-2,-4,-2,-4,-3, 1,-4,-3, 4, 4,-3, 4, 2,-4,-3,-3,-3,-3,-1, 2,-2,-1,-1,-3,-5 },
/* K */ {-1, 0,-3,-1, 1,-4,-2, 0,-3,-3, 6,-3,-2, 0,-1, 2, 3, 0,-1,-3,-3,-1,-2, 1,-5 },
/* L */ {-2,-4,-2,-4,-3, 1,-4,-3, 2, 4,-3, 5, 3,-4,-4,-2,-3,-3,-1, 1, 2,-1,-1,-3,-5 },
/* M */ {-1,-3,-2,-4,-2, 0,-3,-1, 2, 2,-2, 3, 7,-2,-3, 0,-2,-2,-1, 1,-1,-1, 0,-1,-5 },
/* N */ {-1, 5,-2, 2, 0,-4, 0, 1,-3,-4, 0,-4,-2, 7,-2, 0,-1, 1, 0,-3,-4,-1,-2, 0,-5 },
/* P */ {-1,-2,-4,-1,-1,-4,-2,-2,-3,-3,-1,-4,-3,-2,10,-1,-3,-1,-1,-3,-4,-1,-3,-1,-5 },
/* Q */ {-1, 0,-3, 0, 2,-4,-2, 1,-3,-3, 2,-2, 0, 0,-1, 7, 1, 0,-1,-3,-1,-1,-1, 4,-5 },
/* R */ {-2,-1,-4,-2, 0,-3,-3, 0,-4,-3, 3,-3,-2,-1,-3, 1, 7,-1,-1,-3,-3,-1,-1, 0,-5 },
/* S */ { 1, 0,-1, 0,-1,-3, 0,-1,-3,-3, 0,-3,-2, 1,-1, 0,-1, 5, 2,-2,-4,-1,-2, 0,-5 },
/* T */ { 0, 0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,-1, 0,-1,-1,-1, 2, 5, 0,-3,-1,-2,-1,-5 },
/* V */ { 0,-3,-1,-4,-3,-1,-4,-4, 4, 2,-3, 1, 1,-3,-3,-3,-3,-2, 0, 5,-3,-1,-1,-3,-5 },
/* W */ {-3,-5,-5,-5,-3, 1,-3,-3,-3,-2,-3,-2,-1,-4,-4,-1,-3,-4,-3,-3,15,-1, 2,-2,-5 },
/* X */ {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-5 },
/* Y */ {-2,-3,-3,-3,-2, 4,-3, 2,-1,-1,-2,-1, 0,-2,-3,-1,-1,-2,-2,-1, 2,-1, 8,-2,-5 }, 
/* Z */ {-1, 1,-3, 1, 5,-4,-2, 0,-3,-3, 1,-3,-1, 0,-1, 4, 0, 0,-1,-3,-2,-1,-2, 5,-5 },
/* * */ {-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5, 1 }
};

using namespace std;

__constant__ char constQuery[1024];
__constant__ int constSubstitutionMatrix[625];

char convertStringToChar(char character) {
    switch(character)
    {
        case 'A': { return 'A'; }
        case 'R': { return 'R'; }
        case 'N': { return 'N'; }
        case 'D': { return 'D'; }
        case 'C': { return 'C'; }
        case 'Q': { return 'Q'; }
        case 'E': { return 'E'; }
        case 'G': { return 'G'; }
        case 'H': { return 'H'; }
        case 'I': { return 'I'; }
        case 'L': { return 'L'; }
        case 'K': { return 'K'; }
        case 'M': { return 'M'; }
        case 'F': { return 'F'; }
        case 'P': { return 'P'; }
        case 'S': { return 'S'; }
        case 'T': { return 'T'; }
        case 'W': { return 'W'; }
        case 'Y': { return 'Y'; }
        case 'V': { return 'V'; }
        case 'B': { return 'B'; }
        case 'J': { return 'J'; }
        case 'Z': { return 'Z'; }
        case 'X': { return 'X'; }
    }
    return '*';
}

// Kernel function for computing the scoring matrix of a sequence
__global__ void f_scoreSequence(char* subject, float* scoringMatrix, float* maxScoreList, int width, int height, int numSubjects) {

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
				
				int similarityScore = 0;
				int query = (int)constQuery[i-1];
				int subjectValue = (int)subject[width*yIndex + j -1];
				int searchMatrix;
				
				// first determine how to search the constSubstitutionMatrix
				// have both in V-*, 1 in V-Z 1 in P-T (vice versa flip them), 1 in V-Z 1 in A-N (vice versa flip them),
				// both in P-T, 1 in P-T 1 in A-N (vice versa flip them), both in A-N
				// 6 possibilities after collapsing vice versas
				if ((query > 85 || query == 42) && (subjectValue > 85 || subjectValue == 42)) // both in V-*
					searchMatrix = 1;
				else if ((query > 85 || query == 42) && (subjectValue > 79)) // query in V-*, subjectValue in P-T
					searchMatrix = 2;
				else if ((query > 79) && (subjectValue > 85 || subjectValue == 42)) {// query in P-T, subjectValue in V-*
					int temp = query;
					query = subjectValue;
					subjectValue = temp;
					searchMatrix = 2;
				}
				else if (query > 85 || query == 42) // query in V-*, subjectValue in A-N
					searchMatrix = 3;
				else if (subjectValue > 85 || subjectValue == 42) { // query in A-N, subjectValue in V-*
					int temp = query;
					query = subjectValue;
					subjectValue = temp;
					searchMatrix = 3;
				}

				else if (query > 79 && subjectValue > 79) // both in P-T
					searchMatrix = 4;
				else if (query > 79) // query in P-T, subjectValue in A-N
					searchMatrix = 5;
				else if (subjectValue > 79) { // query in A-N, subjectValue in P-T
					int temp = query;
					query = subjectValue;
					subjectValue = temp;
					searchMatrix = 5;
				}

				else // both in A-N
					searchMatrix = 6;

				// based on the searchMatrix value, use a switch case and calculate the similarityScore
				// if value is * (42), result is either -5 (not similar) or 1 (similar)
				switch(searchMatrix)
				{
					case 1: {
						if (query == 42 && subjectValue == 42) similarityScore = 1;
						else if ((query == 42 && subjectValue != 42) || (query != 42 && subjectValue == 42)) similarityScore = -5;
						else similarityScore = constSubstitutionMatrix[((query-67) * 25) + subjectValue-67]; 
						break; }

					case 2: {
						if (query == 42) similarityScore = -5;
						else similarityScore = constSubstitutionMatrix[((query-67) * 25) + subjectValue-66];
						break; }
				
					case 3: {
						if (query == 42) similarityScore = -5;
						else similarityScore = constSubstitutionMatrix[((query-67) * 25) + subjectValue-65];
						break; }
			
					case 4: {
						similarityScore = constSubstitutionMatrix[((query-66) * 25) + subjectValue-66];
						break; }

					case 5: {
						similarityScore = constSubstitutionMatrix[((query-66) * 25) + subjectValue-65];
						break; }
			
					case 6: {
						similarityScore = constSubstitutionMatrix[((query-65) * 25) + subjectValue-65];
						break; }
				}
					
                //int similarityScore = constSubstitutionMatrix[((int)constQuery[i - 1] * 25) + (int)subject[width*yIndex + j - 1]];
                score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + ((i - 1) * (width + 1)) + j - 1] + similarityScore);

                maxScore = max(maxScore, score);

                scoringMatrix[(width + 1)*(height + 1)*yIndex + (i * (width + 1)) + j] = score;
            }
        }
        maxScoreList[yIndex] = maxScore;
    }
}

vector<seqid_score> smith_waterman_cuda_char(FASTAQuery &query, FASTADatabase &db) {
    string querySequenceBuffer = query.get_buffer();
	char querySequence[querySequenceBuffer.length() + 1];
	copy(querySequenceBuffer.begin(), querySequenceBuffer.end(), querySequence);
	querySequence[querySequenceBuffer.length()] = '\0';
    vector<seqid_score> scores;

    char* d_input_subject;
    cudaMallocManaged((void**) &d_input_subject, (db.largestSubjectLength * db.numSubjects) * sizeof(char));
	/*
	memcpy(d_input_subject, subjectSequences[0], ((largestSubjectLength * numSubjects) + 1) * sizeof(char));
	
	for (int i = 0; i < numSubjects; i++) {
		strcat(d_input_subject, subjectSequences[i]);
	}
	*/

	// Set up offsets 
    int grid_y_dim = ceil(db.numSubjects / BLOCK_Y_DIM);
    
    char* d_input_offsets;
    cudaMallocManaged((void**) &d_input_offsets, grid_y_dim * sizeof(char));

    float* d_output_scoring;
    cudaMallocManaged((void**) &d_output_scoring, ((strlen(querySequence) + 1) * (db.largestSubjectLength + 1) * db.numSubjects) * sizeof(float));

	float* d_output_max_score;
    cudaMallocManaged((void**) &d_output_max_score, db.numSubjects * sizeof(float));

	for (int i = 0; i < db.numSubjects; i++) {
        for (int j = 0; j < db.largestSubjectLength; j++) { // Will need to pad here
            if (j < db.subjectSequences[i].sequence.length()) {
                d_input_subject[i*db.largestSubjectLength + j] = convertStringToChar(db.subjectSequences[i].sequence[j]);
            }
            else d_input_subject[i*db.largestSubjectLength + j] = '*';
        }
    }

	cudaMemcpyToSymbol(constQuery, querySequence, sizeof(char)*strlen(querySequence));
	cudaMemcpyToSymbol(constSubstitutionMatrix, blosum50_alpha, sizeof(int)*625);

    // Call GPU
    dim3 block(1, BLOCK_Y_DIM);
    dim3 grid(1, grid_y_dim);
 
    f_scoreSequence<<<grid, block>>>(d_input_subject, d_output_scoring, d_output_max_score, db.largestSubjectLength, strlen(querySequence), db.numSubjects);

    cudaDeviceSynchronize();

    cout << endl << db.subjectSequences[0].sequence << endl;
    /*
    // Print results for 1 subject query
    for (int subject = 0; subject < db.numSubjects; subject++) {
        string seqA = querySequence;
        string seqB = db.subjectSequences[subject].sequence;

        cout << "    ";
        for (int j = 0; j < (seqB.length() + 1); j++) {
            cout << seqB[j] << " ";
        }
        cout << endl;

        for (int i = 0; i < (seqA.length() + 1); i++) {
            if (i != 0) cout << seqA[i - 1] << " ";
            else cout << "  ";
            for (int j = 0; j < (seqB.length() + 1); j++) {
                cout << d_output_scoring[((db.largestSubjectLength + 1) * (querySequence.length() + 1) * subject) + (i * (seqB.length() + 1)) + j] << " ";
            }
            cout << endl;
        }
    }*/
    
    // Print results for 1 subject query
	
	
    for (int subject = 0; subject < db.numSubjects; subject++) {
        scores.push_back(make_pair(db.subjectSequences[subject].id, d_output_max_score[subject])); // change this
    }

    // Free device memory
    cudaFree(d_input_subject);
    cudaFree(d_input_offsets);
    cudaFree(d_output_scoring);
    cudaFree(d_output_max_score);
    cudaDeviceReset();
	
    return scores;
}

