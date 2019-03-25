#include <sys/time.h>
#include <iostream>
#include <algorithm>
#include "boost/program_options.hpp"

#include "FASTAParsers.h"

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
namespace po = boost::program_options;

__constant__ float constQuery[1024];

// Time stamp function
double getTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_usec/1000000 + tv.tv_sec;
}

// Kernel function for computing the scoring matrix of a sequence
__global__ void f_scoreSequence(float* subject, float* scoringMatrix, float* maxScoreList, 
                int width /*largestSubjectLength*/, int height /*querySequence.length()*/, int numSubjects) {
    
    int substitutionMatrix[2] = {SEQ_EQUAL, SEQ_DIFF};

    //register int xIndex = threadIdx.x + blockIdx.x * blockDim.x;
    register int yIndex = threadIdx.y + blockIdx.y * blockDim.y;

    float maxScore = 0;
        if (yIndex < numSubjects) {
        for (int i = 1; i < (height + 1); i++) {
            for (int j = 1; j < (width + 1); j++) {
                float score = 0;

                score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + (i * (width + 1)) + j - 1] - GAP_PENALTY);
                score = max(score, scoringMatrix[(width + 1)*(height + 1)*yIndex + ((i - 1) * (width + 1)) + j] - GAP_PENALTY);

                int similarityScore = 0;

                if (constQuery[i - 1] == subject[width*yIndex + j - 1]) similarityScore = substitutionMatrix[0];
                else similarityScore = substitutionMatrix[1];

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
    
    int substitutionMatrix[2] = {SEQ_EQUAL, SEQ_DIFF};

    //register int xIndex = threadIdx.x + blockIdx.x * blockDim.x;
    register int yIndex = threadIdx.y + blockIdx.y * blockDim.y;
    
    // Use map for different offsets (Change the width)
    int blockOffset = (blockIdx.y * blockDim.y)*(width + 1)*(height + 1);

    float maxScore = 0;
        if (yIndex < numSubjects) {
        for (int i = 1; i < (height + 1); i++) {
            for (int j = 1; j < (width + 1); j++) {
                float score = 0;

                score = max(score, scoringMatrix[blockOffset + (threadIdx.y + ((j - 1) * blockDim.y * (height + 1))) + (blockDim.y * i)] - GAP_PENALTY);
                score = max(score, scoringMatrix[blockOffset + (threadIdx.y + (j * blockDim.y * (height + 1))) + (blockDim.y * (i - 1))] - GAP_PENALTY);

                int similarityScore = 0;

                if (constQuery[i - 1] == subject[threadIdx.y + ((j - 1) * blockDim.y)]) similarityScore = substitutionMatrix[0];
                else similarityScore = substitutionMatrix[1];

                score = max(score, scoringMatrix[blockOffset + (threadIdx.y + ((j - 1) * blockDim.y * (height + 1))) + (blockDim.y * (i - 1))] + similarityScore);

                maxScore = max(maxScore, score);
                
                scoringMatrix[blockOffset + (threadIdx.y + (j * blockDim.y * (height + 1))) + (blockDim.y * i)] = score;
            }
        }
        maxScoreList[yIndex] = maxScore;
    }
}

int main( int argc, char *argv[] ) {
    double time_start = getTimeStamp();

    po::options_description desc("Smith-Waterman CUDA Usage");
    po::variables_map vm;
    
    try {
        desc.add_options()
            ("help", "Display this help message")
            ("query", po::value<std::string>()->required(),"Path to query file (required)")
            ("db", po::value<std::string>()->required(), "Path to database file (required)");

        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if(vm.count("help") || argc <= 1){
            cout << desc;
            return 1;
        }
    } catch (const po::required_option & e) {
        cout << desc;
        return 1;
    }

    std::string querypath = vm["query"].as<std::string>();
    FASTAQuery query(querypath, true);
    cout << "Input buffer:";
    query.print_buffer();
    cout << endl;
    string querySequence = query.get_buffer();

    // Parse database file
    std::string datapath = vm["db"].as<std::string>();
    FASTADatabase db(datapath);


    /*
    for (map<int, vector<subject_sequence> >::iterator it = parsedDB.begin(); it != parsedDB.end(); ++it) {
        cout << it->first 
             << ":"
             << it->second.size()
             << endl;
    } */
    
    cout << "Largest subject: " << db.largestSubjectLength << endl;
    cout << "Num subjects: " << db.numSubjects << endl;
    cout << "Accumulated db length: " << db.subjectLengthSum << endl;

    // alloc memory on GPU
    float* d_input_query = new float[querySequence.length()];
    memset(d_input_query, 0, sizeof(float) * querySequence.length());

    float* d_input_subject;
    cudaMallocManaged((void**) &d_input_subject, (db.largestSubjectLength * db.numSubjects) * sizeof(float));

    float* d_output_scoring;
    cudaMallocManaged((void**) &d_output_scoring, ((querySequence.length() + 1) *
                (db.largestSubjectLength + 1) * db.numSubjects) * sizeof(float));
    
    float* d_output_max_score;
    cudaMallocManaged((void**) &d_output_max_score, db.numSubjects * sizeof(float));

    // Convert string to float representation (can't really use strings on the GPU)
    for (int i = 0; i < querySequence.length();i++) { // Pad to nearest 8 eventually here
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

    for (int i = 0; i < db.numSubjects; i++) {
        for (int j = 0; j < db.largestSubjectLength; j++) { // Will need to pad here
            if (j < db.subjectSequences[i].length()) {
                switch(db.subjectSequences[i][j])
                {
                    case 'A': { d_input_subject[i*db.largestSubjectLength + j] = A;
                                break;
                            }
                    case 'G': { d_input_subject[i*db.largestSubjectLength + j] = G;
                                break;
                            }
                    case 'C': { d_input_subject[i*db.largestSubjectLength + j] = C;
                                break;
                            }
                    case 'T': { d_input_subject[i*db.largestSubjectLength + j] = T;
                                break;
                            }
                }
            }
        }
    }

    cudaMemcpyToSymbol(constQuery, d_input_query, sizeof(float)*querySequence.length());

    int grid_y_dim = ceil(db.numSubjects / 32.0);
    
    // Call GPU
    dim3 block(1, 32);
    dim3 grid(1, grid_y_dim);
    
    f_scoreSequence<<<grid, block>>>(d_input_subject, d_output_scoring, d_output_max_score, db.largestSubjectLength, querySequence.length(), db.numSubjects);

    cudaDeviceSynchronize();

    /*
    // Print results for 1 subject query
    for (int subject = 0; subject < numSubjects; subject++) {
        string seqA = querySequence;
        string seqB = subjectSequences[subject];

        cout << "    ";
        for (int j = 0; j < (seqB.length() + 1); j++) {
            cout << seqB[j] << " ";
        }
        cout << endl;

        for (int i = 0; i < (seqA.length() + 1); i++) {
            if (i != 0) cout << seqA[i - 1] << " ";
            else cout << "  ";
            for (int j = 0; j < (seqB.length() + 1); j++) {
                cout << d_output_scoring[((largestSubjectLength + 1) * (querySequence.length() + 1) * subject) + (i * (seqB.length() + 1)) + j] << " ";
            }
            cout << endl;
        }
    }
    */
    
    // Print results for 1 subject query
    for (int subject = 0; subject < db.numSubjects; subject++) {
        //cout << d_output_max_score[subject] << endl;
    }

    double time_end = getTimeStamp();
    double seconds_elapsed = time_end - time_start;

    cout << std::string(80, '=') << endl;
    cout << "METRICS:" << endl;
    cout << "Query length: " << querySequence.length() << " chars." << endl;
    cout << "Sum of DB length: " << db.subjectLengthSum << " chars." << endl;
    cout << "Time elapsed: " << seconds_elapsed << " seconds." << endl;
    cout << "Performance: " << 1E-9 * (querySequence.length() * db.subjectLengthSum)
            / seconds_elapsed << " GCUPS." << endl;

    delete[] d_input_query;

    // Free device memory
    cudaFree(d_input_query);
    cudaFree(d_input_subject);
    cudaFree(d_output_scoring);
    cudaFree(d_output_max_score);
    cudaDeviceReset();
}
