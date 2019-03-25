#include <sys/time.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <math.h>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include "boost/program_options.hpp"

#define GAP_PENALTY 2
// define affine penalty ?

#define FROM_LEFT 1
#define FROM_TOP 2
#define FROM_TOP_LEFT 3

#define MAX_BLOCK_SIZE 1024
#define MAX_GRID_DIM 65535

#define LENGTH_THRESHOLD 100

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

int blosum50[25][25] = {
//        A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
/* A */ { 5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-2,-1,-1,-3,-1, 1, 0,-3,-2, 0,-2,-2,-1,-1,-5},
/* R */ {-2, 7,-1,-2,-4, 1, 0,-3, 0,-4,-3, 3,-2,-3,-3,-1,-1,-3,-1,-3,-1,-3, 0,-1,-5},
/* N */ {-1,-1, 7, 2,-2, 0, 0, 0, 1,-3,-4, 0,-2,-4,-2, 1, 0,-4,-2,-3, 5,-4, 0,-1,-5},
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
/* * */ {-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5 }
};

using namespace std;
namespace po = boost::program_options;

__constant__ float constQuery[1024];
__constant__ int constSubstitutionMatrix[625];

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

struct subject_sequence {
    int id;
    string sequence;
};

// Time stamp function
double getTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_usec/1000000 + tv.tv_sec;
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
    
    int substitutionMatrix[2] = {3, -3};

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

class ParsedFASTA {
private:
    bool isQuery;
    stringstream header;
    string buffer;
public:
    ParsedFASTA(std::string filepath, bool _isQuery) {
        isQuery = _isQuery;

        ifstream filestream;
        filestream.open(filepath.c_str());

        buffer.reserve(10000); // optimization -- reserve some arbitrary length
        string tmp;
        getline(filestream, tmp); // Skip first line
        while (getline(filestream, tmp)) {
            buffer.append(tmp);
        }
        filestream.close();
    };

    ~ParsedFASTA() {
    };

    void print_buffer() {
        cout << buffer << endl;
    };

    string get_buffer() {
        return buffer;
    };
};

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
    ParsedFASTA query(querypath, true);
    cout << "Input buffer:";
    query.print_buffer();
    cout << endl;
    string querySequence = query.get_buffer();

    // Parse database file
    ifstream databaseFile;
    std::string datapath = vm["db"].as<std::string>();
    databaseFile.open(datapath.c_str());

    int subjectLengthSum = 0;

    string temp;

    // key is sequence length, value is a vector of subject_sequence struct
    map<int, vector<subject_sequence> > parsedDB;

    vector<string> subjectSequences;
    string subjectSequence = "";
    int largestSubjectLength = 0;
    int numSubjects = 0;
    bool isFirst = true;

    subject_sequence tmp;

    int _id = 0;

    while (getline(databaseFile, temp)) {
        
        // This line denotes the start of a sequence
        if (temp[0] == '>') {
            if (!isFirst) {
                if (subjectSequence.length() <= LENGTH_THRESHOLD) {
                    tmp.id = _id++;
                    tmp.sequence = subjectSequence;
                    parsedDB[subjectSequence.length()].push_back(tmp);

                    subjectSequences.push_back(subjectSequence);
                    subjectLengthSum += subjectSequence.length();
                    largestSubjectLength = max(largestSubjectLength, (int)subjectSequence.length());
                    
                    numSubjects++;
                }
            }
            isFirst = false;
            
            //cout << subjectSequence << endl;
            subjectSequence = "";
        }
        else {
            subjectSequence += temp;
        }
        
    }
    // Adding last sequence 
    if (subjectSequence.length() <= LENGTH_THRESHOLD) {
        tmp.id = _id++;
        tmp.sequence = subjectSequence;
        parsedDB[subjectSequence.length()].push_back(tmp);

        subjectSequences.push_back(subjectSequence);
        subjectLengthSum += subjectSequence.length();
        largestSubjectLength = max(largestSubjectLength, (int)subjectSequence.length());
        
        numSubjects++;
    }


    databaseFile.close();

    /*
    for (map<int, vector<subject_sequence> >::iterator it = parsedDB.begin(); it != parsedDB.end(); ++it) {
        cout << it->first 
             << ":"
             << it->second.size()
             << endl;
    } */

    // alloc memory on GPU
    float* d_input_query = new float[querySequence.length()];
    memset(d_input_query, 0, sizeof(float) * querySequence.length());

    float* d_input_subject;
    cudaMallocManaged((void**) &d_input_subject, (largestSubjectLength * numSubjects) * sizeof(float));
    
    // Set up offsets 
    int grid_y_dim = ceil(numSubjects / BLOCK_Y_DIM);
    
    float* d_input_offsets;
    cudaMallocManaged((void**) &d_input_offsets, grid_y_dim * sizeof(float));

    float* d_output_scoring;
    cudaMallocManaged((void**) &d_output_scoring, ((querySequence.length() + 1) * (largestSubjectLength + 1) * numSubjects) * sizeof(float));
    
    float* d_output_max_score;
    cudaMallocManaged((void**) &d_output_max_score, numSubjects * sizeof(float));

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
    for (int i = 0; i < numSubjects; i++) {
        for (int j = 0; j < largestSubjectLength; j++) { // Will need to pad here
            if (j < subjectSequences[i].length()) {
                d_input_subject[i*largestSubjectLength + j] = convertStringToFloat(subjectSequences[i][j]);
            }
            else d_input_subject[i*largestSubjectLength + j] = STAR;
        }
    }

    cudaMemcpyToSymbol(constQuery, d_input_query, sizeof(float)*querySequence.length());
    cudaMemcpyToSymbol(constSubstitutionMatrix, blosum50, sizeof(int)*625);
    
    // Call GPU
    dim3 block(1, BLOCK_Y_DIM);
    dim3 grid(1, grid_y_dim);
    
    f_scoreSequence<<<grid, block>>>(d_input_subject, d_output_scoring, d_output_max_score, largestSubjectLength, querySequence.length(), numSubjects);

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
    for (int subject = 0; subject < numSubjects; subject++) {
        cout << d_output_max_score[subject] << endl;
    }

    double time_end = getTimeStamp();
    double seconds_elapsed = time_end - time_start;

    cout << std::string(80, '=') << endl;
    cout << "METRICS:" << endl;
    cout << "Query length: " << querySequence.length() << " chars." << endl;
    cout << "Num subjects: " << numSubjects << endl;
    cout << "Sum of DB length: " << subjectLengthSum << " chars." << endl;
    cout << "Time elapsed: " << seconds_elapsed << " seconds." << endl;
    cout << "Performance: " << 1E-9 * (querySequence.length() * subjectLengthSum)
            / seconds_elapsed << " GCUPS." << endl;

    delete[] d_input_query;

    // Free device memory
    cudaFree(d_input_query);
    cudaFree(d_input_subject);
    cudaFree(d_input_offsets);
    cudaFree(d_output_scoring);
    cudaFree(d_output_max_score);
    cudaDeviceReset();
}
