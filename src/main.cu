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

#define SEQ_EQUAL 3
#define SEQ_DIFF -3
#define GAP_PENALTY 2
// define affine penalty ?

#define FROM_LEFT 1
#define FROM_TOP 2
#define FROM_TOP_LEFT 3

#define MAX_BLOCK_SIZE 1024
#define MAX_GRID_DIM 65535

#define LENGTH_THRESHOLD 50

#define A 1
#define G 2
#define C 3
#define T 4

using namespace std;
namespace po = boost::program_options;

__constant__ float constQuery[1024];

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
__global__ void f_scoreSequence(float* subject, float* scoringMatrix, float* maxScoreList, int width, int height, int numSubjects) {
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

                // Just index scoring matrix from shared/constant memory in the future
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
        filestream.ignore(numeric_limits<streamsize>::max(), '\n');

        stringstream fasta_stream;
        fasta_stream << filestream.rdbuf();
        buffer.reserve(10000); // optimization -- reserve some arbitrary length
        string tmp;
        while (fasta_stream) {
            fasta_stream >> tmp;
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

    for (map<int, vector<subject_sequence> >::iterator it = parsedDB.begin(); it != parsedDB.end(); ++it) {
        cout << it->first 
             << ":"
             << it->second.size()
             << endl;
    }

    cout << largestSubjectLength << endl;
    cout << numSubjects << endl;
    cout << subjectLengthSum << endl;

    // alloc memory on GPU
    float* d_input_query = new float[querySequence.length()];
    memset(d_input_query, 0, sizeof(float) * querySequence.length());

    float* d_input_subject;
    cudaMallocManaged((void**) &d_input_subject, (largestSubjectLength * numSubjects) * sizeof(float));

    float* d_output_scoring;
    cudaMallocManaged((void**) &d_output_scoring, ((querySequence.length() + 1) * (largestSubjectLength + 1) * numSubjects) * sizeof(float));
    
    float* d_output_max_score;
    cudaMallocManaged((void**) &d_output_max_score, numSubjects * sizeof(float));

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

    for (int i = 0; i < numSubjects; i++) {
        for (int j = 0; j < largestSubjectLength; j++) { // Will need to pad here
            switch(subjectSequences[i][j])
            {
                case 'A': { d_input_subject[i*largestSubjectLength + j] = A;
                            break;
                        }
                case 'G': { d_input_subject[i*largestSubjectLength + j] = G;
                            break;
                        }
                case 'C': { d_input_subject[i*largestSubjectLength + j] = C;
                            break;
                        }
                case 'T': { d_input_subject[i*largestSubjectLength + j] = T;
                            break;
                        }
            }
        }
    }

    cudaMemcpyToSymbol(constQuery, d_input_query, sizeof(float)*querySequence.length());

    int grid_y_dim = ceil(numSubjects / 32.0);
    
    // Call GPU
    dim3 block(1, 32);
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
    cout << "Sum of DB length: " << subjectLengthSum << " chars." << endl;
    cout << "Time elapsed: " << seconds_elapsed << " seconds." << endl;
    cout << "Performance: " << 1E-9 * (querySequence.length() * subjectLengthSum)
            / seconds_elapsed << " GCUPS." << endl;

    // Free device memory
    cudaFree(d_input_query);
    cudaFree(d_input_subject);
    cudaFree(d_output_scoring);
    cudaFree(d_output_max_score);
    cudaDeviceReset();
}
