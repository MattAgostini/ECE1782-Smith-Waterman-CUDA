#ifndef FASTAPARSERS_H
#define FASTAPARSERS_H

#include <string>
#include <iostream>
#include <fstream>

#include <map>
#include <vector>

//#define LENGTH_THRESHOLD 3000
#define TILE_SIZE 8

using namespace std;

struct subject_sequence {
    int id;
    string sequence;
};

static inline int roundUp(int numToRound, int multiple)
{
    if (multiple == 0)
        return numToRound;

    int remainder = numToRound % multiple;
    if (remainder == 0)
        return numToRound;

    return numToRound + multiple - remainder;
}

class FASTAQuery {
private:
    bool isQuery;
    string buffer;
public:
    FASTAQuery(std::string filepath, bool _isQuery) {
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

    ~FASTAQuery() {
    };

    void print_buffer() {
        cout << buffer << endl;
    };

    string get_buffer() {
        return buffer;
    };
};

class FASTADatabase {
public:
    // key is sequence length, value is a vector of subject_sequence struct
    map<int, vector<subject_sequence> > parsedDB;
    int largestSubjectLength;
    int numSubjects;
    int subjectLengthSum;

    FASTADatabase(std::string filepath) {
        ifstream databaseFile;
        databaseFile.open(filepath.c_str());
        
        largestSubjectLength = 0;
        numSubjects = 0;
        subjectLengthSum = 0;
        
        string temp;
        int _id = -1;
        bool isFirst = true;
        
        subject_sequence tmp;
        string subjectSequence = "";

        while (getline(databaseFile, temp)) {
            // This line denotes the start of a sequence
            if (temp[0] == '>') {
                if (!isFirst) {
                    //while (subjectSequence.length() % 8 != 0) // pad to nearest 8
                    //    subjectSequence = subjectSequence + "/";
                    size_t length = subjectSequence.length();
                    size_t rounded = roundUp(length, TILE_SIZE);
                    subjectSequence.append(rounded - length, '/');

                    tmp.id = _id;
                    tmp.sequence = subjectSequence;
                    //if (tmp.sequence.length() % 8 != 0) exit(0); sequence length checkers for pad to 8
                    parsedDB[subjectSequence.length()].push_back(tmp);

                    subjectLengthSum += subjectSequence.length();
                    largestSubjectLength = max(largestSubjectLength, (int)subjectSequence.length());
                    
                    numSubjects++;
                }
                isFirst = false;
                
                //cout << subjectSequence << endl;
                subjectSequence = "";
                _id++;
            }
            else {
                subjectSequence.append(temp);
            }
            
        }
        // Adding last sequence 
        size_t length = subjectSequence.length();
        size_t rounded = roundUp(length, TILE_SIZE);
        subjectSequence.append(rounded - length, '/');
        //while (subjectSequence.length() % 8 != 0) // pad to nearest 8
        //    subjectSequence = subjectSequence + "/";

        tmp.id = _id;
        tmp.sequence = subjectSequence;
                    //if (tmp.sequence.length() % 8 != 0) exit(0); sequence length checkers for pad to 8
        parsedDB[subjectSequence.length()].push_back(tmp);
        
        subjectLengthSum += subjectSequence.length();
        largestSubjectLength = max(largestSubjectLength, (int)subjectSequence.length());
        
        numSubjects++;
        databaseFile.close();
    };

};

#endif /* FASTAPARSERS_H */

