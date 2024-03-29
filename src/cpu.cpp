#include <sys/time.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>

#define SEQ_EQUAL 3
#define SEQ_DIFF -3
#define GAP_PENALTY 2

#define FROM_LEFT 1
#define FROM_TOP 2
#define FROM_TOP_LEFT 3

using namespace std;

int main( int argc, char *argv[] ) {
    // get program arguments
    if (argc != 3) {
        printf("Error: wrong number of args\n");
        exit(1);
    }
    
    int subMatrix[2] = {SEQ_EQUAL, SEQ_DIFF};

    string seqA = argv[1];
    string seqB = argv[2];
    
    int scoringMatrix[(seqA.length() + 1)][(seqB.length() + 1)];
    int tracebackMatrix[(seqA.length() + 1)][(seqB.length() + 1)];
    
    // Initialize Matrix
    for (int i = 0; i < (seqA.length() + 1); i++) {
        for (int j = 0; j < (seqB.length() + 1); j++) {
            scoringMatrix[i][j] = 0;
            tracebackMatrix[i][j] = 0;
        }
    }
    
    // Do the scoring
    int maxScore = 0;
    int maxIndexI = 0;
    int maxIndexJ = 0;
    for (int i = 1; i < (seqA.length() + 1); i++) {
        for (int j = 1; j < (seqB.length() + 1); j++) {
            int score = 0;
            
            if (scoringMatrix[i][j - 1] - GAP_PENALTY > score) {
                score = scoringMatrix[i][j - 1] - GAP_PENALTY;
                tracebackMatrix[i][j] = FROM_LEFT;
            }
            
            if (scoringMatrix[i - 1][j] - GAP_PENALTY > score) {
                score = scoringMatrix[i - 1][j] - GAP_PENALTY;
                tracebackMatrix[i][j] = FROM_TOP;
            }
            
            int similarityScore = 0;
            if (seqA[i - 1] == seqB[j - 1]) similarityScore = subMatrix[0];
            else similarityScore = subMatrix[1];
            
            if (scoringMatrix[i - 1][j - 1] + similarityScore > score) {
                score = scoringMatrix[i - 1][j - 1] + similarityScore;
                tracebackMatrix[i][j] = FROM_TOP_LEFT;
            }
            
            if (score > maxScore) {
                maxScore = score;
                maxIndexI = i;
                maxIndexJ = j;
            }
                    
            scoringMatrix[i][j] = score;
        }
    }
    
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
            cout << scoringMatrix[i][j] << " "; 
        }
        cout << endl;
    }
}
