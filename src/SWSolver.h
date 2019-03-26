#ifndef SWSOLVER_H
#define SWSOLVER_H

#include <vector>
#include "FASTAParsers.h"

typedef std::pair<int, int> seqid_score;

std::vector<seqid_score> smith_waterman_cuda(FASTAQuery &query, FASTADatabase &db);


#endif /* SWSOLVER_H */