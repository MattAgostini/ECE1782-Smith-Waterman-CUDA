#ifndef SWSOLVER_H
#define SWSOLVER_H

#include <vector>
#include "FASTAParsers.h"

typedef std::pair<int, int> seqid_score;

void smith_waterman_cuda(FASTAQuery &query, FASTADatabase &db, std::vector<seqid_score> &result);

#endif /* SWSOLVER_H */
