#ifndef SWSOLVERCHAR_H
#define SWSOLVERCHAR_H

#include <vector>
#include "FASTAParsers.h"

typedef std::pair<int, int> seqid_score;

std::vector<seqid_score> smith_waterman_cuda_char(FASTAQuery &query, FASTADatabase &db);

#endif /* SWSOLVERCHAR_H */
