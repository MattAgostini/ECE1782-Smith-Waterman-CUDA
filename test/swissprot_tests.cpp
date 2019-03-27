#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE QueryComparison
#include <boost/test/unit_test.hpp>
#include <string>
#include <map>
#include <sstream>
#include <iostream>

#include "../src/FASTAParsers.h"
#include "../src/SWSolver.h"

std::map<int, int> parse_golden_results(std::string filepath) {
    ifstream filestream;
    filestream.open(filepath.c_str());
    
    map<int, int> parsed_results;
    
    string tmp;
    
    int idx = 0;
    int score;
    
    while (getline(filestream, tmp)) {
        std::istringstream(tmp) >> score;
        parsed_results[idx++] = score;
    }
    
    filestream.close();
    return parsed_results;
}

BOOST_AUTO_TEST_CASE( P02232 )
{
    FASTAQuery query("data/queries/P02232.fasta", true);
    FASTADatabase db("data/dbs/uniprot_sprot.fasta");
    
    std::vector<seqid_score> result = smith_waterman_cuda(query, db);
    std::map<int, int> reference_results = parse_golden_results("test/reference/P02232.txt");
    
    for (vector<seqid_score>::iterator it = result.begin(); it != result.end(); ++it) {
        BOOST_CHECK( (*it).second == reference_results[(*it).first] );
    }
    
    std::cout << "Number of queries executed in CUDA: " << result.size() << endl;
    //std::cout << reference_results[result[0].first] << std::endl;
}

/*BOOST_AUTO_TEST_CASE( P01008 )
{
    FASTAQuery query("data/queries/P01008.fasta", true);
    FASTADatabase db("data/dbs/uniprot_sprot.fasta");
    
    std::vector<seqid_score> result = smith_waterman_cuda(query, db);
    std::map<int, int> reference_results = parse_golden_results("test/reference/P01008.txt");
    
    for (vector<seqid_score>::iterator it = result.begin(); it != result.end(); ++it) {
        BOOST_CHECK( (*it).second == reference_results[(*it).first] );
    }
    
    std::cout << "Number of queries executed in CUDA: " << result.size() << endl;
    //std::cout << reference_results[result[0].first] << std::endl;
}*/
