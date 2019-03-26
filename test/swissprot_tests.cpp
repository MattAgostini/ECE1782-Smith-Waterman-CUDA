#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
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
    //BOOST_CHECK( add( 2,2 ) == 4 );        // #1 continues on error
    FASTAQuery query("data/queries/P02232.fasta", true);
    FASTADatabase db("data/dbs/uniprot_sprot.fasta");
    
    std::vector<seqid_score> result = smith_waterman_cuda(query, db);
    std::map<int, int> reference_results = parse_golden_results("test/reference/P02232.txt");
    
    //std::cout << result[0].first << ":" << result[0].second << std::endl;
    //query.print_buffer();
    
    for (vector<seqid_score>::iterator it = result.begin(); it != result.end(); ++it) {
        //cout << (*it).first << ":" << (*it).second << endl;
        BOOST_CHECK( (*it).second == reference_results[(*it).first] );
    }
    
    
    
    //std::cout << reference_results[result[0].first] << std::endl;
}