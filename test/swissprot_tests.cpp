#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE QueryComparison
#include <boost/test/unit_test.hpp>
#include <string>
#include <map>
#include <sstream>
#include <iostream>
#include <sys/time.h>

#include "../src/FASTAParsers.h"
#include "../src/SWSolver.h"
#include "../src/SWSolver_char.h"

double getTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_usec/1000000 + tv.tv_sec;
}

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

void run_query_performance(std::string querypath, std::string dbpath) {
    //std::vector<double> times;
    
    double time_start = getTimeStamp();

    FASTAQuery query(querypath, true);
    FASTADatabase db(dbpath);

    
    std::vector<seqid_score> result;
    result.reserve(600000);
    smith_waterman_cuda(query, db, result);

    double time_end = getTimeStamp();
    double seconds_elapsed = time_end - time_start;

    cout << "Query " << querypath << " length " << query.get_buffer().length() << ", Performance: " << 1E-9 * (query.get_buffer().length() * db.subjectLengthSum)
            / seconds_elapsed << " GCUPS, Time: " << seconds_elapsed << endl;
}

void run_query_against_reference(std::string querypath, std::string dbpath, std::string refpath) {
    FASTAQuery query(querypath, true);
    FASTADatabase db(dbpath);
    
    std::vector<seqid_score> result;
    result.reserve(600000);
    smith_waterman_cuda(query, db, result);

    std::map<int, int> reference_results = parse_golden_results(refpath);
    
    for (vector<seqid_score>::iterator it = result.begin(); it != result.end(); ++it) {
        BOOST_CHECK_MESSAGE( (*it).second == reference_results[(*it).first], (*it).first << ": Ours: " << (*it).second << " | Theirs: " << reference_results[(*it).first] );
    }
    
    std::cout << "Number of queries executed in CUDA: " << result.size() << endl;
}

BOOST_AUTO_TEST_SUITE( Comparison );

/*
BOOST_AUTO_TEST_CASE(P02232) { 
    run_query_against_reference(
        "data/queries/P02232.fasta",
        "data/dbs/uniprot_sprot.fasta",
        "test/reference/P02232.txt"
    );
}
*/

BOOST_AUTO_TEST_CASE(P01008) { 
    run_query_against_reference(
        "data/queries/P01008.fasta",
        "data/dbs/uniprot_sprot.fasta",
        "test/reference/P01008.txt"
    );
}

BOOST_AUTO_TEST_SUITE_END();
BOOST_AUTO_TEST_SUITE( Performance );
BOOST_AUTO_TEST_CASE(P02232) {run_query_performance("data/queries/P02232.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P05013) {run_query_performance("data/queries/P05013.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P14942) {run_query_performance("data/queries/P14942.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P07327) {run_query_performance("data/queries/P07327.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P01008) {run_query_performance("data/queries/P01008.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P03435) {run_query_performance("data/queries/P03435.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P42357) {run_query_performance("data/queries/P42357.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P21177) {run_query_performance("data/queries/P21177.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P27895) {run_query_performance("data/queries/P27895.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P07756) {run_query_performance("data/queries/P07756.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P04775) {run_query_performance("data/queries/P04775.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P19096) {run_query_performance("data/queries/P19096.fasta","data/dbs/uniprot_sprot.fasta");}
/*
BOOST_AUTO_TEST_CASE(P28167) {run_query_performance("data/queries/P28167.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P20930) {run_query_performance("data/queries/P20930.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P0C6B8) {run_query_performance("data/queries/P0C6B8.fasta","data/dbs/uniprot_sprot.fasta");}
BOOST_AUTO_TEST_CASE(P08519) {run_query_performance("data/queries/P08519.fasta","data/dbs/uniprot_sprot.fasta");} 
BOOST_AUTO_TEST_CASE(P33450) {run_query_performance("data/queries/P33450.fasta","data/dbs/uniprot_sprot.fasta");}
*/
BOOST_AUTO_TEST_SUITE_END();

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
