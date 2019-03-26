#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
#include <boost/test/unit_test.hpp>
#include <string>

#include "../src/FASTAParsers.h"
#include "../src/SWSolver.h"

int add( int i, int j ) { return i+j; }

vector<seqid_score> parse_golden_results(std::string filepath) {
	
}

BOOST_AUTO_TEST_CASE( my_test )
{
    // seven ways to detect and report the same error:
    BOOST_CHECK( add( 2,2 ) == 4 );        // #1 continues on error

    BOOST_REQUIRE( add( 2,2 ) == 4 );      // #2 throws on error

    if( add( 2,2 ) != 4 )
      BOOST_ERROR( "Ouch..." );            // #3 continues on error

    if( add( 2,2 ) != 4 )
      BOOST_FAIL( "Ouch..." );             // #4 throws on error

    if( add( 2,2 ) != 4 ) throw "Ouch..."; // #5 throws on error

    BOOST_CHECK_MESSAGE( add( 2,2 ) == 4,  // #6 continues on error
                         "add(..) result: " << add( 2,2 ) );

    
	
	FASTAQuery query("../data/queries/P02232.fasta", true);
	FASTADatabase db("data/dbs/uniprot_sprot.fasta");
	
	std::vector<seqid_score> result = smith_waterman_cuda(query, db);
	BOOST_CHECK(!result.empty());
	
	
}