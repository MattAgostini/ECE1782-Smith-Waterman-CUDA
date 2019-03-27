#include <sys/time.h>
#include <vector>

#include "boost/program_options.hpp"

#include "FASTAParsers.h"
#include "SWSolver.h"

using namespace std;
namespace po = boost::program_options;

// Time stamp function
double getTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_usec/1000000 + tv.tv_sec;
}

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
    FASTAQuery query(querypath, true);
    cout << "Input buffer:";
    query.print_buffer();
    cout << endl;
    string querySequence = query.get_buffer();

    // Parse database file
    std::string datapath = vm["db"].as<std::string>();
    FASTADatabase db(datapath);
    
    vector<seqid_score> result = smith_waterman_cuda(query, db);
	
    for (vector<seqid_score>::iterator it = result.begin(); it != result.end(); ++it) {
        cout << (*it).first << ":" << (*it).second << endl;
    }

    double time_end = getTimeStamp();
    double seconds_elapsed = time_end - time_start;

    cout << std::string(80, '=') << endl;
    cout << "METRICS:" << endl;
    cout << "Query length: " << querySequence.length() << " chars." << endl;
    cout << "Num subjects: " << db.numSubjects << endl;
    cout << "Sum of DB length: " << db.subjectLengthSum << " chars." << endl;
    cout << "Time elapsed: " << seconds_elapsed << " seconds." << endl;
    cout << "Performance: " << 1E-9 * (querySequence.length() * db.subjectLengthSum)
            / seconds_elapsed << " GCUPS." << endl;

}

