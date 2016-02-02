#include "mpi.h"
#include <unistd.h>
#include <string>

#include "IOHandler.h"
#include "MCC.h"
#include "Fingerprint.h"
#include "Functions.h"

using namespace std;

IOHandler::IOHandler(int argc, char *argv[])
{
	int c;
	vector<string> fp_types;

	string valid_args = string("-a:b:c:d:ef:hi:k:m:o:q:r:s:t:u:v:w");

	// Standard output initialization, to (somewhat) keep the correct writing order of the threads
	if(setvbuf(stdout,NULL,_IONBF,0))
		ParallelHandler::exit("Error in setvbuf");

	// Set default values
	threads = omp_get_num_procs();
	stop = ParallelHandler::STOP_MAX;
	threshold = vector<double>(1, 0);
	ranking = 10;
	quality = 0;
	templatefiles.clear();
	classfiles.clear();
	inputfiles.clear();
	featurefiles.clear();
	output = OUT_SHORT;
	fusion = ParallelHandler::FUSION_SUM;
	outpath = "output";
	inputclasses = false;

	// getopt variables initialization
	optind = 1;
	opterr = 0;

	while((c = getopt (argc, argv, valid_args.c_str())) != -1)
	{
		switch (c)
		{
			case 'a':
				fp_types = stringSplit(optarg, ",");
				fp_type.resize(fp_types.size());

				for(unsigned int i = 0; i < fp_types.size(); ++i)
					fp_type[i] = getFpType(fp_types[i]);
				break;
			case 'h':
				printSyntax();
				ParallelHandler::exit("", 0);
				break;
			case 'k':
				threads = atoi(optarg);
				break;
			case 'u':
				threshold = stringSplitDouble(optarg, ",");
				break;
			case 'r':
				ranking = atoi(optarg);
				break;
			case 's':
				stop = getStopType(optarg);
				break;
			case 'q':
				quality = atoi(optarg);
				break;
			case 't':
				templatefiles = stringSplit(optarg, ",");
				break;
			case 'c':
				classfiles = stringSplit(optarg, ",");
				break;
			case 'b':
				featurefiles = stringSplit(optarg, ",");
				break;
			case 'i':
				inputfiles = stringSplit(optarg, ",");
				break;
			case 'o':
				outpath = optarg;
				break;
			case 'd':
				classifierpaths = stringSplit(optarg, ",");
				break;
			case 'w':
				inputclasses = true;
				break;
			case 'm':
				output = getOutputType(optarg);
				break;
			case 'f':
				fusion = getFusionType(optarg);
				break;
			case 'e':
				clearCache();
				break;
			case 'v':
				if(strstr(argv[0], "DPDDFF") == 0)
					cerr << "Syntax warning: option -v (DPDDFF variant) does not make sense with an executable other than DPDDFF" << endl;

				variant = getVariant(optarg);
				break;
			case '?':
				if (valid_args.find(optopt) != string::npos)
					cerr << "Syntax error: option -" << (char)optopt << " requires an argument. Using default." << endl;
		}
	}

	// Consider non-valid combinations
	if(templatefiles.empty())
		ParallelHandler::exit("Parameter -t <database_filenames> is mandatory", -1);
	if(inputfiles.empty())
		ParallelHandler::exit("Parameter -i <input_filenames> is mandatory", -1);
	if(templatefiles.size() != inputfiles.size())
		ParallelHandler::exit("Parameters -t and -i must have the same number of files", -1);
	if(classfiles.size() > 1)
		ParallelHandler::exit("Parameter -c must have a single file. Classification for multiple fingers is not yet supported.", -1);
	if(featurefiles.size() != classifierpaths.size())
		ParallelHandler::exit("Parameter -b (features of the input fingerprints) and -d (trained classifiers) must be provided together.", -1);
	if(!classfiles.empty() && classfiles.size() != templatefiles.size())
		ParallelHandler::exit("Parameters -t and -c must have the same number of files", -1);
	if(!classfiles.empty() && strstr(argv[0], "DPDDFF") > 0)
		ParallelHandler::exit("Parameter -c is ignored when using double phase", -1);

	// More default values
	if(fp_type.empty())
		fp_type.resize(inputfiles.size(), ParallelHandler::FP_JIANG);

	inputfd.resize(inputfiles.size());

	for(unsigned int i = 0; i < inputfiles.size(); ++i)
	{
		inputfd[i] = new ifstream(inputfiles[i].c_str());

		if(!inputfd[i]->good())
			ParallelHandler::exit(string("Could not open input file ") + inputfiles[i], -1);
	}
}

IOHandler::~IOHandler()
{
	for(unsigned int i = 0; i < inputfd.size(); ++i)
	{
		inputfd[i]->close();
		delete inputfd[i];
	}

	if(scorefile)
		scorefile.close();

	if(timefile)
		timefile.close();
}

void IOHandler::printSyntax()
{
	cerr << endl;
	cerr << "PARALLEL FINGERPRINT RECOGNITION FRAMEWORK" << endl;
	cerr << "Usage: mpirun -np <number_of_processes> ... <executable> [parameters]" << endl;
	cerr << endl;
	cerr << "number_of_processes should be greater or equal to 2" << endl;
	cerr << endl;
	cerr << "List of framework parameters:" << endl;
	cerr << "\t-h : shows this help message" << endl;
	cerr << "\t-a <algorithm> : matching algorithm to be used. Defaults to Jiang's algorithm." << endl;
	cerr << "\t-k <threads> : number of threads for OpenMP. Defaults to the number of detected processors (which usually is twice the number of physical cores)." << endl;
	cerr << "\t-s <stop> : stop criterion for the search." << endl;
	cerr << "\t-r <ranking> : ranking for the search (only applicable to \"r\" stop criterion). Defaults to 10." << endl;
	cerr << "\t-u <threshold> : threshold for the search. Defaults to 0." << endl;
	cerr << "\t-q <quality> : minimum minutiae quality. Defaults to 0." << endl;
	cerr << "\t-t <template_files> : comma-separated list of files that in turn contain lists of the template fingerprints. Each list should correspond to a finger. No default." << endl;
	cerr << "\t-c <class_files> : comma-separated list of files that contain the classes of the template fingerprints. Each list should correspond to a finger (although the multi-finger case is not implemented yet). If not provided, no classification is performed." << endl;
	cerr << "\t-i <input_files> : comma-separated list of files that in turn contain lists of the input fingerprints. Each list should correspond to a finger. No default." << endl;
	cerr << "\t-b <feature_matrix> : comma-separated list of files. Each of these files must contain the feature matrix of the input fingerprints, in CSV format, where each row corresponds to a fingerprint. No default." << endl;
	cerr << "\t-d <classifier_model> : file (in RData format) in which the trained classifier model is stored. No default." << endl;
	cerr << "\t-o <ouptut_file> : root name for the output files. Defaults to \"output\"." << endl;
	cerr << "\t-m <output_mode> : output mode. Defaults to \"short\"" << endl;
	cerr << endl;
	cerr << "List of parameters for the matching algorithms:" << endl;
	cerr << "\t-F <fusion> : fusion aggregation function. Defaults to using only minutiae, without fusion." << endl;
	cerr << endl;
	cerr << "List of parameters for the MCC matching algorithm:" << endl;
	cerr << "\t-N <ns> : number of cells on each side of the cylinder. Defaults to 8." << endl;
	cerr << "\t-C <consolidation> : consolidation algorithm. Defaults to LSSR." << endl;
	cerr << "\t-B : activates the bit-encoding. Disabled by default." << endl;
	cerr << "\t-H : deactivates the convex-hull filter. Disabled by default." << endl;

}

void IOHandler::readFileName(int finger, string &cinput, const ParallelHandler *parallel_handler)
{
	int string_size;

	if(parallel_handler->isMaster())
	{
		readFileNameLocal(finger, cinput);
		string_size = cinput.size();
	}

	MPI::COMM_WORLD.Bcast(&string_size, 1, MPI::INT, ParallelHandler::MASTER_ID);
	cinput.resize(string_size);
	MPI::COMM_WORLD.Bcast(&cinput[0], string_size, MPI::CHAR, ParallelHandler::MASTER_ID);
	cinput[string_size] = '\0';
}

vector< Matrix<double> > IOHandler::readFeatureMatrices()
{
	vector< Matrix<double> > res(featurefiles.size());

	for(unsigned int i = 0; i < featurefiles.size(); ++i)
		res[i].readFromCSV(featurefiles[i], inputclasses);

	return res;
}

void IOHandler::readFileName(vector<string> &cinput, const ParallelHandler *parallel_handler)
{
	cinput.resize(inputfd.size());
	for(unsigned int i = 0; i < inputfd.size(); ++i)
		readFileName(i, cinput[i], parallel_handler);
}


void IOHandler::readFileNameLocal(int finger, string &cinput)
{
	do {
		getline(*inputfd[finger], cinput);
	} while(*inputfd[finger] && cinput.empty());

	if(!*inputfd[finger])
		cinput = "exit";
}


ParallelHandler::fingerprint_t IOHandler::getFpType(const string &cinput)
{
	if(cinput == "jiang")
		return ParallelHandler::FP_JIANG;
	else if(cinput == "mcc")
		return ParallelHandler::FP_MCC;
	else
	{
		ParallelHandler::exit(string("Syntax error: the <algorithm> argument should be {jiang|mcc} and found ") + string(cinput) + " instead");
		return ParallelHandler::FP_UNKNOWN;
	}
}


ParallelHandler::fusion_t IOHandler::getFusionType(const string &cinput)
{
	if(cinput == "sum")
		return ParallelHandler::FUSION_SUM;
	else if(cinput == "max")
		return ParallelHandler::FUSION_MAX;
	else if(cinput == "min")
		return ParallelHandler::FUSION_MIN;
	else if(cinput == "product")
		return ParallelHandler::FUSION_PROD;
	else
	{
		ParallelHandler::exit(string("Syntax error: the <fusion> argument should be {sum|hybrid|max|min|product} and found ") + string(cinput) + " instead");
		return ParallelHandler::FUSION_SUM;
	}
}


ParallelHandler::stop_t IOHandler::getStopType(const string &cinput)
{
	if(cinput == "y")
		return ParallelHandler::STOP_YES;
	else if(cinput[0] == 'n')
		return ParallelHandler::STOP_NO;
	else if(cinput == "m")
		return ParallelHandler::STOP_MAX;
	else if(cinput == "a")
		return ParallelHandler::STOP_ALL;
	else if(cinput[0] == 'r')
		return ParallelHandler::STOP_RANKING;
	else
	{
		ParallelHandler::exit("Syntax error: the <stop> argument should be {y|n|m|a|r}");
		return ParallelHandler::STOP_MAX;
	}
}


IOHandler::output_t IOHandler::getOutputType(const string &cinput)
{
	if(cinput == "short")
		return OUT_SHORT;
	else if(cinput == "long")
		return OUT_LONG;
	else if(cinput == "files")
		return OUT_BYFILE;
	else
	{
		ParallelHandler::exit("Syntax error: the <output> argument should be {short|long|files}");
		return OUT_SHORT;
	}
}


ParallelHandler::variant_t IOHandler::getVariant(const string &cinput)
{
	if(cinput == "SS")
		return ParallelHandler::DPDDFF_SS;
	else if(cinput == "SD")
		return ParallelHandler::DPDDFF_SD;
	else if(cinput == "DS")
		return ParallelHandler::DPDDFF_DS;
	else if(cinput == "DD")
		return ParallelHandler::DPDDFF_DD;
	else
	{
		ParallelHandler::exit("Syntax error: the <variant> argument should be {SS|SD|DS|DD}");
		return ParallelHandler::DPDDFF_SS;
	}
}


void IOHandler::printInitialOutput()
{
	if(output == IOHandler::OUT_BYFILE)
	{
		scorefile.open("scores.dat");
		timefile.open("times.dat");
		meantime = 0.0;

		if(!scorefile || ! timefile)
			ParallelHandler::exit(string("ERROR: Could not open files for output"), -1);
	}
}


void IOHandler::printIdentificationOutput(const string &cinput, const vector<Score> &matches, double avgmatchingtime, const vector<double> &times, int iteration, float num_candidates, float penetrationrate)
{
	if(output == IOHandler::OUT_SHORT)
	{
		for(vector< Score >::const_iterator i = matches.begin(); i != matches.end(); ++i)
			cout << cinput << "\t" << i->getScore() << "\t" << i->getId() << "\t" << printTabSeparated(times) << "\t" << avgmatchingtime << "\t" << num_candidates << "\t" << penetrationrate << endl;

		if(matches.empty())
			cout << cinput << "\t" <<    -1    << "\t" <<    -1     << "\t" << printTabSeparated(times) << "\t" << avgmatchingtime << "\t" << num_candidates << "\t" << penetrationrate << endl;

	}
	else if(output == IOHandler::OUT_LONG)
	{
		cout << "INICIO" << endl;
		cout << "Huella: " << cinput << endl;

		cout << "Average matching time: " << avgmatchingtime << endl;
		cout << "Penetration rate: " << penetrationrate << endl;

		for(vector< Score >::const_iterator i = matches.begin(); i != matches.end(); ++i)
			if(i->getScore() != -1)
				cout << i->getScore() << " " << i->getId() << endl;

		cout << "Tiempo medio " << getFpType() << " en la huella " << iteration << " : " << times[1] << " segundos" << endl;
		cout << "FIN" << endl;
	}
	else if(output == IOHandler::OUT_BYFILE)
	{
		int bestindex = -1;

		cout << "INICIO" << endl;
		cout << "Huella: " << cinput << endl;

		cout << "Average matching time: " << avgmatchingtime << endl;
		cout << "Penetration rate: " << penetrationrate << endl;

		for(unsigned int i = 0; i < matches.size(); ++i)
		{
			scorefile << matches[i].getScore() << " ";
			if(Fingerprint::better(matches[i].getScore(), matches[bestindex].getScore()))
				bestindex = i;
		}

		scorefile << endl;
		timefile << times[1] << endl;
		meantime += times[1];

		if(matches[bestindex].getScore() != -1)
			cout << matches[bestindex].getScore() << " " << matches[bestindex].getId() << endl;

		cout << "Tiempo medio " << getFpType() << " en la huella " << iteration << " : " << times[1] << " segundos" << endl;
		cout << "FIN" << endl;
	}
}

void IOHandler::printFinalOutput(const ParallelHandler *parallel_handler)
{
	if(output == IOHandler::OUT_LONG)
		cout << "FINFICHERO" << endl;
	else if(output == IOHandler::OUT_BYFILE)
	{
		timefile << "Average: " << meantime/(parallel_handler->getIteration()) << endl;
		scorefile.close();
		timefile.close();
	}
}

int IOHandler::getNumInputFingerprints(int finger) const
{
	int tmp = std::count(istreambuf_iterator<char>(*inputfd[finger]), istreambuf_iterator<char>(), '\n');

	inputfd[finger]->seekg(0);

	return tmp;
}



void IOHandler::printScores(const vector<Score> &matches) const
{
	for(vector<Score>::const_iterator i = matches.begin(); i != matches.end(); ++i)
		cout << i->getScore() << "\t" << i->getId() << endl;
}


void IOHandler::printParameters() const
{
// 	cerr <<
//
// 		vector<ParallelHandler::fingerprint_t> fp_type;
// 	ParallelHandler::fusion_t fusion;
// 	unsigned int threads;
// 	ParallelHandler::stop_t stop;
// 	unsigned int ranking;
// 	vector<double> threshold;
// 	unsigned int quality;
// 	output_t output;
// 	ParallelHandler::variant_t variant;
//
// 	vector<string> templatefiles;
// 	vector<string> inputfiles;
//
// 	vector<ifstream *> inputfd;
//
// 	ofstream scorefile;
// 	ofstream timefile;
// 	string outpath;
// 	double meantime;
}

void IOHandler::clearCache() const
{
	sync();

	std::ofstream ofs("/proc/sys/vm/drop_caches");
	ofs << "3" << std::endl;
}
