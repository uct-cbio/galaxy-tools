/*
 * revcomp - A tool to calculate the reverse compliment of a nucleotide sequence
 * version 1.1
 * originally written by N. Immelman - immelman@science.uct.ac.za
 * modified by G.R. Botha - gerrit@cbio.uct.ac.za
 *
 * Usage:
 * revcomp [-v] [-w] [-i] input_dir/* [-o] output_dir
 *
 * -v : 	turns on verbose output. revcomp will print out the original file names,
 * 	their headers and both the original and reverse complimented sequence.
 *	Default behavior to be silent, except for errors
 * -w :	turns on overwriting. revcomp will overwrite the original sequences with
 *	the reverse complimented sequence. Default behaviour to put the new
 *	sequences in a directory called output_dir, in the current working
 *	directory.
 * -i input_dir/* :	a list of files containing the sequences to be processed. These
 *		must be in FASTA format, and contain only one sequence per file.
 * -o output_dir : directory where reverse complimented sequence will be stored.
 *
 */

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <error.h>
#include <errno.h>
using namespace std;

/*
 * struct fasfile
 * a structure containing all the file information per sequence.
 */
struct fasfile {
	string filename; //The name of the file the original sequence came from
	string header; //The header (sequence name, id, etc)
	string seq; //The original sequence
	string rcseq; //The new reverse, complimented sequence
	fasfile(string fn) //a struct constructor
	:
		filename(fn) {
	}
};
/*
 * usage() - Prints out the usage and then exits
 * Usually when the user has entered an invalid parameter.
 */
void usage() {
	cerr
			<< "revcomp - Calculate the reverse compliment of a nucleotide sequence.\nUsage:\n\trevcomp [-v] [-w] [-i] input_dir [-o] output_dir \nWhere:\n\t-v\t\tverbose output\n\t-w\t\toverwrite the input files with the new sequence (if specified -o flag is ignored) \n\tinput_dir\tcontains a nucleotide sequences in FASTA format (can be specified by wildcards e.g. input_dir/*_R_*.fsa)\n\toutput_dir\tcontains the reverse complemented sequence files\n\t"
			<< endl;
	exit(1);
}

int main(int argc, char** argv) {
	if (argc == 1) //Check if there are no arguments
		usage();
	bool verbose = false, overwrite = false;//some flags set by the parameters
	vector<fasfile *> files; //a vector which stores all the sequences in the struct fasfile
	string input_dir;
	string output_dir;
	bool iflag = false, oflag = false;

	int x = 1;
	while (x < argc) {
		if (argv[x][0] == '-') { //if the argument starts with -, it is not a file name
			switch (argv[x][1]) //check which parameter passed in
			{
				case 'v': //setting verbose output
				{
					verbose = true;
					break;
				}
				case 'w': //setting overwriting
				{
					overwrite = true;
					break;
				}
				case 'i': /* input directory */
				{
					if (!iflag) {
						iflag = true;
						oflag = false;
						break;
					} else {
						usage();
						break;
					}
				}

				case 'o': /* output directory */
				{
					if (!oflag) {
						oflag = true;
						iflag = false;
						break;
					} else {
						usage();
						break;
					}
				}

				default: //invalid option
				{
					usage();
					break;
				}
			}
		} else { //if it doesnt start with a "-", then it must be the name of a file
			if (iflag) {
				input_dir = string(argv[x]);
				files.push_back(new fasfile(input_dir)); //so push it into the vector
			}else if (oflag) {
				output_dir = string(argv[x]);
				oflag = false;
			}else
				usage();
		}
 		x++;
	}

	if (files.size() == 0) //we havent found any files listed in the command line
		usage();
	char * tmp = new char[82]; //a temp variable to handle file reading
	int count; //a variable to tell us how many characters read in at any one time from the file
	for (unsigned int x = 0; x < files.size(); x++) //set up a loop to read in all the files
	{
		ifstream in(files.at(x)->filename.c_str()); //open the stream
		if (!in) //check if the file stream is open, but dont stress if we cant open it
		{
			cerr << "ERROR: Cannot open file " << files.at(x)->filename << endl;
		}
		getline(in, files.at(x)->header); //get the first line, the header
		if (files.at(x)->header.at(0) != '>') //check if it starts with a ">" - which is a fasta file
		{
			cerr << "ERROR: file " << files.at(x)->filename
					<< " does not appear to be a FASTA file" << endl;
			exit(1);
		}

		while (!in.eof() && in.good()) //do this until we reach the eof
		{
			if (in.peek() == '>')//if we find a ">", then there are more than 1 sequences in the file
			{
				cerr
						<< "Error: Cannot handle more than 1 sequences in a file, atm"
						<< endl;
				exit(1);
			}
			if (in.peek() == '\n')
				in.get();
			in.getline(tmp, 81); //read a line. each lien should not be longer than 80 chars (FASTA)

			if (in.bad()) //try find a read error
			{
				cerr << "Read Error" << endl;
				exit(1);
			}
			count = in.gcount(); //find out how much data was read
			files.at(x)->seq.append(tmp, count); //and append it the the previously read sequence
		}
		in.close(); //close the stream
		//in the sequence there should be just characters. i have found some newlines and nulls sneaking in
		unsigned int pos; //used in the next 4 lines

		while ((files.at(x)->seq.find(char(0), 0)) != string::npos) //try find a null character
		{
			pos = files.at(x)->seq.find(char(0), 0);
			files.at(x)->seq.erase(pos, 1); //and erase it from the sequence
		}

		while ((files.at(x)->seq.find(char(11), 0)) != string::npos) //try find a newline character
		{
			pos = files.at(x)->seq.find(char(11), 0);
			files.at(x)->seq.erase(pos, 1);//and erase it from the sequence
		}

	}
	//read files in, now lets revcomp
	//Setup a map, to do a sneaky compliment and to uppercase
	map<char, char> comp;
	comp['A'] = 'T';
	comp['C'] = 'G';
	comp['G'] = 'C';
	comp['T'] = 'A';
	comp['U'] = 'A';
	comp['R'] = 'Y';
	comp['Y'] = 'R';
	comp['K'] = 'M';
	comp['M'] = 'K';
	comp['S'] = 'S';
	comp['W'] = 'W';
	comp['B'] = 'V';
	comp['D'] = 'H';
	comp['H'] = 'D';
	comp['V'] = 'B';
	comp['N'] = 'N';
	comp['-'] = '-';
	comp['a'] = 'T';
	comp['c'] = 'G';
	comp['g'] = 'C';
	comp['t'] = 'A';
	comp['u'] = 'A';
	comp['r'] = 'Y';
	comp['y'] = 'R';
	comp['k'] = 'M';
	comp['m'] = 'K';
	comp['s'] = 'S';
	comp['w'] = 'W';
	comp['b'] = 'V';
	comp['d'] = 'H';
	comp['h'] = 'D';
	comp['v'] = 'B';
	comp['n'] = 'N';
	comp['F'] = 'N'; /* do not know why but 'F' and 'f' is found in some sequences files, for know replace with 'N' */
	comp['f'] = 'N';

	for (unsigned int y = 0; y < files.size(); y++) //a loop to the do revcomp
	{
		for (unsigned int x = 0; x < files.at(y)->seq.length(); x++)
			files.at(y)->rcseq.insert(x,1,comp[files.at(y)->seq.at(files.at(y)->seq.length() - x - 1)]);

		if (verbose){ /* output the new sequences, if verbose */
			cout << "_______________________________" << endl;
			cout <<  files.at(y)->filename << "\t" << files.at(y)->header << endl;
			cout << "Original Sequence:" << endl << files.at(y)->seq << endl;
			cout << "RC sequence:" << endl << files.at(y)->rcseq << endl;
			cout << "_______________________________" << endl;
		}
	}

	if (!overwrite) //if we arent overwriting, then create a directory for output
	{
		/* need to implement remove directory without calling a system command */
		string rm_dir_command = "rm -rf " + output_dir;
		system(rm_dir_command.c_str());

		if (mkdir(output_dir.c_str(), 0777) == -1 && errno != EEXIST) //this creates a directory, and check for errors, apart from if the directory exists
		{
			cerr << "Error: Cannot create directory for output" << endl;
			exit(1); //die
		}
		if (chdir(output_dir.c_str())) //change directory to the newly created (or existing one)
		{
			cerr << "Error: Cannot change to directory for output" << endl;
			exit(1);
		}
	}

	//write out the results
	for (unsigned int x = 0; x < files.size(); x++) {
		ofstream out;

		if((files.at(x)->filename.find_last_of('/',files.at(x)->filename.size()-1)) != string::npos){
			int filename_start_pos = files.at(x)->filename.find_last_of('/',files.at(x)->filename.size()-1) + 1;
			string output_filename = files.at(x)->filename.substr(filename_start_pos, files.at(x)->filename.size()-1);
			out.open((output_filename).c_str());
		}else{
			cerr
					<< "Need to FIX: specify input_dir as absolute path e.g. ./Xh_* not Xh_*" << endl;
			exit(1);
		}

		 //Open the outfile stream
		out << files.at(x)->header; //write the header
		for (unsigned int y = 0; y < files.at(x)->rcseq.length(); y++) //loop to write the sequence out
		{
			if (y % 80 == 0) //cant have lines more than 80 chars
				out.put('\n');
			out.put(files.at(x)->rcseq.at(y));
		}
		out.put('\n');
		out.close();
	}

	return 0; //success
}
