#!/home/mkasa/lcl/bin/cpas

#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

static const size_t BUFFER_SIZE = 16 * 1024 * 1024u;

static void output_read_name(const char* header)
{
	if(*header++ == '\0') return;
	while(*header != '\0' && *header != ' ') cout << *header++;
	cout << '\n';
}

void show_read_names_in_file(const char* fname)
{
	ifstream ist(fname);
	if(!ist) {
		cerr << "Cannot open '" << fname << "'" << endl;
		return;
	}
    char* b = new char[BUFFER_SIZE];
	size_t line_count = 0;
    if(ist.getline(b, BUFFER_SIZE)) {
		++line_count;
		output_read_name(b);
        if(b[0] != '@') { 
			// This should be FASTA
            while(ist.getline(b, BUFFER_SIZE)) {
				++line_count;
                if(b[0] == '>') output_read_name(b);
			}
		} else {
			// This should be FASTQ
        	size_t number_of_nucleotides_in_read = 0;
            while(ist.getline(b, BUFFER_SIZE)) {
				++line_count;
                if(b[0] == '+' && b[1] == '\0') { // EOS
                    long long n = number_of_nucleotides_in_read;
                    while(ist.getline(b, BUFFER_SIZE)) {
                        ++line_count;
                        const size_t number_of_qvchars_in_line = strlen(b);
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    if(ist.peek() != '@' && !ist.eof()) {
                        cerr << "WARNING: bad file format? at line " << line_count << endl;
                    }
                    number_of_nucleotides_in_read = 0;
                    if(!ist.getline(b, BUFFER_SIZE)) break;
                    ++line_count;
					output_read_name(b);
                } else {
                    const size_t number_of_nucleotides_in_line = strlen(b);
                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                }
            }
		}
	}
    delete b;
}

void count_number_of_reads_in_file(const char* fname)
{
    ifstream ist(fname);
    if(!ist) {
        cerr << "Cannot open '" << fname << "'" << endl;
        return;
    }
	cout << fname << flush;
    size_t number_of_sequences = 0;
    size_t number_of_nucleotides = 0;
	size_t line_count = 0;
	size_t min_read_len = -1;
	size_t max_read_len = 0;
	#define UPDATE_MIN_AND_MAX(len) min_read_len = std::min<size_t>(min_read_len, len); max_read_len = std::max<size_t>(max_read_len, len);
    char* b = new char[BUFFER_SIZE];
    if(ist.getline(b, BUFFER_SIZE)) {
		++line_count;
        number_of_sequences++;
        size_t number_of_nucleotides_in_read = 0;
        if(b[0] != '@') { 
            // This should be FASTA
            while(ist.getline(b, BUFFER_SIZE)) {
				++line_count;
                if(b[0] == '>') {
                    number_of_sequences++;
					UPDATE_MIN_AND_MAX(number_of_nucleotides_in_read);
					number_of_nucleotides_in_read = 0;
                } else {
                    number_of_nucleotides += strlen(b);
                }
            }
        } else {
            // This is FASTQ
            while(ist.getline(b, BUFFER_SIZE)) {
				++line_count;
                if(b[0] == '+' && b[1] == '\0') { // EOS
                    long long n = number_of_nucleotides_in_read;
                    while(ist.getline(b, BUFFER_SIZE)) {
                        ++line_count;
                        const size_t number_of_qvchars_in_line = strlen(b);
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    if(ist.peek() != '@' && !ist.eof()) {
                        cerr << "WARNING: bad file format? at line " << line_count << endl;
                    }
					UPDATE_MIN_AND_MAX(number_of_nucleotides_in_read);
                    number_of_nucleotides_in_read = 0;
                    if(!ist.getline(b, BUFFER_SIZE))
                        break;
                    ++line_count;
					++number_of_sequences;
                } else {
                    const size_t number_of_nucleotides_in_line = strlen(b);
                    number_of_nucleotides += number_of_nucleotides_in_line;
                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                }
            }
        }
    }
	cout << '\t' << number_of_sequences << '\t' << number_of_nucleotides << '\t' << (double(number_of_nucleotides) / number_of_sequences);
	cout << '\t' << min_read_len << '\t' << max_read_len << '\n';
    delete b;
	#undef UPDATE_MIN_AND_MAX
}

void do_count(int argc, char** argv)
{
	cout << "FILE\tNUM_READS\tNUM_NUCLS\tAVG_READ_LEN\tMIN_READ_LEN\tMAX_READ_LEN\n";
	for(int i = 2; i < argc; ++i) {
		count_number_of_reads_in_file(argv[i]);
	}
}

void do_name(int argc, char** argv)
{
	for(int i = 2; i < argc; ++i) {
		show_read_names_in_file(argv[i]);
	}
}

void show_usage()
{
    cerr << "Usage: fatt <command> [options...]" << endl;
}

void show_help()
{
    show_usage();
    cerr << "\nCommands:\n";
	cerr << "\tcount\tcount the number of reads/nucleotides\n";
	cerr << "\tname\toutput the names of reads\n";
    cerr << "\thelp\tshow help message\n";
    cerr << "\nType 'fatt help <command>' to show the detail of the command.\n";
}

void dispatchByCommand(const string& commandString, int argc, char** argv)
{
	if(commandString == "count") {
		do_count(argc, argv);
		return;
	}
	if(commandString == "name") {
		do_name(argc, argv);
		return;
	}
    // Help or error.
    if(commandString != "help") {
        cerr << "ERROR: Unknown command '" << commandString << "'" << endl;
    }
    show_help();
}

int main(int argc, char** argv)
{
    if(argc < 2) {
        show_usage();
        return 1;
    }
    const char* commandStr = argv[1];
    dispatchByCommand(commandStr, argc, argv);
    return 0;
}
