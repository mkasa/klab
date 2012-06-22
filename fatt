#!/home/mkasa/lcl/bin/cpas
// -*- mode:C++; c-basic-offset:4; tab-width:4 -*-
//opt: -g
//
// @author Masahiro Kasahara <masahiro.kasahara.ws>
//

#include <iostream>
#include <fstream>
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <getopt.h>
//#include <stackdump.h>
//#include <debug.h>

using namespace std;

static const size_t BUFFER_SIZE = 16 * 1024 * 1024u;

static void output_read_name(const char* header)
{
	if(*header++ == '\0') return;
	while(*header != '\0' && *header != ' ') cout << *header++;
	cout << '\n';
}

void show_read_names_in_file(const char* fname, bool show_name) // if show_name is false, output read length
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
		if(show_name) output_read_name(b);
        size_t number_of_nucleotides_in_read = 0;
        if(b[0] != '@') { 
			// This should be FASTA
            while(ist.getline(b, BUFFER_SIZE)) {
				++line_count;
                if(b[0] == '>') {
                    if(show_name)
                        output_read_name(b);
                    else
                        cout << number_of_nucleotides_in_read << "\n";
                    number_of_nucleotides_in_read = 0;
                } else {
                    number_of_nucleotides_in_read += strlen(b);
                }
			}
            if(!show_name && 0 < number_of_nucleotides_in_read)
                cout << number_of_nucleotides_in_read << "\n";
		} else {
			// This should be FASTQ
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
                    if(!ist.getline(b, BUFFER_SIZE)) break;
                    ++line_count;
					if(show_name)
                        output_read_name(b);
                    else
                        cout << number_of_nucleotides_in_read << "\n";
                    number_of_nucleotides_in_read = 0;
                } else {
                    const size_t number_of_nucleotides_in_line = strlen(b);
                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                }
            }
            if(!show_name && 0 < number_of_nucleotides_in_read)
                cout << number_of_nucleotides_in_read << "\n";
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

static string get_read_name_from_header(const char* header)
{
    if(*header++ == '\0') return "";
    string rv;
    rv.reserve(64);
	while(*header != '\0' && *header != ' ') rv += *header++;
    return rv;
}

static void add_read_name_and_show_error_if_duplicates(map<string, int>& readName2fileIndex, char** argv, char* headerLine, int fileIndex, bool doNotShowFileName)
{
	if(*headerLine++ == '\0') return;
    char* start = headerLine;
    char* end = start;
    while(*end != '\0' && *end != ' ') ++end;
    const char saved_char = *end; *end = '\0';
    string readName = start;
    const map<string, int>::const_iterator cit = readName2fileIndex.find(readName);
    if(cit == readName2fileIndex.end()) {
        readName2fileIndex[readName] = fileIndex;
    } else {
        cerr << start;
        if(!doNotShowFileName)
            cerr << '\t' << argv[cit->second];
        cerr << '\n';
    }
    *end = saved_char;
}

void do_check_same_names(int argc, char** argv)
{
    static struct option long_options[] = {
        {"read", no_argument , 0, 'r'},
        {0, 0, 0, 0} // end of long options
    };
    bool flag_do_not_show_file_name = false;
    while(true) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "", long_options, &option_index);
		if(c == -1) break;
		switch(c) {
		case 0:
			// you can see long_options[option_index].name/flag and optarg (null if no argument).
			break;
		case 'r':
            flag_do_not_show_file_name = true;
			break;
		}
	}
    char* b = new char[BUFFER_SIZE];
    map<string, int> readName2fileIndex;
    for(int findex = optind + 1; findex < argc; ++findex) {
        const char* file_name = argv[findex];
        ifstream ist(file_name);
        if(!ist) {
            cerr << "Cannot open '" << file_name << "'" << endl;
        }
        size_t number_of_sequences = 0;
        size_t number_of_nucleotides = 0;
        size_t line_count = 0;
        if(ist.getline(b, BUFFER_SIZE)) {
            ++line_count;
            number_of_sequences++;
            size_t number_of_nucleotides_in_read = 0;
			add_read_name_and_show_error_if_duplicates(readName2fileIndex, argv, b, findex, flag_do_not_show_file_name);
            if(b[0] != '@') { 
                // This should be FASTA
                while(ist.getline(b, BUFFER_SIZE)) {
                    ++line_count;
                    if(b[0] == '>') {
                        number_of_sequences++;
						add_read_name_and_show_error_if_duplicates(readName2fileIndex, argv, b, findex, flag_do_not_show_file_name);
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
                        number_of_nucleotides_in_read = 0;
                        if(!ist.getline(b, BUFFER_SIZE))
                            break;
                        ++line_count;
                        ++number_of_sequences;
						add_read_name_and_show_error_if_duplicates(readName2fileIndex, argv, b, findex, flag_do_not_show_file_name);
                    } else {
                        const size_t number_of_nucleotides_in_line = strlen(b);
                        number_of_nucleotides += number_of_nucleotides_in_line;
                        number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                    }
                }
            }
        }
    }
    delete b;
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
		show_read_names_in_file(argv[i], true);
	}
}

void do_len(int argc, char** argv)
{
	for(int i = 2; i < argc; ++i) {
		show_read_names_in_file(argv[i], false);
	}
}

void do_extract(int argc, char** argv)
{
    bool flag_reverse_condition = false;
    bool flag_read_from_stdin = false;
    bool flag_output_unique = false;

    static struct option long_options[] = {
        {"reverse", no_argument , 0, 'r'},
        {"seq", required_argument, 0, 's'},
        {"file", required_argument, 0, 'f'},
        {"stdin", no_argument, 0, 'c'},
        {"unique", no_argument, 0, 'u'},
        {0, 0, 0, 0} // end of long options
    };

    set<string> readNamesToTake;
    vector<string> fileInputs;

    while(true) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "", long_options, &option_index);
		if(c == -1) break;
		switch(c) {
		case 0:
			// you can see long_options[option_index].name/flag and optarg (null if no argument).
			break;
		case 'r':
            flag_reverse_condition = true;
			break;
        case 's':
            readNamesToTake.insert(optarg);
            break;
        case 'f':
            fileInputs.push_back(optarg);
            break;
        case 'c':
            flag_read_from_stdin = true;
            break;
        case 'u':
            flag_output_unique = true;
            break;
		}
	}
    if(flag_output_unique) {
        flag_reverse_condition = !flag_reverse_condition;
    }
    if(flag_read_from_stdin) {
        string line;
        while(getline(cin, line)) {
            readNamesToTake.insert(line);
        }
    }
    {
        for(int i = 0; i < fileInputs.size(); ++i) {
            ifstream ist(fileInputs[i].c_str());
            if(!ist) {
                cerr << "ERROR: Cannot open '" << fileInputs[i] << "'" << endl;
                return;
            }
            string line;
            while(getline(ist, line)) {
                readNamesToTake.insert(line);
            }
        }
    }
    char* b = new char[BUFFER_SIZE];
    for(int findex = optind + 1; findex < argc; ++findex) {
        const char* file_name = argv[findex];
        ifstream ist(file_name);
        if(!ist) {
            cerr << "Cannot open '" << file_name << "'" << endl;
        }
        size_t number_of_sequences = 0;
        size_t number_of_nucleotides = 0;
        size_t line_count = 0;
        bool current_read_has_been_taken = false;
        if(ist.getline(b, BUFFER_SIZE)) {
            ++line_count;
            number_of_sequences++;
            size_t number_of_nucleotides_in_read = 0;
            current_read_has_been_taken = (readNamesToTake.count(get_read_name_from_header(b)) != 0) ^ flag_reverse_condition;
            if(current_read_has_been_taken) cout << b << endl;
            if(flag_output_unique) readNamesToTake.insert(get_read_name_from_header(b));
            if(b[0] != '@') { 
                // This should be FASTA
                while(ist.getline(b, BUFFER_SIZE)) {
                    ++line_count;
                    if(b[0] == '>') {
                        number_of_sequences++;
                        current_read_has_been_taken = (readNamesToTake.count(get_read_name_from_header(b)) != 0) ^ flag_reverse_condition;
                        if(current_read_has_been_taken) cout << b << endl;
                        if(flag_output_unique) readNamesToTake.insert(get_read_name_from_header(b));
                        number_of_nucleotides_in_read = 0;
                    } else {
                        if(current_read_has_been_taken) cout << b << endl;
                        number_of_nucleotides += strlen(b);
                    }
                }
            } else {
                // This is FASTQ
                while(ist.getline(b, BUFFER_SIZE)) {
                    ++line_count;
                    if(b[0] == '+' && b[1] == '\0') { // EOS
                        long long n = number_of_nucleotides_in_read;
                        if(current_read_has_been_taken) cout << b << endl;
                        while(ist.getline(b, BUFFER_SIZE)) {
                            ++line_count;
                            const size_t number_of_qvchars_in_line = strlen(b);
                            n -= number_of_qvchars_in_line;
                            if(current_read_has_been_taken) cout << b << endl;
                            if(n <= 0) break;
                        }
                        if(ist.peek() != '@' && !ist.eof()) {
                            cerr << "WARNING: bad file format? at line " << line_count << endl;
                        }
                        number_of_nucleotides_in_read = 0;
                        if(!ist.getline(b, BUFFER_SIZE))
                            break;
                        ++line_count;
                        ++number_of_sequences;
                        current_read_has_been_taken = (readNamesToTake.count(get_read_name_from_header(b)) != 0) ^ flag_reverse_condition;
                        if(current_read_has_been_taken) cout << b << endl;
                        if(flag_output_unique) readNamesToTake.insert(get_read_name_from_header(b));
                    } else {
                        const size_t number_of_nucleotides_in_line = strlen(b);
                        number_of_nucleotides += number_of_nucleotides_in_line;
                        number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                        if(current_read_has_been_taken) cout << b << endl;
                    }
                }
            }
        }
    }
    delete b;
}

void show_usage()
{
    cerr << "Usage: fatt <command> [options...]" << endl;
}

void show_help(const char* subcommand)
{
	const string subcmd = subcommand;
	if(subcmd == "count") {
        cerr << "Usage: fatt count [options...] <FAST(A|Q) files>\n\n";
        cerr << "Currently, no options available." << endl;
        return;
	}
    if(subcmd == "name") {
        cerr << "Usage: fatt name [options...] <FAST(A|Q) files>\n\n";
        cerr << "Currently, no options available." << endl;
        return;
	}
    if(subcmd == "chksamename") {
        cerr << "Usage: fatt chksamename [options...] <FAST(A|Q) files>\n\n";
        cerr << "--read\tOutput only read names. (Omit file names)\n";
        return;
    }
    if(subcmd == "extract") {
        cerr << "Usage: fatt extract [options...] <FAST(A|Q) files>\n\n";
        cerr << "--unique\tOutput only unique reads. Reads with the same read name are removed.\n";
        cerr << "--seq\tSpecify the name of the read to be retrieved. You can specify this option as many times as you wish.\n";
        cerr << "--file\tSpecify a file in which you listed the read names. One line, one read.\n";
        cerr << "--stdin\tRead the list of read names from stdin. It may be useful when you combine with *NIX pipe.\n";
        cerr << "--reverse\tReverse the extracting condition. It is like -v option of grep.\n";
        return;
    }
    if(subcmd == "len") {
        cerr << "Usage: fatt len [options...] <FAST(A|Q) files>\n\n";
        cerr << "Currently, no options available." << endl;
        return;
    }
    if(subcmd == "help") {
        cerr << "Uh? No detailed help for help.\n";
        cerr << "Read the manual, or ask the author.\n";
        return;
    }
    show_usage();
    cerr << "\nCommands:\n";
	cerr << "\tcount\tcount the number of reads/nucleotides\n";
	cerr << "\tname\toutput the names of reads\n";
    cerr << "\tchksamename\toutput the names of reads if the read name is duplicated\n";
	cerr << "\textract\textract a set of reads with condition\n";
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
    if(commandString == "chksamename") {
        do_check_same_names(argc, argv);
        return;
    }
	if(commandString == "extract") {
		do_extract(argc, argv);
        return;
	}
    if(commandString == "len") {
        do_len(argc, argv);
        return;
    }
    // Help or error.
    if(commandString != "help") {
        cerr << "ERROR: Unknown command '" << commandString << "'" << endl;
    }
	if(argc < 3) {
    	show_help("");
	} else {
		show_help(argv[2]);
	}
}

int main(int argc, char** argv)
{
	//GDB_On_SEGV gos(argv[0]);
    if(argc < 2) {
        show_usage();
        return 1;
    }
    const char* commandStr = argv[1];
    dispatchByCommand(commandStr, argc, argv);
    return 0;
}
