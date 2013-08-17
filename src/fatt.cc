// fatt: FASTA/FASTQ manipulation tool.
//
// -*- mode:C++; c-basic-offset:4; tab-width:4 -*-
//
// @author Masahiro Kasahara <masahiro.kasahara.ws>
//

#include <iostream>
#include <fstream>
#include <string.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <map>
#include <set>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <getopt.h>
#include <unistd.h>
#include "sqdb.h"
//#include <stackdump.h>
//#include <debug.h>

using namespace std;

static const size_t BUFFER_SIZE = 16 * 1024 * 1024u;

struct CSVEscape
{
	const char* p;
	CSVEscape(const char* p) : p(p) {}
};

static ostream& operator << (ostream& os, const CSVEscape& c)
{
	for(const char* p = c.p; *p; ++p) {
		if(*p == '"') os << '"';
		os << *p;
	}
	return os;
}

static size_t strlen_without_n(const char * p)
{
    size_t retval = 0;
    for(; *p; p++){ if(!(*p == 'N' || *p == 'n')) retval++; }
    return retval;
}

static void output_read_name(const char* header)
{
	if(*header++ == '\0') return;
	while(*header != '\0' && *header != ' ') cout << *header++;
	cout << '\n';
}

static string get_index_file_name(const char* fastq_file_name)
{
    string index_file_name = fastq_file_name;
    index_file_name += ".index";
    return index_file_name;
}

static bool doesIndexExist(const char* fastq_file_name)
{
    return access(get_index_file_name(fastq_file_name).c_str(), F_OK) == 0;
}

static bool is_file_fastq(const char* fastq_file_name)
{
    ifstream ist(fastq_file_name);
    if(!ist) return false;
    string line;
    if(!getline(ist, line)) return false;
    if(line.empty() || line[0] != '@') return false;
    return true;
}

void calculate_n50_statistics(const char* fname,
                              vector<size_t>& length_of_scaffolds_wgap,
                              vector<size_t>& length_of_scaffolds_wogap,
                              vector<size_t>& length_of_contigs)
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
        size_t length_as_scaffolds_wgap = 0;
        size_t length_as_scaffolds_wogap = 0;
        size_t length_as_contigs = 0;
        if(b[0] != '@') { 
            if(b[0] != '>') { // NOTE: this check must be copied to other functions.
                cerr << fname << " does not look like FASTA/FASTQ!\n";
                return;
            }
			// This should be FASTA
            while(ist.getline(b, BUFFER_SIZE)) {
				++line_count;
                if(b[0] == '>') {
                    length_of_scaffolds_wgap.push_back(length_as_scaffolds_wgap);
                    length_of_scaffolds_wogap.push_back(length_as_scaffolds_wogap);
                    length_as_scaffolds_wgap = 0;
                    length_as_scaffolds_wogap = 0;
                    if(0 < length_as_contigs) {
                        length_of_contigs.push_back(length_as_contigs);
                        length_as_contigs = 0;
                    }
                } else {
                    length_as_scaffolds_wgap += strlen(b);
                    length_as_scaffolds_wogap += strlen_without_n(b);
                    for(const char* p = b; *p != '\0'; p++) {
                        if(*p == 'N' || *p == 'n') {
                            if(0 < length_as_contigs) {
                                length_of_contigs.push_back(length_as_contigs);
                                length_as_contigs = 0;
                            }
                        } else {
                            length_as_contigs++;
                        }
                    }
                }
			}
            length_of_scaffolds_wgap.push_back(length_as_scaffolds_wgap);
            length_of_scaffolds_wogap.push_back(length_as_scaffolds_wogap);
            if(0 < length_as_contigs) {
                length_of_contigs.push_back(length_as_contigs);
            }
		} else {
			// This should be FASTQ
            while(ist.getline(b, BUFFER_SIZE)) {
				++line_count;
                if(b[0] == '+' && b[1] == '\0') { // EOS
                    long long n = length_as_scaffolds_wgap;
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
                    length_of_scaffolds_wgap.push_back(length_as_scaffolds_wgap);
                    length_of_scaffolds_wogap.push_back(length_as_scaffolds_wogap);
                    length_as_scaffolds_wgap = 0;
                    length_as_scaffolds_wogap = 0;
                    if(0 < length_as_contigs) {
                        length_of_contigs.push_back(length_as_contigs);
                        length_as_contigs = 0;
                    }
                } else {
                    const size_t number_of_nucleotides_in_line = strlen(b);
                    length_as_scaffolds_wgap += number_of_nucleotides_in_line;
                    length_as_scaffolds_wogap += strlen_without_n(b);
                    for(const char* p = b; *p != '\0'; p++) {
                        if(*p == 'N' || *p == 'n') {
                            if(0 < length_as_contigs) {
                                length_of_contigs.push_back(length_as_contigs);
                                length_as_contigs = 0;
                            }
                        } else {
                            length_as_contigs++;
                        }
                    }
                }
            }
            length_of_scaffolds_wgap.push_back(length_as_scaffolds_wgap);
            length_of_scaffolds_wogap.push_back(length_as_scaffolds_wogap);
            if(0 < length_as_contigs) {
                length_of_contigs.push_back(length_as_contigs);
            }
		}
	}
    delete b;
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
            continue;
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

struct DeleteOnFailure {
	const string filename;
	bool hasToDelete;
	DeleteOnFailure(const string& filename) : filename(filename), hasToDelete(true) {}
	void doNotDelete() { hasToDelete = false; }
	~DeleteOnFailure() {
		if(hasToDelete) {
			unlink(filename.c_str());
		}
	}
};

void create_index(const char* fname)
{
    const string index_file_name = get_index_file_name(fname);
    if(doesIndexExist(fname)) {
        cerr << "'" << index_file_name << "' already exists!" << endl;
        return;
    }
    cerr << "Creating index" << flush;
    try {
		DeleteOnFailure dof(index_file_name);
    	sqdb::Db db(index_file_name.c_str());
    	db.MakeItFasterAndDangerous();
    	db.Do("create table seqpos(name text primary key, pos integer, readindex integer)");
    	db.Do("begin");
    	cerr << "." << flush;
        sqdb::Statement stmt = db.Query("insert into seqpos values(?, ?, ?)");
        ifstream ist(fname);
        if(!ist) {
            cerr << "Cannot open '" << fname << "'" << endl;
            return;
        }
        char* b = new char[BUFFER_SIZE];
        size_t line_count = 0;
		long long sequence_count = 0;
        iostream::pos_type last_pos = ist.tellg();
        if(ist.getline(b, BUFFER_SIZE)) {
            ++line_count;
            #define INSERT_NAME_INTO_TABLE() { stmt.Bind(1, get_read_name_from_header(b)); stmt.Bind(2, static_cast<long long>(last_pos)); stmt.Bind(3, sequence_count); stmt.Next(); }
            INSERT_NAME_INTO_TABLE();
			++sequence_count;
            size_t number_of_nucleotides_in_read = 0;
            if(b[0] != '@') { 
                // This should be FASTA
                last_pos = ist.tellg();
                while(ist.getline(b, BUFFER_SIZE)) {
                    ++line_count;
                    if(b[0] == '>') {
                        INSERT_NAME_INTO_TABLE();
						++sequence_count;
                        number_of_nucleotides_in_read = 0;
                    } else {
                        number_of_nucleotides_in_read += strlen(b);
                    }
                }
            } else {
                // This should be FASTQ
                last_pos = ist.tellg();
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
                        last_pos = ist.tellg();
                        if(!ist.getline(b, BUFFER_SIZE)) break;
                        ++line_count;
                        INSERT_NAME_INTO_TABLE();
						++sequence_count;
                        number_of_nucleotides_in_read = 0;
                    } else {
                        const size_t number_of_nucleotides_in_line = strlen(b);
                        number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                    }
                }
            }
            #undef INSERT_NAME_INTO_TABLE
        }
        delete b;
    	db.Do("end");
    	cerr << "." << flush;
    	db.Do("create index seqpos_name_index on seqpos(name);");
    	cerr << "." << flush;
    	db.Do("create index read_index_index on seqpos(readindex);");
    	cerr << endl;
		dof.doNotDelete();
    } catch(size_t line_num) {
        cerr << "DB Creation Error. (insert) at line " << line_num << endl;
    } catch(const sqdb::Exception& e) {
		cerr << "DB Creation Error. " << e.GetErrorMsg() << endl;
	}
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

void print_n50(vector<size_t>& lengths, const bool flag_html, const bool flag_json, const bool is_contig = false)
{
    sort(lengths.rbegin(), lengths.rend());
    const size_t total_length = accumulate(lengths.begin(), lengths.end(), 0ull);
    size_t sequence_index = (size_t)-1;
    size_t n50_length = 0;
    size_t max_length = 0;
    size_t min_length = 0;
    size_t avg_length = 0;
    if(!lengths.empty()) {
        size_t sum = 0;
        const size_t n50_total_length = (total_length + 1ull) / 2;
        for(sequence_index = 0; sequence_index < lengths.size(); sequence_index++) {
            sum += lengths[sequence_index];
            if(n50_total_length <= sum) break;
        }
        max_length = lengths.front();
        n50_length = lengths[sequence_index];
        min_length = lengths.back();
        avg_length = (total_length + lengths.size() / 2) / lengths.size();
    }
    if(flag_html) {
        cout << "<tr><td></td><td>size (bp)</td><td>number</td></tr>\n";
        cout << "<tr><td>max</td><td>" << max_length << "</td><td>1</td></tr>\n";
        cout << "<tr><td>N50</td><td>" << n50_length << "</td><td>" << (sequence_index + 1) << "</td></tr>\n";
        cout << "<tr><td>min</td><td>" << min_length << "</td><td>" << lengths.size() << "</td></tr>\n";
        cout << "<tr><td>avg</td><td>" << avg_length << "</td><td></td></tr>\n";
    } else if(flag_json) {
        cout << "{\"total_length\": " << total_length;
        cout << ",\"max_length\": " << max_length;
        cout << ",\"n50\": " << n50_length;
        cout << ",\"n50num\": " << (sequence_index + 1);
        cout << ",\"avg\": " << avg_length;
        cout << ",\"min\": " << min_length;
        cout << "}";
    } else {
        cout << "Total # of bases = " << total_length << "\n";
        cout << "Max size = " << max_length << "\n";
        const string scaffold_str = is_contig ? "contig" : "scaffold";
        cout << "N50 " << scaffold_str << " size = " << n50_length << "\n";
        cout << "N50 " << scaffold_str << " # = " << (sequence_index + 1) << "\n";
        cout << "Avg size = " << avg_length << "\n";
        cout << "Min size = " << min_length << "\n";
    }
}

void do_stat(int argc, char** argv)
{
	bool flag_html = false;
    bool flag_json = false;
    static struct option long_options[] = {
        {"html", no_argument, 0, 'h'},
        {"json", no_argument, 0, 'j'},
        {0, 0, 0, 0} // end of long options
    };
    while(true) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "", long_options, &option_index);
		if(c == -1) break;
		switch(c) {
		case 0:
			// you can see long_options[option_index].name/flag and optarg (null if no argument).
			break;
		case 'h':
            flag_html = true;
			break;
        case 'j':
            flag_json = true;
            break;
        }
	}
    if(flag_html && flag_json) {
        cerr << "ERROR: You can use either --html or --json\n";
        return;
    }
    vector<size_t> length_of_scaffolds_wgap;
    vector<size_t> length_of_scaffolds_wogap;
    vector<size_t> length_of_contigs;
	for(int i = optind + 1; i < argc; ++i) {
		calculate_n50_statistics(argv[i], length_of_scaffolds_wgap, length_of_scaffolds_wogap, length_of_contigs);
	}
    if(!(length_of_scaffolds_wgap.size() == length_of_scaffolds_wogap.size())) {
        cerr << "Assertion failed. Maybe you found a bug! Please report to the author.\n";
        return;
    }
    if(flag_html) {
        cout << "<table border=\"2\" bgcolor=\"#ffffff\">\n";
        cout << "<tr><th colspan=\"3\" bgcolor=\"#fdfdd4\">Scaffold (w/gap) statistics</th></tr>\n";
    } else if(flag_json) {
        cout << "{\"scaffold_wgap\": ";
    } else {
        cout << "Scaffold (w/gap) statistics\n";
    }
    print_n50(length_of_scaffolds_wgap, flag_html, flag_json, false);
    if(flag_html) {
        cout << "<tr><th colspan=\"3\" bgcolor=\"#fdfdd4\">Scaffold (wo/gap) statistics</th></tr>\n";
    } else if(flag_json) {
        cout << ",\"scaffold_wogap\": ";
    } else {
        cout << "\nScaffold (wo/gap) statistics\n";
    }
    print_n50(length_of_scaffolds_wogap, flag_html, flag_json, false);
    if(flag_html) {
        cout << "<tr><th colspan=\"3\" bgcolor=\"#fdfdd4\">Contig statistics</th></tr>\n";
    } else if(flag_json) {
        cout << ",\"contig\": ";
    } else {
        cout << "\nContig statistics\n";
    }
    print_n50(length_of_contigs, flag_html, flag_json, true);
    if(flag_html) {
        cout << "</table>\n";
    } else if(flag_json) {
        cout << "}\n";
    }
}

void do_index(int argc, char** argv)
{
	for(int i = 2; i < argc; ++i) {
		create_index(argv[i]);
	}
}

void do_extract(int argc, char** argv)
{
    bool flag_reverse_condition = false;
    bool flag_read_from_stdin = false;
    bool flag_output_unique = false;
	bool flag_noindex = false;
	bool flag_index = false;
	long long param_start = -1;
	long long param_end = -1;
	long long param_num = -1;

    static struct option long_options[] = {
        {"reverse", no_argument , 0, 'r'},
        {"seq", required_argument, 0, 's'},
        {"file", required_argument, 0, 'f'},
        {"stdin", no_argument, 0, 'c'},
        {"unique", no_argument, 0, 'u'},
		{"index", no_argument, 0, 'i'},
		{"noindex", no_argument, 0, 'n'},
    	{"start", required_argument, 0, 'a'},
    	{"end", required_argument, 0, 'e'},
    	{"num", required_argument, 0, 'q'},
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
		case 'n':
			flag_noindex = true;
			break;
		case 'i':
			flag_index = true;
			break;
		case 'a':
			param_start = atoll(optarg);
			break;
		case 'e':
			param_end = atoll(optarg);
			break;
		case 'q':
			param_num = atoll(optarg);
			break;
		}
	}
	if(flag_index && flag_noindex) {
		cerr << "ERROR: do not specify --index and --noindex at once" << endl;
		return;
	}
	if(param_end != -1 && param_num != -1) {
		cerr << "ERROR: you cannot specify --end and --num at once" << endl;
		return;
	}
    if(param_num != -1) {
        if(param_start == -1) {
            param_start = 0;
            param_end = param_num;
        } else {
            param_end = param_start + param_num;
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
    if(!readNamesToTake.empty() && param_start != -1) {
        cerr << "ERROR: you can either select the range or the sequence names, but not both." << endl;
        return;
    }
    char* b = new char[BUFFER_SIZE];
    for(int findex = optind + 1; findex < argc; ++findex) {
        const char* file_name = argv[findex];
		if(flag_index && !doesIndexExist(file_name)) {
			create_index(file_name);
		}
		const bool use_index = (flag_index || (!flag_noindex && doesIndexExist(file_name))) && !flag_reverse_condition;
        ifstream ist(file_name);
        if(!ist) {
            cerr << "Cannot open '" << file_name << "'" << endl;
            continue;
        }
        if(use_index) {
            const bool is_fastq = is_file_fastq(file_name);
            const string index_file_name = get_index_file_name(file_name);
            try {
                sqdb::Db db(index_file_name.c_str());
                if(param_start == -1) {
                    sqdb::Statement stmt = db.Query("select pos from seqpos where name=?");
                    for(set<string>::const_iterator it = readNamesToTake.begin(); it != readNamesToTake.end(); ++it) {
                        const string& read_name = *it;
                        stmt.Bind(1, read_name);
                        if(stmt.Next()) {
                            const long long pos = stmt.GetField(0);
                            ist.seekg(pos);
                            if(ist.fail() || !ist.getline(b, BUFFER_SIZE)) {
                                cerr << "WARNING: " << read_name << " is missing in the file. Maybe the index is old?\n";
                                continue;
                            }
                            cout << b << "\n";
                            if(is_fastq) {
                                size_t number_of_nucleotides_in_read = 0;
                                while(ist.getline(b, BUFFER_SIZE)) {
                                    if(b[0] == '+' && b[1] == '\0') { // EOS
                                        long long n = number_of_nucleotides_in_read;
                                        cout << b << endl;
                                        while(ist.getline(b, BUFFER_SIZE)) {
                                            const size_t number_of_qvchars_in_line = strlen(b);
                                            n -= number_of_qvchars_in_line;
                                            cout << b << endl;
                                            if(n <= 0) break;
                                        }
                                        break;
                                    } else {
                                        const size_t number_of_nucleotides_in_line = strlen(b);
                                        number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                                        cout << b << endl;
                                    }
                                }
                            } else {
                                while(ist.getline(b, BUFFER_SIZE)) {
                                    if(b[0] == '>') break;
                                    cout << b << "\n";
                                }
                            }
                        } else {
                            cerr << "WARNING: " << read_name << " was not found.\n";
                        }
                    }
                } else {
                    sqdb::Statement stmt = db.Query("select pos from seqpos limit 1 offset ?");
                    long long sequence_index = param_start;
                    stmt.Bind(1, param_start);
                    if(stmt.Next()) {
                        const long long pos = stmt.GetField(0);
                        ist.seekg(pos);
                        if(ist.fail() || !ist.getline(b, BUFFER_SIZE)) {
                            cerr << "WARNING: Cannot seek to that far. Maybe the index is old?\n";
                            continue;
                        }
                        cout << b << "\n";
                        if(is_fastq) {
                            size_t number_of_nucleotides_in_read = 0;
                            while(ist.getline(b, BUFFER_SIZE)) {
                                if(b[0] == '+' && b[1] == '\0') { // EOS
                                    long long n = number_of_nucleotides_in_read;
                                    cout << b << endl;
                                    while(ist.getline(b, BUFFER_SIZE)) {
                                        const size_t number_of_qvchars_in_line = strlen(b);
                                        n -= number_of_qvchars_in_line;
                                        cout << b << endl;
                                        if(n <= 0) break;
                                    }
                                    sequence_index++;
                                    if(param_end <= sequence_index)
                                        break;
                                    const long long line_start_pos = ist.tellg();
                                    if(!ist.getline(b, BUFFER_SIZE)) {
                                        cerr << "WARNING: reached the end of file.\n";
                                        return;
                                    }
                                    if(b[0] != '@') {
                                        cerr << "ERROR: bad file format. The line does not start with '@' at pos " << line_start_pos << endl;
                                        return;
                                    }
                                    cout << b << "\n";
									number_of_nucleotides_in_read = 0;
                                } else {
                                    const size_t number_of_nucleotides_in_line = strlen(b);
                                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                                    cout << b << endl;
                                }
                            }
                        } else {
                            while(ist.getline(b, BUFFER_SIZE)) {
                                if(b[0] == '>') {
                                    sequence_index++;
                                    if(param_end <= sequence_index)
                                        break;
                                }
                                cout << b << "\n";
                            }
                        }
                    } else {
                        cerr << "WARNING: the start index (" << param_start << ") is larger than the number of sequences in the file." << endl;
                    }
                }
            } catch (const sqdb::Exception& e) {
                cerr << "ERROR: db error. " << e.GetErrorMsg() << endl;
                return;
            }
        } else {
            size_t number_of_sequences = 0;
            size_t number_of_nucleotides = 0;
            size_t line_count = 0;
            bool current_read_has_been_taken = false;
            if(ist.getline(b, BUFFER_SIZE)) {
                ++line_count;
                number_of_sequences++;
                size_t number_of_nucleotides_in_read = 0;
                if(param_start == -1) {
                    current_read_has_been_taken = (readNamesToTake.count(get_read_name_from_header(b)) != 0) ^ flag_reverse_condition;
                } else {
                    current_read_has_been_taken = param_start + 1 <= number_of_sequences && number_of_sequences <= param_end;
                    // NOTE: the latter condition never hold, if I properly implemented.
                }
                if(current_read_has_been_taken) cout << b << endl;
                if(flag_output_unique) readNamesToTake.insert(get_read_name_from_header(b));
                if(b[0] != '@') { 
                    // This should be FASTA
                    while(ist.getline(b, BUFFER_SIZE)) {
                        ++line_count;
                        if(b[0] == '>') {
                            number_of_sequences++;
                            if(param_start == -1) {
                                current_read_has_been_taken = (readNamesToTake.count(get_read_name_from_header(b)) != 0) ^ flag_reverse_condition;
                            } else {
                                current_read_has_been_taken = param_start + 1 <= number_of_sequences && number_of_sequences <= param_end;
                                // NOTE: the latter condition never hold, if I properly implemented.
                            }
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
                            if(param_start == -1) {
                                current_read_has_been_taken = (readNamesToTake.count(get_read_name_from_header(b)) != 0) ^ flag_reverse_condition;
                            } else {
                                current_read_has_been_taken = param_start + 1 <= number_of_sequences && number_of_sequences <= param_end;
                                // NOTE: the latter condition never hold, if I properly implemented.
                            }
                            if(param_end < number_of_sequences) // NOTE: param_end is 0-origin, number_of_sequences is 1-origin.
                                break;
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
    }
    delete b;
}

void do_guess_qv_type(int argc, char** argv)
{
    char* b = new char[BUFFER_SIZE];

    for(int findex = 2; findex < argc; ++findex) {
        const char* file_name = argv[findex];
        ifstream ist(file_name);
        if(!ist) {
            cerr << "Cannot open '" << file_name << "'" << endl;
            continue;
        }
        size_t histogram[256];
        for(int i = 0; i < 256; ++i) histogram[i] = 0;
        size_t number_of_sequences = 0;
        size_t number_of_nucleotides = 0;
        size_t line_count = 0;
        if(ist.getline(b, BUFFER_SIZE)) {
            ++line_count;
            number_of_sequences++;
            size_t number_of_nucleotides_in_read = 0;
            if(b[0] != '@') { 
                cerr << "ERROR: the input file '" << file_name << "' does not seem to be a FASTQ file at line " << line_count << endl;
                return;
            }
            while(ist.getline(b, BUFFER_SIZE)) {
                ++line_count;
                if(b[0] == '+' && b[1] == '\0') { // EOS
                    long long n = number_of_nucleotides_in_read;
                    while(ist.getline(b, BUFFER_SIZE)) {
                        ++line_count;
                        const size_t number_of_qvchars_in_line = strlen(b);
                        for(unsigned char* p = reinterpret_cast<unsigned char*>(b); *p; ++p) histogram[*p]++;
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
                } else {
                    const size_t number_of_nucleotides_in_line = strlen(b);
                    number_of_nucleotides += number_of_nucleotides_in_line;
                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                }
            }
        }
        {
            size_t numBadQVchars = 0;
            for(int i = 0; i < 32; ++i)    numBadQVchars += histogram[i];
            for(int i = 127; i < 256; ++i) numBadQVchars += histogram[i];
            if(0 < numBadQVchars) {
                cout << file_name << '\t' << "bad\tthe QV strings have " << numBadQVchars << " characters that do not look like QV chars.\n";
                continue;
            }
        }
        {
            size_t numSangerOnlyQVchars = 0;
            for(int i = 33; i <= 58; ++i)  numSangerOnlyQVchars += histogram[i];
            if(0 < numSangerOnlyQVchars) {
                cout << file_name << '\t' << "sanger\tIt must be Sanger FASTQ or Illumina 1.8+.\n";
                continue;
            }
        }
        {
            size_t numSolexaOnlyQVchars = 0;
            for(int i = 59; i <= 63; ++i)  numSolexaOnlyQVchars += histogram[i];
            if(0 < numSolexaOnlyQVchars) {
                cout << file_name << '\t' << "solexa\tIt must be Solexa FASTQ.\n";
                continue;
            }
        }
        if(0 < histogram[64] + histogram[65]) {
            cout << file_name << '\t' << "illumina13\tIt must be Illumina 1.3+\n";
            continue;
        }
        cout << file_name << '\t' << "illumina15\tIt looks like Illumina 1.5+\n";
    }
    delete b;
}

void to_csv(const char* file_name, bool does_not_output_header, bool output_in_tsv)
{
    ifstream ist(file_name);
    if(!ist) {
        cerr << "Cannot open '" << file_name << "'" << endl;
        return;
    }
    char* b = new char[BUFFER_SIZE];
    size_t number_of_sequences = 0;
    size_t number_of_nucleotides = 0;
    size_t line_count = 0;
    if(ist.getline(b, BUFFER_SIZE)) {
        ++line_count;
        number_of_sequences++;
        size_t number_of_nucleotides_in_read = 0;
        if(b[0] != '@') { 
            // This should be FASTA
            if(!does_not_output_header) {
                if(output_in_tsv) {
                    cout << "id\tdesc\tseq\n";
                } else {
                    cout << "id,desc,seq\n";
                }
            }
            #define OUTPUT_HEADER() {                                           \
                const string& idstr = get_read_name_from_header(b);\
                cout << idstr << (output_in_tsv ? "\t" : ",\"");\
                const char* descp = b + 1 + idstr.size();\
                if(*descp != '\0') descp++;\
                if(output_in_tsv) cout << descp; else cout << CSVEscape(descp);\
                cout << (output_in_tsv ? "\t" : "\",");\
            }
            OUTPUT_HEADER();
            while(ist.getline(b, BUFFER_SIZE)) {
                ++line_count;
                if(b[0] == '>') {
                    number_of_sequences++;
                    number_of_nucleotides_in_read = 0;
                    cout << '\n';
                    OUTPUT_HEADER();
                } else {
                    cout << b;
                    number_of_nucleotides += strlen(b);
                }
            }
            cout << '\n';
        } else {
            // This is FASTQ
            if(!does_not_output_header) {
                if(output_in_tsv) {
                    cout << "id\tdesc\tseq\tqv\n";
                } else {
                    cout << "id,desc,seq,qv\n";
                }
            }
            OUTPUT_HEADER();
            while(ist.getline(b, BUFFER_SIZE)) {
                ++line_count;
                if(b[0] == '+' && b[1] == '\0') { // EOS
                    cout << (output_in_tsv ? '\t' : ',') << '"';
                    long long n = number_of_nucleotides_in_read;
                    while(ist.getline(b, BUFFER_SIZE)) {
                        ++line_count;
                        const size_t number_of_qvchars_in_line = strlen(b);
						if(output_in_tsv)
                        	cout << b;
						else
							cout << CSVEscape(b);
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    if(ist.peek() != '@' && !ist.eof()) {
                        cerr << "WARNING: bad file format? at line " << line_count << endl;
                    }
                    number_of_nucleotides_in_read = 0;
                    if(!ist.getline(b, BUFFER_SIZE))
                        break;
                    cout << '"' << '\n';
                    OUTPUT_HEADER();
                    ++line_count;
                    ++number_of_sequences;
                } else {
                    const size_t number_of_nucleotides_in_line = strlen(b);
                    number_of_nucleotides += number_of_nucleotides_in_line;
                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                    cout << b;
                }
            }
            cout << '"' << '\n';
            #undef OUTPUT_HEADER
        }
    }
    delete b;
}

void fold_fastx(const char* file_name, int length_of_line, bool is_folding)
{
    ifstream ist(file_name);
    if(!ist) {
        cerr << "Cannot open '" << file_name << "'" << endl;
        return;
    }
    char* b = new char[BUFFER_SIZE];
    size_t line_count = 0;
    if(ist.getline(b, BUFFER_SIZE)) {
        ++line_count;
        size_t number_of_nucleotides_in_output_line = 0;
        cout << b << '\n';
        if(b[0] != '@') { 
            // This should be FASTA
            while(ist.getline(b, BUFFER_SIZE)) {
                ++line_count;
                if(b[0] == '>') {
                    if(0 < number_of_nucleotides_in_output_line) {
                        cout << '\n';
                        number_of_nucleotides_in_output_line = 0;
                    }
                    cout << b << '\n';
                } else {
                    const int number_of_chars_in_line = strlen(b);
                    if(is_folding) {
                        int off = 0;
                        while(off < number_of_chars_in_line) {
                            int s = number_of_chars_in_line - off;
                            if(length_of_line - number_of_nucleotides_in_output_line <= s) s = length_of_line - number_of_nucleotides_in_output_line;
                            for(int i = 0; i < s; ++i) cout << b[off + i];
                            off += s;
                            number_of_nucleotides_in_output_line += s;
                            if(length_of_line <= number_of_nucleotides_in_output_line) {
                                cout << '\n';
                                number_of_nucleotides_in_output_line = 0;
                            }
                        }
                    } else {
                        cout << b;
                        number_of_nucleotides_in_output_line += number_of_chars_in_line;
                    }
                }
            }
            if(0 < number_of_nucleotides_in_output_line)
                cout << '\n';
        } else {
            size_t number_of_nucleotides_in_read = 0;
            // This is FASTQ
            while(ist.getline(b, BUFFER_SIZE)) {
                ++line_count;
                if(b[0] == '+' && b[1] == '\0') { // EOS
                    if(0 < number_of_nucleotides_in_output_line) {
                        cout << '\n';
                        number_of_nucleotides_in_output_line = 0;
                    }
					cout << "+\n";
                    long long n = number_of_nucleotides_in_read;
                    while(ist.getline(b, BUFFER_SIZE)) {
                        ++line_count;
                        const size_t number_of_qvchars_in_line = strlen(b);
                        if(is_folding) {
                            int off = 0;
                            while(off < number_of_qvchars_in_line) {
                                int s = number_of_qvchars_in_line - off;
                                if(length_of_line - number_of_nucleotides_in_output_line <= s) s = length_of_line - number_of_nucleotides_in_output_line;
                                for(int i = 0; i < s; ++i) cout << b[off + i];
                                off += s;
                                number_of_nucleotides_in_output_line += s;
                                if(length_of_line <= number_of_nucleotides_in_output_line) {
                                    cout << '\n';
                                    number_of_nucleotides_in_output_line = 0;
                                }
                            }
                        } else {
                            cout << b;
                            number_of_nucleotides_in_output_line += number_of_qvchars_in_line;
                        }
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    if(ist.peek() != '@' && !ist.eof()) {
                        cerr << "WARNING: bad file format? at line " << line_count << endl;
                    }
                    number_of_nucleotides_in_read = 0;
                    if(!ist.getline(b, BUFFER_SIZE))
                        break;
                    if(0 < number_of_nucleotides_in_output_line) {
                        cout << '\n';
                        number_of_nucleotides_in_output_line = 0;
                    }
                    cout << b << '\n';
                    ++line_count;
                } else {
                    const int number_of_chars_in_line = strlen(b);
                    number_of_nucleotides_in_read += number_of_chars_in_line;
                    if(is_folding) {
                        int off = 0;
                        while(off < number_of_chars_in_line) {
                            int s = number_of_chars_in_line - off;
                            if(length_of_line - number_of_nucleotides_in_output_line <= s) s = length_of_line - number_of_nucleotides_in_output_line;
                            for(int i = 0; i < s; ++i) cout << b[off + i];
                            off += s;
                            number_of_nucleotides_in_output_line += s;
                            if(length_of_line <= number_of_nucleotides_in_output_line) {
                                cout << '\n';
                                number_of_nucleotides_in_output_line = 0;
                            }
                        }
                    } else {
                        cout << b;
                        number_of_nucleotides_in_output_line += number_of_chars_in_line;
                    }
                }
            }
            if(0 < number_of_nucleotides_in_output_line)
                cout << '\n';
        }
    }
    delete b;
}

void fastq_to_fasta(const char* file_name)
{
    ifstream ist(file_name);
    if(!ist) {
        cerr << "Cannot open '" << file_name << "'" << endl;
        return;
    }
    char* b = new char[BUFFER_SIZE];
    size_t line_count = 0;
    if(ist.getline(b, BUFFER_SIZE)) {
        ++line_count;
        size_t number_of_nucleotides_in_output_line = 0;
        if(b[0] != '@') { 
            if(b[0] == '>') {
                cerr << "Input is already FASTA." << endl;
            } else {
                cerr << "Input does not look like FASTA/FASTQ." << endl;
            }
            delete b;
            return;
        }
        // This is FASTQ
        b[0] = '>';
        cout << b << '\n';
        size_t number_of_nucleotides_in_read = 0;
        while(ist.getline(b, BUFFER_SIZE)) {
            ++line_count;
            if(b[0] == '+' && b[1] == '\0') { // EOS
                long long n = number_of_nucleotides_in_read;
                while(ist.getline(b, BUFFER_SIZE)) {
                    ++line_count;
                    const size_t number_of_qvchars_in_line = strlen(b);
                    number_of_nucleotides_in_output_line += number_of_qvchars_in_line;
                    n -= number_of_qvchars_in_line;
                    if(n <= 0) break;
                }
                if(ist.peek() != '@' && !ist.eof()) {
                    cerr << "WARNING: bad file format? at line " << line_count << endl;
                    delete b;
                    return;
                }
                number_of_nucleotides_in_read = 0;
                if(!ist.getline(b, BUFFER_SIZE))
                    break;
                b[0] = '>';
                cout << b << '\n';
                ++line_count;
            } else {
                const int number_of_chars_in_line = strlen(b);
                number_of_nucleotides_in_read += number_of_chars_in_line;
                cout << b << '\n';
            }
        }
    }
    delete b;
}

void clean_nucleotide_line(char* buffer, bool process_n, int char_change_into)
{
    for(char* p = buffer; *p != '\0'; p++) {
        const char c = toupper(*p);
        if(c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            if(process_n || c != 'N') {
                if(char_change_into != -1) {
                    *p = char_change_into;
                } else {
                    const char t = rand() / ((RAND_MAX / 2 + 1) / 2);
                    *p = "ACGTT"[t];
                }
            }
        }
    }
}

void clean_fastx(const char* file_name, bool flag_process_n, int char_change_into/* -1 means random*/)
{
    ifstream ist(file_name);
    if(!ist) {
        cerr << "Cannot open '" << file_name << "'" << endl;
        return;
    }
    srand(time(NULL));
    char* b = new char[BUFFER_SIZE];
    size_t line_count = 0;
    if(ist.getline(b, BUFFER_SIZE)) {
        ++line_count;
        cout << b << '\n';
        if(b[0] != '@') { 
            // This should be FASTA
            while(ist.getline(b, BUFFER_SIZE)) {
                ++line_count;
                if(b[0] != '>') { clean_nucleotide_line(b, flag_process_n, char_change_into); }
                cout << b << '\n';
            }
        } else {
            size_t number_of_nucleotides_in_read = 0;
            // This is FASTQ
            while(ist.getline(b, BUFFER_SIZE)) {
                ++line_count;
                if(b[0] == '+' && b[1] == '\0') { // EOS
					cout << "+\n";
                    long long n = number_of_nucleotides_in_read;
                    while(ist.getline(b, BUFFER_SIZE)) {
                        ++line_count;
                        const size_t number_of_qvchars_in_line = strlen(b);
                        cout << b << '\n';
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    if(ist.peek() != '@' && !ist.eof()) {
                        cerr << "WARNING: bad file format? at line " << line_count << endl;
                    }
                    number_of_nucleotides_in_read = 0;
                    if(!ist.getline(b, BUFFER_SIZE))
                        break;
                    cout << b << '\n';
                    ++line_count;
                } else {
                    const int number_of_chars_in_line = strlen(b);
                    number_of_nucleotides_in_read += number_of_chars_in_line;
                    clean_nucleotide_line(b, flag_process_n, char_change_into);
                    cout << b << '\n';
                }
            }
        }
    }
    delete b;
}

void do_to_csv(int argc, char** argv)
{
    static struct option long_options[] = {
        {"noheader", no_argument , 0, 'n'},
        {"tsv", no_argument, 0, 't'},
        {0, 0, 0, 0} // end of long options
    };
    bool flag_no_header = false;
    bool flag_output_in_tsv = false;

    while(true) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "", long_options, &option_index);
		if(c == -1) break;
		switch(c) {
		case 0:
			// you can see long_options[option_index].name/flag and optarg (null if no argument).
			break;
		case 'n':
            flag_no_header = true;
			break;
        case 't':
            flag_output_in_tsv = true;
            break;
		}
	}
    for(int i = optind + 1; i < argc; ++i) {
        to_csv(argv[i], flag_no_header, flag_output_in_tsv);
    }
}

void do_fold(int argc, char** argv)
{
    static struct option long_options[] = {
        {"len", required_argument, 0, 'l'},
        {0, 0, 0, 0} // end of long options
    };
    int length_of_line = 70;

    while(true) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "", long_options, &option_index);
		if(c == -1) break;
		switch(c) {
		case 0:
			// you can see long_options[option_index].name/flag and optarg (null if no argument).
			break;
		case 'l':
            length_of_line = atoi(optarg);
			break;
		}
	}
    for(int i = optind + 1; i < argc; ++i) {
        fold_fastx(argv[i], length_of_line, true);
    }
}

void do_unfold(int argc, char** argv)
{
    for(int i = 2; i < argc; ++i) {
        fold_fastx(argv[i], 0, false);
    }
}

void do_tofasta(int argc, char** argv)
{
    for(int i = 2; i < argc; ++i) {
        fastq_to_fasta(argv[i]);
    }
}

void do_clean(int argc, char** argv)
{
    bool flag_change_into_a = false;
    bool flag_change_into_c = false;
    bool flag_change_into_g = false;
    bool flag_change_into_t = false;
    bool flag_change_into_n = false;
    bool flag_process_n = false;
    bool flag_random_change = false;
    static struct option long_options[] = {
        {"processn", no_argument, 0, 'p'},
        {"a", no_argument, 0, 'a'},
        {"c", no_argument, 0, 'c'},
        {"g", no_argument, 0, 'g'},
        {"t", no_argument, 0, 't'},
        {"n", no_argument, 0, 'n'},
        {"random", no_argument, 0, 'r'},
        {0, 0, 0, 0} // end of long options
    };

    while(true) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "", long_options, &option_index);
		if(c == -1) break;
		switch(c) {
		case 0:
			// you can see long_options[option_index].name/flag and optarg (null if no argument).
			break;
		case 'a':
            flag_change_into_a = true;
			break;
		case 'c':
            flag_change_into_c = true;
			break;
		case 'g':
            flag_change_into_g = true;
			break;
		case 't':
            flag_change_into_t = true;
			break;
		case 'n':
            flag_change_into_n = true;
			break;
        case 'p':
            flag_process_n = true;
            break;
        case 'r':
            flag_random_change = true;
            break;
		}
	}
    const int num_change_flags = int(flag_change_into_a) + int(flag_change_into_c) + int(flag_change_into_g) + int(flag_change_into_t) + int(flag_change_into_n) + int(flag_random_change);
    if(1 < num_change_flags) {
        cerr << "ERROR: --a, --c, --g, --t, --n, --random are exclusive." << endl;
        exit(1);
    }
    int char_change_into = -1; // -1 means a random char
    if(flag_change_into_a) char_change_into = 'A';
    if(flag_change_into_c) char_change_into = 'C';
    if(flag_change_into_g) char_change_into = 'G';
    if(flag_change_into_t) char_change_into = 'T';
    if(flag_change_into_n) char_change_into = 'N';
    for(int i = optind + 1; i < argc; ++i) {
        clean_fastx(argv[i], flag_process_n, char_change_into);
    }
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
        cerr << "Currently, no options available.\n\n";
        cerr << "It counts the number of the sequences in each given file.\n";
        return;
	}
    if(subcmd == "name") {
        cerr << "Usage: fatt name [options...] <FAST(A|Q) files>\n\n";
        cerr << "Currently, no options available.\n\n";
        cerr << "It outputs the name of the sequences in each given file.\n";
        return;
	}
    if(subcmd == "chksamename") {
        cerr << "Usage: fatt chksamename [options...] <FAST(A|Q) files>\n\n";
        cerr << "--read\tOutput only read names. (Omit file names)\n\n";
        cerr << "It outputs the name of the sequences that appear more than\n";
        cerr << "once in the given files. Note that it looks only names.\n";
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
        cerr << "Currently, no options available.\n\n";
        cerr << "It outputs the length of the sequences in given files.\n";
        return;
    }
    if(subcmd == "stat") {
        cerr << "Usage: fatt stat [options...] <FAST(A|Q) files>\n\n";
        cerr << "--html\tOutput in HTML format.\n";
        cerr << "--json\tOutput in JSON format.\n";
        return;
    }
    if(subcmd == "index") {
        cerr << "Usage: fatt index [options...] <FAST(A|Q) files>\n\n";
        cerr << "Currently, no options available.\n\n";
        cerr << "It creates an index on the name of the sequences in each given file.\n";
        cerr << "Subsequent access may get faster if the file is very large and you\n";
        cerr << "retrieve only a few sequences.\n";
        return;
    }
    if(subcmd == "guessqvtype") {
        cerr << "Usage: fatt guessqvtype [options...] <FAST(A|Q) files>\n\n";
        cerr << "Currently, no options available.\n\n";
        return;
    }
    if(subcmd == "tocsv") {
        cerr << "Usage: fatt tocsv [options...] <FAST(A|Q) files>\n\n";
        cerr << "--noheader\tSuppress header output.\n";
        cerr << "--tsv\tUse TSV instead of CSV.\n";
        return;
    }
    if(subcmd == "fold") {
        cerr << "Usage: fatt fold [options...] <FAST(A|Q) files>\n\n";
        cerr << "--len=n\tFold lines at n characters. n is 70 by default.\n";
        return;
    }
    if(subcmd == "unfold") {
        cerr << "Usage: fatt unfold [options...] <FAST(A|Q) files>\n\n";
        cerr << "Currently, no options available.\n\n";
        return;
    }
    if(subcmd == "tofasta") {
        cerr << "Usage: fatt tofasta [options...] <FASTQ> files>\n\n";
        cerr << "Currently, no options available.\n\n";
        return;
    }
    if(subcmd == "clean") {
        cerr << "Usage: fatt clean [options...] <FAST(A|Q) files>\n\n";
        cerr << "--processn\tBy default, [^ACGTNacgtn] will be changed. With --processn, [^ACGTacgt] will be changed\n";
        cerr << "--a\tChange into 'A'\n";
        cerr << "--c\tChange into 'C'\n";
        cerr << "--g\tChange into 'G'\n";
        cerr << "--t\tChange into 'T'\n";
        cerr << "--n\tChange into 'N'\n";
        cerr << "--random\tChange into A/C/G/T randomly\n";
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
	cerr << "\tlen\toutput the lengths of reads\n";
    cerr << "\tstat\tshow the statistics of input sequences\n";
    cerr << "\tindex\tcreate an index on read names\n";
    cerr << "\tclean\tconvert non-ACGT(N) characters to ACGT\n";
    cerr << "\tguessqvtype\tguess the type of FASTQ (Sanger/Illumina1.3/Illumina1.5/...)\n";
    cerr << "\ttocsv\tconvert sequences into CSV format\n";
    cerr << "\tfold\tfold sequences\n";
    cerr << "\tunfold\tunfold sequences\n";
    cerr << "\ttofasta\tconvert a FASTQ file into a FASTA file\n";
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
    if(commandString == "stat") {
        do_stat(argc, argv);
        return;
    }
    if(commandString == "index") {
        do_index(argc, argv);
        return;
    }
    if(commandString == "guessqvtype") {
        do_guess_qv_type(argc, argv);
        return;
    }
    if(commandString == "tocsv") {
        do_to_csv(argc, argv);
        return;
    }
    if(commandString == "fold") {
        do_fold(argc, argv);
        return;
    }
    if(commandString == "unfold") {
        do_unfold(argc, argv);
        return;
    }
    if(commandString == "tofasta") {
        do_tofasta(argc, argv);
        return;
    }
    if(commandString == "clean") {
        do_clean(argc, argv);
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
