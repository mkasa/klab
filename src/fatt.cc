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
#include <cctype>
#include <iomanip>
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

#ifndef VERSION_STRING
#define VERSION_STRING ""
#endif
void show_version()
{
    cout << "fatt version " << VERSION_STRING << endl;
}

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

class CoutBuffering
{
    vector<char> buffer;
public:
    CoutBuffering(bool) {}
    void setBufferSize(size_t size = 256 * 1024u * 1024u) {
        if(isatty(fileno(stdin))) return; // No buffering if TTY.
        buffer.resize(size);
        cout.rdbuf()->pubsetbuf(&*buffer.begin(), buffer.size());
        std::ios_base::sync_with_stdio(false);
    }
    CoutBuffering(size_t size = 256 * 1024u * 1024u) {
        setBufferSize(size);
    }
};

class FileLineBufferWithAutoExpansion
{
    ifstream ist;
    vector<char> bufferForIFStream;
    string fileName;
    static const size_t INITIAL_BUFFER_SIZE = 8 * 1024u;
    static const size_t STREAM_BUFFER_SIZE = 16 * 1024u * 1024u;
    size_t currentBufferSize;
    size_t bufferOffsetToBeFill;
    size_t line_count; ///< 1-origin
    vector<char> headerID;
    bool isFASTAMode;
    bool isFASTQMode;

public:
    char* b;

private:
    void expandBufferDouble() {
        const size_t nextBufferSize = currentBufferSize * 2u;
        char* newb = new char[nextBufferSize];
        memcpy(newb, b, bufferOffsetToBeFill);
        delete[] b;
        b = newb;
        currentBufferSize = nextBufferSize;
        isFASTAMode = false;
        isFASTQMode = false;
    }

public:
    FileLineBufferWithAutoExpansion() {
        b = new char[INITIAL_BUFFER_SIZE];
        currentBufferSize = INITIAL_BUFFER_SIZE;
        line_count = 0; // Just for safety
        headerID.reserve(INITIAL_BUFFER_SIZE);
        bufferForIFStream.resize(STREAM_BUFFER_SIZE);
    }
    ~FileLineBufferWithAutoExpansion() {
        delete[] b;
    }
    bool open(const char* file_name) {
        ist.rdbuf()->pubsetbuf(&*bufferForIFStream.begin(), bufferForIFStream.size());
        ist.open(file_name, ios::binary);
        line_count = 0;
        fileName = file_name;
        return ist;
    }
    void close() {
        ist.close();
    }
    bool getline() {
        bufferOffsetToBeFill = 0u;
        do {
            if(ist.getline(b + bufferOffsetToBeFill, currentBufferSize - bufferOffsetToBeFill)) {
                line_count++;
                const bool isFirstLine = line_count == 1;
                if(isFirstLine) {
                    if(looksLikeFASTQHeader()) {
                        isFASTQMode = true;
                        registerHeaderLine();
                    } else if(looksLikeFASTAHeader()) {
                        isFASTAMode = true;
                    } else {
                        cerr << fileName << " does not look like either of FASTA/FASTQ!\n";
                        exit(1);
                    }
                }
                return true;
            }
            if(ist.eof()) return false;
            bufferOffsetToBeFill += ist.gcount();
            expandBufferDouble();
            ist.clear();
        } while(true);
    }
    bool looksLikeFASTQHeader() const { return b[0] == '@'; }
    bool looksLikeFASTAHeader() const { return b[0] == '>'; }
    bool notFollowedByHeaderOrEOF() {
        if(ist.peek() == '@') return false;
        return !ist.eof();
    }
    void expectHeaderOfEOF() {
        if(notFollowedByHeaderOrEOF()) {
            cerr << "WARNING: Bad file format at line " << getLineCount() << ".\n";
            cerr << "         We expected a new sequence (starting with '@') or EOF,\n";
            cerr << "         but there is not. Check the file (" << fileName << ") first.\n";
        }
    }
    size_t getLineCount() const { return line_count; }
    bool fail() { return ist.fail(); }
    size_t tellg() { return ist.tellg(); }
    size_t len() { return strlen(b); }
    void seekg(size_t offset) { ist.seekg(offset); }
    void registerHeaderLine() {
        const size_t len_with_gt = len();
        if(len_with_gt == 0u) return;
        const size_t l = len_with_gt - 1u;
        headerID.resize(l);
        for(size_t i = 1; i < len_with_gt; i++) {
            if(b[i] != ' ') {
                headerID[i - 1u] = b[i];
            } else {
                headerID.resize(i - 1);
                break;
            }
        }
    }
    bool looksLikeFASTQSeparator() const {
        if(b[0] != '+') return false;
        if(b[1] == '\0') return true;
        for(size_t i = 0; i < headerID.size(); ++i) {
            const char next_b = b[i + 1u];
            if(next_b == '\0' || next_b == ' ' || next_b != headerID[i]) return false;
        }
        const char stopChar = b[headerID.size() + 1u];
        return stopChar == '\0' || stopChar == ' ';
    }
};

void calculate_n50_statistics(const char* fname,
                              vector<size_t>& length_of_scaffolds_wgap,
                              vector<size_t>& length_of_scaffolds_wogap,
                              vector<size_t>& length_of_contigs)
{
    FileLineBufferWithAutoExpansion f;
    if(!f.open(fname)) {
        cerr << "Cannot open '" << fname << "'" << endl;
        return;
    }
    if(f.getline()) {
        size_t length_as_scaffolds_wgap = 0;
        size_t length_as_scaffolds_wogap = 0;
        size_t length_as_contigs = 0;
        if(!f.looksLikeFASTQHeader()) { 
			// This should be FASTA
            while(f.getline()) {
                if(f.looksLikeFASTAHeader()) {
                    length_of_scaffolds_wgap.push_back(length_as_scaffolds_wgap);
                    length_of_scaffolds_wogap.push_back(length_as_scaffolds_wogap);
                    length_as_scaffolds_wgap = 0;
                    length_as_scaffolds_wogap = 0;
                    if(0 < length_as_contigs) {
                        length_of_contigs.push_back(length_as_contigs);
                        length_as_contigs = 0;
                    }
                } else {
                    length_as_scaffolds_wgap += f.len();
                    length_as_scaffolds_wogap += strlen_without_n(f.b);
                    for(const char* p = f.b; *p != '\0'; p++) {
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
            while(f.getline()) {
                if(f.looksLikeFASTQSeparator()) { // EOS
                    long long n = length_as_scaffolds_wgap;
                    while(f.getline()) {
                        const size_t number_of_qvchars_in_line = f.len();
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    f.expectHeaderOfEOF();
                    if(!f.getline()) break;
                    f.registerHeaderLine();
                    length_of_scaffolds_wgap.push_back(length_as_scaffolds_wgap);
                    length_of_scaffolds_wogap.push_back(length_as_scaffolds_wogap);
                    length_as_scaffolds_wgap = 0;
                    length_as_scaffolds_wogap = 0;
                    if(0 < length_as_contigs) {
                        length_of_contigs.push_back(length_as_contigs);
                        length_as_contigs = 0;
                    }
                } else {
                    const size_t number_of_nucleotides_in_line = f.len();
                    length_as_scaffolds_wgap += number_of_nucleotides_in_line;
                    length_as_scaffolds_wogap += strlen_without_n(f.b);
                    for(const char* p = f.b; *p != '\0'; p++) {
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
}

void show_read_names_in_file(const char* fname, bool show_name) // if show_name is false, output read length
{
    FileLineBufferWithAutoExpansion f;
	if(!f.open(fname)) {
		cerr << "Cannot open '" << fname << "'" << endl;
		return;
	}
    if(f.getline()) {
		if(show_name) output_read_name(f.b);
        size_t number_of_nucleotides_in_read = 0;
        if(!f.looksLikeFASTQHeader()) { 
            while(f.getline()) {
                if(f.looksLikeFASTAHeader()) {
                    if(show_name)
                        output_read_name(f.b);
                    else
                        cout << number_of_nucleotides_in_read << "\n";
                    number_of_nucleotides_in_read = 0;
                } else {
                    number_of_nucleotides_in_read += f.len();
                }
			}
            if(!show_name && 0 < number_of_nucleotides_in_read)
                cout << number_of_nucleotides_in_read << "\n";
		} else {
            while(f.getline()) {
                if(f.looksLikeFASTQSeparator()) {
                    long long n = number_of_nucleotides_in_read;
                    while(f.getline()) {
                        const size_t number_of_qvchars_in_line = f.len();
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    f.expectHeaderOfEOF();
                    if(!f.getline()) break;
                    f.registerHeaderLine();
					if(show_name)
                        output_read_name(f.b);
                    else
                        cout << number_of_nucleotides_in_read << "\n";
                    number_of_nucleotides_in_read = 0;
                } else {
                    const size_t number_of_nucleotides_in_line = f.len();
                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                }
            }
            if(!show_name && 0 < number_of_nucleotides_in_read)
                cout << number_of_nucleotides_in_read << "\n";
		}
	}
}

void count_number_of_reads_in_file(const char* fname)
{
    FileLineBufferWithAutoExpansion f;
    if(!f.open(fname)) {
        cerr << "Cannot open '" << fname << "'" << endl;
        return;
    }
	cout << fname << flush;
    size_t number_of_sequences = 0;
    size_t number_of_nucleotides = 0;
	size_t min_read_len = -1;
	size_t max_read_len = 0;
	#define UPDATE_MIN_AND_MAX(len) min_read_len = std::min<size_t>(min_read_len, len); max_read_len = std::max<size_t>(max_read_len, len);
    if(f.getline()) {
        number_of_sequences++;
        size_t number_of_nucleotides_in_read = 0;
        if(!f.looksLikeFASTQHeader()) { 
            while(f.getline()) {
                if(f.looksLikeFASTAHeader()) {
                    number_of_sequences++;
					UPDATE_MIN_AND_MAX(number_of_nucleotides_in_read);
					number_of_nucleotides_in_read = 0;
                } else {
                    number_of_nucleotides += f.len();
                }
            }
        } else {
            while(f.getline()) {
                if(f.looksLikeFASTQSeparator()) {
                    long long n = number_of_nucleotides_in_read;
                    while(f.getline()) {
                        const size_t number_of_qvchars_in_line = f.len();
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    f.expectHeaderOfEOF();
					UPDATE_MIN_AND_MAX(number_of_nucleotides_in_read);
                    number_of_nucleotides_in_read = 0;
                    if(!f.getline()) break;
                    f.registerHeaderLine();
					++number_of_sequences;
                } else {
                    const size_t number_of_nucleotides_in_line = f.len();
                    number_of_nucleotides += number_of_nucleotides_in_line;
                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                }
            }
        }
    }
	cout << '\t' << number_of_sequences << '\t' << number_of_nucleotides << '\t' << (double(number_of_nucleotides) / number_of_sequences);
	cout << '\t' << min_read_len << '\t' << max_read_len << '\n';
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
    FileLineBufferWithAutoExpansion f;
    map<string, int> readName2fileIndex;
    for(int findex = optind + 1; findex < argc; ++findex) {
        const char* file_name = argv[findex];
        if(!f.open(file_name)) {
            cerr << "Cannot open '" << file_name << "'" << endl;
            continue;
        }
        size_t number_of_sequences = 0;
        size_t number_of_nucleotides = 0;
        if(f.getline()) {
            number_of_sequences++;
            size_t number_of_nucleotides_in_read = 0;
			add_read_name_and_show_error_if_duplicates(readName2fileIndex, argv, f.b, findex, flag_do_not_show_file_name);
            if(!f.looksLikeFASTQHeader()) { 
                while(f.getline()) {
                    if(f.looksLikeFASTAHeader()) {
                        number_of_sequences++;
						add_read_name_and_show_error_if_duplicates(readName2fileIndex, argv, f.b, findex, flag_do_not_show_file_name);
                        number_of_nucleotides_in_read = 0;
                    } else {
                        number_of_nucleotides += f.len();
                    }
                }
            } else {
                while(f.getline()) {
                    if(f.looksLikeFASTQSeparator()) {
                        long long n = number_of_nucleotides_in_read;
                        while(f.getline()) {
                            const size_t number_of_qvchars_in_line = f.len();
                            n -= number_of_qvchars_in_line;
                            if(n <= 0) break;
                        }
                        f.expectHeaderOfEOF();
                        number_of_nucleotides_in_read = 0;
                        if(!f.getline()) break;
                        f.registerHeaderLine();
                        ++number_of_sequences;
						add_read_name_and_show_error_if_duplicates(readName2fileIndex, argv, f.b, findex, flag_do_not_show_file_name);
                    } else {
                        const size_t number_of_nucleotides_in_line = f.len();
                        number_of_nucleotides += number_of_nucleotides_in_line;
                        number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                    }
                }
            }
        }
    }
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
        FileLineBufferWithAutoExpansion f;
        if(!f.open(fname)) {
            cerr << "Cannot open '" << fname << "'" << endl;
            return;
        }
		long long sequence_count = 0;
        size_t last_pos = f.tellg();
        if(f.getline()) {
            #define INSERT_NAME_INTO_TABLE() { stmt.Bind(1, get_read_name_from_header(f.b)); stmt.Bind(2, static_cast<long long>(last_pos)); stmt.Bind(3, sequence_count); stmt.Next(); }
            INSERT_NAME_INTO_TABLE();
			++sequence_count;
            size_t number_of_nucleotides_in_read = 0;
            if(!f.looksLikeFASTQHeader()) { 
                last_pos = f.tellg();
                while(f.getline()) {
                    if(f.looksLikeFASTAHeader()) {
                        INSERT_NAME_INTO_TABLE();
						++sequence_count;
                        number_of_nucleotides_in_read = 0;
                    } else {
                        number_of_nucleotides_in_read += f.len();
                    }
                }
            } else {
                last_pos = f.tellg();
                while(f.getline()) {
                    if(f.looksLikeFASTQSeparator()) {
                        long long n = number_of_nucleotides_in_read;
                        while(f.getline()) {
                            const size_t number_of_qvchars_in_line = f.len();
                            n -= number_of_qvchars_in_line;
                            if(n <= 0) break;
                        }
                        f.expectHeaderOfEOF();
                        last_pos = f.tellg();
                        if(!f.getline()) break;
                        f.registerHeaderLine();
                        INSERT_NAME_INTO_TABLE();
						++sequence_count;
                        number_of_nucleotides_in_read = 0;
                    } else {
                        const size_t number_of_nucleotides_in_line = f.len();
                        number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                    }
                }
            }
            #undef INSERT_NAME_INTO_TABLE
        }
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
    size_t n50_sequence_index = 0;
    size_t n50_length = 0;
    size_t n70_sequence_index = 0;
    size_t n70_length = 0;
    size_t n80_sequence_index = 0;
    size_t n80_length = 0;
    size_t n90_sequence_index = 0;
    size_t n90_length = 0;
    size_t max_length = 0;
    size_t min_length = 0;
    size_t avg_length = 0;
    if(total_length > 0){
        size_t sequence_index;
        size_t sum = 0;
        const size_t n50_total_length = (size_t)((total_length + 1ull) * 0.5);
        for(sequence_index = 0; sum < n50_total_length && sequence_index < lengths.size(); sequence_index++) sum += lengths[sequence_index];
        n50_length = lengths[sequence_index - 1];
        n50_sequence_index = sequence_index;
        const size_t n70_total_length = (size_t)((total_length + 1ull) * 0.7); 
        for(; sum < n70_total_length && sequence_index < lengths.size(); sequence_index++) sum += lengths[sequence_index];
        n70_length = lengths[sequence_index - 1];
        n70_sequence_index = sequence_index;
        const size_t n80_total_length = (size_t)((total_length + 1ull) * 0.8); 
        for(; sum < n80_total_length && sequence_index < lengths.size(); sequence_index++) sum += lengths[sequence_index];
        n80_length = lengths[sequence_index - 1];
        n80_sequence_index = sequence_index;
        const size_t n90_total_length = (size_t)((total_length + 1ull) * 0.9);
        for(; sum < n90_total_length && sequence_index < lengths.size(); sequence_index++) sum += lengths[sequence_index];
        n90_length = lengths[sequence_index - 1];
        n90_sequence_index = sequence_index;

        min_length = lengths.back();
        max_length = lengths.front();
        avg_length = (total_length + lengths.size() / 2) / lengths.size();
    }
    if(flag_html) {
        cout << "<tr><td></td><td>size (bp)</td><td>number</td></tr>\n";
        cout << "<tr><td>max</td><td>" << max_length << "</td><td>1</td></tr>\n";
        cout << "<tr><td>N50</td><td>" << n50_length << "</td><td>" << (n50_sequence_index) << "</td></tr>\n";
        cout << "<tr><td>N70</td><td>" << n70_length << "</td><td>" << (n70_sequence_index) << "</td></tr>\n";
        cout << "<tr><td>N80</td><td>" << n80_length << "</td><td>" << (n80_sequence_index) << "</td></tr>\n";
        cout << "<tr><td>N90</td><td>" << n90_length << "</td><td>" << (n90_sequence_index) << "</td></tr>\n";
        cout << "<tr><td>min</td><td>" << min_length << "</td><td>" << lengths.size() << "</td></tr>\n";
        cout << "<tr><td>avg</td><td>" << avg_length << "</td><td></td></tr>\n";
    } else if(flag_json) {
        cout << "{\"total_length\": " << total_length;
        cout << ",\"count\": " << lengths.size();
        cout << ",\"max_length\": " << max_length;
        cout << ",\"n50\": " << n50_length << ",\"n50num\": " << (n50_sequence_index);
        cout << ",\"n70\": " << n70_length << ",\"n70num\": " << (n70_sequence_index);
        cout << ",\"n80\": " << n80_length << ",\"n80num\": " << (n80_sequence_index);
        cout << ",\"n90\": " << n90_length << ",\"n90num\": " << (n90_sequence_index);
        cout << ",\"avg\": " << avg_length;
        cout << ",\"min\": " << min_length;
        cout << "}";
    } else {
        cout << "Total # of bases = " << total_length << "\n";
        const string scaffold_str = is_contig ? "contig" : "scaffold";
        cout << "# of " << scaffold_str << "s = " << lengths.size() << "\n";
        cout << "Max size = " << max_length << " (# = 1)\n";
        cout << "N50 " << scaffold_str << " size = " << n50_length << " (# = " << (n50_sequence_index) << ")\n";
        cout << "N70 " << scaffold_str << " size = " << n70_length << " (# = " << (n70_sequence_index) << ")\n";
        cout << "N80 " << scaffold_str << " size = " << n80_length << " (# = " << (n80_sequence_index) << ")\n";
        cout << "N90 " << scaffold_str << " size = " << n90_length << " (# = " << (n90_sequence_index) << ")\n";
        cout << "Min size = " << min_length << "\n";
        cout << "Total " << scaffold_str << " # = " << lengths.size() << "\n";
        cout << "Avg size = " << avg_length << "\n";
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
    for(int findex = optind + 1; findex < argc; ++findex) {
        const char* file_name = argv[findex];
		if(flag_index && !doesIndexExist(file_name)) {
			create_index(file_name);
		}
		const bool use_index = (flag_index || (!flag_noindex && doesIndexExist(file_name))) && !flag_reverse_condition;
        FileLineBufferWithAutoExpansion f;
        if(!f.open(file_name)) {
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
                            f.seekg(pos);
                            if(f.fail() || !f.getline()) {
                                cerr << "WARNING: " << read_name << " is missing in the file. Maybe the index is old?\n";
                                continue;
                            }
                            cout << f.b << "\n";
                            if(is_fastq) {
                                f.registerHeaderLine();
                                size_t number_of_nucleotides_in_read = 0;
                                while(f.getline()) {
                                    if(f.looksLikeFASTQSeparator()) {
                                        long long n = number_of_nucleotides_in_read;
                                        cout << f.b << endl;
                                        while(f.getline()) {
                                            const size_t number_of_qvchars_in_line = f.len();
                                            n -= number_of_qvchars_in_line;
                                            cout << f.b << endl;
                                            if(n <= 0) break;
                                        }
                                        break;
                                    } else {
                                        const size_t number_of_nucleotides_in_line = f.len();
                                        number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                                        cout << f.b << endl;
                                    }
                                }
                            } else {
                                while(f.getline()) {
                                    if(f.looksLikeFASTAHeader()) break;
                                    cout << f.b << "\n";
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
                        f.seekg(pos);
                        if(f.fail() || !f.getline()) {
                            cerr << "WARNING: Cannot seek to that far. Maybe the index is old?\n";
                            continue;
                        }
                        cout << f.b << "\n";
                        if(is_fastq) {
                            f.registerHeaderLine();
                            size_t number_of_nucleotides_in_read = 0;
                            while(f.getline()) {
                                if(f.looksLikeFASTQSeparator()) {
                                    long long n = number_of_nucleotides_in_read;
                                    cout << f.b << endl;
                                    while(f.getline()) {
                                        const size_t number_of_qvchars_in_line = f.len();
                                        n -= number_of_qvchars_in_line;
                                        cout << f.b << endl;
                                        if(n <= 0) break;
                                    }
                                    sequence_index++;
                                    if(param_end <= sequence_index) break;
                                    const long long line_start_pos = f.tellg();
                                    if(!f.getline()) {
                                        cerr << "WARNING: reached the end of file.\n";
                                        return;
                                    }
                                    if(!f.looksLikeFASTQHeader()) {
                                        cerr << "ERROR: bad file format. The line does not start with '@' at line " << f.getLineCount() << " (pos " << line_start_pos << ")" << endl;
                                        return;
                                    }
                                    cout << f.b << "\n";
									number_of_nucleotides_in_read = 0;
                                } else {
                                    const size_t number_of_nucleotides_in_line = f.len();
                                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                                    cout << f.b << endl;
                                }
                            }
                        } else {
                            while(f.getline()) {
                                if(f.looksLikeFASTAHeader()) {
                                    sequence_index++;
                                    if(param_end <= sequence_index) break;
                                }
                                cout << f.b << "\n";
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
            bool current_read_has_been_taken = false;
            if(f.getline()) {
                number_of_sequences++;
                size_t number_of_nucleotides_in_read = 0;
                if(param_start == -1) {
                    current_read_has_been_taken = (readNamesToTake.count(get_read_name_from_header(f.b)) != 0) ^ flag_reverse_condition;
                } else {
                    current_read_has_been_taken = param_start + 1 <= number_of_sequences && number_of_sequences <= param_end;
                    // NOTE: the latter condition never hold, if I properly implemented.
                }
                if(current_read_has_been_taken) cout << f.b << endl;
                if(flag_output_unique) readNamesToTake.insert(get_read_name_from_header(f.b));
                if(!f.looksLikeFASTQHeader()) { 
                    while(f.getline()) {
                        if(f.looksLikeFASTAHeader()) {
                            number_of_sequences++;
                            if(param_start == -1) {
                                current_read_has_been_taken = (readNamesToTake.count(get_read_name_from_header(f.b)) != 0) ^ flag_reverse_condition;
                            } else {
                                current_read_has_been_taken = param_start + 1 <= number_of_sequences && number_of_sequences <= param_end;
                                // NOTE: the latter condition never hold, if I properly implemented.
                            }
                            if(current_read_has_been_taken) cout << f.b << endl;
                            if(flag_output_unique) readNamesToTake.insert(get_read_name_from_header(f.b));
                            number_of_nucleotides_in_read = 0;
                        } else {
                            if(current_read_has_been_taken) cout << f.b << endl;
                            number_of_nucleotides += f.len();
                        }
                    }
                } else {
                    while(f.getline()) {
                        if(f.looksLikeFASTQSeparator()) {
                            long long n = number_of_nucleotides_in_read;
                            if(current_read_has_been_taken) cout << f.b << endl;
                            while(f.getline()) {
                                const size_t number_of_qvchars_in_line = f.len();
                                n -= number_of_qvchars_in_line;
                                if(current_read_has_been_taken) cout << f.b << endl;
                                if(n <= 0) break;
                            }
                            f.expectHeaderOfEOF();
                            number_of_nucleotides_in_read = 0;
                            if(!f.getline()) break;
                            f.registerHeaderLine();
                            ++number_of_sequences;
                            if(param_start == -1) {
                                current_read_has_been_taken = (readNamesToTake.count(get_read_name_from_header(f.b)) != 0) ^ flag_reverse_condition;
                            } else {
                                current_read_has_been_taken = param_start + 1 <= number_of_sequences && number_of_sequences <= param_end;
                                // NOTE: the latter condition never hold, if I properly implemented.
                            }
                            if(param_end < number_of_sequences) // NOTE: param_end is 0-origin, number_of_sequences is 1-origin.
                                break;
                            if(current_read_has_been_taken) cout << f.b << endl;
                            if(flag_output_unique) readNamesToTake.insert(get_read_name_from_header(f.b));
                        } else {
                            const size_t number_of_nucleotides_in_line = f.len();
                            number_of_nucleotides += number_of_nucleotides_in_line;
                            number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                            if(current_read_has_been_taken) cout << f.b << endl;
                        }
                    }
                }
            }
        }
    }
}

void do_convert_qv_type(int argc, char** argv)
{
    int param_from_base = 64;
    int param_to_base = 33;
    int param_out_qv_min = 0;
    int param_out_qv_max = 40;
    int param_in_qv_min = -5;
    int param_in_qv_max = 99;

    static struct option long_options[] = {
        {"fromsanger",   no_argument, 0, 's'},
        {"fromsolexa",   no_argument, 0, 'x'},
        {"fromillumina", no_argument, 0, 'i'},
        {"toillumina13", no_argument, 0, '3'},
        {"toillumina15", no_argument, 0, '5'},
        {"toillumina18", no_argument, 0, '8'},
        {"tosanger",     no_argument, 0, 'q'},
    	{"frombase", required_argument, 0, 'f'},
    	{"tobase", required_argument, 0, 't'},
    	{"min", required_argument, 0, 'a'},
    	{"max", required_argument, 0, 'b'},
    	{"inmin", required_argument, 0, 'c'},
    	{"inmax", required_argument, 0, 'd'},
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
            param_out_qv_min = atoi(optarg);
			break;
        case 'b':
            param_out_qv_max = atoi(optarg);
            break;
		case 'c':
            param_in_qv_min = atoi(optarg);
			break;
        case 'd':
            param_in_qv_max = atoi(optarg);
            break;
        case 'f':
            param_from_base = atoi(optarg);
            break;
        case 't':
            param_to_base = atoi(optarg);
            break;
        case 's':
            param_from_base = 33;
            break;
        case 'i':
            param_from_base = 64;
            break;
        case 'x':
            param_from_base = 64;
            param_in_qv_min = -5;
            param_in_qv_max = 40;
            break;
        case '3':
            param_to_base = 64;
            param_out_qv_min = 0;
            param_out_qv_max = 40;
            break;
        case '5':
            param_to_base = 64;
            param_out_qv_min = 3;
            param_out_qv_max = 40;
            break;
        case '8':
            param_to_base = 33;
            param_out_qv_min = 0;
            param_out_qv_max = 40;
            break;
        case 'q':
            param_to_base = 33;
            param_out_qv_min = 0;
            param_out_qv_max = 40;
            break;
		}
	}
	if(param_in_qv_max < param_in_qv_min) {
		cerr << "ERROR: Input QV range [" << param_in_qv_min << ", " << param_in_qv_max << ") is empty." << endl;
		return;
	}
	if(param_out_qv_max < param_out_qv_min) {
		cerr << "ERROR: Output QV range [" << param_out_qv_min << ", " << param_out_qv_max << ") is empty." << endl;
		return;
	}
    cerr << "QV [" << param_in_qv_min << ", " << param_in_qv_max << "):base" << param_from_base << " ==> "
            "QV [" << param_out_qv_min << ", " << param_out_qv_max << "):base" << param_to_base << endl;
    for(int findex = optind + 1; findex < argc; ++findex) {
        const char* file_name = argv[findex];
        FileLineBufferWithAutoExpansion f;
        if(!f.open(file_name)) {
            cerr << "Cannot open '" << file_name << "'" << endl;
            continue;
        }
        if(f.getline()) {
            size_t number_of_nucleotides_in_read = 0;
            if(!f.looksLikeFASTQHeader()) { 
                cerr << "ERROR: the input file '" << file_name << "' does not seem to be a FASTQ file at line " << f.getLineCount() << endl;
                return;
            }
            cout << f.b << "\n";
            while(f.getline()) {
                if(f.looksLikeFASTQSeparator()) {
                    cout << f.b << "\n";
                    long long n = number_of_nucleotides_in_read;
                    while(f.getline()) {
                        const size_t number_of_qvchars_in_line = f.len();
                        for(unsigned char* p = reinterpret_cast<unsigned char*>(f.b); *p; ++p) {
                            int qv = *p - param_from_base;
                            if(qv < param_in_qv_min || param_in_qv_max < qv) {
                                cerr << "ERROR: the input file '" << file_name << "' contains an invalid QV (" << qv << "; chr = '" << *p << "'; ord = '" << int(*p) << "') at line " << f.getLineCount() << endl;
                                return;
                            }
                            if(param_out_qv_max < qv) qv = param_out_qv_max;
                            if(qv < param_out_qv_min) qv = param_out_qv_min;
                            *p = param_to_base + qv;
                        }
                        cout << f.b << '\n';
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    f.expectHeaderOfEOF();
                    number_of_nucleotides_in_read = 0;
                    if(!f.getline()) break;
                    f.registerHeaderLine();
                    cout << f.b << "\n";
                } else {
                    const size_t number_of_nucleotides_in_line = f.len();
                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                    cout << f.b << "\n";
                }
            }
        }
    }
}

void do_guess_qv_type(int argc, char** argv)
{
    for(int findex = 2; findex < argc; ++findex) {
        const char* file_name = argv[findex];
        FileLineBufferWithAutoExpansion f;
        if(!f.open(file_name)) {
            cerr << "Cannot open '" << file_name << "'" << endl;
            continue;
        }
        size_t histogram[256];
        for(int i = 0; i < 256; ++i) histogram[i] = 0;
        size_t number_of_sequences = 0;
        size_t number_of_nucleotides = 0;
        if(f.getline()) {
            number_of_sequences++;
            size_t number_of_nucleotides_in_read = 0;
            if(!f.looksLikeFASTQHeader()) { 
                cerr << "ERROR: the input file '" << file_name << "' does not seem to be a FASTQ file at line " << f.getLineCount() << endl;
                return;
            }
            while(f.getline()) {
                if(f.looksLikeFASTQSeparator()) {
                    long long n = number_of_nucleotides_in_read;
                    while(f.getline()) {
                        const size_t number_of_qvchars_in_line = f.len();
                        for(unsigned char* p = reinterpret_cast<unsigned char*>(f.b); *p; ++p) histogram[*p]++;
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    f.expectHeaderOfEOF();
                    number_of_nucleotides_in_read = 0;
                    if(!f.getline()) break;
                    f.registerHeaderLine();
                    ++number_of_sequences;
                } else {
                    const size_t number_of_nucleotides_in_line = f.len();
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
}

void to_csv(const char* file_name, bool does_not_output_header, bool output_in_tsv)
{
    FileLineBufferWithAutoExpansion f;
    if(!f.open(file_name)) {
        cerr << "Cannot open '" << file_name << "'" << endl;
        return;
    }
    size_t number_of_sequences = 0;
    size_t number_of_nucleotides = 0;
    if(f.getline()) {
        number_of_sequences++;
        size_t number_of_nucleotides_in_read = 0;
        if(!f.looksLikeFASTQHeader()) { 
            // This should be FASTA
            if(!does_not_output_header) {
                if(output_in_tsv) {
                    cout << "id\tdesc\tseq\n";
                } else {
                    cout << "id,desc,seq\n";
                }
            }
            #define OUTPUT_HEADER() {                                           \
                const string& idstr = get_read_name_from_header(f.b);\
                cout << idstr << (output_in_tsv ? "\t" : ",\"");\
                const char* descp = f.b + 1 + idstr.size();\
                if(*descp != '\0') descp++;\
                if(output_in_tsv) cout << descp; else cout << CSVEscape(descp);\
                cout << (output_in_tsv ? "\t" : "\",");\
            }
            OUTPUT_HEADER();
            while(f.getline()) {
                if(f.looksLikeFASTAHeader()) {
                    number_of_sequences++;
                    number_of_nucleotides_in_read = 0;
                    cout << '\n';
                    OUTPUT_HEADER();
                } else {
                    cout << f.b;
                    number_of_nucleotides += f.len();
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
            while(f.getline()) {
                if(f.looksLikeFASTQSeparator()) {
                    cout << (output_in_tsv ? '\t' : ',') << '"';
                    long long n = number_of_nucleotides_in_read;
                    while(f.getline()) {
                        const size_t number_of_qvchars_in_line = f.len();
						if(output_in_tsv)
                        	cout << f.b;
						else
							cout << CSVEscape(f.b);
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    f.expectHeaderOfEOF();
                    number_of_nucleotides_in_read = 0;
                    if(!f.getline()) break;
                    f.registerHeaderLine();
                    cout << '"' << '\n';
                    OUTPUT_HEADER();
                    ++number_of_sequences;
                } else {
                    const size_t number_of_nucleotides_in_line = f.len();
                    number_of_nucleotides += number_of_nucleotides_in_line;
                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                    cout << f.b;
                }
            }
            cout << '"' << '\n';
            #undef OUTPUT_HEADER
        }
    }
}

void fold_fastx(const char* file_name, int length_of_line, bool is_folding)
{
    FileLineBufferWithAutoExpansion f;
    if(!f.open(file_name)) {
        cerr << "Cannot open '" << file_name << "'" << endl;
        return;
    }
    if(f.getline()) {
        size_t number_of_nucleotides_in_output_line = 0;
        cout << f.b << '\n';
        if(!f.looksLikeFASTQHeader()) { 
            // This should be FASTA
            while(f.getline()) {
                if(f.looksLikeFASTAHeader()) {
                    if(0 < number_of_nucleotides_in_output_line) {
                        cout << '\n';
                        number_of_nucleotides_in_output_line = 0;
                    }
                    cout << f.b << '\n';
                } else {
                    const int number_of_chars_in_line = f.len();
                    if(is_folding) {
                        int off = 0;
                        while(off < number_of_chars_in_line) {
                            int s = number_of_chars_in_line - off;
                            if(length_of_line - number_of_nucleotides_in_output_line <= s) s = length_of_line - number_of_nucleotides_in_output_line;
                            for(int i = 0; i < s; ++i) cout << f.b[off + i];
                            off += s;
                            number_of_nucleotides_in_output_line += s;
                            if(length_of_line <= number_of_nucleotides_in_output_line) {
                                cout << '\n';
                                number_of_nucleotides_in_output_line = 0;
                            }
                        }
                    } else {
                        cout << f.b;
                        number_of_nucleotides_in_output_line += number_of_chars_in_line;
                    }
                }
            }
            if(0 < number_of_nucleotides_in_output_line)
                cout << '\n';
        } else {
            size_t number_of_nucleotides_in_read = 0;
            while(f.getline()) {
                if(f.looksLikeFASTQSeparator()) {
                    if(0 < number_of_nucleotides_in_output_line) {
                        cout << '\n';
                        number_of_nucleotides_in_output_line = 0;
                    }
					cout << "+\n";
                    long long n = number_of_nucleotides_in_read;
                    while(f.getline()) {
                        const size_t number_of_qvchars_in_line = f.len();
                        if(is_folding) {
                            int off = 0;
                            while(off < number_of_qvchars_in_line) {
                                int s = number_of_qvchars_in_line - off;
                                if(length_of_line - number_of_nucleotides_in_output_line <= s) s = length_of_line - number_of_nucleotides_in_output_line;
                                for(int i = 0; i < s; ++i) cout << f.b[off + i];
                                off += s;
                                number_of_nucleotides_in_output_line += s;
                                if(length_of_line <= number_of_nucleotides_in_output_line) {
                                    cout << '\n';
                                    number_of_nucleotides_in_output_line = 0;
                                }
                            }
                        } else {
                            cout << f.b;
                            number_of_nucleotides_in_output_line += number_of_qvchars_in_line;
                        }
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    f.expectHeaderOfEOF();
                    number_of_nucleotides_in_read = 0;
                    if(!f.getline()) break;
                    f.registerHeaderLine();
                    if(0 < number_of_nucleotides_in_output_line) {
                        cout << '\n';
                        number_of_nucleotides_in_output_line = 0;
                    }
                    cout << f.b << '\n';
                } else {
                    const int number_of_chars_in_line = f.len();
                    number_of_nucleotides_in_read += number_of_chars_in_line;
                    if(is_folding) {
                        int off = 0;
                        while(off < number_of_chars_in_line) {
                            int s = number_of_chars_in_line - off;
                            if(length_of_line - number_of_nucleotides_in_output_line <= s) s = length_of_line - number_of_nucleotides_in_output_line;
                            for(int i = 0; i < s; ++i) cout << f.b[off + i];
                            off += s;
                            number_of_nucleotides_in_output_line += s;
                            if(length_of_line <= number_of_nucleotides_in_output_line) {
                                cout << '\n';
                                number_of_nucleotides_in_output_line = 0;
                            }
                        }
                    } else {
                        cout << f.b;
                        number_of_nucleotides_in_output_line += number_of_chars_in_line;
                    }
                }
            }
            if(0 < number_of_nucleotides_in_output_line) cout << '\n';
        }
    }
}

void fastq_to_fasta(const char* file_name)
{
    FileLineBufferWithAutoExpansion f;
    if(!f.open(file_name)) {
        cerr << "Cannot open '" << file_name << "'" << endl;
        return;
    }
    if(f.getline()) {
        size_t number_of_nucleotides_in_output_line = 0;
        if(!f.looksLikeFASTQHeader()) { 
            if(f.looksLikeFASTAHeader()) {
                cerr << "Input is already FASTA." << endl;
            } else {
                cerr << "Input does not look like FASTA/FASTQ." << endl;
            }
            return;
        }
        // This is FASTQ
        f.b[0] = '>';
        cout << f.b << '\n';
        size_t number_of_nucleotides_in_read = 0;
        while(f.getline()) {
            if(f.looksLikeFASTQSeparator()) {
                long long n = number_of_nucleotides_in_read;
                while(f.getline()) {
                    const size_t number_of_qvchars_in_line = f.len();
                    number_of_nucleotides_in_output_line += number_of_qvchars_in_line;
                    n -= number_of_qvchars_in_line;
                    if(n <= 0) break;
                }
                if(f.notFollowedByHeaderOrEOF()) {
                    cerr << "WARNING: bad file format? at line " << f.getLineCount() << "." << endl;
                    return;
                }
                number_of_nucleotides_in_read = 0;
                if(!f.getline()) break;
                f.registerHeaderLine();
                f.b[0] = '>';
                cout << f.b << '\n';
            } else {
                const int number_of_chars_in_line = f.len();
                number_of_nucleotides_in_read += number_of_chars_in_line;
                cout << f.b << '\n';
            }
        }
    }
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
    FileLineBufferWithAutoExpansion f;
    if(!f.open(file_name)) {
        cerr << "Cannot open '" << file_name << "'" << endl;
        return;
    }
    srand(time(NULL));
    if(f.getline()) {
        cout << f.b << '\n';
        if(!f.looksLikeFASTQHeader()) { 
            // This should be FASTA
            while(f.getline()) {
                if(!f.looksLikeFASTAHeader()) { clean_nucleotide_line(f.b, flag_process_n, char_change_into); }
                cout << f.b << '\n';
            }
        } else {
            size_t number_of_nucleotides_in_read = 0;
            while(f.getline()) {
                if(f.looksLikeFASTQSeparator()) {
					cout << "+\n";
                    long long n = number_of_nucleotides_in_read;
                    while(f.getline()) {
                        const size_t number_of_qvchars_in_line = f.len();
                        cout << f.b << '\n';
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    f.expectHeaderOfEOF();
                    number_of_nucleotides_in_read = 0;
                    if(!f.getline()) break;
                    f.registerHeaderLine();
                    cout << f.b << '\n';
                } else {
                    const int number_of_chars_in_line = f.len();
                    number_of_nucleotides_in_read += number_of_chars_in_line;
                    clean_nucleotide_line(f.b, flag_process_n, char_change_into);
                    cout << f.b << '\n';
                }
            }
        }
    }
}

static char asterisk_if_nulchar(int i)
{
    if(i != 0) return i;
    return '*';
}

void investigate_composition(const char* file_name, bool ignore_case, bool flag_only_monomer, bool flag_only_bimer, bool flag_only_trimer, bool flag_dapi_check)
{
    FileLineBufferWithAutoExpansion f;
    if(!f.open(file_name)) {
        cerr << "Cannot open '" << file_name << "'" << endl;
        return;
    }
    // We will investigate 1-mer to 3-mer composition of the given input file.
    // NOTE: the variables below are static only for avoiding stack overflow.
    const size_t N = 256;
    static size_t freq_1_mer[N];
    static size_t freq_2_mer[N][N];
    static size_t freq_3_mer[N][N][N];
    memset(freq_1_mer, 0, sizeof(freq_1_mer));
    memset(freq_2_mer, 0, sizeof(freq_2_mer));
    memset(freq_3_mer, 0, sizeof(freq_3_mer));
    if(f.getline()) {
        const size_t LOOK_BEHIND_SIZE = 2;
        char previousCharacters[LOOK_BEHIND_SIZE];
        memset(previousCharacters, 0, sizeof(previousCharacters));
        if(f.looksLikeFASTQHeader()) {
            size_t number_of_nucleotides_in_read = 0;
            while(f.getline()) {
                if(f.looksLikeFASTQSeparator()) {
                    long long n = number_of_nucleotides_in_read;
                    while(f.getline()) {
                        const size_t number_of_qvchars_in_line = f.len();
                        n -= number_of_qvchars_in_line;
                        if(n <= 0) break;
                    }
                    {
                        const char c = 0;
                        freq_1_mer[c]++;
                        freq_2_mer[previousCharacters[0]][c]++;
                        freq_3_mer[previousCharacters[1]][previousCharacters[0]][c]++;
                        for(size_t j = 0; j < sizeof(previousCharacters) / sizeof(char); j++) previousCharacters[j] = 0;
                    }
                    f.expectHeaderOfEOF();
                    number_of_nucleotides_in_read = 0;
                    if(!f.getline()) break;
                    f.registerHeaderLine();
                } else {
                    const int number_of_chars_in_line = f.len();
                    number_of_nucleotides_in_read += number_of_chars_in_line;
                    for(size_t i = 0; i < number_of_chars_in_line; i++) {
                        const char c = ignore_case ? toupper(f.b[i]) : f.b[i];
                        freq_1_mer[c]++;
                        freq_2_mer[previousCharacters[0]][c]++;
                        freq_3_mer[previousCharacters[1]][previousCharacters[0]][c]++;
                        for(size_t j = sizeof(previousCharacters) / sizeof(char) - 1u; 0 < j; j--) previousCharacters[j] = previousCharacters[j - 1];
                        previousCharacters[0] = c;
                    }
                }
            }
        } else {
            while(f.getline()) {
                if(f.looksLikeFASTAHeader()) {
                    {
                        const char c = 0;
                        freq_1_mer[c]++;
                        freq_2_mer[previousCharacters[0]][c]++;
                        freq_3_mer[previousCharacters[1]][previousCharacters[0]][c]++;
                        for(size_t j = 0; j < sizeof(previousCharacters) / sizeof(char); j++) previousCharacters[j] = 0;
                    }
                } else {
                    const size_t number_of_chars_in_line = f.len();
                    for(size_t i = 0; i < number_of_chars_in_line; i++) {
                        const char c = ignore_case ? toupper(f.b[i]) : f.b[i];
                        freq_1_mer[c]++;
                        freq_2_mer[previousCharacters[0]][c]++;
                        freq_3_mer[previousCharacters[1]][previousCharacters[0]][c]++;
                        for(size_t j = sizeof(previousCharacters) / sizeof(char) - 1u; 0 < j; j--) previousCharacters[j] = previousCharacters[j - 1];
                        previousCharacters[0] = c;
                    }
                }
            }
        }
    }
    // show the results
    const size_t total_n_nmers = accumulate(freq_1_mer, freq_1_mer + sizeof(freq_1_mer) / sizeof(size_t), 0llu);
    cout << "Total # n-mers\n\t" << total_n_nmers << "\n";
    cout.setf(ios_base::fixed, ios_base::floatfield);
    if(!flag_only_trimer && !flag_only_bimer && !flag_dapi_check) {
        cout << "1-mer stats\n";
        for(size_t i = 0; i < N; ++i) {
            if(freq_1_mer[i] != 0) cout << "\t" << asterisk_if_nulchar(i) << "\t" << freq_1_mer[i] << "\t" << (double(freq_1_mer[i]) / total_n_nmers) << "\n";
        }
    }
    if(!flag_only_trimer && !flag_only_monomer && !flag_dapi_check) {
        cout << "2-mer stats\n";
        for(size_t i = 0; i < N; ++i) {
            for(size_t j = 0; j < N; ++j) {
                if(freq_2_mer[i][j] != 0) cout << "\t" << asterisk_if_nulchar(i) << asterisk_if_nulchar(j) << "\t" << freq_2_mer[i][j] << "\t" << (double(freq_2_mer[i][j]) / total_n_nmers) << "\n";
            }
        }
    }
    if(!flag_only_monomer && !flag_only_bimer && !flag_dapi_check) {
        cout << "3-mer stats\n";
        for(size_t i = 0; i < N; ++i) {
            for(size_t j = 0; j < N; ++j) {
                for(size_t k = 0; k < N; ++k) {
                    if(freq_3_mer[i][j][k] != 0) cout << "\t" << asterisk_if_nulchar(i) << asterisk_if_nulchar(j) << asterisk_if_nulchar(k) << "\t" << freq_3_mer[i][j][k] << "\t" << (double(freq_3_mer[i][j][k]) / total_n_nmers) << "\n";
                }
            }
        }
    }
    if(flag_dapi_check) {
        cout << "Stats for DAPI intensity\n";
        const size_t sum_at = freq_1_mer['A'] + freq_1_mer['a'] + freq_1_mer['T'] + freq_1_mer['t'];
        cout << "\tA/T\t" << sum_at << "\t" << (double(sum_at) / total_n_nmers) << "\n";
        for(size_t i = 0; i < N; ++i) {
            if(toupper(i) != 'A' && toupper(i) != 'T') continue;
            for(size_t j = 0; j < N; ++j) {
                if(toupper(j) != 'A' && toupper(j) != 'T') continue;
                for(size_t k = 0; k < N; ++k) {
                    if(toupper(k) != 'A' && toupper(k) != 'T') continue;
                    if(isupper(i) && isupper(j) && isupper(k)) continue;
                    freq_3_mer[toupper(i)][toupper(j)][toupper(k)] += freq_3_mer[i][j][k];
                }
            }
        }
        const size_t freq_AAA = freq_3_mer['A']['A']['A'] + freq_3_mer['T']['T']['T'];
        const size_t freq_AAT = freq_3_mer['A']['A']['T'] + freq_3_mer['A']['T']['T'];
        const size_t freq_ATA = freq_3_mer['A']['T']['A'] + freq_3_mer['T']['A']['T'];
        const size_t freq_TTA = freq_3_mer['T']['T']['A'] + freq_3_mer['T']['A']['A'];
        cout << "\tAAA\t" << freq_AAA << "\t" << (double(freq_AAA) / total_n_nmers) << "\n";
        cout << "\tAAT\t" << freq_AAT << "\t" << (double(freq_AAT) / total_n_nmers) << "\n";
        cout << "\tATA\t" << freq_ATA << "\t" << (double(freq_ATA) / total_n_nmers) << "\n";
        cout << "\tTTA\t" << freq_TTA << "\t" << (double(freq_TTA) / total_n_nmers) << "\n";
    }
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

void do_composition(int argc, char** argv)
{
    bool flag_ignore_case  = false;
    bool flag_only_monomer = false;
    bool flag_only_bimer   = false;
    bool flag_only_trimer  = false;
    bool flag_dapi_check   = false;
    static struct option long_options[] = {
        {"ignorecase", no_argument, 0, 'i'},
        {"monomer", no_argument, 0, '1'},
        {"bimer", no_argument, 0, '2'},
        {"trimer", no_argument, 0, '3'},
        {"dapicheck", no_argument, 0, 'd'},
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
		case 'i':
            flag_ignore_case = true;
			break;
        case '1':
            flag_only_monomer = true;
            break;
        case '2':
            flag_only_bimer = true;
            break;
        case '3':
            flag_only_trimer = true;
            break;
        case 'd':
            flag_dapi_check = true;
            break;
        }
    }
    for(int i = optind + 1; i < argc; ++i) {
        investigate_composition(argv[i], flag_ignore_case, flag_only_monomer, flag_only_bimer, flag_only_trimer, flag_dapi_check);
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
        cerr << "--start\tSpecify the start index of reads to be output. 0-based, inclusive.\n";
        cerr << "--end\tSpecify the end index of reads to be output. 0-based, exclusive.\n";
        cerr << "--num\tSpecify the number of reads to be output.\n";
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
    if(subcmd == "convertqv") {
        cerr << "Usage: fatt convertqv [options...] <FASTQ files>\n\n";
        cerr << "Normal options:\n";
        cerr << "\t--fromsolexa\tInput is Solexa FASTQ\n";
        cerr << "\t--fromsanger\tInput is Sanger FASTQ\n";
        cerr << "\t--fromillumina\tInput is Illumina FASTQ (of any types other than Illumina 1.8+, which is essentially Sanger FASTQ)\n";
        cerr << "\t--toillumina13\tOutput is Illumina 1.3\n";
        cerr << "\t--toillumina15\tOutput is Illumina 1.5\n";
        cerr << "\t--toillumina18\tOutput is Illumina 1.8\n";
        cerr << "\t--tosanger\tOutput is Sanger FASTQ (default)\n";
        cerr << "\nCustom options:\n";
        cerr << "\t--frombase\tSet the base of input (the ord of the character that corresponds to QV = 0)\n";
        cerr << "\t--tobase\tSet the base of output (the ord of the character that corresponds to QV = 0)\n";
        cerr << "\t--min d\tSet the minimum QV for output. When a QV is lower than d, it will be set to d.\n";
        cerr << "\t--max d\tSet the maximum QV for output. When a QV is higher than d, it will be set to d.\n";
        cerr << "\t--inmin d\tSet the minimum QV for input. When a QV in input is lower than d, print an error.\n";
        cerr << "\t--inmax d\tSet the maximum QV for input. When a QV in input is higher than d, print an error.\n";
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
    if(subcmd == "compisition") {
        cerr << "Usage: fatt composition [options...] <FAST(A|Q) files>\n\n";
        cerr << "--ignorecase\tIgnore case ('A' and 'a' will be considered as identical)\n";
        cerr << "--monomer\tShow only monomers\n";
        cerr << "--bimer\tShow only bimers\n";
        cerr << "--trimer\tShow only trimers\n";
        cerr << "--dapicheck\tShow DAPI-staining related stats\n";
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
    cerr << "\tconvertqv\tconvert into a different type of FASTQ (Sanger/Illumina1.3/Illumina1.5/...)\n";
    cerr << "\ttocsv\tconvert sequences into CSV format\n";
    cerr << "\tfold\tfold sequences\n";
    cerr << "\tunfold\tunfold sequences\n";
    cerr << "\ttofasta\tconvert a FASTQ file into a FASTA file\n";
    cerr << "\thelp\tshow help message\n";
    cerr << "\nType 'fatt help <command>' to show the detail of the command.\n";
}

void dispatchByCommand(const string& commandString, int argc, char** argv)
{
    if(commandString == "--version") {
        show_version();
        return;
    }
    CoutBuffering cb;
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
    if(commandString == "convertqv") {
        do_convert_qv_type(argc, argv);
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
    if(commandString == "composition") {
        do_composition(argc, argv);
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
