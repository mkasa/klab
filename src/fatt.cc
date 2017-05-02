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
#include <regex>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
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

inline bool isN(char c)
{
    return c == 'N' || c == 'n';
}

static ostream& operator << (ostream& os, const CSVEscape& c)
{
	for(const char* p = c.p; *p; ++p) {
		if(*p == '"') os << '"';
		os << *p;
	}
	return os;
}

static ostream& operator << (ostream& os, const vector<char>& v)
{
    for(size_t i = 0; i < v.size(); i++) {
        os << v[i];
    }
	return os;
}

static string sep_comma(size_t val)
{
    char buf[32];
    sprintf(buf, "%lu", val);
    string retval;
    size_t offset = (3 - strlen(buf) % 3) % 3;
    for(char* p = buf; *p; p++) {
        if(p != buf && (p - buf + offset) % 3 == 0) retval += ',';
        retval += *p;
    }
    return retval;
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
}

static string get_index_file_name(const char* fastq_file_name)
{
    string index_file_name = fastq_file_name;
    index_file_name += ".index";
    return index_file_name;
}

static bool index_older_than_file(const string file_name, const string index_file_name)
{
    struct stat s, s2;
    {
        int retf = stat(file_name.c_str(), &s);
        if(retf != 0) return false;
    }
    {
        int retf = stat(index_file_name.c_str(), &s2);
        if(retf != 0) return false;
    }
    return s.st_mtime > s2.st_mtime;
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

class FileBuffering
{
    vector<char> buffer;
public:
    FileBuffering(size_t size = 256 * 1024u * 1024u) {
        buffer.resize(size);
    }
    void connectToStream(ostream& os) {
        os.rdbuf()->pubsetbuf(&*buffer.begin(), buffer.size());
        std::ios_base::sync_with_stdio(false);
    }
};

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
    bool is_first_open;
    vector<char> bufferForIFStream;
    string fileName;
    static const size_t INITIAL_BUFFER_SIZE = 8 * 1024u;
    static const size_t STREAM_BUFFER_SIZE = 16 * 1024u * 1024u;
    size_t currentBufferSize;
    size_t bufferOffsetToBeFill;
    size_t line_count; ///< 1-origin
    off_t off_count;
    vector<char> headerID;
    vector<char> headerIDwithDesc;
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
    bool isDirectory(const char* fname) {
        struct stat s;
        const int ret = stat(fname, &s);
        if(ret != 0) return false;
        return S_ISDIR(s.st_mode);
    }

public:
    FileLineBufferWithAutoExpansion() {
        b = new char[INITIAL_BUFFER_SIZE];
        currentBufferSize = INITIAL_BUFFER_SIZE;
        line_count = 0; // Just for safety
        off_count = 0;
        headerID.reserve(INITIAL_BUFFER_SIZE);
        bufferForIFStream.resize(STREAM_BUFFER_SIZE);
        is_first_open = true;
    }
    ~FileLineBufferWithAutoExpansion() {
        if(!is_first_open) close();
        delete[] b;
    }
    bool open(const char* file_name) {
        if(isDirectory(file_name)) { return false; }
        if(is_first_open) {
            is_first_open = false;
            ist.rdbuf()->pubsetbuf(&*bufferForIFStream.begin(), bufferForIFStream.size());
        } else {
            close();
        }
        ist.open(file_name, ios::binary);
        line_count = 0;
        off_count = 0;
        fileName = file_name;
        return true;
    }
    void close() {
        ist.close();
    }
    bool getline() {
        bufferOffsetToBeFill = 0u;
        do {
            if(ist.getline(b + bufferOffsetToBeFill, currentBufferSize - bufferOffsetToBeFill)) {
                line_count++;
                off_count += strlen(b) + 1; // for the delimiter
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
//    size_t tellg() { return ist.tellg(); }
    off_t get_offset() {return off_count; }
    size_t len() { return strlen(b); }
    void seekg(off_t offset) {ist.clear(); ist.seekg(offset); } // clear eofbit before seeking
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
    void registerHeaderLineWithDesc() {
        registerHeaderLine();
        headerIDwithDesc.assign(b, b + len());
    }
    string getSequenceName() const {
        return string(headerID.begin(), headerID.end());
    }
    string getSequenceDescription() const {
        const vector<char>::const_iterator p = find(headerIDwithDesc.begin(), headerIDwithDesc.end(), ' ');
        const vector<char>::const_iterator st = p == headerIDwithDesc.end() ? headerIDwithDesc.end() : p + 1;
        return string(st, headerIDwithDesc.end());
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

// returns true if succeeded.
bool calculate_n50_statistics(const char* fname,
                              vector<size_t>& length_of_scaffolds_wgap,
                              vector<size_t>& length_of_scaffolds_wogap,
                              vector<size_t>& length_of_contigs)
{
    FileLineBufferWithAutoExpansion f;
    if(!f.open(fname)) {
        cerr << "Cannot open '" << fname << "'" << endl;
        return false;
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
    return true;
}

void show_read_names_in_file(const char* fname, bool show_name, bool show_length)
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
                    if(show_length) {
                        if(show_name) cout << "\t";
                        cout << number_of_nucleotides_in_read;
                    }
                    cout << "\n";
                    if(show_name) output_read_name(f.b);
                    number_of_nucleotides_in_read = 0;
                } else {
                    number_of_nucleotides_in_read += f.len();
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
                    if(!f.getline()) break;
                    f.registerHeaderLine();
                    if(show_length) {
                        if(show_name) cout << "\t";
                        cout << number_of_nucleotides_in_read;
                    }
                    cout << "\n";
                    if(show_name) output_read_name(f.b);
                    number_of_nucleotides_in_read = 0;
                } else {
                    const size_t number_of_nucleotides_in_line = f.len();
                    number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                }
            }
		}
        if(show_length) {
            if(show_name) cout << "\t";
            cout << number_of_nucleotides_in_read;
        }
        cout << "\n";
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

void create_index(const char* fname, bool flag_force)
{
    const string index_file_name = get_index_file_name(fname);
    if(doesIndexExist(fname)) {
        cerr << "'" << index_file_name << "' already exists!" << endl;
        if(flag_force) {
            cerr << "However, --force flag is given, so we remove it first." << endl;
            if(unlink(index_file_name.c_str()) != 0) {
                cerr << "Could not delete '" << index_file_name << "'. Abort." << endl;
                return;
            }
        } else {
            return;
        }
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
        off_t last_pos = f.get_offset();
        if(f.getline()) {
            #define INSERT_NAME_INTO_TABLE() { stmt.Bind(1, get_read_name_from_header(f.b)); stmt.Bind(2, static_cast<long long>(last_pos)); stmt.Bind(3, sequence_count); stmt.Next(); }
            INSERT_NAME_INTO_TABLE();
			++sequence_count;
            size_t number_of_nucleotides_in_read = 0;
            if(!f.looksLikeFASTQHeader()) { 
                last_pos = f.get_offset();
                while(f.getline()) {
                    if(f.looksLikeFASTAHeader()) {
                        INSERT_NAME_INTO_TABLE();
						++sequence_count;
                        number_of_nucleotides_in_read = 0;
                    } else {
                        number_of_nucleotides_in_read += f.len();
                    }
                    last_pos = f.get_offset();
                }
            } else {
                last_pos = f.get_offset();
                while(f.getline()) {
                    if(f.looksLikeFASTQSeparator()) {
                        long long n = number_of_nucleotides_in_read;
                        while(f.getline()) {
                            const size_t number_of_qvchars_in_line = f.len();
                            n -= number_of_qvchars_in_line;
                            if(n <= 0) break;
                        }
                        f.expectHeaderOfEOF();
                        last_pos = f.get_offset();
                        if(!f.getline()) break;
                        f.registerHeaderLine();
                        INSERT_NAME_INTO_TABLE();
						++sequence_count;
                        number_of_nucleotides_in_read = 0;
                    } else {
                        const size_t number_of_nucleotides_in_line = f.len();
                        number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                    }
                    last_pos = f.get_offset();
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
		show_read_names_in_file(argv[i], true, false);
	}
}

void do_len(int argc, char** argv)
{
	bool flag_output_name = false;
    static struct option long_options[] = {
        {"name", no_argument, 0, 'n'},
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
        case 'n':
            flag_output_name = true;
            break;
        }
	}
	for(int i = optind + 1; i < argc; ++i) {
		show_read_names_in_file(argv[i], flag_output_name, true);
	}
}

void print_n50(vector<size_t>& lengths, const bool flag_html, const bool flag_json, const char * contig_or_scaff = "scaffold")
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
        cout << "<tr><td>max</td><td>" << sep_comma(max_length) << "</td><td>1</td></tr>\n";
        cout << "<tr><td>N50</td><td>" << sep_comma(n50_length) << "</td><td>" << sep_comma(n50_sequence_index) << "</td></tr>\n";
        cout << "<tr><td>N70</td><td>" << sep_comma(n70_length) << "</td><td>" << sep_comma(n70_sequence_index) << "</td></tr>\n";
        cout << "<tr><td>N80</td><td>" << sep_comma(n80_length) << "</td><td>" << sep_comma(n80_sequence_index) << "</td></tr>\n";
        cout << "<tr><td>N90</td><td>" << sep_comma(n90_length) << "</td><td>" << sep_comma(n90_sequence_index) << "</td></tr>\n";
        cout << "<tr><td>min</td><td>" << sep_comma(min_length) << "</td><td>" << sep_comma(lengths.size()) << "</td></tr>\n";
        cout << "<tr><td>avg</td><td>" << sep_comma(avg_length) << "</td><td></td></tr>\n";
        cout << "<tr><td>total</td><td>" << sep_comma(total_length) << "</td><td></td></tr>\n";
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
        cout << "Total # of bases = " << sep_comma(total_length) << "\n";
        cout << "# of " << contig_or_scaff << "s = " << sep_comma(lengths.size()) << "\n";
        cout << "Max size = " << sep_comma(max_length) << " (# = 1)\n";
        cout << "N50 " << contig_or_scaff << " size = " << sep_comma(n50_length) << " (# = " << sep_comma(n50_sequence_index) << ")\n";
        cout << "N70 " << contig_or_scaff << " size = " << sep_comma(n70_length) << " (# = " << sep_comma(n70_sequence_index) << ")\n";
        cout << "N80 " << contig_or_scaff << " size = " << sep_comma(n80_length) << " (# = " << sep_comma(n80_sequence_index) << ")\n";
        cout << "N90 " << contig_or_scaff << " size = " << sep_comma(n90_length) << " (# = " << sep_comma(n90_sequence_index) << ")\n";
        cout << "Min size = " << sep_comma(min_length) << "\n";
        cout << "Total " << contig_or_scaff << " # = " << sep_comma(lengths.size()) << "\n";
        cout << "Avg size = " << sep_comma(avg_length) << "\n";
    }
}

void do_stat(int argc, char** argv)
{
    bool flag_html = false;
    bool flag_json = false;
    bool flag_contig = false;
    bool flag_scaffold = false;
    bool flag_all = true;
    static struct option long_options[] = {
        {"html", no_argument, 0, 'h'},
        {"json", no_argument, 0, 'j'},
        {"contig", no_argument, 0, 'c'},
        {"scaffold", no_argument, 0, 's'},
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
            case 'c':
                flag_contig = true;
                flag_all = false;
                break;
            case 's':
                flag_scaffold = true;
                flag_all = false;
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
    int number_of_successfully_processed_files = 0;
    for(int i = optind + 1; i < argc; ++i) {
    	if(calculate_n50_statistics(argv[i], length_of_scaffolds_wgap, length_of_scaffolds_wogap, length_of_contigs)) {
            number_of_successfully_processed_files++;
        }
    }
    if(number_of_successfully_processed_files <= 0) {
        cerr << "ERROR: All file(s) could not be opened." << endl;
        return;
    }
    if(!(length_of_scaffolds_wgap.size() == length_of_scaffolds_wogap.size())) {
        cerr << "Assertion failed. Maybe you found a bug! Please report to the author.\n";
        return;
    }
    if(flag_all || flag_scaffold) {
        if(flag_html) {
            cout << "<table border=\"2\" bgcolor=\"#ffffff\">\n";
            cout << "<tr><th colspan=\"3\" bgcolor=\"#fdfdd4\">Scaffold (w/gap) statistics</th></tr>\n";
        } else if(flag_json) {
            cout << "{\"scaffold_wgap\": ";
        } else {
            cout << "Scaffold (w/gap) statistics\n";
        }
        print_n50(length_of_scaffolds_wgap, flag_html, flag_json, "scaffold");
    }
    if(flag_all) {
        if(flag_html) {
            cout << "<tr><th colspan=\"3\" bgcolor=\"#fdfdd4\">Scaffold (wo/gap) statistics</th></tr>\n";
        } else if(flag_json) {
            cout << ",\"scaffold_wogap\": ";
        } else {
            cout << "\nScaffold (wo/gap) statistics\n";
        }
        print_n50(length_of_scaffolds_wogap, flag_html, flag_json, "scaffold");
    }
    if(flag_all || flag_contig) {
        if(flag_html) {
            cout << "<tr><th colspan=\"3\" bgcolor=\"#fdfdd4\">Contig statistics</th></tr>\n";
        } else if(flag_json) {
            cout << ",\"contig\": ";
        } else {
            cout << "\nContig statistics\n";
        }
        print_n50(length_of_contigs, flag_html, flag_json, "contig");
    }
    if(flag_html) {
        cout << "</table>\n";
    } else if(flag_json) {
        cout << "}\n";
    }
}

void do_split(int argc, char** argv)
{
    bool flag_num = false;
    bool flag_max = false;
    bool flag_exclude_n = false;
    bool flag_return_number_of_partitions_by_errcode = false;
    long long param_specified_num = -1;
    long long param_specified_max = -1;
    string status_file_name;
    string output_file_prefix;
    static struct option long_options[] = {
        {"num", required_argument, 0, 'n'},
        {"max", required_argument, 0, 'm'},
        {"excn", no_argument, 0, 'e'},
        {"prefix", required_argument, 0, 'p'},
        {"retstat", no_argument, 0, 's'},
        {"filestat", required_argument, 0, 'f'},
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
		case 'n':
            flag_num = true;
            param_specified_num = atoll(optarg);
			break;
        case 'm':
            flag_max = true;
            param_specified_max = atoll(optarg);
            break;
        case 'e':
            flag_exclude_n = true;
            break;
        case 'p':
            output_file_prefix = optarg;
            break;
        case 's':
            flag_return_number_of_partitions_by_errcode = true;
            break;
        case 'f':
            status_file_name = optarg;
            break;
        }
	}
    if(output_file_prefix.empty()) {
        if(optind + 1 < argc) {
            output_file_prefix = argv[optind + 1];
        } else {
            cerr << "ERROR: No input files." << endl;
            return;
        }
    }
    if(flag_max && flag_num) {
        cerr << "ERROR: You can use either --max or --num, but not both.\n";
        return;
    }
    int out_file_index = 0;
    long bases_per_file = -1;
    if(flag_max) {
        bases_per_file = param_specified_max;
    } else if(flag_num) {
        vector<size_t> length_of_scaffolds_wgap;
        vector<size_t> length_of_scaffolds_wogap;
        vector<size_t> length_of_contigs;
        for(int i = optind + 1; i < argc; ++i) {
            cerr << "Counting the number of bases ('" << argv[i] << "')\r" << flush;
            calculate_n50_statistics(argv[i], length_of_scaffolds_wgap, length_of_scaffolds_wogap, length_of_contigs);
        }
        cerr << "\n";
        if(!(length_of_scaffolds_wgap.size() == length_of_scaffolds_wogap.size())) {
            cerr << "Assertion failed. Maybe you have found a bug! Please report to the author.\n";
            return;
        }
        const long long total_bases = flag_exclude_n ?
            accumulate(length_of_scaffolds_wogap.begin(), length_of_scaffolds_wogap.end(), 0ll):
            accumulate(length_of_scaffolds_wgap.begin(),  length_of_scaffolds_wgap.end(),  0ll);
        const long long bases_per_file = (total_bases + param_specified_num - 1) / param_specified_num;
        cerr << "Total " << total_bases << " bases (" << (flag_exclude_n ? "wo/ gaps" : "w/ gaps") << ") ";
        cerr << bases_per_file << " bases per file\n" << flush;
    } else { /* never come here */ cerr << "ERROR: Please report to the author." << endl; exit(-1); }
    {
        size_t number_of_nucleotides_in_output_file = bases_per_file + 1;
        ofstream ost;
        for(int findex = optind + 1; findex < argc; ++findex) {
            const char* file_name = argv[findex];
            FileLineBufferWithAutoExpansion f;
            if(!f.open(file_name)) {
                cerr << "Cannot open '" << file_name << "'" << endl;
                continue;
            }
#define OPEN_NEXT_FILE_IF_NEEDED() if(bases_per_file <= number_of_nucleotides_in_output_file) { \
                                       ost.close(); \
                                       number_of_nucleotides_in_output_file = 0; \
                                       out_file_index++; \
                                       char buf[16]; sprintf(buf, "%d", out_file_index); \
                                       const string output_file_name = output_file_prefix + "." + buf; \
                                       ost.open(output_file_name.c_str()); \
                                       if(!ost) { cerr << "ERROR: cannot open an output file '" << output_file_name << "'" << endl; return; } \
                                   }
            if(f.getline()) {
                OPEN_NEXT_FILE_IF_NEEDED();
                ost << f.b << "\n";
                if(!f.looksLikeFASTQHeader()) { 
                    while(f.getline()) {
                        if(f.looksLikeFASTAHeader()) {
                            OPEN_NEXT_FILE_IF_NEEDED();
                            ost << f.b << "\n";
                        } else {
                            ost << f.b << "\n";
                            number_of_nucleotides_in_output_file += f.len();
                            if(flag_exclude_n)
                                number_of_nucleotides_in_output_file -= count_if(f.b, f.b + f.len(), isN);
                        }
                    }
                } else {
                    long long number_of_nucleotides_in_read = 0;
                    while(f.getline()) {
                        if(f.looksLikeFASTQSeparator()) {
                            long long n = number_of_nucleotides_in_read;
                            ost << f.b << "\n";
                            while(f.getline()) {
                                const size_t number_of_qvchars_in_line = f.len();
                                n -= number_of_qvchars_in_line;
                                ost << f.b << "\n";
                                if(n <= 0) break;
                            }
                            f.expectHeaderOfEOF();
                            number_of_nucleotides_in_read = 0;
                            if(!f.getline()) break;
                            f.registerHeaderLine();
                            OPEN_NEXT_FILE_IF_NEEDED();
                            ost << f.b << "\n";
                        } else {
                            const size_t number_of_nucleotides_in_line = f.len();
                            number_of_nucleotides_in_output_file += number_of_nucleotides_in_line;
                            if(flag_exclude_n)
                                number_of_nucleotides_in_output_file -= count_if(f.b, f.b + f.len(), isN);
                            number_of_nucleotides_in_read += number_of_nucleotides_in_line;
                            ost << f.b << "\n";
                        }
                    }
                }
            }
#undef OPEN_NEXT_FILE_IF_NEEDED
        }
    }
    if(out_file_index <= 0) {
        cerr << "No output.\n";
    } else {
        cerr << out_file_index << " files output.\n";
    }
    if(!status_file_name.empty()) {
        ofstream ost(status_file_name.c_str());
        if(!ost) {
            cerr << "ERROR: cannot open a status file '" << status_file_name << "'" << endl;
        } else {
            ost << out_file_index << endl;
        }
    }
    if(flag_return_number_of_partitions_by_errcode)
        exit(out_file_index < 100 ? out_file_index : 100);
    exit(0);
}

void do_index(int argc, char** argv)
{
    bool flag_force = false;
    static struct option long_options[] = {
        {"force", no_argument , 0, 'f'},
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
		case 'f':
            flag_force = true;
			break;
        }
	}
	for(int i = optind + 1; i < argc; ++i) {
		create_index(argv[i], flag_force);
	}
}

static bool check_read_conditions(long long param_start, const char* head_inf, set <string> readNamesToTake, regex re, bool flag_reverse_condition, long long number_of_sequences, long long param_end)
{
    bool current_read_has_been_taken;
    if(param_start == -1) {
        string read_name = get_read_name_from_header(head_inf);
        if(!readNamesToTake.empty()) {
            current_read_has_been_taken = (readNamesToTake.count(read_name) != 0) ^ flag_reverse_condition;
        } else {
            current_read_has_been_taken = regex_search(head_inf, re) ^ flag_reverse_condition;
        }
   } else {
      current_read_has_been_taken = param_start + 1 <= number_of_sequences && number_of_sequences <= param_end;
   }
   return current_read_has_been_taken;
}

void do_extract(int argc, char** argv)
{
    bool flag_reverse_condition = false;
    bool flag_read_from_stdin = false;
    bool flag_output_unique = false;
	bool flag_noindex = false;
	bool flag_index = false;
    bool flag_force = false;
	long long param_start = -1;
	long long param_end = -1;
	long long param_num = -1;

    static struct option long_options[] = {
        {"reverse", no_argument , 0, 'r'},
        {"force", no_argument , 0, 'F'},
        {"seq", required_argument, 0, 's'},
        {"regex", required_argument, 0, 'x'}, //x is chosen because r and e was used for other options
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
    std::regex re;
    vector<string> fileInputs;

    while(true) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "", long_options, &option_index);
		if(c == -1) break;
		switch(c) {
		case 0:
			// you can see long_options[option_index].name/flag and optarg (null if no argument).
			break;
        case 'F':
            flag_force = true;
            break;
		case 'r':
            flag_reverse_condition = true;
			break;
        case 's':
            readNamesToTake.insert(optarg);
            break;
        case 'x':
            re = std::regex(optarg);
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
			create_index(file_name, flag_force);
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
            if(index_older_than_file(file_name, index_file_name)){
                cerr << "Warning: index file " << index_file_name << 
                      " is older than " << file_name << endl;
            }
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
                                    const long long line_start_pos = f.get_offset();
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
                current_read_has_been_taken = check_read_conditions(param_start, f.b, readNamesToTake, re, flag_reverse_condition, number_of_sequences, param_end);
                if(current_read_has_been_taken) cout << f.b << endl;
                if(flag_output_unique) readNamesToTake.insert(get_read_name_from_header(f.b));
                if(!f.looksLikeFASTQHeader()) { 
                    while(f.getline()) {
                        if(f.looksLikeFASTAHeader()) {
                            number_of_sequences++;
                            current_read_has_been_taken = check_read_conditions(param_start, f.b, readNamesToTake, re, flag_reverse_condition, number_of_sequences, param_end);
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
                            current_read_has_been_taken = check_read_conditions(param_start, f.b, readNamesToTake, re, flag_reverse_condition, number_of_sequences, param_end);
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

void investigate_composition(const char* file_name, bool ignore_case, bool flag_only_monomer, bool flag_only_bimer, bool flag_only_trimer, bool flag_dapi_check, bool flag_count_ends)
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
        freq_2_mer[previousCharacters[0]]['\0']++;
        freq_3_mer[previousCharacters[1]][previousCharacters[0]]['\0']++;
        freq_3_mer[previousCharacters[0]]['\0']['\0']++;
    }
    // if you do not need the counts from both ends,
    if(!flag_count_ends) {
        freq_1_mer[0] = 0;
        for(size_t i = 0; i < N; ++i) {
            freq_2_mer[i][0] = 0;
            freq_2_mer[0][i] = 0;
        }
        for(size_t i = 0; i < N; ++i) for(size_t j = 0; j < N; ++j) {
            freq_3_mer[0][i][j] = 0;
            freq_3_mer[i][0][j] = 0;
            freq_3_mer[i][j][0] = 0;
        }
    }
    // show the results
    const size_t total_n_nmers = accumulate(freq_1_mer, freq_1_mer + sizeof(freq_1_mer) / sizeof(size_t), 0llu);
    const size_t total_n2_nmers = accumulate(reinterpret_cast<size_t*>(freq_2_mer), reinterpret_cast<size_t*>(freq_2_mer) + sizeof(freq_2_mer) / sizeof(size_t), 0llu);
    const size_t total_n3_nmers = accumulate(reinterpret_cast<size_t*>(freq_3_mer), reinterpret_cast<size_t*>(freq_3_mer) + sizeof(freq_3_mer) / sizeof(size_t), 0llu);
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
                if(freq_2_mer[i][j] != 0) cout << "\t" << asterisk_if_nulchar(i) << asterisk_if_nulchar(j) << "\t" << freq_2_mer[i][j] << "\t" << (double(freq_2_mer[i][j]) / total_n2_nmers) << "\n";
            }
        }
    }
    if(!flag_only_monomer && !flag_only_bimer && !flag_dapi_check) {
        cout << "3-mer stats\n";
        for(size_t i = 0; i < N; ++i) {
            for(size_t j = 0; j < N; ++j) {
                for(size_t k = 0; k < N; ++k) {
                    if(freq_3_mer[i][j][k] != 0) cout << "\t" << asterisk_if_nulchar(i) << asterisk_if_nulchar(j) << asterisk_if_nulchar(k) << "\t" << freq_3_mer[i][j][k] << "\t" << (double(freq_3_mer[i][j][k]) / total_n3_nmers) << "\n";
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
    bool flag_count_ends   = false;
    static struct option long_options[] = {
        {"ignorecase", no_argument, 0, 'i'},
        {"monomer", no_argument, 0, '1'},
        {"bimer", no_argument, 0, '2'},
        {"trimer", no_argument, 0, '3'},
        {"dapicheck", no_argument, 0, 'd'},
        {"countends", no_argument, 0, 'c'},
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
        case 'c':
            flag_count_ends = true;
            break;
        }
    }
    for(int i = optind + 1; i < argc; ++i) {
        investigate_composition(argv[i], flag_ignore_case, flag_only_monomer, flag_only_bimer, flag_only_trimer, flag_dapi_check, flag_count_ends);
    }
}

struct Sequence {
    size_t file_position;
    string name;
    string description;
    vector<char> sequence;
    vector<char> qv;
    Sequence() {}
    Sequence(size_t file_position, string name, string description, const vector<char>& sequence) : file_position(file_position), name(name), description(description), sequence(sequence) {}
    Sequence(size_t file_position, string name, string description, const vector<char>& sequence, const vector<char>& qv) : file_position(file_position), name(name), description(description), sequence(sequence), qv(qv) {}
    inline bool operator < (const Sequence& b) const {
        return name < b.name;
    }
};

char complement_char(char c)
{
//      'atgcrymkdhvbswn'
//      'tacgyrkmhdbvswn'
// a compiler might optimize to a table lookup...
  switch(c){
    case 'a': return 't';
    case 'c': return 'g';
    case 'g': return 'c';
    case 't': return 'a';
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';

    case 'n': return 'n';
    case 'N': return 'N';

    case 'r': return 'y';
    case 'y': return 'r';
    case 'm': return 'k';
    case 'k': return 'm';
    case 'd': return 'h';
    case 'h': return 'd';
    case 'v': return 'b';
    case 'b': return 'v';
    case 's': return 's';
    case 'w': return 'w';

    case 'R': return 'Y';
    case 'Y': return 'R';
    case 'M': return 'K';
    case 'K': return 'M';
    case 'D': return 'H';
    case 'H': return 'D';
    case 'V': return 'B';
    case 'B': return 'V';
    case 'S': return 'S';
    case 'W': return 'W';
  }
  return c;
}
class GenomeEditScript {
    bool is_fastq;
    bool has_file_type_determined;
    bool is_verbose;
    typedef map<string, Sequence> SequenceMap;
    SequenceMap sequences;

    static const size_t FOLD_WITH_THIS_SIZE;

    string trim(const string& s) {
        string::size_type f = s.find_first_not_of(" \t");
        if(f == string::npos) f = 0;
        string::size_type l = s.find_last_not_of(" \t");
        if(l == string::npos) l = s.size();
        return s.substr(f, l - f + 1);
    }
    vector<string> split(const string& s) {
        vector<string> retval;
        string::size_type i = 0;
        while(i < s.size()) {
            const bool this_term_is_quoted = s[i] == '"';
            if(this_term_is_quoted) {
                string::size_type j = i + 1;
                while(j < s.size()) {
                    if(s[j] == '"') break;
                    if(s[j] == '\\') {
                        j++;
                        if(s.size() <= j) break;
                    }
                    j++;
                }
                if(s.size() <= j) {
                    cerr << "ERROR: unmatched quote '\"'\nTHIS LINE: " << s << endl;
                    exit(2);
                }
                retval.push_back(s.substr(i + 1, j - i - 1));
                string::size_type k = s.find_first_not_of(" \t", j + 1);
                if(k == string::npos) break; // I think this never happens because s is already trimmed.
                i = k;
            } else {
                string::size_type j = s.find_first_of(" \t", i);
                if(j == string::npos) {
                    retval.push_back(s.substr(i));
                    break;
                }
                retval.push_back(s.substr(i, j - i));
                if(s.size() <= j + 1) break;
                string::size_type k = s.find_first_not_of(" \t", j + 1);
                if(k == string::npos) break; // I think this never happens because s is already trimmed.
                i = k;
            }
        }
        return retval;
    }
public:
    bool loadEntireSeq(const char* seq_file_name) {
        FileLineBufferWithAutoExpansion f;
        if(!f.open(seq_file_name)) {
            cerr << "Cannot open '" << seq_file_name << "'" << endl;
            return false;
        }
        string current_read_name_and_description;
        if(f.getline()) {
            f.registerHeaderLineWithDesc();
            vector<char> current_sequence;
            current_sequence.reserve(64 * 1024);
            if(!f.looksLikeFASTQHeader()) { 
                // This should be FASTA
                if(has_file_type_determined) {
                    if(is_fastq) {
                        cerr << "ERROR: The file '" << seq_file_name << "' looks like a FASTA file, but the previous file(s) is in FASTQ format.\n";
                        cerr << "       You cannot mix both the formats.\n";
                        return false;
                    }
                } else {
                    is_fastq = false;
                }
                while(f.getline()) {
                    if(f.looksLikeFASTAHeader()) {
                        const string sequence_name = f.getSequenceName();
                        sequences[sequence_name] = Sequence(f.getLineCount(), sequence_name, f.getSequenceDescription(), current_sequence);
                        current_sequence.resize(0);
                        f.registerHeaderLineWithDesc();
                    } else {
                        const int number_of_chars_in_line = f.len();
                        current_sequence.insert(current_sequence.end(), f.b, f.b + f.len());
                    }
                }
                const string sequence_name = f.getSequenceName();
                sequences[sequence_name] = Sequence(f.getLineCount(), sequence_name, f.getSequenceDescription(), current_sequence);
            } else {
                if(has_file_type_determined) {
                    if(!is_fastq) {
                        cerr << "ERROR: The file '" << seq_file_name << "' looks like a FASTQ file, but the previous file(s) is in FASTA format.\n";
                        cerr << "       You cannot mix both the formats.\n";
                        return false;
                    }
                } else {
                    is_fastq = true;
                }
                vector<char> current_qv;
                current_qv.reserve(64 * 1024);
                while(f.getline()) {
                    if(f.looksLikeFASTQSeparator()) {
                        long long n = current_sequence.size();
                        while(f.getline()) {
                            const size_t number_of_qvchars_in_line = f.len();
                            current_qv.insert(current_qv.end(), f.b, f.b + number_of_qvchars_in_line);
                            n -= number_of_qvchars_in_line;
                            if(n <= 0) break;
                        }
                        f.expectHeaderOfEOF();
                        {
                            const string sequence_name = f.getSequenceName();
                            sequences[sequence_name] = Sequence(f.getLineCount(), sequence_name, f.getSequenceDescription(), current_sequence, current_qv);
                            current_sequence.resize(0);
                            current_qv.resize(0);
                        }
                        if(!f.getline()) break;
                        f.registerHeaderLineWithDesc();
                    } else {
                        const int number_of_chars_in_line = f.len();
                        current_sequence.insert(current_sequence.end(), f.b, f.b + number_of_chars_in_line);
                    }
                }
            }
        }
        return true;
    }
    ostream& output_with_fold(ostream& os, const vector<char>& s) {
        size_t cursor = 0;
        while(cursor < s.size()) {
            const size_t len = std::min<size_t>(FOLD_WITH_THIS_SIZE, s.size() - cursor);
            for(size_t i = 0; i < len; i++) os << s[cursor + i];
            os << '\n';
            cursor += len;
        }
        return os;
    }
    typedef pair<size_t, const Sequence *> OrderPreservedSequence;
    bool saveEntireSeq(const char* sequence_file_name) {
        ofstream ost(sequence_file_name);
        if(!ost) {
            cerr << "Could not open file '" << sequence_file_name << "'" << endl;
            return false;
        }
        vector<OrderPreservedSequence> ops;
        for(SequenceMap::const_iterator cit = sequences.begin(); cit != sequences.end(); ++cit) {
            ops.push_back(OrderPreservedSequence(cit->second.file_position, &(cit->second)));
        }
        sort(ops.begin(), ops.end());
        if(is_fastq) {
            for(size_t i = 0; i < ops.size(); i++) {
                const Sequence& s = *(ops[i].second);
                ost << '@' << s.name;
                if(!s.description.empty()) ost << ' ' << s.description;
                ost << '\n' << s.sequence << "\n+\n" << s.qv << '\n';
            }
        } else {
            for(size_t i = 0; i < ops.size(); i++) {
                const Sequence& s = *(ops[i].second);
                ost << '>' << s.name;
                if(!s.description.empty()) ost << ' ' << s.description;
                ost << '\n';
                output_with_fold(ost, s.sequence);
            }
        }
        return true;
    }
    bool saveOneSeq(const string& sequence_file_name, const string& sequence_name) {
        if(sequences.count(sequence_name) == 0) {
            cerr << "Could not find a sequence '" << sequence_name << "'\n"; return false;
        }
        ofstream ost(sequence_file_name.c_str());
        if(!ost) {
            cerr << "Could not open file '" << sequence_file_name << "'" << endl;
            return false;
        }
        const Sequence& s = sequences[sequence_name];
        if(is_fastq) {
            ost << '@' << s.name;
            if(!s.description.empty()) ost << ' ' << s.description;
            ost << '\n' << s.sequence << "\n+\n" << s.qv << '\n';
        } else {
            ost << '>' << s.name;
            if(!s.description.empty()) ost << ' ' << s.description;
            ost << '\n';
            ost << s.sequence << endl;
        }
        return true;
    }
    bool loadOneSeq(const string& sequence_file_name, const string& sequence_name) {
        if(sequences.count(sequence_name)) {
            cerr << "There is already a sequence '" << sequence_name << "'\n"; return false;
        }
        if(!doesIndexExist(sequence_file_name.c_str())) {
            cerr << "You have to create the index of '" << sequence_file_name << "' first.\n"; return false;
        }
        FileLineBufferWithAutoExpansion f;
        if(!f.open(sequence_file_name.c_str())) {
            cerr << "Could not open file '" << sequence_file_name << "'" << endl; return false;
        }
        const bool is_fastq = is_file_fastq(sequence_file_name.c_str());
        if(has_file_type_determined) {
            if(this->is_fastq != is_fastq) {
                cerr << "You cannot mix both FASTA/FASTQ files." << endl; return false;
            }
        } else {
            this->is_fastq = is_fastq;
            has_file_type_determined = true;
        }
        if(is_verbose) { cerr << "MODE: " << (is_fastq ? "fastq" : "fasta") << endl; }
        const string index_file_name = get_index_file_name(sequence_file_name.c_str());
        try {
            sqdb::Db db(index_file_name.c_str());
            sqdb::Statement stmt = db.Query("select pos from seqpos where name=?");
            stmt.Bind(1, sequence_name);
            if(stmt.Next()) {
                const long long pos = stmt.GetField(0);
                f.seekg(pos);
                if(is_verbose) { cerr << "SEEK to " << pos << endl; }
                if(f.fail()) {
                    cerr << "'" << sequence_name << "' is missing in the file. The file is too small. Maybe the index is old?\n"; return false;
                }
                if(!f.getline()) {
                    cerr << "'" << sequence_name << "' is missing in the file. The header is not as expected. Maybe the index is old?\n"; return false;
                }
                f.registerHeaderLineWithDesc();
                vector<char> current_sequence, current_qv;
                if(is_verbose) { cerr << "HEADER: " << f.b << endl; }
                if(is_fastq) {
                    size_t number_of_nucleotides_in_read = 0;
                    while(f.getline()) {
                        if(f.looksLikeFASTQSeparator()) {
                            long long n = current_sequence.size();
                            while(f.getline()) {
                                const size_t number_of_qvchars_in_line = f.len();
                                current_qv.insert(current_qv.end(), f.b, f.b + number_of_qvchars_in_line);
                                n -= number_of_qvchars_in_line;
                                if(n <= 0) break;
                            }
                            sequences[f.getSequenceName()] = Sequence(f.getLineCount(), f.getSequenceName(), f.getSequenceDescription(), current_sequence, current_qv);
                            break;
                        } else {
                            current_sequence.insert(current_sequence.end(), f.b, f.b + f.len());
                        }
                    }
                } else {
                    while(f.getline()) {
                        if(f.looksLikeFASTAHeader()) break;
                        current_sequence.insert(current_sequence.end(), f.b, f.b + f.len());
                    }
                    sequences[f.getSequenceName()] = Sequence(f.getLineCount(), f.getSequenceName(), f.getSequenceDescription(), current_sequence);
                }
            } else {
                cerr << "'" << sequence_name << "' was not found.\n"; return false;
            }
        } catch (const sqdb::Exception& e) {
            cerr << "ERROR: database error. " << e.GetErrorMsg() << endl;
            return false;
        }
        return true;
    }
    bool renameSequence(const string& old_name, const string& new_name) {
        if(sequences.count(old_name) == 0) {
            cerr << "Could not find a sequence '" << old_name << "', which you tried to rename into '" << new_name << "'\n"; return false;
        }
        if(sequences.count(new_name)) {
            cerr << "There is already a sequence '" << new_name << "', to which you tried to rename a sequence '" << old_name << "'\n"; return false;
        }
        sequences[new_name] = sequences[old_name];
        sequences[new_name].name = new_name;
        sequences.erase(sequences.find(old_name));
        return true;
    }
    bool deleteSequence(const string& old_name) {
        if(sequences.count(old_name) == 0) {
            cerr << "Could not find a sequence '" << old_name << "', which you tried to delete\n"; return false;
        }
        sequences.erase(old_name);
        return true;
    }
    bool duplicateSequence(const string& old_name, const string& new_name) {
        if(sequences.count(old_name) == 0) {
            cerr << "Could not find a sequence '" << old_name << "', which you tried to duplicate to '" << new_name << "'\n"; return false;
        }
        if(sequences.count(new_name)) {
            cerr << "There is already a sequence '" << new_name << "', to which you tried to duplicate a sequence '" << old_name << "'\n"; return false;
        }
        sequences[new_name] = sequences[old_name];
        sequences[new_name].name = new_name;
        return true;
    }
    bool complementSequence(const string& old_name, const string& new_name) {
        if(sequences.count(old_name) == 0) {
            cerr << "Could not find a sequence '" << old_name << "', which you tried to duplicate to '" << new_name << "'\n"; return false;
        }
        if(sequences.count(new_name)) {
            cerr << "There is already a sequence '" << new_name << "', to which you tried to duplicate a sequence '" << old_name << "'\n"; return false;
        }
        size_t length = sequences[old_name].sequence.size();
        sequences[new_name].name = new_name;
        sequences[new_name].description = sequences[old_name].description;
        sequences[new_name].sequence = vector<char>(length);
        for(int i = 0; i < length; ++i) {
          char c = complement_char(sequences[old_name].sequence[i]);
          sequences[new_name].sequence[length-1-i] = c;
        }
        sequences[new_name].qv = sequences[old_name].qv;
        reverse(sequences[new_name].qv.begin(),sequences[new_name].qv.end());
        sequences.erase(old_name);
        return true;
    }
    bool joinSequence(const string& left_sequence_name, const string& right_sequence_name, const string& new_sequence_name) {
        if(sequences.count(left_sequence_name) == 0) {
            cout << "Could not find a sequence '" << left_sequence_name << "', which you tried to join to '" << new_sequence_name << "'\n"; return false;
        }
        if(sequences.count(right_sequence_name) == 0) {
            cout << "Could not find a sequence '" << right_sequence_name << "', which you tried to join to '" << new_sequence_name << "'\n"; return false;
        }
        if(sequences.count(new_sequence_name)) {
            cout << "There is already a sequence '" << new_sequence_name << "', to which you tried to join sequences ('" << left_sequence_name << "', '" << right_sequence_name << "')\n"; return false;
        }
        Sequence& left_seq = sequences[left_sequence_name];
        Sequence& right_seq = sequences[right_sequence_name];
        Sequence& new_seq = sequences[new_sequence_name];
        new_seq.name = new_sequence_name;
        new_seq.sequence.resize(left_seq.sequence.size() + right_seq.sequence.size());
        copy(left_seq.sequence.begin(), left_seq.sequence.end(), new_seq.sequence.begin());
        copy(right_seq.sequence.begin(), right_seq.sequence.end(), new_seq.sequence.begin() + left_seq.sequence.size());
        if(!left_seq.qv.empty() && !right_seq.qv.empty()) {
            new_seq.qv.resize(left_seq.qv.size() + right_seq.qv.size());
            copy(left_seq.qv.begin(), left_seq.qv.end(), new_seq.qv.begin());
            copy(right_seq.qv.begin(), right_seq.qv.end(), new_seq.qv.begin() + left_seq.qv.size());
        }
        sequences.erase(left_sequence_name);
        sequences.erase(right_sequence_name);
        return true;
    }
    bool setSequenceDescription(const vector<string>& subargs) {
        const string& seq_name = subargs.front();
        if(sequences.count(seq_name) == 0) {
            cerr << "Could not find a sequence '" << seq_name << "'\n"; return false;
        }
        string desc;
        for(size_t i = 1; i < subargs.size(); i++) {
            if(1 < i) desc += ' ';
            desc += subargs[i];
        }
        sequences[seq_name].description = desc;
        return true;
    }
    bool splitSequence(const string& seq_name, const string& split_position_str, const string& left_seq_name, const string& right_seq_name) {
        if(sequences.count(seq_name) == 0) {
            cerr << "Could not find a sequence '" << seq_name << "'\n"; return false;
        }
        if(sequences.count(left_seq_name)) {
            cerr << "There is already a sequence '" << left_seq_name << "'\n"; return false;
        }
        if(sequences.count(right_seq_name)) {
            cerr << "There is already a sequence '" << right_seq_name << "'\n"; return false;
        }
        const long long pos = atoll(split_position_str.c_str());
        Sequence& orig_seq = sequences[seq_name];
        if(pos < 0 || orig_seq.sequence.size() < pos) {
            cerr << "Sequence '" << seq_name << "' was " << orig_seq.sequence.size() << " bp in length, but you tried to split it at " << pos << " bp, which is out of range" << endl; return false;
        }
        Sequence& left_seq = sequences[left_seq_name];
        Sequence& right_seq = sequences[right_seq_name];
        left_seq.file_position = orig_seq.file_position;
        right_seq.file_position = orig_seq.file_position + 1;
        left_seq.name = left_seq_name;
        right_seq.name = right_seq_name;
        left_seq.sequence.assign(orig_seq.sequence.begin(), orig_seq.sequence.begin() + pos);
        right_seq.sequence.assign(orig_seq.sequence.begin() + pos, orig_seq.sequence.end());
        if(!orig_seq.qv.empty()) {
            left_seq.qv.assign(orig_seq.qv.begin(), orig_seq.qv.begin() + pos);
            right_seq.qv.assign(orig_seq.qv.begin() + pos, orig_seq.qv.end());
        }
        sequences.erase(seq_name);
        return true;
    }
    void printSequence(const vector<string>& subargs) {
        const string& seq_name = subargs[0];
        if(sequences.count(seq_name) == 0) {
            cerr << "Could not find a sequence '" << seq_name << "'\n"; return;
        }
        const Sequence& s = sequences[seq_name];
        const long long start_pos = subargs.size() <= 1 ? 0ull : atoll(subargs[1].c_str());
        const long long end_pos = subargs.size() <= 2 ? s.sequence.size() : atoll(subargs[2].c_str());
        if(start_pos < 0 || s.sequence.size() < start_pos) {
            cerr << "Start position (" << start_pos << " bp; 0-origin, inclusive) is out of range. The length of the sequence '" << seq_name << "' is " << s.sequence.size() << endl; return;
        }
        if(end_pos < 0 || s.sequence.size() < end_pos) {
            cerr << "End position (" << end_pos << " bp; 0-origin, exclusive) is out of range. The length of the sequence '" << seq_name << "' is " << s.sequence.size() << endl; return;
        }
        cout << (is_fastq ? '@' : '>') << seq_name;
        if(0 < start_pos || end_pos < s.sequence.size()) cout << ' ' << (start_pos + 1) << ':' << end_pos;
        cout << '\n';
        if(is_fastq) {
            for(size_t i = start_pos; i < end_pos; i++) cout << s.sequence[i];
            cout << "\n+\n";
            for(size_t i = start_pos; i < end_pos; i++) cout << s.qv[i];
            cout << '\n';
        } else {
            for(size_t i = start_pos; i < end_pos; ) {
                size_t len = min<size_t>(FOLD_WITH_THIS_SIZE, end_pos - i); 
                for(size_t j = 0; j < len; j++) cout << s.sequence[i + j];
                cout << '\n';
                i += len;
            }
        }
    }
    bool trimSequence(int direction, const string& seq_name, const string& amount) {
        if(sequences.count(seq_name) == 0) {
            cerr << "Could not find a sequence '" << seq_name << "'\n"; return false;
        }
        if(direction != 5 && direction != 3) {
            cerr << "ERROR: logic error (trimSequence)\n"; return false;
        }
        long long amount_int = atoll(amount.c_str());
        if(amount_int == 0) {
            cerr << "WARNING: You tried to trim 0 bases from '" << seq_name << "'\n";
        }
        if(is_verbose) {
            cerr << "TRIM '" << seq_name << "' at " << direction << "'-end by " << amount_int << " bp" << endl;
        }
        vector<char>& seq = sequences[seq_name].sequence;
        vector<char>& qv = sequences[seq_name].qv;
        if(seq.size() <= amount_int) {
            cerr << "WARNING: You tried to trim '" << seq_name << "' by " << amount_int << "bp, but its size is only " << seq.size() << " bp.\n";
            seq.resize(0);
            sequences[seq_name].qv.resize(0);
            return true;
        }
        if(direction == 5) {
            seq.erase(seq.begin(), seq.begin() + amount_int);
            if(!qv.empty())
                qv.erase(qv.begin(), qv.begin() + amount_int);
        } else {
            seq.resize(seq.size() - amount_int);
            if(!qv.empty())
                qv.resize(qv.size() - amount_int);
        }
        return true;
    }
    void executeScript(const char* script_file_name) {
        ifstream ist(script_file_name);
        if(!ist) {
            cerr << "ERROR: Cannot open '" << script_file_name << "'\n" << endl;
            return;
        }
        size_t script_line_number = 0;
        string line;
        while(getline(ist, line)) {
            script_line_number++;
            line = trim(line);
            if(line.empty() || line[0] == '#') continue;
            vector<string> args = split(line);
            if(args.empty()) continue;
            const string& cmd = args.front();
            vector<string> subargs(args.begin() + 1, args.end());
            if(is_verbose) {
                cerr << "COMMAND: " << cmd << " SUBARGS: [";
                for(size_t i = 0; i < subargs.size(); i++) { if(i) cerr << ", "; cerr << '"' << subargs[i] << '"'; }
                cerr << "]" << endl;
            }
            if(cmd == "loadall") {
                if(subargs.size() != 1) {
                    cerr << "ERROR: # of the arguments is invalid.\n";
                    cerr << "usage: loadall <file name (FASTA or FASTQ)>" << endl;
                    return;
                }
                if(!loadEntireSeq(subargs.front().c_str())) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "saveall") {
                if(subargs.size() != 1) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: saveall <file name>\n";
                    cerr << "The file type (either FASTA or FASTQ) is automatically determined by the input file" << endl;
                    return;
                }
                if(!saveEntireSeq(subargs.front().c_str())) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "loadone") {
                if(subargs.size() != 2) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: loadone <file name> <sequence name>\n";
                    return;
                }
                if(!loadOneSeq(subargs[0], subargs[1])) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "saveone") {
                if(subargs.size() != 2) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: saveone <file name> <sequence name>\n";
                    return;
                }
                if(!saveOneSeq(subargs[0], subargs[1])) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "rename") {
                if(subargs.size() != 2) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: rename <old name> <new name>\n";
                    return;
                }
                if(!renameSequence(subargs[0], subargs[1])) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "setdesc") {
                if(subargs.size() < 1) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: setdesc <sequence name> <descriptions (spaces are allowed)>\n";
                    cerr << "Note that successive space characters are compressed into one.\n";
                    return;
                }
                if(!setSequenceDescription(subargs)) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "trim5") {
                if(subargs.size() != 2) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: trim5 <sequence name> <trim size (in bp)>\n";
                    return;
                }
                if(!trimSequence(5, subargs[0], subargs[1])) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "trim3") {
                if(subargs.size() != 2) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: trim3 <sequence name> <trim size (in bp)>\n";
                    return;
                }
                if(!trimSequence(3, subargs[0], subargs[1])) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "print") {
                if(subargs.size() < 1 || 3 < subargs.size()) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: print <sequence name> [<start position (0-origin, inclusive)> <end position (0-origin, exclusive)]\n";
                    return;
                }
                printSequence(subargs);
            } else if(cmd == "split") {
                if(subargs.size() != 4) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: split <sequence name> <split position (0-origin)> <left sequence name> <right sequence name>\n";
                    cerr << "The base at the split position goes to the right sequence.\n" << endl;
                    return;
                }
                if(!splitSequence(subargs[0], subargs[1], subargs[2], subargs[3])) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "dupseq") {
                if(subargs.size() != 2) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: dupseq <sequence name> <new sequence name>\n";
                    return;
                }
                if(!duplicateSequence(subargs[0], subargs[1])) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "delete") {
                if(subargs.size() != 1) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: delete <sequence name>\n";
                    return;
                }
                if(!deleteSequence(subargs[0])) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "complement") {
                if(subargs.size() != 2) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: complement <sequence name> <new sequence name>\n";
                    return;
                }
                if(!complementSequence(subargs[0], subargs[1])) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else if(cmd == "join") {
                if(subargs.size() != 3) {
                    cerr << "ERROR: # of the argument is invalid.\n";
                    cerr << "usage: join <left sequence> <right sequence>\n";
                    return;
                }
                if(!joinSequence(subargs[0], subargs[1], subargs[2])) {
                    cerr << "ERROR: an error occurred." << endl;
                    exit(2);
                }
            } else {
                cerr << "ERROR: unknown command '" << cmd << "' at line " << script_line_number << endl;
                return;
            }
        }
    }
    GenomeEditScript(char** argv, int start, int end, bool verbose) : is_verbose(verbose) {
        is_fastq = false;
        has_file_type_determined = false;
        for(int i = start + 1; i < end; ++i) {
            if(!loadEntireSeq(argv[i])) {
                cerr << "ERROR: an error occurred." << endl;
                exit(2);
            }
        }
        executeScript(argv[start]);
    }
};

void do_edit(int argc, char** argv)
{
    bool flag_verbose = false;
    static struct option long_options[] = {
        {"verbose", no_argument, 0, 'v'},
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
		case 'v':
            flag_verbose = true;
			break;
        }
    }
    GenomeEditScript ges(argv, optind + 1, argc, flag_verbose);
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
        cerr << "--force\tForce on error.\n";
        return;
    }
    if(subcmd == "len") {
        cerr << "Usage: fatt len [options...] <FAST(A|Q) files>\n\n";
        cerr << "--name\tAdd the name of the sequences in the second column.\n\n";
        cerr << "It outputs the length of the sequences in given files.\n";
        return;
    }
    if(subcmd == "stat") {
        cerr << "Usage: fatt stat [options...] <FAST(A|Q) files>\n\n";
        cerr << "--html\tOutput in HTML format.\n";
        cerr << "--json\tOutput in JSON format.\n";
        cerr << "--contig\tOutput contig statistics.\n";
        cerr << "--scaffold\tOutput statistics of scaffold with gaps.\n";
        cerr << "If neither of --contig nor --scaffold is specified, statistics of scaffold with gaps, scaffold without gaps, and contigs are reported.\n";
        return;
    }
    if(subcmd == "index") {
        cerr << "Usage: fatt index [options...] <FAST(A|Q) files>\n\n";
        cerr << "--force\tRemove an existing index if any.\n\n";
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
    if(subcmd == "composition") {
        cerr << "Usage: fatt composition [options...] <FAST(A|Q) files>\n\n";
        cerr << "--ignorecase\tIgnore case ('A' and 'a' will be considered as identical)\n";
        cerr << "--monomer\tShow only monomers\n";
        cerr << "--bimer\tShow only bimers\n";
        cerr << "--trimer\tShow only trimers\n";
        cerr << "--dapicheck\tShow DAPI-staining related stats\n";
        return;
    }
    if(subcmd == "edit") {
        cerr << "Usage: fatt edit [options...] <edit script> [FAST(A|Q) files]\n\n";
        cerr << "No options available.\n";
        cerr << "Genome Edit Script (GES) specification:\n";
        cerr << "\tYou can put at most one command in a line; you need n lines for n commands.\n";
        cerr << "\tBlank lines and lines starting with '#' will be ignored.\n";
        cerr << "\tA line with a command looks like this:\n";
        cerr << "\t\tCommand arg1 arg2 arg3 ...\n";
        cerr << "\tThe list of available commands are the following:\n";
        cerr << "\t\tloadall\tload entire sequences in a file (arg1) into memory\n";
        cerr << "\t\tsaveall\tsave entire sequences in memory into a file (arg1)\n";
        cerr << "\t\tloadone\tload a specified sequence in a file (arg1) with name arg2 into memory (index is used when available)\n";
        cerr << "\t\tsaveone\tsave a specified sequence (arg2) into a file (arg1)\n";
        cerr << "\t\trename\trename a sequence (arg1) into arg2\n";
        cerr << "\t\tdelete\tdelete a sequence (arg1)\n";
        cerr << "\t\tsetdesc\tset a description (arg2) to a specified sequence (arg1)\n";
        cerr << "\t\ttrim5\ttrim the 5'-end of a specified sequence (arg1) in memory by arg2 bp\n";
        cerr << "\t\ttrim3\ttrim the 3'-end of a specified sequence (arg1) in memory by arg2 bp\n";
        cerr << "\t\tprint\tprint a specified sequence (arg1) in memory; range [arg2, arg3) is optional\n";
        cerr << "\t\tsplit\tsplit a specified sequence (arg1) at position (arg2; 0-origin; the base at arg2 belongs to the latter fragment) into arg3 and arg4\n";
        cerr << "\t\tdupseq\tduplicate a specified sequence (arg1) and name it arg2\n";
        cerr << "\t\tcomplement\tconstruct reverse complement of a specified sequence (arg1) and name it arg2\n";
        cerr << "\t\tjoin\tjoin two specified sequences (arg1, arg2) into one (arg3)\n";
        cerr << "\nFiles given in the command line will be loaded by loadall command before executing the edit script.\n";
        return;
    }
    if(subcmd == "split") {
        cerr << "Usage: fatt split --num=n <FAST(A|Q) file>\n";
        cerr << "                split the file into n files.\n";
        cerr << "       fatt split --max=n <FAST(A|Q) file>\n";
        cerr << "                split the file into files of around n bytes\n";
        cerr << "Split the input files into multiple files.\n\n";
        cerr << "--prefix=name\tSpecify the prefix of output file name. If not specified, it will be the first input file\n";
        cerr << "\n";
        cerr << "'fatt split --num=3 huge.fastq' will split huge.fastq into 3 files.\n";
        cerr << "fatt counts the number of bases in huge.fastq in the first phase.\n";
        cerr << "Then, fatt copies sequences in the input to a chunk until the chunk has more than (total/3) bases.\n";
        cerr << "The last chunk is usually smaller than the other chunks. Note that the number of output chunks\n";
        cerr << "might be smaller than the specified number in some cases (e.g., --num=10 for 7 sequences).\n";
        cerr << "The output chunks will be huge.fastq.1, huge.fastq.2, ..., and so on.\n";
        cerr << "You can specify the prefix for output files by giving --prefix option.\n";
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
    cerr << "\tcomposition\tcalculate the 1-, 2-, 3-mer composition.\n";
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
    cerr << "\tedit\tedit sequences by DSL (domain-specific language)\n";
    cerr << "\tsplit\tsplit sequences into multiple files\n";
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
    if(commandString == "edit") {
        do_edit(argc, argv);
        return;
    }
    if(commandString == "split") {
        do_split(argc, argv);
        return; // NOTE: do_split never returns.
    }
    // Help or error.
    if(commandString != "help") {
        cerr << "ERROR: Unknown command '" << commandString << "'" << endl;
    }
	if(commandString == "help" && 3 <= argc) {
		show_help(argv[2]);
	} else {
    	show_help("");
	}
}

const size_t GenomeEditScript::FOLD_WITH_THIS_SIZE = 70u;
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
