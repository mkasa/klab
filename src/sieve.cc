
// -*- mode:C++; c-basic-offset:2; tab-width:2 -*-
// @author Masahiro Kasahara (masahiro@kasahara.ws)
//
#include <iostream>
#include <fstream>
#include <string>
#include <cerrno>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <getopt.h>
#include <sys/time.h>
#include <unistd.h>

using namespace std;

static int flag_debug = 0;

static void print_usage(ostream& os)
{
	os << "usage: sieve -p 30 <file(s)>     Take 30% of the lines in the file(s).\n";
	os << "       sieve -c 40 <file(s)>     Take 40 lines out of the lines in the file(s).\n";
	os << "\n";
	os << "options:\n";
	os << "       -p <percent>              Take <percent>% of the lines (0 <= percent <= 100).\n";
	os << "       -c <count>                Take exactly <count> lines (or every line if the\n";
	os << "                                 input has fewer lines than that).\n";
	os << "       --debug                   Show progress messages on stderr.\n";
	os << "\n";
	os << "Lines are chosen uniformly at random; every line has the same chance of\n";
	os << "being taken. The default (no -p/-c) is to take 50% of the lines.\n";
}

// Seed drand48() with something that differs between two invocations started
// within the same second (time(NULL) alone does not).
static void seed_random()
{
	struct timeval tv;
	if(gettimeofday(&tv, NULL) != 0) {
		tv.tv_sec  = static_cast<long>(time(NULL));
		tv.tv_usec = 0;
	}
	// Unsigned arithmetic: wrap-around is well defined, signed overflow is not.
	unsigned long seed = static_cast<unsigned long>(tv.tv_sec);
	seed ^= static_cast<unsigned long>(tv.tv_usec) * 1000003UL;
	seed ^= static_cast<unsigned long>(getpid()) * 7919UL;
	seed ^= static_cast<unsigned long>(reinterpret_cast<size_t>(&tv)); // stack address (ASLR)
	srand48(static_cast<long>(seed & 0x7fffffffUL));
}

int do_sieve(int argc, char *argv[], int optind, long long param_percent, long long param_lines)
{
	long long lines = 0;
	string tmp;
	// count lines
	for(int i = optind; i < argc; i++) {
		if(flag_debug) {
			cerr << "Opening " << argv[i] << endl;
		}
		// NOTE: ifstream, not fstream. The default fstream mode is
		// ios::in|ios::out, which needs write permission on the file even
		// though we only ever read it.
		ifstream fin(argv[i]);
		if(!fin) {
			cerr << "Could not open '" << argv[i] << "'" << endl;
			return 1;
		}
		while(getline(fin, tmp)) lines++;
		if(flag_debug) {
			cerr << lines << " lines so far" << endl;
		}
	}
	// NOTE: computed in long long; 'lines * percent' overflows a 32bit int
	// for inputs of more than ~21M lines, which used to yield a negative
	// take_lines and thus no output at all.
	long long take_lines = lines / 2;
	if(param_percent != -1) { take_lines = lines * param_percent / 100; }
	if(param_lines   != -1) { take_lines = param_lines; if(lines < take_lines) take_lines = lines; }
	if(take_lines < 0) take_lines = 0;
	if(flag_debug) {
		cerr << "We will take " << take_lines << " lines" << endl;
	}
	// Selection sampling (Knuth, TAOCP vol. 2, algorithm 3.4.2 S): while
	// walking over the input, take the current line with probability
	// (number of lines still to take) / (number of lines still to look at).
	// This picks each of the C(lines, take_lines) subsets with equal
	// probability, so every single line is taken with probability
	// take_lines/lines.
	seed_random();
	long long remaining_lines = lines;
	long long remaining_take  = take_lines;
	for(int i = optind; i < argc; i++) {
		if(flag_debug) {
			cerr << "Opening " << argv[i] << endl;
		}
		ifstream fin(argv[i]);
		if(!fin) {
			cerr << "Could not open '" << argv[i] << "'" << endl;
			return 1;
		}
		while(getline(fin, tmp)) {
			if(remaining_lines <= 0) break; // the input grew since the counting pass
			if(0 < remaining_take && drand48() * remaining_lines < static_cast<double>(remaining_take)) {
				cout << tmp << "\n";
				remaining_take--;
			}
			remaining_lines--;
		}
	}
	// Do not exit with a success status when the output could not be written
	// (full disk, closed pipe, ...); a truncated result would look complete.
	cout.flush();
	if(!cout) {
		cerr << "Error while writing the output" << endl;
		return 1;
	}
	return 0;
}

// Parses an integer option argument strictly (atoi() cannot tell "0" from "abc").
static bool parse_int_option(const char* arg, char optchar, long long& result)
{
	if(arg == NULL || *arg == '\0') {
		cerr << "-" << optchar << " requires a number" << endl;
		return false;
	}
	char* endp = NULL;
	errno = 0;
	const long value = strtol(arg, &endp, 10);
	if(errno != 0 || endp == arg || *endp != '\0') {
		cerr << "'" << arg << "' is not a valid number for -" << optchar << endl;
		return false;
	}
	result = static_cast<long long>(value);
	return true;
}

int main(int argc, char *argv[]) {
	int c;
	long long param_percent = -1;
	long long param_lines   = -1;
	while(1) {
		static struct option long_options[] = {
			{"debug", no_argument       , &flag_debug,   1 /*value to set*/}, // note that int_flag must be static
			//{"opt",      optional_argument,         0, 'o' /*equiv. short flag*/},
			//{"file",     required_argument,         0, 'f' /*equiv. short flag*/},
			{0, 0, 0, 0} // end of long options
		};
		int option_index = 0;
		c = getopt_long(argc, argv, "c:p:", long_options, &option_index);
		if(c == -1) break;
		switch(c) {
		case 0:
			// you can see long_options[option_index].name/flag and optarg (null if no argument).
			break;
		case 'c':
			if(!parse_int_option(optarg, 'c', param_lines)) return 1;
			if(param_lines < 0) {
				cerr << "the number of lines given to -c must not be negative" << endl;
				return 1;
			}
			break;
		case 'p':
			if(!parse_int_option(optarg, 'p', param_percent)) return 1;
			if(!(0 <= param_percent && param_percent <= 100)) {
				cerr << "percentage must be within 0 to 100" << endl;
				return 1;
			}
			break;
		case '?':
		default:
			// getopt_long() has already complained about the offending option.
			print_usage(cerr);
			return 1;
		}
	}
	if(argc < optind + 1 /* # of non-option arguments */) {
		print_usage(cerr);
		return 1;
	}
	if(param_percent != -1 && param_lines != -1) {
		cerr << "You cannot specify -p & -c at once\n" << flush;
		return 1;
	}
	return do_sieve(argc, argv, optind, param_percent, param_lines);
}

/*
=pod

=head1 NAME

sieve - extract lines randomly

=head1 SYNOPSIS

    sieve -p 20 file1 file2 ...
    sieve -c 30 file1

=head1 DESCRIPTION

B<sieve> reads the given text files (which are treated as a single
concatenated stream of lines) and writes a uniformly random sample of their
lines to the standard output, in the order in which they appear in the input.

The input files are read twice, once to count the lines and once to emit the
sample, so B<sieve> cannot read from a pipe.

=head1 OPTIONS

=over 4

=item B<-p> I<percent>

Take I<percent>% of the lines. I<percent> must be between 0 and 100.

=item B<-c> I<count>

Take exactly I<count> lines. If the input has fewer than I<count> lines, every
line is taken.

=item B<--debug>

Show progress messages on the standard error.

=back

If neither B<-p> nor B<-c> is given, half of the lines are taken.

=cut
 */
