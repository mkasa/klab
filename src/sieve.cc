
// -*- mode:C++; c-basic-offset:2; tab-width:2 -*-
// @author Masahiro Kasahara (masahiro@kasahara.ws)
//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <set>
#include <getopt.h>

using namespace std;

static int flag_debug = 0;

void do_sieve(int argc, char *argv[], int optind, int param_percent, int param_lines)
{
	int lines = 0;
	string tmp;
	// count lines
	for(int i = optind; i < argc; i++) {
		if(flag_debug) {
			cerr << "Opening " << argv[i] << endl;
		}
		fstream fin(argv[i]);
		if(!fin) {
			cerr << "Could not open '" << argv[i] << "'" << endl;
			exit(1);
		}
		while(getline(fin, tmp)) lines++;
		if(flag_debug) {
			cerr << lines << " lines so far" << endl;
		}
	}
	int take_lines = lines / 2;
	if(param_percent != -1) { take_lines = lines * param_percent / 100; }
	if(param_lines   != -1) { take_lines = param_lines; if(lines < take_lines) take_lines = lines; }
	if(flag_debug) {
		cerr << "We will take " << take_lines << " lines" << endl;
		cerr << "Shuffling..." << endl;
	}
	// create set
	set<int> is_this_line_will_taken;
	for(int i = 0; i < take_lines; i++) is_this_line_will_taken.insert(i);
	srand48(time(NULL));
	for(int i = 0; i < lines; i++) {
		if(!is_this_line_will_taken.count(i)) continue;
		const int target = drand48() * lines;
		if(is_this_line_will_taken.count(target)) continue;
		is_this_line_will_taken.insert(target);
		is_this_line_will_taken.erase(i);
	}
	if(flag_debug) {
		cerr << "output" << endl;
	}
	int cursor = 0;
	for(int i = optind; i < argc; i++) {
		if(flag_debug) {
			cerr << "Opening " << argv[i] << endl;
		}
		fstream fin(argv[i]);
		if(!fin) {
			cerr << "Could not open '" << argv[i] << "'" << endl;
			exit(1);
		}
		while(getline(fin, tmp)) {
			if(is_this_line_will_taken.count(cursor))
				cout << tmp << "\n";
			cursor++;
		}
		if(flag_debug) {
			cerr << lines << " lines so far" << endl;
		}
	}
}

int main(int argc, char *argv[]) {
	int c;
	int param_percent = -1;
	int param_lines   = -1;
	while(1) {
		static struct option long_options[] = {
			{"debug", no_argument      , &flag_debug, 1   /*value to set*/}, // note that int_flag must be static
			//{"sflag",    no_argument      ,         0, 's' /*equiv. short flag*/},
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
			param_lines = atoi(optarg);
			break;
		case 'p':
			param_percent = atoi(optarg);
			if(!(0 <= param_percent && param_percent <= 100)) {
				cerr << "percentage must be within 0 to 100" << endl;
				return 1;
			}
			break;
		}
	}
	if(argc < optind + 1 /* # of non-option arguments */) {
		cerr << "usage: sieve -p 30 <file(s)>     Take 30% of the lines in the file.\n";
		cerr << "       sieve -c 40 <file(s)>     Take 40 lines of the lines in the file.\n" << flush;
		return 1;
	}
	if(param_percent != -1 && param_lines != -1) {
		cerr << "You cannot specify -p & -c at once\n" << flush;
		return 1;
	}
	do_sieve(argc, argv, optind, param_percent, param_lines);
}

/*
=pod

=head1 NAME

sieve - extract lines randomly

=head1 SYNOPSIS

    sieve -p 20 file1 file2 ...
    sieve -c 30 file1

=head1 DESCRIPTION


=cut
 */
