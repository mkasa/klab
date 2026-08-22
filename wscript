
# -*- python -*-
APPNAME = 'kalab'
VERSION = '1.57.0'

def options(opt):
    opt.load(['compiler_c', 'compiler_cxx', 'python', 'perl'])
    opt.add_option('--enable-perl', action = 'store_true', default = False, help = 'enable Perl modules')

def configure(conf):
    conf.load(['compiler_c', 'compiler_cxx', 'python', 'perl'])
    # --enable-perl is given at CONFIGURE time, but it is needed again at
    # BUILD/INSTALL time. conf.options is only populated from the command line
    # of the command being run, so the flag has to be remembered in conf.env.
    conf.env.ENABLE_PERL = bool(conf.options.enable_perl)
    if conf.env.ENABLE_PERL:
        conf.check_perl_version((5,6,0))
        conf.check_perl_ext_devel()
    conf.check_python_version((3,8,0))
    conf.check_cxx(lib='z', header_name='zlib.h', uselib_store='ZLIB', mandatory=True)
    # Apple clang (and older gcc) default to C++98 when no -std= is given,
    # which rejects the C++11 <random> code in src/sieve.cc.
    conf.env.append_unique('CXXFLAGS', ['-std=c++17', '-O2', '-DVERSION_STRING="' + VERSION + '"'])
    # src/sqlite3.c is compiled as C, so it needs its own optimization flags;
    # without these the SQLite amalgamation is built with no optimization at all.
    conf.env.append_unique('CFLAGS', ['-O2'])
    # INCLUDES is a list: appending a plain string would extend it one
    # character at a time.
    conf.env.INCLUDES += ['.']

def build(bld):
    from waflib import Utils
    bld(features = 'cxx cxxprogram', source = 'src/sieve.cc', target = 'sieve')
    # pthread/dl are only required by the SQLite amalgamation, so they are
    # linked into fatt alone rather than into every target.
    bld(features = 'cxx c cxxprogram', source = ['src/fatt.cc', 'src/sqlite3.c', 'src/sqdb.cc'], target = 'fatt', use = 'ZLIB', lib = ['pthread', 'dl'])
    executables = ['convertsequence', 'fixshebang', 'icc-color', 'gcc-color',
                   'mydaemon', 'rep', 'sha_scanp', 'gfwhich', 'json2csv', 'csv2html', 'plotr',
                   'ispcr', 'headtail', 'recompressbyxz', 'split_paf', 'reduce_genome_feature', 'cco',
                   'imgcat2', 'pbjellysummary2json', 'quastreport2json']
    bld.install_files('${PREFIX}/bin', ['script/' + x for x in executables], chmod=Utils.O755)
    if bld.env.ENABLE_PERL:
        bld.install_files('${ARCHDIR_PERL}', ['script/BLASTM8Parse.pm'])
