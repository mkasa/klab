
# -*- python -*-
APPNAME = 'kalab'
VERSION = '1.37'

def options(opt):
    opt.load(['compiler_c', 'compiler_cxx', 'python', 'perl'])
    opt.add_option('--enable-perl', action = 'store_true', default = False, help = 'enable Perl modules')

def configure(conf):
    conf.load(['compiler_c', 'compiler_cxx', 'python', 'perl'])
    if conf.options.enable_perl:
        conf.check_perl_version((5,6,0))
        conf.check_perl_ext_devel()
    conf.check_python_version((2,4,2))
    conf.env.append_unique('CXXFLAGS', ['-O2', '-DVERSION_STRING=' + VERSION])
    conf.env.INCLUDES += '.'
    conf.env.LIB += ['pthread', 'dl']

def build(bld):
    from waflib import Utils
    bld(features = 'cxx cxxprogram', source = 'src/sieve.cc', target = 'sieve')
    bld(features = 'cxx c cxxprogram', source = ['src/fatt.cc', 'src/sqlite3.c', 'src/sqdb.cc'], target = 'fatt')
    executables = ['convertsequence', 'fixshebang', 'icc-color', 'gcc-color',
                  'mydaemon', 'rep', 'sq', 'mddoc', 'gmddoc', 'sha_scan', 'gfwhich', 'json2csv', 'csv2html', 'csv2md',
                   'ods2xls', 'ods2xlsx', 'pbjellysummary2json', 'ispcr', 'headtail', 'recompressbyxz', 'taw', 'quastreport2json']
    # bld.install_files('${PREFIX}/bin', ['script/' + x for x in executables], chmod=0755)
    bld.install_files('${PREFIX}/bin', ['script/' + x for x in executables], chmod=Utils.O755)
    if bld.options.enable_perl:
        bld.install_files('${ARCHDIR_PERL}', ['script/BLASTM8Parse.pm'])
