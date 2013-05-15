
# -*- python -*-
APPNAME = 'kalab'
VERSION = '1.14'

def options(opt):
    opt.load(['compiler_c', 'compiler_cxx', 'python'])

def configure(conf):
    conf.load(['compiler_c', 'compiler_cxx', 'python'])
    conf.check_python_version((2,4,2))
    conf.env.append_unique('CXXFLAGS', ['-O2'])
    conf.env.INCLUDES += '.'
    conf.env.LIB += ['pthread', 'dl']

def build(bld):
    bld(features = 'cxx cxxprogram', source = 'src/sieve.cc', target = 'sieve')
    bld(features = 'cxx c cxxprogram', source = ['src/fatt.cc', 'src/sqlite3.c', 'src/sqdb.cc'], target = 'fatt')
    executables = ['convertsequence', 'fixshebang', 'icc-color', 'gcc-color',
                   'mydaemon', 'rep', 'sql', 'mddoc', 'sha_scan', 'gfwhich', 'json2csv',
                   'ods2xls', 'ods2xlsx']
    bld.install_files('${PREFIX}/bin', ['script/' + x for x in executables], chmod=0755)
