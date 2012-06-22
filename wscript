
# -*- python -*-
APPNAME = 'kalab'
VERSION = '1.01'

def options(opt):
    opt.load('compiler_cxx python')

def configure(conf):
    conf.load('compiler_cxx python')
    conf.check_python_version((2,4,2))
    conf.env.append_unique('CXXFLAGS', ['-O2'])
    conf.env.INCLUDES += '.'
    conf.env.LIB += ['pthread']

def build(bld):
    bld(features = 'cxx cprogram', source = 'src/fatt.cc', target = 'fatt')
    bld(features = 'cxx cprogram', source = 'src/sieve.cc', target = 'sieve')
    executables = ['convertsequence', 'fixshebang', 'icc-color', 'mydaemon', 'rep', 'sql']
    bld.install_files('${PREFIX}/bin', ['script/' + x for x in executables], chmod=0755)
