
# -*- python -*-
APPNAME = 'kalab'
VERSION = '1.00'

def options(opt):
    opt.load('compiler_cxx python')

def configure(conf):
    conf.load('compiler_cxx python')
    conf.check_python_version((2,4,2))
    conf.env.append_unique('CXXFLAGS', ['-O2'])
    conf.env.INCLUDES += '.'
    conf.env.LIB += ['pthread']

def build(bld):
    bld(features = 'cxx cprogram', source = 'fatt.cc', target = 'fatt')
    executables = ['convertsequence', 'fixshebang', 'icc-color']
    bld.install_files('${PREFIX}/scripts', ['script/' + x for x in executables], chmod=0755)
