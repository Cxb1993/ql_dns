# scons build file for building MRIDSS outside of xcode

import os

# Compiler flags, libraries etc.
env = Environment(CXX = '/usr/local/hdf5/bin/h5c++',HDF5_CC = 'mpic++',CCFLAGS = '-Wall -O3')
print "CXX is:", env['CXX']
env.Append(CPPPATH = ['/usr/local/include/'])
env.Append(LIBPATH = ['/usr/local/lib/'])
env.Append(LIBS = ['fftw3'])

# Source files

base = 'QL_DNS/'
models = base+'Models/'
aux = base+'Auxiliary/'
integs = base+'Integrators/'
s = '   '

basesrc = Split(base + 'main.cpp'    +s+    base + 'solution.cpp')
modsrc = Split(models+'ConstantDamping.cpp')
intsrc = Split(integs+'Euler.cpp'   +s+     integs+'EulerCN.cpp')
auxsrc = Split(aux+'Input_parameters.cpp'     +s+     aux+'MPIdata.cpp'    +s+    aux+'fftwPlans.cpp'    +s+    aux+'Kdata.cpp')

# Compile
env.Program('mridss_prog',basesrc+modsrc+intsrc+auxsrc)



