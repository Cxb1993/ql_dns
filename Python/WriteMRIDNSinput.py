#!/usr/bin/python

import os,shutil,math,sys
from WriteMRIDNS_Aux import *

MRIdir='../QL_DNS/'

aspect_ratio_list = [[1.0,2.0,1.0],[0.5,2.0,1.0],[0.25,2.0,1.0],[1.0,1.0,1.0],[0.5,1.0,1.0],[0.25,1.0,1.0],[0.5,0.5,1.0],[0.25,0.5,1.0]]

Rm = 16000.
Pm=2.
noise = 2.

for AR in aspect_ratio_list:
    filename='FullQL_Pm2Noise2Rm' + n2s(Rm) + '_L' + n2s(AR[0]) + n2s(AR[1]) + n2s(AR[2])

    class RD: # Passed to function to write input file
        equation_type='MHD_FullQlin'
        L=AR
        N=[128*AR[0], 64*AR[1], 128*AR[2]]
        q=1.5
        nu=Pm / Rm
        eta=1. / Rm
        f_noise = noise
        dt=-0.1
        CFL = 1.5
        time_interval = [0,500.]
        tv_save= 0.2
        full_save=time_interval[1] + 1.
        save_enQ=1
        save_amQ=1
        save_dissQ=1
        save_ReyQ=0
        save_mfQ=1
        init_By = -0.01
        remapQ=1
        QuasiLinearQ=1
        StartFromSavedQ=0


    DELETE_OLD_FILES = 1 # Wether or not to delete old .DSSinput files

    if DELETE_OLD_FILES: # If desired, move .DSSinput files to .OldInputFiles
        file_list=os.listdir(MRIdir)
        for fname in file_list:
            if fname.find('.DSSinput') != -1:
                shutil.move(MRIdir+fname, MRIdir+'.OldInputFiles/'+fname)
                print('Moved file ' + fname+' to .OldInputFiles/')

    write_input_file(MRIdir + filename + '.DSSinput', RD)

    print '<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>'
    print 'Done writing '+ filename + '.DSSinput'
    print '<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>'

    # Run mpiexec
    print ''
    os.system('mpiexec -np 8 ../mridns_prog '+ filename)
    print 'Done case '+filename



