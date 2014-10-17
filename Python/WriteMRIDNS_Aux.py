def n2s( num ): # num2str replacing decimals
    if abs(num - round(num)) < 1e-6:
        num = int(num)
    outstr = str(num)
    outstr=outstr.replace('.','')
    return outstr


def write_input_file(filename, RD):
    file = open(filename,'w')
    # Write contents
    file.write('Input file for QL_DNS - Generated automatically by WriteMRIDNSinput.py\n\n')
    file.write('Format: variable_ = var (can have arbitrary number of spaces around = )\n\n')
    file.write('// Simulation type \n')
    
    file.write('equations_to_use_ = %s \n' % RD.equation_type)
    file.write('// Box parameters\n')
    file.write('lx_ = %15.15f, ly_ = %15.15f, lz_ = %15.15f \n'%(RD.L[0] ,RD.L[1] ,RD.L[2]))
    file.write('nx_ = %d,  ny_ = %d,  nz_ = %d \n\n' %(RD.N[0], RD.N[1], RD.N[2]))
    file.write('// Time parameters\n')
    file.write('dt_ = %15.15f\n' %RD.dt)
    file.write('CFL_ = %15.15f\n' %RD.CFL)
    file.write('t_initial_ = %15.15f,   t_final_ = %15.15f,\n\n'%(RD.time_interval[0],RD.time_interval[1]))
    file.write('// Saving\n')
    file.write('tv_save_ = %15.15f,   // Save energy, momentum etc. time-step\n'%RD.tv_save)
    file.write('full_save_ = %15.15f,   // Save full solution\n'%RD.full_save)
    file.write('save_energy? =%d  save_angular_mom? =%d     save_dissipation?=%d, // On-off flags\n'%(RD.save_enQ,RD.save_amQ,RD.save_dissQ))
    file.write('save_reynolds_stress?= %d \n'%RD.save_ReyQ)
    file.write('mean_field_save? = %d  // Save mean fields\n\n'%RD.save_mfQ)
    file.write('// General physical parameters\n')
    file.write('q_ = %15.15f \n'%RD.q)
    file.write('nu_ = %15.15f,  eta_ = %15.15f,   // Dissipation parameters \n'%(RD.nu,RD.eta))
    file.write('f_noise_ = %15.15f \n\n'%RD.f_noise)
    file.write('// Initial conditions\n')
    file.write('initial_By_ = %15.15f \n'%RD.init_By)
    file.write('// Flags for simulation\n')
    file.write('Remap? = %d \n'%RD.remapQ)
    file.write('QuasiLinear? = %d \n\n'%RD.QuasiLinearQ)
    file.write('// Restart simulation\n')
    file.write('start_from_saved? = %d \n'%RD.StartFromSavedQ)
    
    
    file.close()

