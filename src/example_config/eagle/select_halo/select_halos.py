"""
reads from a subfind file the location and masses of fof haloes
"""

import h5py_wrap as h5
import numpy as np
import pylab as plt
import random
import os


##################################################################
#
#
#  Charm school, just modify things in here
#  First, produce config files for ibin = 11
#  Second, produce config files for your groups mass bin 
#
ibin = 11
#
#
##################################################################



# discover dark matter particle mass 
#--------------------------------------
data_dir = '/disks/strw5/galtay/Charm-Data'
fname = data_dir + '/snapshot/snap_012_z003p017.0.hdf5'
header = h5.raa( fname, 'Header' )
dm_par_mass = header['MassTable'][1] * 1.0e10 # in Msun/h



# set group file directory 
#--------------------------------------
grp_dir_name = '/disks/strw5/galtay/Charm-Data/group'
grp_base = 'subfind_tab_012_z003p017'

fname      = grp_dir_name + '/' + grp_base + '.0.hdf5'

header = h5.raa( fname, 'Header' )
totngroups = h5.ra(fname,'Header','TotNgroups')
numfiles   = h5.ra(fname,'Header','NumFilesPerSnapshot')

r200  = np.zeros(0, dtype=np.float32)
m200  = np.zeros(0, dtype=np.float32)
coord = np.zeros( (0,3), dtype=np.float32 )


# read in r200 and m200 from group file
#--------------------------------------
for infile in range(numfiles):
    fname   = grp_dir_name + '/' + grp_base + '.' + str(infile) + '.hdf5'

    r200_l  = h5.rd(fname,'FOF/Group_R_TopHat200')
    m200_l  = h5.rd(fname,'FOF/Group_M_TopHat200')
    coord_l = h5.rd(fname,'FOF/GroupCentreOfPotential')
    
    r200  = np.concatenate( (r200, r200_l) )
    m200  = np.concatenate( (m200, m200_l) )
    coord = np.concatenate( (coord, coord_l))


# put m200 in Msun/h
#--------------------------------------
m200 = m200 * 1.0e10

# mass function 
#-----------------------------------------
m_limit = dm_par_mass * 64
igood = np.where( m200 >= m_limit )

logM200 = np.log10( m200[igood] )

fig = plt.figure( figsize=(10,10) )
(cnt, bins_edge, patchs) = plt.hist( logM200, bins=15, log=True )
plt.close('all')
bins_lo = bins_edge[0:-1]
bins_hi = bins_edge[1:]
bins_cen = (bins_lo + bins_hi ) / 2








n_max_halos_to_choose = 5

c1 = m200 > 10**bins_lo[ibin]
c2 = m200 <= 10**bins_hi[ibin]
ihalos_in_bin = np.where( c1 & c2 )[0]

nhalos = min( n_max_halos_to_choose, ihalos_in_bin.size )
print 'choosing ' + str(nhalos) + ' halos from mass bin ' + str(ibin)


# for now, just choose the most massive halos in the bin
#--------------------------------------------------------
ihalo_list = []
for i in range(nhalos):
    ihalo_list.append( ihalos_in_bin[ihalos_in_bin.size-1-i] )



for ihalo in ihalo_list:

    # report 
    #----------------------------------------------------
    print 'writing config for halo ' + str(ihalo)

    # create urchin config file
    #----------------------------------------------------
    fname = 'urchin-fof' + str(ihalo) + '.config'
    cmnd = 'cp templates/urchin-template.config ' + fname
    os.system( cmnd )

    outdir  = '../../Planck1/Urchin/output/fof' + str(ihalo)
    cmnd = 'mkdir -p ' + outdir
    os.system( cmnd )

    with open( fname, 'a') as config_file:
        option = '\nDoSelection: T \n'
        option += 'SelectRadius: ' + str(1.5*r200[ihalo]) + '\n'
        option += 'SelectXcoord: ' + str(coord[ihalo,0]) + '\n'
        option += 'SelectYcoord: ' + str(coord[ihalo,1]) + '\n'
        option += 'SelectZcoord: ' + str(coord[ihalo,2]) + '\n'
        option += 'OutputDir: ./output/fof' + str(ihalo) + '\n'
        config_file.write( option )

    cmnd = 'cp ' + fname + ' ../../Planck1/Urchin/urchin_config'
    print cmnd
    os.system( cmnd )

    cmnd = 'mv ' + fname + ' urchin_config/'
    print cmnd
    os.system( cmnd )

    # create sph2grid config files
    #----------------------------------------------------
    fname = 'grid-fof' + str(ihalo) + '-H1.config'
    cmnd = 'cp templates/grid-template.config ' + fname
    os.system( cmnd )

    urchin_fname = 'output/fof' + str(ihalo) + '/snap_999.hdf5'

    with open( fname, 'a') as config_file:
        option =  '\n[Owls]\n'
        option += 'owls_fname = ' + urchin_fname + '\n'
        option += 'owls_eos_var_name = OnEquationOfState \n'
        option += 'owls_ism_temp = 1.0e4 \n'
        option += '\n[Urchin]\n'
        option += 'urchin_fname = ' + urchin_fname + '\n'
        option += 'urchin_selection = 1 \n'
        option += '\n[Grid]\n'
        option += 'xy_len = 0.2 \n'
        option += 'z_depth = 0.2 \n'
        option += 'n_cells = 1024 \n'
        option += 'quantity = h1 \n'
        option += 'periodic = 0 \n'
        option += 'x_center = ' + str(coord[ihalo,0]) + '\n'
        option += 'y_center = ' + str(coord[ihalo,1]) + '\n'
        option += 'z_center = ' + str(coord[ihalo,2]) + '\n'
        option += '\n[Output]\n'
        option += 'output_dir = ./output/fof' + str(ihalo) + '\n'
        option += 'output_base = grid_999_H1 \n\n'
        config_file.write( option )

    cmnd = 'cp ' + fname + ' ../../Planck1/Urchin/grid_H1_config'
    print cmnd
    os.system( cmnd )

    cmnd = 'mv ' + fname + ' grid_H1_config/'
    print cmnd
    os.system( cmnd )



    fname = 'grid-fof' + str(ihalo) + '-H.config'
    cmnd = 'cp templates/grid-template.config ' + fname
    os.system( cmnd )

    urchin_fname = 'output/fof' + str(ihalo) + '/snap_999.hdf5'

    with open( fname, 'a') as config_file:
        option =  '\n[Owls]\n'
        option += 'owls_fname = ' + urchin_fname + '\n'
        option += 'owls_eos_var_name = OnEquationOfState \n'
        option += 'owls_ism_temp = 1.0e4 \n'
        option += '\n[Urchin]\n'
        option += 'urchin_fname = ' + urchin_fname + '\n'
        option += 'urchin_selection = 1 \n'
        option += '\n[Grid]\n'
        option += 'xy_len = 0.2 \n'
        option += 'z_depth = 0.2 \n'
        option += 'n_cells = 1024 \n'
        option += 'quantity = H \n'
        option += 'periodic = 0 \n'
        option += 'x_center = ' + str(coord[ihalo,0]) + '\n'
        option += 'y_center = ' + str(coord[ihalo,1]) + '\n'
        option += 'z_center = ' + str(coord[ihalo,2]) + '\n'
        option += '\n[Output]\n'
        option += 'output_dir = ./output/fof' + str(ihalo) + '\n'
        option += 'output_base = grid_999_H \n\n'
        config_file.write( option )

    cmnd = 'cp ' + fname + ' ../../Planck1/Urchin/grid_H_config'
    print cmnd
    os.system( cmnd )

    cmnd = 'mv ' + fname + ' grid_H_config/'
    print cmnd
    os.system( cmnd )

    print
