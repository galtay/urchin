import os
import sys
import h5py
import h5py_wrap as h5

import numpy as np



def create_skeleton( fname0 ):

    """ Given the name of any file in an Eagle snapshot, this function 
    creates a set of skeleton output files which contain fields for gas
    particle IDs and HI number density. """ 

    if not os.path.isfile( fname0 ):
        raise IOError( fname0 + ' is not a file.' )


    cmnd = 'rm -rf example_skeleton'
    os.system( cmnd )
    cmnd = 'mkdir example_skeleton'
    os.system( cmnd )

    n_files = h5.ra( fname0, '/Header', 'NumFilesPerSnapshot' )

    for ifile in range( n_files ):

        in_file = fname0.split('.')[0] + '.' + str(ifile) + '.hdf5'
        sk_file = 'example_skeleton/out.' + str(ifile) + '.hdf5'

        if os.path.exists( sk_file ):
            raise IOError( 'output file already exists.' )

        print 'working on file: ', in_file
        print 'sk file: ', sk_file

        head = h5.raa( in_file, 'Header' )
        ngas = head['NumPart_ThisFile'][0]

        h5_in = h5py.File( in_file, 'r' )
        h5_sk = h5py.File( sk_file, 'w' )

        # use copy method to copy attribute groups
        #-------------------------------------------------------------
        groups = ['Config', 'Constants', 'HashTable', 'Header', 
                  'Parameters/ChemicalElements', 'RuntimePars', 'Units' ]

        for grp in groups:
            h5_sk.copy( h5_in[grp], h5_sk['/'] )


        # copy over gas IDs
        #-------------------------------------------------------------
        h5_sk.create_group( 'PartType0' )
        dset = 'PartType0/ParticleIDs'
        h5_sk.copy( h5_in[dset], dset )


        # make dummy arrays for Urchin particle data
        #-------------------------------------------------------------
        dum = np.ones( ngas, dtype=np.float32 ) * -1


        dset_name = 'PartType0/HydrogenOneFraction'
        h5_sk.create_dataset( dset_name, data=dum )
        attrs = {'CGSConversionFactor': 1.0,
                 'h-scale-exponent': 0.0, 
                 'aexp-scale-exponent': 0.0,
                 'VarDescription': 'Hydrogen neutral fraction = n_HI / (n_HI + n_HII)'}

        for k,v in attrs.items():        
            h5_sk[dset_name].attrs[k] = v

        h5_in.close()
        h5_sk.close()



if __name__ == '__main__':

    if len(sys.argv) != 2:
        txt = 'create_skeleton takes a single file name as its only argument.'
        raise SyntaxError( txt )
    create_skeleton( sys.argv[1] )


