#!/bin/env python

import read_eagle

infile = '/gpfs/data/Eagle/mfTestRuns/EAGLE_RUNS/L0100N1440/ANARCHY_CSP_048/COOL_OWLS_LIM_1/AGN_HALO_11p60_SNe_T0_5p5_x2Fe_x2SNIA/data/snipshot_065_z009p266/snip_065_z009p266.11.hdf5'

snap = read_eagle.EagleSnapshot("./test_data/snapshot_020/snap_020.0.hdf5")

snap.select_region(4.0, 5.0, 2.0, 3.0, 3.0, 4.0)

pos = snap.read_dataset(0, "Coordinates")
ids = snap.read_dataset(0, "ParticleIDs")

snap.close()

for (id, x, y, z) in zip(ids, pos[:,0], pos[:,1], pos[:,2]):
    print id, x, y, z
