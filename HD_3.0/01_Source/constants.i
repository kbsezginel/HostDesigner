c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     constants.i -- atom labels, masses, and van der Waals radii, 
c     bond lengths, dihedral angles, and tolerances
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c     xmass        atomic masses
c     vdwrad       van der Waals radii
c     elneg        electronegativity of the atom
c     covrad       covalent radius for bonding to a single bonded atom
c     atomname     periodic table label 
c     rmstol       minimum difference in rmsd values to reject a hit
c     moitol       maximum proportional difference in all individual
c                  axis'' moment of inertia to accept a hit
c     shapethresh  max (rmsd/avelen) of guest metal coordination sphere
c     metshape     Shape of guest metal coordination sphere
c     enrot        energy multiplier for rotatable bonds
c     PI           the ratio of the circumference to the radius of a circle
c     HDCONPATH    path to the constants file and default library file

c     THESE VARIABLES ARE USED IN THE BUMPCHECKS, all distances in angstroms:
c     all are stored as the square of the input value
c     scalevdw12   used to scale the sum of radii 
c     scalevdw22   used to scale guest - attached atom distances
c     HHpass2      value used to check H--H contacts 

c     THESE VARIABLES ARE USED IN OVERLAY:
c     toldist      maximum difference between host and guest distances
c     tolphi       maximum difference between host and guest dihedrals
      


      character*2  atomname  
      character*4  metshape
      character*120 HDCONPATH

      real    xmass, elneg, vdwrad, covrad
      real    toldist, tolphi
      real    rmstol, scalevdw12, scalevdw22
      real    moitol, shapethresh, HHpass2
      real    rotval, roten, enrot
      real    PI

      integer   natoms, nrotval
      
      common /tolerance/ scalevdw12,scalevdw22,HHpass2,toldist,
     &         tolphi, rmstol, moitol, shapethresh, enrot
 
      common /atprops/ xmass(MAXAT), elneg(MAXAT), vdwrad(MAXAT), 
     &         covrad(MAXAT), natoms
     
      common /atomchar/ atomname(MAXAT), metshape, HDCONPATH
     
      common /torint/ nrotval(16,4,16,4)
      
      common /torreal/ rotval(16,4,16,4,6),roten(16,4,16,4,6), PI
