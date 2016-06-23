c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     rotateme.i -- allowed values for rotation about the bonds between
c     host A and the link, bond associated with A, and between host B 
c     molecule C, bond associated with B. 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c     nrotvala     number of dihedrals for bond A
c     rotvala      array of dihedral angle values
c     nrotvalb     number of dihedrals for bond B
c     rotvalb      array of dihedral angle values

      real  rotvala, rotvalb
      
      integer      nrotvala, nrotvalb
      

      common /rotvalue1/ rotvala(MAXTORVAL), rotvalb(MAXTORVAL)
      
      common /rotvalue2/ nrotvala, nrotvalb
     
