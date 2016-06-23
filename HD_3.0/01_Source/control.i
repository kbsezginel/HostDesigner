c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     control.i -- variables read in by readcontrol
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c dowhat :  string specifying link or over
c linklib:  location of the link library, optionally with path.
c numkeep:  number of structures to save in storage array.
c numview:  number of structures to write to output
c nrota:  number of rotatable bonds in the hosta fragment
c nrotb:  number of rotatable bonds in the hostb fragment
c hosta:  name of the file describing the first fragment in linker,
c hostb:  name of the file describing the second linker fragment.
c nochiral:  excludes all chiral links. Default=false    
c noprochiral:  excludes all chiral and prochiral links. Default=false.
c noasym:  excludes all asymmetric links. Default=false
c testdrive:  print out geometry drives. Default=false
c xyz:  logical that turns on simple xyz output format. Default=false
c notype:  logical that turns off reading and writing of atoms types.  
c maxrmsd:  maximum amount rmsd may be and still be stored. Default=100 
c tightness:  real used in torsion rejections in overlay. Defaults=1.0

      character*120 linklib
      character*4  dowhat
      character*120 hosta,hostb
      character*20 archive
      logical testdrive, xyz, notype
      logical nochiral, noprochiral, noasym
      integer numkeep, numview
      real maxrmsd, tightness

      common /controlchar/ linklib,dowhat,hosta,hostb,archive
 
      common /controllg/ testdrive, xyz, nochiral, noprochiral,
     &   noasym, notype

      common /controlint/ numkeep, numview

      common /controlreal/ maxrmsd, tightness
