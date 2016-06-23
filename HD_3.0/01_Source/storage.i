c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     storage.i -- molecular information for hits stored in ram
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c     stonatom     total number of atoms
c     stoatlab     atom label
c     stoguest     tag to identify guest atoms
c     stox         x-coordinates 
c     stoy         y-coordinates
c     stoz         z-coordinates
c     stotype      MM3 atom type
c     ston12       total number of atoms directly bonded to each atom
c     stoi12       serial numbers of atoms directly bonded to each atom
c     stormsd      rmsd values
c     pointer      pointer
c     stoix        first moment of inertia
c     stoiy        second moment of inertia
c     stoiz        third moment of inertia
c     stoen        energy estimate (sum of link and bond formed energies)
c     storotab     sum of number of rotatable bonds in hosta and hostb

      character*2  stoatlab

      real stox, stoy, stoz, stormsd, stoix, stoiy, stoiz, stoen
      
      integer      stonatom, stotype, stoguest
      integer      ston12, stoi12, pointer, stoii, storotab
      

      common /stomol/ stox(MAXHITS,MAXATOMD), stoy(MAXHITS,MAXATOMD), 
     &        stoz(MAXHITS,MAXATOMD), stormsd(MAXHITS), stoen(MAXHITS),
     &        stoix(MAXHITS),stoiy(MAXHITS),stoiz(MAXHITS)
     
      common /stoatlab/ stoatlab(MAXHITS,MAXATOMD)
    
      common /stomol2/ stonatom(MAXHITS), stotype(MAXHITS,MAXATOMD),
     &         stoguest(MAXHITS,MAXATOMD), pointer(MAXHITS), 
     &         stoii(MAXHITS), storotab(MAXHITS)
                    
      common /stoconn12/ ston12(MAXHITS,MAXATOMD), 
     &         stoi12(MAXHITS,MAXVAL,MAXATOMD)
      
