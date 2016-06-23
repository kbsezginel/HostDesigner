c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      params.i  --  fixed parameters
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c     MAXVAL       maximum valence

c     MAXATOMA     maximum number of atoms in host A
c     MAXATTACHA   maximum number of host A attachment sites
c     MAXDRIVEA    maximum number of host A drives

c     MAXATOMB     maximum number of atoms in host B
c     MAXATTACHB   maximum number of host B attachment sites
c     MAXDRIVEB    maximum number of host B drives

c     MAXATOMG     maximum number of atoms in link
c     MAXLINKS     maximum number of links

c     the molecule c is made by merging host A with link
c     MAXATOMC     maximum number of atoms in c

c     the molecule d is made by merging c with host B
c     MAXATOMD     maximum number of atoms in d

c     MAXGST       maximum number of atoms in guest

c     MAXAT        maximum number of atoms in CONSTANTS
c     MAXTORVAL    maximum number of possible torsions per bond 

c     MAXHITS      maximum number of D molecules to be stored in ram 

c     MAXSYMATTACHA  maximum number symmetrical attachment sites


      integer      MAXVAL
      integer      MAXATOMA, MAXATTACHA, MAXDRIVEA
      integer      MAXATOMB, MAXATTACHB, MAXDRIVEB
      integer      MAXATOMG, MAXLINKS
      integer      MAXATOMC, MAXATOMD
      integer      MAXGST
      integer      MAXAT, MAXTORVAL
      integer      MAXSYMATTACHA      
      integer      MAXHITS, MAXRCP

      parameter (MAXVAL = 10)

      parameter (MAXATOMA = 200)
      parameter (MAXATTACHA = 20)
      parameter (MAXDRIVEA = 2000)
      
      parameter (MAXATOMB = 200)
      parameter (MAXATTACHB = 10)
      parameter (MAXDRIVEB = 200)

      parameter (MAXATOMG = 35)
      parameter (MAXLINKS = 15000)

      parameter (MAXSYMATTACHA = 5)

      parameter (MAXATOMC=MAXATOMA+(MAXSYMATTACHA+1)*MAXATOMG)
      
      parameter (MAXATOMD = MAXATOMC + MAXATOMB)
     
      parameter (MAXGST = 50)

      parameter (MAXAT = 150)
      parameter (MAXTORVAL = 6)
      
      parameter (MAXRCP = 8*MAXTORVAL*MAXTORVAL*MAXATTACHB)

      parameter (MAXHITS = 1000)
