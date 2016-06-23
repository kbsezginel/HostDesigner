c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     moleculec.i -- parameters for molecule C formed by bonding Host A 
c                    and link
c
c     molecule C has one attachment site, originally site b on the link.
c     molecule C has one rotatable bond, where we connected host A with 
c     the link.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


c     nc           total number of atoms
c     atnumc       atomic number
c     atlabc       atom label
c     atwtc        atomic weight
c     iguestc      tag to identify guest atoms
c     xc           x-coordinates 
c     yc           y-coordinates
c     zc           z-coordinates
c     radc         van der Waals radius
c     idihedc      originally itwista
c     alt_idihedc  alternate value for idihedc
c     jdihedc      originally ibondera
c     kdihedc      originally ibonderga
c     ldihedc      originally itwistb
c     itypec       MM3 atom type
c     looserc      serial number of atom to be replaced by host B
c     ibonderc     serial number of atom that forms bond with host B
c     torrefc      serial number of atom to measure new bond dihedral from
c     alt_torrefc  alternate serial number of atom to measure new dihedral
c     itwistc      serial number of atom used to define dihedral
c     n12c         total number of atoms directly bonded to each atom
c     i12c         serial numbers of atoms directly bonded to each atom
c     listc        connectivity values for each atom, used in bump D
c     listc0       connectivity values for each atom, used in bump C
c     symlooserc   serial number of symmetrical atom, used in overlay 
c     symibonderc  serial number of symmetrical atom that forms bond with
c                  the symmetrical atom to be replaced by the link  
c     ibonderaca   updated serial number of ibondera(jj) 
c     ibonderacb   updated serial number of ibondera(jj+1)
c     torrefaca    updated serial number of torrefa(jj)
c     torrefacb    updated serial number of torrefa(jj+1)

c     ibondergca   updated serial number of ibonderga
c     ibondergcb   updated serial number of ibondergb
c     torrefgca    updated serial number of torrefga
c     torrefgcb    updated serial number of torrefgb

 
      character*2  atlabc
      
      real  xc, yc, zc, radc, atwtc
      
      integer      nc, itypec, iguestc, atnumc
      integer      idihedc, jdihedc, kdihedc, ldihedc
      integer      looserc, ibonderc, torrefc, itwistc
      integer      alt_torrefc, alt_idihedc
      integer      n12c, i12c, listc, listc0
      integer      symlooserc,symibonderc
      integer      ibonderaca, ibonderacb, torrefaca, torrefacb
      integer      ibondergca, ibondergcb,torrefgca, torrefgcb   
   

      common /molc1/ xc(MAXATOMC), yc(MAXATOMC), zc(MAXATOMC), 
     &         radc(MAXATOMC), atwtc(MAXATOMC)

      common /molcchar/ atlabc(MAXATOMC)
     
      common /molc2/ nc, atnumc(MAXATOMC), itypec(MAXATOMC), 
     &         iguestc(MAXATOMC), idihedc, jdihedc, kdihedc, ldihedc, 
     &         looserc, ibonderc, torrefc, alt_torrefc, itwistc,
     &         alt_idihedc, symlooserc(MAXSYMATTACHA, MAXATTACHA),
     &         symibonderc(MAXSYMATTACHA, MAXATTACHA),
     &         ibonderaca, ibonderacb, torrefaca, torrefacb,
     &         ibondergca, ibondergcb,torrefgca, torrefgcb
                    
      common /conn12c/ n12c(MAXATOMC), i12c(MAXVAL,MAXATOMC), 
     &       listc(MAXATOMC), listc0(MAXATOMC)
