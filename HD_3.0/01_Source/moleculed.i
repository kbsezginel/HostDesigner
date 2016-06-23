c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     moleculed.i -- parameters for molecule D formed by bonding Host B
c                    and molecule C
c
c     molecule D has one rotatable bond, where we connected host B with 
c     molecule C.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c     nd           total number of atoms
c     atnumd       atomic number
c     atwtd        atomic weight
c     atlabd       atom label
c     iguestd      tag to identify guest atoms
c     xd           x-coordinates 
c     yd           y-coordinates
c     zd           z-coordinates
c     radd         van der Waals radius
c     idihedd      originally itwista
c     jdihedd      originally ibondera
c     kdihedd      originally ibonderga
c     ldihedd      originally itwistb
c     alt_ldihedd  alternate value
c     ityped       MM3 atom type
c     n12d         total number of atoms directly bonded to each atom
c     i12d         serial numbers of atoms directly bonded to each atom
c     listd        connectivity list
    
      character*2  atlabd
      
      real  xd, yd, zd, radd, atwtd
      
      integer      nd, ityped, iguestd, atnumd
      integer      idihedd, jdihedd, kdihedd, ldihedd
      integer      alt_ldihedd
      integer      n12d, i12d, listd

      common /mold1/ xd(MAXATOMD), yd(MAXATOMD), zd(MAXATOMD), 
     &         radd(MAXATOMD), atwtd(MAXATOMD)
     
      common /moldchar/ atlabd(MAXATOMD)
     
      common /mold2/ nd, atnumd(MAXATOMD),ityped(MAXATOMD), 
     &         iguestd(MAXATOMD), idihedd, jdihedd, kdihedd,
     &         ldihedd, alt_ldihedd      
                    
      common /conn12d/ n12d(MAXATOMD), i12d(MAXVAL,MAXATOMD),
     &         listd(MAXATOMD)
      
