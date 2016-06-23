c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'getrmsd' -- reports best rmsd for overlap of guest on host A with
c     guest on host B
c
c     passed variables:
c     integer   start1     serial number for start of first guest
c     integer   n          number of atoms
c     real      x,y,z      coordinates
c     real      rmsd       root mean square deviation for overlap
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine getrmsd(start1,n,x,y,z,rmsd)
      
      implicit none

      include 'params.i'
      include 'mguest.i'

      real dx,dy,dz,rmsd
      real x(MAXATOMD),y(MAXATOMD),z(MAXATOMD), tally     
      
      integer i, j, n, start1, start2
      
      start2 = n-lguest
      rmsd=1.0e10
      
      do i=1,nequiva
        tally=0.0e0
        do j=1,lguest
          dx = x(start1+j) - x(start2+mapeqa(i,j))
          dy = y(start1+j) - y(start2+mapeqa(i,j))
          dz = z(start1+j) - z(start2+mapeqa(i,j))
          tally = tally + (dx * dx + dy * dy + dz * dz)
        enddo
        rmsd=min(tally,rmsd)
      enddo

      rmsd=sqrt(rmsd/lguest)

      return
      end
