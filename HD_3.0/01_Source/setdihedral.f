c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'setdihedral' -- sets a dihedral to specified value    
c                                                              
c Only the atoms numbered between nstart and nend are rotated. This 
c routine assumes that the bond defining the dihedral angle is 
c oriented along the z-axis, thus only the x and y coordinates are
c changed. 
c
c passed variables:
c integer nstart       pointer to first atom of rotating group
c integer nend         pointer to the last atom of rotating group
c real    x,y,z        coordinates
c real    ang          dihedral angle we want (in degrees)
c real    ph           current dihedral angle (in radians)
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine setdihedral(nstart,nend,x,y,ang,ph)
      
      implicit none
      
      include 'params.i'
      include 'constants.i'

      real x(MAXATOMD), y(MAXATOMD)
      real xi, yi, z11, z22, z21, z12 
      real CONV
      real ang, ph, dx, dy
      
      integer i, nstart, nend

      CONV = 180.0e0 / PI

c-----------------------------------------------------------------------
c Find the angle to rotate about the Z axis
c-----------------------------------------------------------------------

      ang = ang / CONV

      ph = ang - ph 

c-----------------------------------------------------------------------
c Rotate about the Z axis, Z coord stays the same
c-----------------------------------------------------------------------

      z11 = cos(ph)
      z22 = z11
      z21 = sin(ph)
      z12 = -z21

      do i = nstart, nend
         xi = x(i)
         yi = y(i)
         dx = xi * z11 + yi * z12 
         dy = xi * z21 + yi * z22 
         x(i) = dx
         y(i) = dy
      end do

      ph = ang

      return
      end
