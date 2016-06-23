c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c
c   The terminology and variable names come from Schaum's Outlines
c   Computer Graphics by Plastock and Kalley, p121.  Adapted from
c   original alignz.f written by John Nicholas 5/11/00.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'alignzpass' -- reorients a molecule such that atom(ii) is moved
c     to the origin and atom(jj) is aligned with either the positive
c     (idir = 1) or negative (idir = -1) z axis.
c
c     passed variables:
c     integer n         number of atoms
c     real    x,y,z     input coordinates
c     real    xo,yo,zo  output coordinates
c     integer ii        atom to be set at the origin
c     integer jj        terminal atom
c     integer idir      +1 = positive z-axis, -1 = negative z-axis
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine alignzpass(n, x, y, z, ii, jj, idir, xo, yo, zo)
      
      implicit none
      
      include 'params.i'

      real  x(MAXATOMA), y(MAXATOMA), z(MAXATOMA)
      real xo(MAXATOMA),yo(MAXATOMA),zo(MAXATOMA)
      real vmag, a, b, c, ab, ac, gamma
      real a11, a12, a13
      real a21, a22, a23
      real a31, a32, a33
      real xi, yi, zi
      real ZERO
      real SMALL

      integer i, ii, jj, n, idir, count

      parameter (ZERO = 0.0e0)
      parameter (SMALL = 1.0e-8)    

      if (abs(idir) .ne. 1) then
         write (6,'('' Error: Can only call this routine ''
     &      ''with idir = +/- 1'')')
         return  
      end if

c-----------------------------------------------------------------------
c Remember x(ii), as it will be the origin in the output xo,yo,zo:
c-----------------------------------------------------------------------

      xi = x(ii)
      yi = y(ii)
      zi = z(ii)

c-----------------------------------------------------------------------
c Assume the atom jj is the terminal atom bonded to atom ii
c Set up rotation matrix to put atom jj on the z axis 
c Positive and negative z axis is controlled by multiplying the 
c components of the vector by idir (+/- 1.0)
c All coordinates are with respect to the future origin (ii)
c
c If gamma is lt small, it means that the vector between ii and jj is
c on the x-axis (both y(jj) and z(jj) are essentially zero).  Will result
c in division by zero in matrix elements.  It is fixed by slightly
c nudging the y and z coordinates of jj (by 0.0001 angstrom)
c-----------------------------------------------------------------------

    5 continue

      a = (x(jj)-xi) * real(idir)
      b = (y(jj)-yi) * real(idir)
      c = (z(jj)-zi) * real(idir)
      ab = a * b
      ac = a * c
      gamma = sqrt(b * b + c * c)
      vmag = sqrt(a * a + b * b + c * c)

      if (gamma .lt. SMALL) then
         y(jj) = y(jj) + 0.0001 
         z(jj) = z(jj) + 0.0001
         count = count + 1
         if(count.gt.10) then
            print *,'Error in alignzpass'
            return
         end if
         goto 5
      end if

c-----------------------------------------------------------------------
c Compute elements of the rotation matrix
c-----------------------------------------------------------------------

      a11 = gamma / vmag
      a12 = -ab / (gamma * vmag)
      a13 = -ac / (gamma * vmag)
      a21 = ZERO 
      a22 = c / gamma
      a23 = -b / gamma
      a31 = a / vmag
      a32 = b / vmag
      a33 = c / vmag

c-----------------------------------------------------------------------
c Do the rotation and translation
c-----------------------------------------------------------------------

      do i = 1, n
         xo(i) = (x(i)-xi) * a11 + (y(i)-yi) * a12 + (z(i)-zi) * a13
         yo(i) = (x(i)-xi) * a21 + (y(i)-yi) * a22 + (z(i)-zi) * a23
         zo(i) = (x(i)-xi) * a31 + (y(i)-yi) * a32 + (z(i)-zi) * a33
      end do

      return
      end
