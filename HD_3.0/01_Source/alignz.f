c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c
c   The terminology and variable names come from Schaum's Outlines
c   Computer Graphics by Plastock and Kalley, p121.  Written by John
c   Nicholas 5/11/00.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'alignz' -- reorients a molecule such that atom(ii) is moved to 
c     the origin and atom(jj) is aligned with either the positive 
c     (idir=1) or negative (idir = -1) z axis.
c
c     passed variables:
c     integer n         number of atoms
c     real    x,y,z     coordinates
c     integer ii        atom to be set at the origin
c     integer jj        atom to be aligned with z axis
c     integer idir      +1 = positive z-axis, -1 = negative z-axis
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine alignz(n, x, y, z, ii, jj, idir)

      implicit none

      real x(*), y(*), z(*)
      real vmag, a, b, c, ab, ac, gamma
      real a11, a12, a13
      real a21, a22, a23
      real a31, a32, a33
      real xi, yi, zi, dx, dy, dz
      real ZERO
      real SMALL

      integer i, ii, jj, n, idir, count

      parameter (ZERO = 0.0e0)
      parameter (SMALL = 1.0e-8)

      if (abs(idir) .ne. 1) then
         write (6,'('' Error: Can only call this routine '',
     &      ''with idir = +/- 1'')')
         return
      end if

c-----------------------------------------------------------------------
c Translate x(ii) to the origin
c-----------------------------------------------------------------------

      xi = x(ii)
      yi = y(ii)
      zi = z(ii)

      do i = 1, n
         x(i) = x(i) - xi
         y(i) = y(i) - yi
         z(i) = z(i) - zi
      end do

c-----------------------------------------------------------------------
c Set up rotation matrix to put atom jj on the z axis
c Positive and negative z axis is controlled by multiplying the
c components of the vector by idir (+/- 1.0)
c
c If gamma is lt small, it means that the vector between ii and jj is
c on the x-axis (both y(jj) and z(jj) are essentially zero).  Will result
c on the x-axis (both y and z are essentially zero).  This will result
c in division by zero in matrix elements.  It is fixed by slightly
c nudging the y and z coordinates of jj (by 0.0001 angstrom)
c-----------------------------------------------------------------------

    5 continue

      a = x(jj) * real(idir)
      b = y(jj) * real(idir)
      c = z(jj) * real(idir)
      ab = a * b
      ac = a * c
      gamma = sqrt(b * b + c * c)
      vmag = sqrt(a * a + b * b + c * c)

      if (gamma .lt. SMALL) then
         y(jj) = y(jj) + 0.0001
         z(jj) = z(jj) + 0.0001
         count = count + 1
         if(count.gt.10) then
            print *,'Error in alignz'
            return
         end if
         goto 5
      end if

c-----------------------------------------------------------------------
c Computer elements of the rotation matrix
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
c Do the rotation
c-----------------------------------------------------------------------

      do i = 1, n
         xi = x(i)
         yi = y(i)
         zi = z(i)
         dx = xi * a11 + yi * a12 + zi * a13
         dy = xi * a21 + yi * a22 + zi * a23
         dz = xi * a31 + yi * a32 + zi * a33
         x(i) = dx
         y(i) = dy
         z(i) = dz
      end do

      return
      end
