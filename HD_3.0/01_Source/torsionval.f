c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'torsionval' -- computes dihedral angle for ii-jj-kk-ll
c     returns the dihedral angle in radians
c 
c     Multiply ang by 180.0e0/pi to get an answer ranging from 0 to 360
c     degrees
c
c     variables passed:
c     real           x, y, z
c     integer        ii, jj, kk, ll
c     real           ang
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine torsionval(x, y, z, ii, jj, kk, ll, ang)
     
      implicit none 
     
      include 'params.i'
      include 'constants.i'

      integer ii, jj, kk, ll

      real x(MAXATOMD), y(MAXATOMD), z(MAXATOMD)
      real xt12, yt12, zt12
      real xt32, yt32, zt32
      real xt34, yt34, zt34
      real cx12x32, cy12x32, cz12x32
      real cx34x32, cy34x32, cz34x32
      real denom1, denom2
      real d12x32, d34x32, sgn, dprod12
      real cosphi, phiang, ang
      real SMALL, ONE

      parameter (SMALL = 1.0e-12)    
      parameter (ONE = 1.0e0)

c-----------------------------------------------------------------------
c First find the torsion angle
c Compute the 3 vectors in the torsion angle.
c-----------------------------------------------------------------------

      xt12 = x(ii) - x(jj)
      yt12 = y(ii) - y(jj)
      zt12 = z(ii) - z(jj)
      xt32 = x(kk) - x(jj)
      yt32 = y(kk) - y(jj)
      zt32 = z(kk) - z(jj)
      xt34 = x(kk) - x(ll)
      yt34 = y(kk) - y(ll)
      zt34 = z(kk) - z(ll)

c-----------------------------------------------------------------------
c Compute the cross product of vectors 1 x 2 and 3 x 2.
c-----------------------------------------------------------------------

      cx12x32 = yt12 * zt32 - yt32 * zt12
      cy12x32 = zt12 * xt32 - zt32 * xt12
      cz12x32 = xt12 * yt32 - xt32 * yt12
      cx34x32 = yt34 * zt32 - yt32 * zt34
      cy34x32 = zt34 * xt32 - zt32 * xt34
      cz34x32 = xt34 * yt32 - xt32 * yt34

c-----------------------------------------------------------------------
c Save some vector lengths and inner products.
c-----------------------------------------------------------------------

      denom1 = max(small,cx12x32 * cx12x32 +
     &                  cy12x32 * cy12x32 + cz12x32 * cz12x32)
      denom2 = max(small,cx34x32 * cx34x32 +
     &                  cy34x32 * cy34x32 + cz34x32 * cz34x32)
      d12x32 = ONE / denom1
      d34x32 = ONE / denom2

c-----------------------------------------------------------------------
c Save the sign of the angle.
c-----------------------------------------------------------------------

      sgn = - xt12 * cx34x32 - yt12 * cy34x32 - zt12 * cz34x32

c-----------------------------------------------------------------------
c Compute the torsional angle.
c-----------------------------------------------------------------------

      dprod12  = cx12x32 * cx34x32 +
     &           cy12x32 * cy34x32 +
     &           cz12x32 * cz34x32
      cosphi   = dprod12 * sqrt(d12x32 * d34x32)
      cosphi   = max(-ONE , cosphi)
      cosphi   = min( ONE , cosphi)
      phiang   = acos(cosphi)

c-----------------------------------------------------------------------
c Adjust phi so that -180 < phi < 180, ph is the currect torsion angle
c-----------------------------------------------------------------------

      ang = PI - sign(phiang,sgn)

      return
      end
