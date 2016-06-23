c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'angleval' -- computes bond angle for ii-jj-kk
c
c     passed variables:
c     real         x, y, z       coordinates
c     integer      ii, jj, kk    serial numbers of atoms in angle
c     real         ang           bond angle in degrees
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine angleval(x, y, z, ii, jj, kk, ang)
      
      implicit none

      include 'params.i'
      include 'constants.i'
      
      integer ii, jj, kk

      real x(MAXATOMD), y(MAXATOMD), z(MAXATOMD)
      real dxji, dyji, dzji, dxjk, dyjk, dzjk
      real dotji, dotjk, dotik
      real cosa, ang
      
      real ONE
      real CONV

      parameter (ONE = 1.0e0)
      CONV = 180.0e0 / PI

      dxji = x(jj) - x(ii)
      dyji = y(jj) - y(ii)
      dzji = z(jj) - z(ii)
      dxjk = x(jj) - x(kk)
      dyjk = y(jj) - y(kk)
      dzjk = z(jj) - z(kk)
      dotji = dxji * dxji + dyji * dyji + dzji * dzji
      dotjk = dxjk * dxjk + dyjk * dyjk + dzjk * dzjk
            
      dotik = ONE / sqrt(dotji * dotjk)
      cosa = (dxji * dxjk + dyji * dyjk + dzji * dzjk) * dotik
            
      cosa = max(-ONE ,cosa)
      cosa = min( ONE ,cosa)
      ang = acos(cosa) * CONV           

      return
      end
