c                Used with permission from Jay Ponder:
c
c          ###################################################
c          ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c          ##          see www.dasher.wustl/tinker          ##
c          ##              All Rights Reserved              ##
c          ###################################################
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'inertia' -- computes the principal moments of inertia
c
c     passed variables:
c     integer  n           number of atoms
c     real     x, y, z     coordinates
c     integer  atwt        atom atomic weights
c     integer  ix,iy,iz    moments of inertia        
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine inertia (n, x, y, z, atwt, ix, iy, iz)
     
      implicit none

      include 'params.i'

      integer i, n
      
      real x(MAXATOMD), y(MAXATOMD), z(MAXATOMD)
      real atwt(MAXATOMD)
      real total, ix, iy, iz
      real xcm, ycm, zcm
      real xx, xy, xz, yy, yz, zz
      real xterm, yterm, zterm
      real moment(3), vec(3,3), work1(3), work2(3)
      real tensor(3,3)

c-----------------------------------------------------------------------
c     compute the position of the center of mass
c-----------------------------------------------------------------------

      total = 0.0e0
      xcm = 0.0e0
      ycm = 0.0e0
      zcm = 0.0e0
      
      do i = 1, n
         total = total + atwt(i)
         xcm = xcm + x(i) * atwt(i)
         ycm = ycm + y(i) * atwt(i)
         zcm = zcm + z(i) * atwt(i)
      end do
      
      xcm = xcm / total
      ycm = ycm / total
      zcm = zcm / total
      
c-----------------------------------------------------------------------
c     compute and then diagonalize the inertial tensor
c-----------------------------------------------------------------------

      xx = 0.0e0
      xy = 0.0e0
      xz = 0.0e0
      yy = 0.0e0
      yz = 0.0e0
      zz = 0.0e0
      
      do i = 1, n
      
         xterm = x(i) - xcm
         yterm = y(i) - ycm
         zterm = z(i) - zcm
         xx = xx + xterm * xterm * atwt(i)
         xy = xy + xterm * yterm * atwt(i)
         xz = xz + xterm * zterm * atwt(i)
         yy = yy + yterm * yterm * atwt(i)
         yz = yz + yterm * zterm * atwt(i)
         zz = zz + zterm * zterm * atwt(i)
         
      end do
      
      
      tensor(1,1) = yy + zz
      tensor(2,1) = -xy
      tensor(3,1) = -xz
      tensor(1,2) = -xy
      tensor(2,2) = xx + zz
      tensor(3,2) = -yz
      tensor(1,3) = -xz
      tensor(2,3) = -yz
      tensor(3,3) = xx + yy
      
      
      call jacobi (3,3,tensor,moment,vec,work1,work2)
      
      
      ix = moment(1)
      iy = moment(2)
      iz = moment(3)
      
      return
      end
