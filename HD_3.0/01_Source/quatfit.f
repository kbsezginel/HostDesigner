c                Used with permission from Jay Ponder:
c
c          ###################################################
c          ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c          ##          see www.dasher.wustl/tinker          ##
c          ##              All Rights Reserved              ##
c          ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine quatfit  --  quaternion superposition of coords  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "quatfit" uses a quaternion-based method to achieve the best
c     fit superposition of two sets of coordinates
c
c     literature reference:
c
c     S. J. Kearsley, "An Algorithm for the Simultaneous Superposition
c     of a Structural Series", Journal of Computational Chemistry,
c     11, 1187-1192 (1990)
c
c     adapted from an original program written by David J. Heisterberg,
c     Ohio Supercomputer Center, Columbus, OH
c
      subroutine quatfit (x1,y1,z1,n2,x2,y2,z2,nfit,ifit)

      implicit none

      integer i,i1,i2,n2
      real xrot,yrot,zrot
      real xxyx,xxyy,xxyz,xyyx,xyyy
      real xyyz,xzyx,xzyy,xzyz
      real rot(3,3),temp1(4),temp2(4)
      real q(4),d(4),c(4,4),v(4,4)
      real x1(*),y1(*),z1(*)
      real x2(*),y2(*),z2(*)

      integer nfit,ifit(2,*)

c-----------------------------------------------------------------------
c     build the upper triangle of the quadratic form matrix
c-----------------------------------------------------------------------

      xxyx = 0.0e0
      xxyy = 0.0e0
      xxyz = 0.0e0
      xyyx = 0.0e0
      xyyy = 0.0e0
      xyyz = 0.0e0
      xzyx = 0.0e0
      xzyy = 0.0e0
      xzyz = 0.0e0
      do i = 1, nfit
         i1 = ifit(1,i)
         i2 = ifit(2,i)
         xxyx = xxyx + x1(i1) * x2(i2)
         xxyy = xxyy + y1(i1) * x2(i2)
         xxyz = xxyz + z1(i1) * x2(i2)
         xyyx = xyyx + x1(i1) * y2(i2)
         xyyy = xyyy + y1(i1) * y2(i2)
         xyyz = xyyz + z1(i1) * y2(i2)
         xzyx = xzyx + x1(i1) * z2(i2)
         xzyy = xzyy + y1(i1) * z2(i2)
         xzyz = xzyz + z1(i1) * z2(i2)
      end do
      c(1,1) = xxyx + xyyy + xzyz
      c(1,2) = xzyy - xyyz
      c(2,2) = xxyx - xyyy - xzyz
      c(1,3) = xxyz - xzyx
      c(2,3) = xxyy + xyyx
      c(3,3) = xyyy - xzyz - xxyx
      c(1,4) = xyyx - xxyy
      c(2,4) = xzyx + xxyz
      c(3,4) = xyyz + xzyy
      c(4,4) = xzyz - xxyx - xyyy

c-----------------------------------------------------------------------
c     diagonalize the quadratic form matrix
c-----------------------------------------------------------------------

      call jacobi (4,4,c,d,v,temp1,temp2)

c-----------------------------------------------------------------------
c     extract the desired quaternion
c-----------------------------------------------------------------------

      q(1) = v(1,4)
      q(2) = v(2,4)
      q(3) = v(3,4)
      q(4) = v(4,4)

c-----------------------------------------------------------------------
c     assemble rotation matrix that superimposes the molecules
c-----------------------------------------------------------------------

      rot(1,1) = q(1)*q(1) + q(2)*q(2) - q(3)*q(3) - q(4)*q(4)
      rot(2,1) = 2.0e0 * (q(2) * q(3) - q(1) * q(4))
      rot(3,1) = 2.0e0 * (q(2) * q(4) + q(1) * q(3))
      rot(1,2) = 2.0e0 * (q(3) * q(2) + q(1) * q(4))
      rot(2,2) = q(1)*q(1) - q(2)*q(2) + q(3)*q(3) - q(4)*q(4)
      rot(3,2) = 2.0e0 * (q(3) * q(4) - q(1) * q(2))
      rot(1,3) = 2.0e0 * (q(4) * q(2) - q(1) * q(3))
      rot(2,3) = 2.0e0 * (q(4) * q(3) + q(1) * q(2))
      rot(3,3) = q(1)*q(1) - q(2)*q(2) - q(3)*q(3) + q(4)*q(4)

c-----------------------------------------------------------------------
c     rotate second molecule to best fit with first molecule
c-----------------------------------------------------------------------

      do i = 1, n2
         xrot = x2(i) * rot(1,1) + y2(i) * rot(1,2) + z2(i) * rot(1,3)
         yrot = x2(i) * rot(2,1) + y2(i) * rot(2,2) + z2(i) * rot(2,3)
         zrot = x2(i) * rot(3,1) + y2(i) * rot(3,2) + z2(i) * rot(3,3)
         x2(i) = xrot
         y2(i) = yrot
         z2(i) = zrot
      end do

      return
      end
