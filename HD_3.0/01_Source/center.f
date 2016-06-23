c                Used with permission from Jay Ponder:
c
c          ###################################################
c          ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c          ##          see www.dasher.wustl/tinker          ##
c          ##              All Rights Reserved              ##
c          ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine center  --  superimpose structure centroids  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "center" moves the weighted centroid of each coordinate
c     set to the origin during least squares superposition
c
c
      subroutine center (n1,x1,y1,z1,n2,x2,y2,z2,xmid,ymid,zmid,
     &                   nfit,ifit)
      implicit none
      integer i,k,n1,n2
      real xmid,ymid,zmid,norm
      real x1(*),y1(*),z1(*)
      real x2(*),y2(*),z2(*) 

      integer nfit,ifit(2,*)

c-----------------------------------------------------------------------
c find the weighted centroid of the second
c structure and translate it to the origin
c-----------------------------------------------------------------------

      xmid = 0.0e0
      ymid = 0.0e0
      zmid = 0.0e0
      norm = real(nfit)
      do i = 1, nfit
         k = ifit(2,i)
         xmid = xmid + x2(k) 
         ymid = ymid + y2(k) 
         zmid = zmid + z2(k) 
      end do

      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, n2
         x2(i) = x2(i) - xmid
         y2(i) = y2(i) - ymid
         z2(i) = z2(i) - zmid
      end do

c-----------------------------------------------------------------------
c now repeat for the first structure, note
c that this centroid position gets returned
c-----------------------------------------------------------------------

      xmid = 0.0e0
      ymid = 0.0e0
      zmid = 0.0e0
      norm = real(nfit)
      do i = 1, nfit
         k = ifit(1,i)
         xmid = xmid + x1(k) 
         ymid = ymid + y1(k) 
         zmid = zmid + z1(k) 
      end do

      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, n1
         x1(i) = x1(i) - xmid
         y1(i) = y1(i) - ymid
         z1(i) = z1(i) - zmid
      end do

      return
      end

