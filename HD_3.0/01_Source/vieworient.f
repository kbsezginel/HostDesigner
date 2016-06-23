c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c    'vieworient' - reorients structures for easier viewing of the 
c    linkage that was added 
c
c   passed variables:
c     integer n         total number of atoms
c     integer n1        last atom in the hosta part
c     integer n2        last atom in the link part
c       real  x, y, z   coordinates 
c     integer n12, i12  connectivity arrays
c     logical dowhat    link or over
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine vieworient(n, n1, n2, x, y, z, n12, i12, dowhat)
  
      implicit none
      include 'params.i'

      integer n, n1, n2, i, j
      integer numat, one
      integer n12(MAXATOMD), i12(MAXVAL,MAXATOMD)

      real x(MAXATOMD), y(MAXATOMD), z(MAXATOMD)
      real xsum, ysum, zsum, torsion
      real xi, yi, zi, zero

      character*4 dowhat


c-----------------------------------------------------------------------
c If we are pushing maximum atom number, return coordinates unchanged
c-----------------------------------------------------------------------

      if(n.gt.MAXATOMD-6) then
         return
      endif

c-----------------------------------------------------------------------
c if we have the null link, then return coordinates unchanged
c-----------------------------------------------------------------------

      if(n1.eq.n2) then
         return
      endif

c-----------------------------------------------------------------------
c With LINKER places the center of coordinates for hosta, hostb,
c and the link all on the xy-plane and rotates the structure so that
c the hosta and hostb are parallel to the x-axis.  With OVERLAY,
c place the center of coordinates for link and the attachment points
c on the host all on the xy-plane and rotate structure so that the
c attachment points are parallel to the x-axis.
c-----------------------------------------------------------------------

      if(dowhat.eq.'LINK') then

c-----------------------------------------------------------------------
c get the center of dimensions for point #1 (hosta part) index n+1
c-----------------------------------------------------------------------
      
         numat = 0
         xsum = 0
         ysum = 0
         zsum = 0
         do i = 1, n1
            xsum = xsum + x(i)
            ysum = ysum + y(i)
            zsum = zsum + z(i)
            numat = numat + 1
         enddo
      
         x(n+1) = xsum/numat
         y(n+1) = ysum/numat
         z(n+1) = zsum/numat
      
c-----------------------------------------------------------------------
c get the center of dimensions for point #2 (link part) index n+2
c-----------------------------------------------------------------------

         numat = 0
         xsum = 0
         ysum = 0
         zsum = 0
         do i = n1+1, n2
            xsum = xsum + x(i)
            ysum = ysum + y(i)
            zsum = zsum + z(i)
            numat = numat + 1
         enddo
         x(n+2) = xsum/numat
         y(n+2) = ysum/numat
         z(n+2) = zsum/numat
         xi = x(n+2)
         yi = y(n+2)
         zi = z(n+2)
      
c-----------------------------------------------------------------------
c get the center of dimensions for point #3 (hostb part) index n+3
c-----------------------------------------------------------------------

         numat = 0
         xsum = 0
         ysum = 0
         zsum = 0
         do i = n2+1, n
            xsum = xsum + x(i)
            ysum = ysum + y(i)
            zsum = zsum + z(i)
            numat = numat + 1
         enddo
         x(n+3) = xsum/numat
         y(n+3) = ysum/numat
         z(n+3) = zsum/numat
   
      else

c-----------------------------------------------------------------------
c get the first host attachment atom point #1, index n+1
c get the second host attachment atom point #3, index n+3
c-----------------------------------------------------------------------

         numat = 0
         do i = n1+1, n
            do j = 1, n12(i)
               if(i12(j,i).le.n1) then
                  numat = numat + 1
                  if(numat.eq.1) then
                     x(n+1) = x(i12(j,i))
                     y(n+1) = y(i12(j,i))
                     z(n+1) = z(i12(j,i))
                  elseif(numat.eq.2) then
                     x(n+3) = x(i12(j,i))
                     y(n+3) = y(i12(j,i))
                     z(n+3) = z(i12(j,i))
                  elseif(numat.ge.3) then
                     write(6,*) 'Subroutine vieworient:'
                     write(6,*) 'Too many attachments to host'
                     stop
                  endif
               endif
            enddo
         enddo

c-----------------------------------------------------------------------
c get the center of dimensions for point #2, (link part), index n+2
c-----------------------------------------------------------------------

         numat = 0
         xsum = 0
         ysum = 0
         zsum = 0
         do i = n1+1, n
            xsum = xsum + x(i)
            ysum = ysum + y(i)
            zsum = zsum + z(i)
            numat = numat + 1
         enddo
         x(n+2) = xsum/numat
         y(n+2) = ysum/numat
         z(n+2) = zsum/numat
         xi = x(n+2)
         yi = y(n+2)
         zi = z(n+2)

      endif

c-----------------------------------------------------------------------
c mv structure so that link center is at the origin
c-----------------------------------------------------------------------

      do i = 1, n+3
         x(i) = x(i) - xi
         y(i) = y(i) - yi
         z(i) = z(i) - zi
      enddo
      
c-----------------------------------------------------------------------
c calculate the normal vector to the plane defined by the centers
c-----------------------------------------------------------------------

      x(n+4) = y(n+1)*z(n+3) - z(n+1)*y(n+3)
      y(n+4) = z(n+1)*x(n+3) - x(n+1)*z(n+3)
      z(n+4) = x(n+1)*y(n+3) - y(n+1)*x(n+3)

c-----------------------------------------------------------------------
c align the normal vector to the negative z-axis, this should put
c all three center points into the xy plane
c-----------------------------------------------------------------------

      call alignz(n+4,x,y,z,n+2,n+4,-1)

c-----------------------------------------------------------------------
c midpoint between hosta and hostb - index n+5
c-----------------------------------------------------------------------

      x(n+5) = (x(n+1) + x(n+3))/2
      y(n+5) = (y(n+1) + y(n+3))/2
      z(n+5) = (z(n+1) + z(n+3))/2

c-----------------------------------------------------------------------
c make a dummy point
c-----------------------------------------------------------------------

      x(n+6) = 0.0
      y(n+6) = -1.0
      z(n+6) = 1.0

c-----------------------------------------------------------------------
c rotate the molecule about the z-axis so that the dihedral angle 
c defined above is 0, the midpoint should now be on the negative y-axis
c-----------------------------------------------------------------------

      call torsionval(x,y,z,n+6,n+4,n+2,n+5,torsion)
      one=1
      zero=0.0e0
      call setdihedral(one, n+5, x, y, zero, torsion)
      
c-----------------------------------------------------------------------
c molecule should now be oriented with the center of the hosts and 
c and the link all in the xy plane.  The link should be in the +y
c hemisphere and the two hosts should be in the -y hemisphere, 
c parallel to the x-axis
c-----------------------------------------------------------------------

      return
      end      
