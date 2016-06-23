c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c  'drive' -- moves atoms a single step of a drive
c
c  passed variables:
c  character  dtype      type of drive (B=Bond,A=Angle,D=Dihedral)
c  character  dstep      length of drive step (in Angstroms or radians)
c  character  dwhich     array specifying which atoms are to be moved
c  integer    n          number of atoms
c  real       x,y,z      coords
c  integer    di,dj,dk   serial numbers of atoms in drive
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine drive(dtype,dstep,dwhich,n,x,y,z,di,dj,dk)
      
      implicit none
      include 'params.i'
      
      real dstep, x(MAXATOMA),y(MAXATOMA),z(MAXATOMA),xj0,yj0,zj0,signr
      real a,b,c,gamma,vl,xi,yi,zi,dx,dy,dz
      real a11,a12,a13,a21,a22,a23,a31,a32,a33,cost,sint,msint

      integer n,di,dj,dk,i

      character*1 dtype

      logical dwhich(MAXATOMA)
      
c-----------------------------------------------------------------------
c     distance driver
c-----------------------------------------------------------------------

      if(dtype.eq.'B') then

        call alignz(n,x,y,z,dj,di,-1)
        if(dwhich(di)) then 
          do i=1,n
            if(dwhich(i)) z(i)=z(i)-dstep
          enddo
        else
          do i=1,n
            if(dwhich(i)) z(i)=z(i)+dstep
          enddo
        endif

c-----------------------------------------------------------------------
c     angle driver
c-----------------------------------------------------------------------
      
      elseif(dtype.eq.'A') then

c-----------------------------------------------------------------------
c     translate dj to the origin
c-----------------------------------------------------------------------

        xj0=x(dj)
        yj0=y(dj)
        zj0=z(dj)
        do i=1,n
          x(i)=x(i)-xj0
          y(i)=y(i)-yj0
          z(i)=z(i)-zj0
        enddo

c-----------------------------------------------------------------------
c     find the cross product of di-dj and dj-dk vectors
c     these vectors are equal to the di and dk coordinates with dj at 0
c-----------------------------------------------------------------------

        a = y(di) * z(dk) - y(dk) * z(di)
        b = z(di) * x(dk) - z(dk) * x(di)
        c = x(di) * y(dk) - x(dk) * y(di)

c-----------------------------------------------------------------------
c     now align the molecule with that cross-product vector
c-----------------------------------------------------------------------

        gamma = sqrt(b*b+c*c)
        vl = sqrt(a*a+b*b+c*c)
        if(vl.lt.(1.0e-5)) then
          write(6,*)'Angle drive failed- atoms are linear!'
        elseif(gamma.gt.(1.0e-5)) then
          a11 = gamma / vl
          a12 = -a*b / (gamma * vl)
          a13 = -a*c / (gamma * vl)
          a21 = 0.0e0
          a22 = c / gamma
          a23 = -b / gamma
          a31 = a / vl
          a32 = b / vl
          a33 = c / vl
          do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            dx = xi * a11 + yi * a12 + zi * a13
            dy = yi * a22 + zi * a23
            dz = xi * a31 + yi * a32 + zi * a33
            x(i) = dx
            y(i) = dy
            z(i) = dz
          end do
        endif
          
c-----------------------------------------------------------------------
c      Now the angle is driven
c      first, determine which rotation direction makes the angle smaller
c-----------------------------------------------------------------------
        
        gamma=x(dk)*y(di)-x(di)*y(dk)
        if((gamma.lt.0.).EQV.dwhich(di)) then
          signr =-1.0e0
        else
          signr = 1.0e0
        endif
        cost = cos(dstep*signr)
        sint = sin(dstep*signr)
        msint = -sint
        
        do i = 1,n
          if(dwhich(i)) then
            xi = x(i)
            yi = y(i)
            dx = xi * cost + yi * msint
            dy = xi * sint + yi * cost
            x(i) = dx
            y(i) = dy
          endif
        end do
      
c-----------------------------------------------------------------------
c     torsion driver
c-----------------------------------------------------------------------

      else

        if(dwhich(di))      call alignz(n,x,y,z,dj,di,1)
        if(.not.dwhich(di)) call alignz(n,x,y,z,dj,di,-1)
        cost = cos(dstep)
        sint = sin(dstep)
        msint = -sint        
        do i = 1,n
          if(dwhich(i)) then
            xi = x(i)
            yi = y(i)
            dx = xi * cost + yi * msint
            dy = xi * sint + yi * cost
            x(i) = dx
            y(i) = dy
          endif
        end do
      endif
      
      return
      end
