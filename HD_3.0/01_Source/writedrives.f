c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'writedrives' -- writes a series of structures to view the effect
c     of driving degrees of freedome in the input fragments
c
c     passed variables:
c     logical samehost     is true if HostA = HostB     
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine writedrives(samehost)
      
      implicit none
      include 'params.i'
      include 'hosta.i'
      include 'hostb.i'
      include 'control.i'
      include 'constants.i'
      
      integer i,j,k,jj
      real vl,theta,x0,y0
 
      
      logical samehost
     
      open (unit = 9, file = archive(1:index(archive,' ') - 1)//
     &    '_testa.hdo', status = 'unknown')


      i=1

      do j=1,ndrivea

c-----------------------------------------------------------------------     
c     align atoms consistently - atom 1 at origin, 2 on z axis
c-----------------------------------------------------------------------     

        call alignz(na,xa(1,j,1),ya(1,j,1),za(1,j,1),1,2,1)

        k=2
 1254   k=k+1
        vl=(xa(k,j,1)*xa(k,j,1)+ya(k,j,1)*ya(k,j,1))

        if(vl.lt.(1.0e-5)) then
          if(k.eq.na) goto 1258
          goto 1254
        endif

        if(abs(ya(k,j,1)).lt.(1.0e-06)) then
          theta=PI/2
        else 
          theta=-atan(-xa(k,j,1)/ya(k,j,1))
        endif

        if((xa(k,j,1)*sin(theta)+ya(k,j,1)*cos(theta)).lt.0.0e0)
     &    theta=theta+PI

        do k=3,na
          x0=xa(k,j,1)
          y0=ya(k,j,1)
          xa(k,j,1)=x0*cos(theta)-y0*sin(theta)
          ya(k,j,1)=x0*sin(theta)+y0*cos(theta)
        enddo

 1258   write(9,'(i4)') na
        write(9,'(a79)')driveinfoa(j)

        if(xyz) then

           do k=1,na
              write(9,'(a2,2x,3F10.5)') 
     &        atlaba(k),xa(k,j,i),ya(k,j,i),za(k,j,i)
           enddo

        else
           if(notype) then
              do k=1,na
                 write(9,'(A2,2X,3F10.5,10(1X,I3))')
     &           atlaba(k),xa(k,j,i),ya(k,j,i),za(k,j,i),
     &           (i12a(jj,k),jj=1,n12a(k))
              enddo
           else
              do k=1,na
                 write(9,'(A2,2X,3F10.5,11(1X,I3))')
     &           atlaba(k),xa(k,j,i),ya(k,j,i),za(k,j,i),
     &           itypea(k),(i12a(jj,k),jj=1,n12a(k))
              enddo
           endif
        endif

      enddo

      close(9)
      
      if(.not.samehost) then
        open (unit = 10, file = archive(1:index(archive,' ') - 1)//
     &    '_testb.hdo', status = 'unknown')

        do j=1,ndriveb
        call alignz(nb,xb(1,j,1),yb(1,j,1),zb(1,j,1),1,2,1)

        k=2
 1255   k=k+1
        vl=(xb(k,j,1)*xb(k,j,1)+yb(k,j,1)*yb(k,j,1))
        if(vl.lt.(1.0e-5)) then
          if(k.eq.nb) goto 1259
          goto 1255
        endif
        if(abs(yb(k,j,1)).lt.(1.0e-06)) then
          theta=PI/2
        else 
          theta=-atan(-xb(k,j,1)/yb(k,j,1))
        endif
        if((xb(k,j,1)*sin(theta)+yb(k,j,1)*cos(theta)).lt.0.0e0)
     &    theta=theta+PI
        do k=3,nb
          x0=xb(k,j,1)
          y0=yb(k,j,1)
          xb(k,j,1)=x0*cos(theta)-y0*sin(theta)
          yb(k,j,1)=x0*sin(theta)+y0*cos(theta)
        enddo

 1259   write(10,'(i4)') nb
        write(10,'(a79)')driveinfob(j)

        if(xyz) then

           do k=1,nb
             write(10,'(a2,2x,3F10.5)') 
     &       atlabb(k),xb(k,j,i),yb(k,j,i),zb(k,j,i)
           enddo

        else
           if(notype) then
              do k=1,nb
                write(10,'(A2,2X,3F10.5,10(1X,I3))') 
     &          atlabb(k),xb(k,j,i),yb(k,j,i),zb(k,j,i),
     &          (i12b(jj,k),jj=1,n12b(k))
              enddo
           else
              do k=1,nb
                write(10,'(A2,2X,3F10.5,11(1X,I3))') 
     &          atlabb(k),xb(k,j,i),yb(k,j,i),zb(k,j,i),
     &          itypeb(k),(i12b(jj,k),jj=1,n12b(k))
              enddo
           endif
        endif

      enddo

      close(10)

      endif
        
      return
      end
