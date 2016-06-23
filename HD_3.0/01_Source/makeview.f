c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'makeview' -- writes a series of molecules into a viewing
c     file, in rmds rank order.  
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine makeview
      
      implicit none
      include 'params.i'
      include 'storage.i'
      include 'linkage.i'
      include 'hosta.i'
      include 'control.i'
      
      integer i,j,ip, ii, jj, length
      integer nlink, n, n1, n2
      integer n12(MAXATOMD),i12(MAXVAL,MAXATOMD)
 
      real x(MAXATOMD),y(MAXATOMD),z(MAXATOMD)
      
      character*80 name
      character*83 name2
     
      length = 20
      do while (archive(length:length).eq.' ')
         length = length-1
      enddo

      open (unit = 9, file = 
     &    archive(1:length)//'_1.hdo', status = 'unknown')
      
      do j = 1, numview
      
         ip = pointer(j)
         ii = stoii(ip)
       
         write(9,'(i4)') stonatom(ip)
         
         write(name,100) j, stormsd(j), nrot(ii),
     &         stoen(ip), linkname(ii), linknum(ii)
         do i=1,77
            if (name(i:i).eq.' ') then
               name(i:i) = '_'
            endif
         enddo

         write(9,'(a)') name
            
         if(nsymattacha.gt.0) goto 50

c-----------------------------------------------------------------------     
c orient the molecules before we write them
c-----------------------------------------------------------------------     

         if(dowhat.eq.'LINK') then
            nlink = -1
            do i = 1, MAXLINKS
               if(linknum(ii).eq.linknum(i)) then
                  nlink = ng(i)
               endif
            enddo
            if(nlink.eq.-1) then
               write(*,*)'Subroutine makeview: '
               write(*,*)'nlink undefined'
               stop
            endif
            if(nlink.ne.0) then
               nlink = nlink - 2
            endif

            n = stonatom(ip)
            n1 = na - 1
            n2 = na - 1 + nlink

         else

            n = stonatom(ip)
            n1 = na - 2
            n2 = n

         endif    

         do i = 1, n
            x(i) = stox(ip,i)
            y(i) = stoy(ip,i)
            z(i) = stoz(ip,i)
            n12(i) = ston12(ip,i)
            do jj = 1, n12(i)
               i12(jj,i) = stoi12(ip,jj,i)
            enddo
         enddo

         call vieworient(n,n1,n2,x,y,z,n12,i12,dowhat)
    
         do i = 1, n
            stox(ip,i) = x(i)
            stoy(ip,i) = y(i)
            stoz(ip,i) = z(i)
         enddo

c-----------------------------------------------------------------------     
c write the molecules out         
c-----------------------------------------------------------------------     

   50    continue    

         if(xyz) then

            do i = 1, stonatom(ip)
               write(9,'(a2,2x,3F10.5)') stoatlab(ip,i), 
     &         stox(ip,i),stoy(ip,i),stoz(ip,i) 
            enddo

         else
            if(notype) then
               do i = 1, stonatom(ip)
                  write(9,101) stoatlab(ip,i),stox(ip,i),
     &            stoy(ip,i),stoz(ip,i),
     &            (stoi12(ip,jj,i), jj = 1,ston12(ip,i))
               enddo
            else
               do i = 1, stonatom(ip)
                  write(9,102) stoatlab(ip,i),stox(ip,i),
     &            stoy(ip,i),stoz(ip,i),stotype(ip,i),
     &            (stoi12(ip,jj,i), jj = 1,ston12(ip,i))
               enddo
            endif
         endif

      enddo   
         
  100 format(I5,':RMSD=',F5.3,',nrot=',I1,',E=',F6.3,' '
     &       ,A37,'(',I8,')')

  101 format(A2,2X,3F10.5,10(1X,I3))

  102 format(A2,2X,3F10.5,11(1X,I3))
          
      close(9)
         
      return
      end
