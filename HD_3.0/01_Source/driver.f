c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c  'driver' -- varies distances, angles, and dihedrals, calls the 
c  drive, which moves each step
c
c  passed variables:
c  character  dtype      drive type (B=Bond,A=Angle,D=Dihedral,+=sym)
c  character  dstep      length of drive step (in Angstroms or degrees)
c  character  dwhich     array specifying which atoms are to be moved
c  integer    n          number of atoms
c  real       x,y,z      coords
c  integer    di,dj,dk   serial numbers of atoms in drive
c  integer    numsym     0 for each drive if there is no symmetry,
c                        1 for the first symmetrical drive, 2 for the
c                        second, and so on.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine driver(ndrives,dtype,dstart,dstep,nstep,dwhat,n,xg0,yg0
     &  ,zg0,di,dj,dk,n12,i12,nattach,ibonder,ilooser,drivinga,drivenum,
     &  numsym)
     
      implicit none
      include 'params.i'
      include 'mguest.i'
      include 'hosta.i'
      include 'hostb.i'
      include 'constants.i'
     
      integer n,di(10),dj(10),dk(10),i,j,ndrives,intd(MAXATOMA)
      integer n12(MAXATOMA),i12(MAXVAL,MAXATOMA)
      integer point(10), depth, nstep(10)
      integer numsym(10), nn, mm 
      integer ibonder(MAXATTACHA),ilooser(MAXATTACHA),nattach
      integer drivenum

      real dstart(10),dstep(10)
      real xg0(MAXATOMA),yg0(MAXATOMA),zg0(MAXATOMA)
      real x(MAXATOMA),y(MAXATOMA),z(MAXATOMA),dback(10)

      character*1 dtype(10)
      character*1 dwhat(10)

      logical dwhich(MAXATOMA,10),drivinga
      
c-----------------------------------------------------------------------
c     partition the molecules based on the rules from input
c-----------------------------------------------------------------------

      do i=1,ndrives

         if(dwhat(i).eq.'I') then

            call bondsplit(n,n12,i12,di(i),dj(i),intd)

            do j=1,n
               if(intd(j).eq.2) then
                  dwhich(j,i)=.true.
               else
                  dwhich(j,i)=.false.
               endif
            enddo

         elseif(dwhat(i).eq.'G') then

            do j=1,n-lguest
               dwhich(j,i)=.false.
            enddo
            do j=n-lguest+1,n
               dwhich(j,i)=.true.
            enddo

c-----------------------------------------------------------------------
c     error checking; ensure one terminus is in each partition
c-----------------------------------------------------------------------

            if((dtype(i).eq.'B').or.(dtype(i).eq.'D')) then
               if(dwhich(di(i),i).EQV.dwhich(dj(i),i)) then
                  write(6,100)
                  stop
               endif
            elseif(dtype(i).eq.'A') then
               if(dwhich(di(i),i).EQV.dwhich(dk(i),i)) then
                  write(6,100)
                  stop
               endif
            endif
  100       format
     &('guest-partitioned drives must have exactly one guest terminus!')

         else

            write(6,*) 'partition rule of drive',i,' not recognized!'
            stop

         endif

c-----------------------------------------------------------------------
c     convert degrees to radians
c-----------------------------------------------------------------------

         if((dtype(i).eq.'A').or.(dtype(i).eq.'D')) then      
            dstep(i) =  dstep(i)*PI/180.0e0
            dstart(i)= dstart(i)*PI/180.0e0
         endif

      enddo
        
c-----------------------------------------------------------------------
c     setup coordinates
c-----------------------------------------------------------------------

      do i=1,n
         x(i)=xg0(i)
         y(i)=yg0(i)
         z(i)=zg0(i)
      enddo

c-----------------------------------------------------------------------
c     initialize structure to the starting point of each drive
c-----------------------------------------------------------------------

      do i=1,ndrives
         call drive(dtype(i),dstart(i),dwhich(1,i)
     &   ,n,x,y,z,di(i),dj(i),dk(i))
         point(i)=1
         dback(i)=-dstep(i)*real(nstep(i)-1)
      enddo

      depth=ndrives
      drivenum=0

c-----------------------------------------------------------------------
c     top of recursive loop; save the structure, drive the deepest loop
c-----------------------------------------------------------------------

    1 continue

      drivenum=drivenum+1

      if (drivinga) then

         do j=1,nattach
            call alignzpass(n, x, y, z, ibonder(j), ilooser(j), 1,
     &      xa(1,drivenum,j), ya(1,drivenum,j), za(1,drivenum,j))
         enddo

         do i=1,nstep(depth)-1
             do nn = 0, numsym(depth)  
               call drive(dtype(depth-nn),dstep(depth-nn),
     &         dwhich(1,depth-nn), n, x, y, z, di(depth-nn),
     &         dj(depth-nn),dk(depth-nn))
             enddo
             drivenum=drivenum+1
             
             do j=1,nattach
                call alignzpass(n, x, y, z, ibonder(j), ilooser(j), 1,
     &          xa(1,drivenum,j), ya(1,drivenum,j), za(1,drivenum,j))
             enddo
         enddo

      else 

         if(numsym(depth).ne.0) then
           write(6,*) 'symmetric drive does not apply to hostb'
           stop
         endif 

         do j=1,nattach
            call alignzpass(n, x, y, z, ibonder(j), ilooser(j), -1,
     &      xb(1,drivenum,j), yb(1,drivenum,j), zb(1,drivenum,j))
         enddo

         do i=1,nstep(depth)-1
            call drive(dtype(depth),dstep(depth),dwhich(1,depth)
     &      ,n,x,y,z,di(depth),dj(depth),dk(depth))
            drivenum=drivenum+1
            do j=1,nattach
               call alignzpass(n, x, y, z, ibonder(j), ilooser(j), -1,
     &         xb(1,drivenum,j), yb(1,drivenum,j), zb(1,drivenum,j))
            enddo
         enddo

      endif
      
c-----------------------------------------------------------------------
c  reset deepest drive to starting position 
c-----------------------------------------------------------------------

      mm = numsym(depth)
      do nn = 0, mm
        call drive(dtype(depth-nn),dback(depth-nn),dwhich(1,depth-nn)
     &       ,n,x,y,z,di(depth-nn),dj(depth-nn),dk(depth-nn))
      enddo
      depth=depth-1-mm

c-----------------------------------------------------------------------
c     exit if done
c-----------------------------------------------------------------------

    2 if(depth.eq.0) goto 999

c-----------------------------------------------------------------------
c     normal driving (all loops but deepest)
c-----------------------------------------------------------------------

      if(point(depth).lt.nstep(depth)) then
         do nn = 0, numsym(depth)   
           call drive(dtype(depth-nn),dstep(depth-nn),
     &          dwhich(1,depth-nn), n, x, y, z, di(depth-nn),
     &          dj(depth-nn),dk(depth-nn))
         enddo
         point(depth)=point(depth)+1
         depth=ndrives
         goto 1
      else
         mm = numsym(depth)
         do nn = 0, mm
           call drive(dtype(depth-nn),dback(depth-nn),
     &          dwhich(1,depth-nn), n, x, y, z, di(depth-nn),
     &          dj(depth-nn),dk(depth-nn))
        enddo
         point(depth)=1
         depth=depth-1-mm
         goto 2
      endif

c-----------------------------------------------------------------------
c     get out
c-----------------------------------------------------------------------

  999 continue
      return
      end
