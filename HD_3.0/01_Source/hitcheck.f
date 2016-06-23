c----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subroutine hitcheck -- evaluates and stores new hits
c
c     hitcheck performs three functions:
c     1. It evaluates each hit for uniqueness by comparing the RMSD of 
c        the guest molecule distance and values of the three principal 
c        moments of inertia to those in the storage array
c     2. If the uniqueness check fails and there is a matching hit:
c      a.  If the RMSD of the current hit is less than that of the stored
c          matching hit, the current hit is written over the stored hit, 
c          resorting the hit array to reflect the new RMSD rank order
c      b.  in either case, the subroutine returns to linker
c     3. If the hit passes the checks, hitcheck stores the hit into the
c        array of stored hits- if the storage array is MAXHITS in size,
c        it is stored in place of the previous highest RMSD hit.
c-----------------------------------------------------------------------     

      subroutine hitcheck(n,x,y,z,atwt,itype,iguest,atlab,ii,
     &    n12,i12,rmsd,numhits,mac,energy,nrotab) 

      implicit none
      include 'params.i'
      include 'constants.i'
      include 'storage.i'
      include 'control.i'
      
      real x(MAXATOMD), y(MAXATOMD), z(MAXATOMD)
      real atwt(MAXATOMD)
      real ix, iy, iz, rmsd, rmsdiff, energy
      character*2 atlab(MAXATOMD)
      integer n, ii, i, j, k, i12(MAXVAL,MAXATOMD),ip, mac
      integer itype(MAXATOMD),iguest(MAXATOMD),n12(MAXATOMD)
      integer numhits
      integer nrotab

c-----------------------------------------------------------------------
c    ix,iy, and iz are the three pricipal moments of inertia of the hit.
c    n,x,y,z,atwt,itype,iguest,atlab,n12,i12 are all the analogously named
c    variables from linker or overlay and are passed in and out unchanged.
c    rmsd is distance between the guest on host A and that on host B.
c    numhits is the running tally of hits found for the run.
c    ii is the number of the current link; the ii'th link examined so far.
c-----------------------------------------------------------------------

      if (numhits.lt.numkeep) then

        call inertia(n,x,y,z,atwt,ix,iy,iz)

c-----------------------------------------------------------------------
c       FIRST CASE - storage array is not yet full
c       First, check for degeneracy
c-----------------------------------------------------------------------

        do ip = 1, numhits
          if(ii.eq.stoii(ip)) then
          rmsdiff = abs(rmsd - stormsd(ip))
            if(rmsdiff.le.rmstol) then
              if((abs((ix-stoix(ip))/stoix(ip))).le.moitol)then
                if((abs((iy-stoiy(ip))/stoiy(ip))).le.moitol)then
                  if((abs((iz-stoiz(ip))/stoiz(ip))).le.moitol)then
                    if(rmsd.lt.stormsd(ip)) then
                      stoix(ip)=ix
                      stoiy(ip)=iy
                      stoiz(ip)=iz
                      stormsd(ip)=rmsd
                      do i = 1, n
                        stox(ip,i) = x(i)
                        stoy(ip,i) = y(i)
                        stoz(ip,i) = z(i)
                        stotype(ip,i) = itype(i)
                        stoatlab(ip,i) = atlab(i)
                        ston12(ip,i) = n12(i)
                        do j = 1, n12(i)
                          stoi12(ip, j, i) = i12(j,i)
                        end do
                      enddo
                    endif
                    goto 400
                  endif
                endif
              endif

            endif
          endif
        end do

c-----------------------------------------------------------------------
c       IF WE GET THIS FAR, THEN WE HAVE A KEEPER - just add it to the
c       list at pointer position = numhits           
c-----------------------------------------------------------------------

        numhits = numhits + 1
        ip = numhits
        pointer(ip) = numhits
        stormsd(ip) = rmsd
        stoii(ip)  = ii
        stoix(ip)=ix
        stoiy(ip)=iy
        stoiz(ip)=iz
        stonatom(ip) = n
        stoen(ip) = energy
        storotab(ip) = nrotab
        do i = 1, n
          stox(ip,i) = x(i)
          stoy(ip,i) = y(i)
          stoz(ip,i) = z(i)
          stotype(ip,i) = itype(i)
          stoguest(ip,i) = iguest(i)
          stoatlab(ip,i) = atlab(i)
          ston12(ip,i) = n12(i)
          do j = 1, n12(i)
            stoi12(ip, j, i) = i12(j,i)
          end do
        end do
        if(numhits.eq.numkeep) call sortreal(numkeep, stormsd, pointer)
        
c-----------------------------------------------------------------------
c     SECOND CASE - storage array is full
c
c     first check to see if rmsd exceeds current maximum which will
c     always be located in the bin (MAXHITS)
c-----------------------------------------------------------------------

      ELSE

        call inertia(n,x,y,z,atwt,ix,iy,iz)

c-----------------------------------------------------------------------
c       Next, check for degeneracy
c-----------------------------------------------------------------------

        do i = 1, numkeep
          ip = pointer(i)

          if(ii.eq.stoii(ip)) then
            rmsdiff = abs(rmsd - stormsd(i))
            if(rmsdiff.le.rmstol) then

              if((abs((ix-stoix(ip))/stoix(ip))).le.moitol)then
                if((abs((iy-stoiy(ip))/stoiy(ip))).le.moitol)then
                  if((abs((iz-stoiz(ip))/stoiz(ip))).le.moitol)then
                    if(rmsd.lt.stormsd(i)) then
                      stoix(ip)=ix
                      stoiy(ip)=iy
                      stoiz(ip)=iz
                      stormsd(i)=rmsd
                      do k = 1, n
                        stox(ip,k) = x(k)
                        stoy(ip,k) = y(k)
                        stoz(ip,k) = z(k)
                        stotype(ip,k) = itype(k)
                        stoatlab(ip,k) = atlab(k)
                        ston12(ip,k) = n12(k)
                        do j = 1, n12(k)
                          stoi12(ip, j, k) = i12(j,k)
                        end do
                      enddo
                      call sortreal(numkeep, stormsd, pointer)
                    endif
                    goto 400
                  endif
                endif
              endif
           
            endif 
          endif
        end do

c-----------------------------------------------------------------------
c       IF WE GET THIS FAR, THEN WE HAVE A KEEPER - now the storage 
c       position will be in the register pointer(MAXHITS).           
c-----------------------------------------------------------------------

        ip = pointer(numkeep)
        stormsd(numkeep) = rmsd
        stoii(ip)  = ii
        stoix(ip)=ix
        stoiy(ip)=iy
        stoiz(ip)=iz
        stonatom(ip) = n
        stoen(ip) = energy
        storotab(ip) = nrotab
        do i = 1, n
          stox(ip,i) = x(i)
          stoy(ip,i) = y(i)
          stoz(ip,i) = z(i)
          stotype(ip,i) = itype(i)
          stoguest(ip,i) = iguest(i)
          stoatlab(ip,i) = atlab(i)
          ston12(ip,i) = n12(i)
          do j = 1, n12(i)
            stoi12(ip,j,i) = i12(j,i)
          end do
        end do
           
        call sortreal(numkeep, stormsd, pointer)
        maxrmsd = stormsd(numkeep)
   
      endif
      
 400  return
      end
