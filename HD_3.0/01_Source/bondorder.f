c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'bondorder' -- determine bond orders needed to prepare MM
c     input files for the MMFF94 forcefield
c
c     currently this routine handles only C,H,N,O molecules.  Other
c     atom types will always be given a single bond
c
c     also the routine has been modified to correctly identify the
c     bond orders for nitrate, perchlorate, and sulfate anions
c
c     It must be further modified to handle other functional groups
c
c     passed variables:
c     integer  n        number of atoms
c     real     x,y,z    coordinates
c     integer  type     MM atom type
c     integer  n12      number of attached atoms
c     integer  i12      serial numbers of attached atoms
c     integer  bo12     bond order
c     char*2   atlab    atom label
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine bondorder(n,x,y,z,type,n12,i12,bo12,atlab)

      implicit none

      include 'params.i'

      integer i,j, k, jj, kk, ibi, ibj
      integer n, type(MAXATOMD),  n12(MAXATOMD)
      integer i12(MAXVAL,MAXATOMD), bo12(MAXVAL,MAXATOMD)

      real    x(MAXATOMD), y(MAXATOMD), z(MAXATOMD)
      real    xx,yy,zz,disij

      character*2  atlab(MAXATOMD)

c-----------------------------------------------------------------------
c     step 1: fill the bo12 array with ones, i.e., all connected atoms
c     have a bond order of one.
c-----------------------------------------------------------------------

      do i = 1, n
         do j = 1, n12(i)
            bo12(j,i)=1
         enddo
      enddo

c-----------------------------------------------------------------------
c     step 2: find all double and triple bonds - search only types for
c     carbon, oxygen, nitrogen
c-----------------------------------------------------------------------

      do i = 1, n
    
       if(atlab(i).eq.'C '.or.atlab(i).eq.'N '.or.atlab(i).eq.'O '.or.
     & atlab(i).eq.'P '.or.atlab(i).eq.'S '.or.atlab(i).eq.'Cl') then

       do jj = 1, n12(i)

       if(atlab(jj).eq.'C '.or.atlab(jj).eq.'N '.or.atlab(jj).eq.'O '
     & .or.atlab(jj).eq.'P '.or.atlab(jj).eq.'S '.or.atlab(jj).eq.'Cl') 
     & then

          j = i12(jj,i)
       
          xx = x(i) - x(j)
          yy = y(i) - y(j)
          zz = z(i) - z(j)
          disij=sqrt(xx*xx+yy*yy+zz*zz)
     
c-----------------------------------------------------------------------
c     C=O bond
c-----------------------------------------------------------------------

          if((atlab(i).eq.'C '.and.atlab(j).eq.'O '.or.atlab(i).eq.'O '
     &    .and.atlab(j).eq.'C ').and.disij.le.1.30E0) then

             if(i.lt.j) then

                bo12(jj,i) = bo12(jj,i) + 1

                do k = 1, n12(j)
                   if(i12(k,j).eq.i) then
                      bo12(k,j) = bo12(k,j) + 1
                   endif
                enddo
          
             endif

c-----------------------------------------------------------------------
c     alkyne
c-----------------------------------------------------------------------

          elseif(atlab(i).eq.'C '.and.atlab(j).eq.'C '.and.disij.le.
     &    1.25E0) then

             if(i.lt.j) then
 
                bo12(jj,i) = bo12(jj,i) + 2

                do k = 1, n12(j)
                   if(i12(k,j).eq.i) then
                      bo12(k,j) = bo12(k,j) + 2
                   endif
                enddo

             endif
             
c-----------------------------------------------------------------------
c     C=N
c-----------------------------------------------------------------------

          elseif((atlab(i).eq.'C '.and.atlab(j).eq.'N '.or.atlab(i)
     &    .eq.'N '.and.atlab(j).eq.'C ').and.disij.le.1.37E0) then

             ibi = 0
             do k = 1, n12(i)
                if(bo12(k,i).ne.9) then
                    ibi = ibi + bo12(k,i)
                endif
             enddo

             ibj = 0
             do k = 1, n12(j)
                if(bo12(k,j).ne.9) then
                   ibj = ibj + bo12(k,j)
                endif
             enddo

             if(atlab(i).eq.'N ') then

                if(ibi.le.2.and.ibj.le.3) then

                   if(i.lt.j) then
                      bo12(jj,i) = bo12(jj,i) + 1
                      do k = 1, n12(j)
                         if(i12(k,j).eq.i) then
                            bo12(k,j) = bo12(k,j) + 1
                         endif
                      enddo
                   endif
                endif

             else

                if(ibi.le.3.and.ibj.le.2) then
                   if(i.lt.j) then
                      bo12(jj,i) = bo12(jj,i) + 1
                      do k = 1, n12(j)
                         if(i12(k,j).eq.i) then
                            bo12(k,j) = bo12(k,j) + 1
                         endif
                      enddo
                   endif
                endif
             endif
           
c-----------------------------------------------------------------------
c     nitrile
c-----------------------------------------------------------------------

          elseif((atlab(i).eq.'C '.and.atlab(j).eq.'N '.or.atlab(i)
     &    .eq.'N '.and.atlab(j).eq.'C ').and.disij.le.1.18E0) then

            if(i.lt.j) then

               bo12(jj,i) = bo12(jj,i) + 2

               do k = 1, n12(j)
                  if(i12(k,j).eq.i) then
                     bo12(k,j) = bo12(k,j) + 2
                  endif
               enddo

            endif

c-----------------------------------------------------------------------
c     alkenes and alkynes
c-----------------------------------------------------------------------

          elseif(atlab(i).eq.'C '.and.atlab(j).eq.'C '.and.
     &    disij.le.1.423E0) then

             ibi = 0
             do k = 1, n12(i)
                if(bo12(k,i).ne.9) then
                    ibi = ibi + bo12(k,i)
                endif
             enddo

             ibj = 0
             do k = 1, n12(j)
                if(bo12(k,j).ne.9) then
                   ibj = ibj + bo12(k,j)
                endif
             enddo

             if(ibi.le.3.and.ibj.le.3) then

                if(i.lt.j) then
 
                    bo12(jj,i) = bo12(jj,i) + 1

                    do k = 1, n12(j)
                       if(i12(k,j).eq.i) then
                          bo12(k,j) = bo12(k,j) + 1
                       endif
                    enddo

                endif

             endif

c-----------------------------------------------------------------------
c     nitrate
c-----------------------------------------------------------------------

          elseif(atlab(i).eq.'O '.and.atlab(j).eq.'N '.or.atlab(i)
     &    .eq.'N '.and.atlab(j).eq.'O ') then

             if(atlab(i).eq.'N ') then
                 kk = i
             else
                 kk = j
             endif

             ibi = 0
             do k = 1, n12(kk)
                if(atlab(i12(k,kk)).eq.'O ') then
                   ibi = ibi + 1
                endif
             enddo

             if(ibi.eq.3) then

                ibi = 0
                do k = 1, n12(i)
                   if(bo12(k,i).ne.9) then
                      ibi = ibi + bo12(k,i)
                   endif
                enddo

                ibj = 0
                do k = 1, n12(j)
                   if(bo12(k,j).ne.9) then
                      ibj = ibj + bo12(k,j)
                   endif
                enddo

                if(atlab(i).eq.'N ') then

                   if(ibi.le.4.and.ibj.le.1) then

                      if(i.lt.j) then
                         bo12(jj,i) = bo12(jj,i) + 1
                         do k = 1, n12(j)
                            if(i12(k,j).eq.i) then
                               bo12(k,j) = bo12(k,j) + 1
                            endif
                         enddo
                      endif
                   endif

                else

                   if(ibi.le.1.and.ibj.le.4) then
                      if(i.lt.j) then
                         bo12(jj,i) = bo12(jj,i) + 1
                         do k = 1, n12(j)
                            if(i12(k,j).eq.i) then
                               bo12(k,j) = bo12(k,j) + 1
                            endif
                         enddo
                      endif
                   endif
                endif

             endif 
      
c-----------------------------------------------------------------------
c     perchlorate
c-----------------------------------------------------------------------

          elseif(atlab(i).eq.'O '.and.atlab(j).eq.'Cl'.or.atlab(i)
     &    .eq.'Cl'.and.atlab(j).eq.'O ') then

             if(atlab(i).eq.'Cl') then
                 kk = i
             else
                 kk = j
             endif

             ibi = 0
             do k = 1, n12(kk)
                if(atlab(i12(k,kk)).eq.'O ') then
                   ibi = ibi + 1
                endif
             enddo

             if(ibi.eq.4) then

                ibi = 0
                do k = 1, n12(i)
                   if(bo12(k,i).ne.9) then
                      ibi = ibi + bo12(k,i)
                   endif
                enddo

                ibj = 0
                do k = 1, n12(j)
                   if(bo12(k,j).ne.9) then
                      ibj = ibj + bo12(k,j)
                   endif
                enddo

                if(atlab(i).eq.'Cl') then

                   if(ibi.le.6.and.ibj.le.1) then

                      if(i.lt.j) then
                         bo12(jj,i) = bo12(jj,i) + 1
                         do k = 1, n12(j)
                            if(i12(k,j).eq.i) then
                               bo12(k,j) = bo12(k,j) + 1
                            endif
                         enddo
                      endif
                   endif

                else

                   if(ibi.le.1.and.ibj.le.6) then
                      if(i.lt.j) then
                         bo12(jj,i) = bo12(jj,i) + 1
                         do k = 1, n12(j)
                            if(i12(k,j).eq.i) then
                               bo12(k,j) = bo12(k,j) + 1
                            endif
                         enddo
                      endif
                   endif
                endif

             endif 
              
c-----------------------------------------------------------------------
c     sulfate
c-----------------------------------------------------------------------

          elseif(atlab(i).eq.'O '.and.atlab(j).eq.'S '.or.atlab(i)
     &    .eq.'S '.and.atlab(j).eq.'O ') then

             if(atlab(i).eq.'S ') then
                 kk = i
             else
                 kk = j
             endif

             ibi = 0
             do k = 1, n12(kk)
                if(atlab(i12(k,kk)).eq.'O ') then
                   ibi = ibi + 1
                endif
             enddo

             if(ibi.eq.4) then

                ibi = 0
                do k = 1, n12(i)
                   if(bo12(k,i).ne.9) then
                      ibi = ibi + bo12(k,i)
                   endif
                enddo

                ibj = 0
                do k = 1, n12(j)
                   if(bo12(k,j).ne.9) then
                      ibj = ibj + bo12(k,j)
                   endif
                enddo

                if(atlab(i).eq.'S ') then

                   if(ibi.le.5.and.ibj.le.1) then

                      if(i.lt.j) then
                         bo12(jj,i) = bo12(jj,i) + 1
                         do k = 1, n12(j)
                            if(i12(k,j).eq.i) then
                               bo12(k,j) = bo12(k,j) + 1
                            endif
                         enddo
                      endif
                   endif

                else

                   if(ibi.le.1.and.ibj.le.5) then
                      if(i.lt.j) then
                         bo12(jj,i) = bo12(jj,i) + 1
                         do k = 1, n12(j)
                            if(i12(k,j).eq.i) then
                               bo12(k,j) = bo12(k,j) + 1
                            endif
                         enddo
                      endif
                   endif
                endif

             endif 
             

              
c-----------------------------------------------------------------------
c     add new cases here
c-----------------------------------------------------------------------

          endif
          
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

       endif
       enddo
       endif
      enddo

      return
      end
