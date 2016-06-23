c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'bumpcheckd' -- checks molecule D for close contacts between 
c     molecule C and host B
c
c     passed variables:
c     logical pass        true if no bumps
c     integer jk,mn       index variables for loops on drivea and driveb
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------     

      subroutine bumpcheckd(pass,jk,mn)

      implicit none

      include 'params.i'
      include 'constants.i'
      include 'hosta.i'
      include 'moleculec.i'
      include 'moleculed.i'
      include 'mguest.i'

      real dx, dy, dz, r2, rpass2

      integer i, j, m, jk, mn
 
      logical pass
 
c-----------------------------------------------------------------------
c     enter the bumpchecking loop.  
c-----------------------------------------------------------------------     

      do i = nc, nd

         do j = 1, nc - 1
            if(listd(i)+listd(j).gt.4) then

               dx = xd(i) - xd(j)
               dy = yd(i) - yd(j)
               dz = zd(i) - zd(j)

               r2 = dx * dx + dy * dy + dz * dz

               if((atlabd(i).eq.'H ').and.(atlabd(j).eq.'H ')) then

                 rpass2 = HHpass2

               else

                 rpass2 = scalevdw12*(radd(i)+radd(j))*(radd(i)+radd(j))

               endif

c-----------------------------------------------------------------------
c     need to allow for the guest on one part coming close to the 
c     attached atom of the guest on the other part.
c----------------------------------------------------------------------- 

               if(i.gt.nd-lguest) then

                 if((j.ge.na-lguest).and.(j.lt.na)) then

                   rpass2 = 0.0e0

                 else

                  do m = 1, gstn12b(i-nd+lguest)
                     if (j.eq.gsti12b(m,i-nd+lguest))
     &                 rpass2=gstbmpb(m,mn,i-nd+lguest)
                  end do

                endif

               else if((j.ge.na-lguest).and.(j.lt.na)) then

                  do m = 1, gstn12a(j-na+lguest+1)
                     if (i-nc+1.eq.gsti12a(m,j-na+lguest+1))
     &                 rpass2=gstbmpa(m,jk,j-na+lguest+1)
                  end do

               end if

c-----------------------------------------------------------------------
c     do the test
c-----------------------------------------------------------------------     

               if(r2.lt.rpass2) then

                  pass = .false. 
                  return

               endif

            end if

         end do

      end do

      pass = .true.

      return
      end
