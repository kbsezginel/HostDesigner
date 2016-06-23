c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'bumpcheckover' -- checks molecule C for close contacts between
c     Host A and link
c
c     passed variables:
c     logical pass        true if no bumps
c     integer ii          link fragment rank in the library
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------     

      subroutine bumpcheckover(pass,ii)
      
      implicit none
      
      include 'params.i'
      include 'constants.i'
      include 'moleculec.i'
      include 'linkage.i'
      
      real dx, dy, dz, r2, rpass2
      integer i, j, ii
      logical pass
      
      do i =1, nc-ng(ii)+2 

        do j = nc-ng(ii)+3, nc
  
            if((listc(i)+listc(j).gt.4).and.(listc0(i)+listc0(j)
     &          .lt.-4)) then
            
               dx = xc(i) - xc(j)
               dy = yc(i) - yc(j)
               dz = zc(i) - zc(j)
               
               r2 = dx*dx + dy*dy + dz*dz

               rpass2 = scalevdw12*(radc(i)+radc(j))*(radc(i)+radc(j))

               if((atlabc(i).eq.'H ').and.(atlabc(j).eq.'H ')) then
                   rpass2 = HHpass2
               end if
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
