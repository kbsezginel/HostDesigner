c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'connlistg' -- fills an array with connectivity values:
c
c            if binding atom, then value is 1
c            if 1,2 connected, then value is 2
c            if 1,3 connected, then value is 3
c            else the value will be 4
c
c     passed variables:
c     integer       n         number of atoms
c     integer       n12       number of attached atoms
c     integer       i12       serial numbers of attached atoms
c     integer       ibond     serial number of bonding atom
c     integer       list      connectivity values
c     integer       intsign   sign of the connectivity value
c
c     This routine is used in place of connlist in order to allow
c     a change in the size of the i12g matrix
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine connlistg(n, n12, i12, ibond, list, intsign)
      
      implicit none

      include 'params.i'
            
      integer  n12(MAXATOMG),i12(MAXVAL,MAXATOMG),ibond
      integer  list(MAXATOMG)
      integer  i, j, k, n, intsign
      integer  jatt, katt
      
      do i = 1, n
          list(i) = intsign * 4
      end do
      
      do j = 1, n12(ibond)
      
         jatt = i12(j,ibond)
      
         do k = 1, n12(jatt)
             katt = i12(k,jatt)
            list(katt) = intsign * 3
         end do
                         
      end do
     
      do j = 1, n12(ibond)
         jatt = i12(j,ibond)
         list(jatt) = intsign * 2
      end do
               
      list(ibond) = intsign * 1
      
      return
      end
