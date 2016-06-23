c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'sortinteger' -- heapsort of integer array
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine sortinteger (n, list)
      implicit none
      integer i,j,k,n,counter
      integer list(*),lists

      k = n/2 + 1
      counter = n
      dowhile (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(counter)
            list(counter) = list(1)
            counter = counter - 1
            if (counter .le. 1) then
               list(1) = lists
               return
            end if
         end if
         i = k
         j = k + k
         dowhile (j .le. counter)
            if (j .lt. counter) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = counter + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
