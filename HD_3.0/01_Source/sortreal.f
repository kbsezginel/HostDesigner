c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'sortreal' -- heap sort of real array with pointers
c
c     passed variables:
c     integer n          number of entries in array
c     real    list       array of real numbers to be sorted
c     integer key        array of pointers to original position in list
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine sortreal(n, list, key)
      
      implicit none
      
      integer i, j, k, n, counter
      integer key(*), keys
      real  list(*), lists

      k = n/2 + 1
      counter = n
      dowhile (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
            keys = key(k)
         else
            lists = list(counter)
            keys = key(counter)
            list(counter) = list(1)
            key(counter) = key(1)
            counter = counter - 1
            if (counter .le. 1) then
               list(1) = lists
               key(1) = keys
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
               key(i) = key(j)
               i = j
               j = j + j
            else
               j = counter + 1
            end if
         end do
         list(i) = lists
         key(i) = keys
      end do

      return
      end

