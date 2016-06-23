c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'sortarrays' - a sorting subroutine
c-----------------------------------------------------------------------   
c-----------------------------------------------------------------------   
      
      subroutine sortarrays(start,n,list,barr,key)

      implicit none
      integer start,n,i,j
      integer k, key(*)
      real  a, list(*)
      real  b, barr(*)
      real small
      
      small = 0.0001e0

      do j = start+1, start+n-1
         a = list(j)
         k = key(j)
         b = barr(j)
         do i = j-1,start,-1
            if(list(i).lt.a) goto 10
            if(abs(list(i)-a).lt.small) goto 10 
            list(i+1)=list(i)
            key(i+1)=key(i)
            barr(i+1)=barr(i)
         enddo
         i=start-1
  10     list(i+1)=a
         key(i+1)=k
         barr(i+1)=b
      enddo
      return
      END 
