c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'guestsym' -- determines symmetry properties of guests and a map
c      of the ways that guest1 can be overlaid upon guest2
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine guestsym

      implicit none
      include 'params.i'
      include 'hosta.i'
      include 'hostb.i'
      include 'mguest.i'

      integer i,j,map(MAXGST),mmap(MAXGST),point(MAXGST)
      integer depth,pa,pb
      logical sdist,sdist2
      real d1(MAXGST,MAXGST),d2(MAXGST,MAXGST), ddist,dx,dy,dz

c-----------------------------------------------------------------------     
c     set up interatomic distance matrix for both guests.
c     In each case, the matrix includes only the guest atoms,
c     so the input coordinates must be offset to include only guest atoms.
c-----------------------------------------------------------------------     

      pa=na-lguest
      pb=nb-lguest
      
      if(lguest.eq.1) then
        nequiva=1
        nequivb=1
        mapeqa(1,1)=1
        mapeqb(1,1)=1
        return
      endif 

      if(lguest.eq.2) then
         if(atlaba(na-1).ne.atlaba(na)) then
           nequiva=1
           nequivb=1
           if(atlaba(na-1).eq.atlabb(nb-1)) then
              mapeqa(1,1)=1
              mapeqa(1,2)=2
              mapeqb(1,1)=1
              mapeqb(1,2)=2
           else
              mapeqa(1,1)=1
              mapeqa(1,2)=2
              mapeqb(1,1)=2
              mapeqb(1,2)=1
           endif
           return
         endif
      endif 

      do i=1,lguest
        do j=i+1,lguest
          dx=(xa(pa+i,1,1)-xa(pa+j,1,1))
          dy=(ya(pa+i,1,1)-ya(pa+j,1,1))
          dz=(za(pa+i,1,1)-za(pa+j,1,1))
          d1(i,j)= sqrt(dx*dx+dy*dy+dz*dz)
          d1(j,i)= d1(i,j)
          dx=(xb(pb+i,1,1)-xb(pb+j,1,1))
          dy=(yb(pb+i,1,1)-yb(pb+j,1,1))
          dz=(zb(pb+i,1,1)-zb(pb+j,1,1))
          d2(i,j)= sqrt(dx*dx+dy*dy+dz*dz)
          d2(j,i)= d2(i,j)
        enddo
      enddo

c-----------------------------------------------------------------------     
c initialize counters for recursive loop
c-----------------------------------------------------------------------     

      depth=1
      do i=1,lguest
        map(i)=i
        mmap(i)=i
        point(i)=1
      enddo
      nequiva=0

c-----------------------------------------------------------------------
c    Top of recursive symmetry checking loop
c    This was written recursively so that the loop can go arbitrarily
c    deep (a C60 fullerene was one of the test cases, requiring 60 
c    nested loops.)
c
c    Each symmetry equivalent set of positions is saved as a mapeq entry.
c    'map' refers to a N-long set of serial numbers which correspond
c    to an ordering of serial numbers of the input atoms.
c
c    'mmap' is the metamap, or the map of 'map'.  It is another N-long
c    set of serial numbers indicating where each integer in 'map' is.
c
c    The basic structure is that the code is examining the distance 
c    between pairs of atoms (there are (N*N-N)/2 distinct pairs) for 
c    similarity.  If all pairwise distances are within tolerance, the 
c    map of atoms is saved.  Note that this does not distinguish between
c    mirror images; this is deliberate as the code currently can mirror
c    the hosts internally.  The incorrect mirror image will simply overlap
c    poorly in the sections of code where mapeq is used.
c
c    'depth' is the index number of the point in 'map' being checked.
c    'point' is the value to which map(depth) is to be set.
c    A single run through the loop checks whether setting map(depth) to
c    be equal to point will result in pairwise distances within tolerance
c    to each of the previously determined points (which are those < depth)
c    if so, it passes down to the depth+1 to try further, otherwise 
c    it tries the next atom at this depth by adding one to point.
c    It also can add to the mapeq if the depth reaches N-1 or back up 
c    in depth if all possibilities at that depth have been examined.
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------     
c     make sure that the point under examination has not yet been used
c-----------------------------------------------------------------------     

 10   if(mmap(point(depth)).lt.depth) goto 20
      
c-----------------------------------------------------------------------     
c     swap the current map(depth) with the one equal to 'point',
c     and adjust the metamap to point to the correct map locations
c-----------------------------------------------------------------------     

      mmap(map(depth))=mmap(point(depth))
      map(mmap(point(depth)))=map(depth)
      mmap(point(depth))=depth
      map(depth)=point(depth)
      
c-----------------------------------------------------------------------     
c     check distances between the new index atom and the older ones
c-----------------------------------------------------------------------     

      sdist=.true.
      do i=1,depth-1
        ddist=abs(d1(i,depth)-d2(map(i),map(depth)))
        if(ddist.gt.(0.1e0)) sdist=.false.
      enddo
      
      if(sdist) then

c-----------------------------------------------------------------------     
c     if we're checking the second-to-last element, and it passed the
c     check, go ahead and check the final, depth+1, element.
c     if that passes, save the whole current map to mapeq and
c     then proceed through 20, back up to smaller depth.
c-----------------------------------------------------------------------     

        if(depth.eq.lguest-1) then
          sdist2=.true.
          do i=1,depth
            ddist=abs(d1(i,lguest)-d2(map(i),map(lguest)))
            if(ddist.gt.(0.1e0)) sdist2=.false.
          enddo
          if(sdist2) then
            nequiva=nequiva+1
            if(nequiva.gt.1000) then
              write(6,*) 
     &        'exceeded 1000 symmetry equivalent guest positions!'
              goto 30
            endif
            do i=1,lguest
              mapeqa(nequiva,i)=map(i)
            enddo
          endif
          
c-----------------------------------------------------------------------     
c     upon a normal success, check deeper starting with point=1.
c-----------------------------------------------------------------------     

        else
          depth=depth+1
          point(depth)=1
          goto 10
        endif
      endif
          
c-----------------------------------------------------------------------     
c     upon failure or competion and saving of the map, either try the
c     next value of point or, if all points at this depth have been
c     tried, go up a level in depth and try the next point there.  
c     Exit the routine after all points at depth 1 (and therefore, all 
c     branches below depth 1) have been tried.
c-----------------------------------------------------------------------     

 20   if(point(depth).eq.lguest) then
        if(depth.eq.1) goto 30
        depth=depth-1
        goto 20
      else
        point(depth)=point(depth)+1
        goto 10
      endif

 30   continue
 
c-----------------------------------------------------------------------
c  Do it again in the reverse direction (needed for bump check scaling)
c-----------------------------------------------------------------------

      depth=1
      do i=1,lguest
        map(i)=i
        mmap(i)=i
        point(i)=1
      enddo
      nequivb=0
      
 110  if(mmap(point(depth)).lt.depth) goto 120
      mmap(map(depth))=mmap(point(depth))
      map(mmap(point(depth)))=map(depth)
      mmap(point(depth))=depth
      map(depth)=point(depth)
      sdist=.true.
      do i=1,depth-1
        ddist=abs(d2(i,depth)-d1(map(i),map(depth)))
        if(ddist.gt.(0.1e0)) sdist=.false.
      enddo
      if(sdist) then
        if(depth.eq.lguest-1) then
          sdist2=.true.
          do i=1,depth
            ddist=abs(d2(i,lguest)-d1(map(i),map(lguest)))
            if(ddist.gt.(0.1e0)) sdist2=.false.
          enddo
          if(sdist2) then
            nequivb=nequivb+1
            if(nequivb.gt.1000) then
              write(6,*) 
     &        'exceeded 1000 symmetry equivalent guest positions!'
              goto 130
            endif
            do i=1,lguest
              mapeqb(nequivb,i)=map(i)
            enddo
          endif
        else
          depth=depth+1
          point(depth)=1
          goto 110
        endif
      endif
 120  if(point(depth).eq.lguest) then
        if(depth.eq.1) goto 130
        depth=depth-1
        goto 120
      else
        point(depth)=point(depth)+1
        goto 110
      endif
 130  continue

      if(nequiva.ne.nequivb) then
        write(6,*)
     &'WARNING: Superimposition of guest A onto B results in different'
        write(6,*)
     &'numbers of symmetry equivalent positions than does B onto A!'
        write(6,*)'Symmetrization of the guests is recommended.'
      endif
      
      return
      end
