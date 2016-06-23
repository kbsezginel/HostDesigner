c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'bondsplit' -- determines atoms on either side of a bond
c
c     passed variables:
c     integer    n          number of atoms
c     integer    n12        number of attached atoms
c     integer    i12        serial numbers of attached atoms
c     integer    cut1       serial number of first atom in bond
c     integer    cut2       serial number of second atom in bond
c     integer    partition  matrix by serial number (0=unassigned),
c                           1=with cut1, 2 with cut2)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine bondsplit(n,n12,i12,cut1,cut2,partition)

      implicit none
      include 'params.i'

      integer  n,n12(MAXATOMA),i12(MAXVAL,MAXATOMA)
      integer  cut1,cut2,partition(MAXATOMA)
      integer  i,bond(MAXATOMA),height,atom(MAXATOMA)

c-----------------------------------------------------------------------
c Check that the two atoms numbers are between 1 and n inclusive
c-----------------------------------------------------------------------

      if((cut1.lt.1).or.(cut1.gt.n).or.(cut2.lt.1).or.(cut2.gt.n)) then
        write(6,*)'atom label out of bounds in bond-based drive'
        stop
      endif

c-----------------------------------------------------------------------
c Check that the two atoms to be cut are actually connected
c-----------------------------------------------------------------------

      do i=1,n12(cut1)
        if(i12(i,cut1).eq.cut2) goto 1
      enddo
      goto 9999

 1    do i=1,n12(cut2)
        if(i12(i,cut2).eq.cut1) goto 2
      enddo

 9999 write(6,*)
     &'WARNING: atoms to be driven are not bound to each other.'

c-----------------------------------------------------------------------
c Initialize the routine
c-----------------------------------------------------------------------

 2    do i=1,n
        partition(i)=0
      enddo
      partition(cut1)=1
      partition(cut2)=2

      height=1
      do i=1,n12(cut1)
        bond(1)=1

c-----------------------------------------------------------------------
c skip the bond to the other atom of the bond to be cut
c-----------------------------------------------------------------------

        if(i12(i,cut1).ne.cut2) then
          atom(1)=i12(i,cut1)
 10       partition(atom(height))=1
          if(partition(i12(bond(height),atom(height))).eq.0) then

c-----------------------------------------------------------------------
c  If the connection being checked has only one attachment point (e.g. H)
c  label it and proceed
c-----------------------------------------------------------------------

            if(n12(i12(bond(height),atom(height))).eq.1) then
              partition(i12(bond(height),atom(height)))=1

c-----------------------------------------------------------------------
c  Or, the connection being checked needs to be checked further- go up the
c  tree to it and start checking its branches
c-----------------------------------------------------------------------

            else
              bond(height+1)=1
              atom(height+1)=i12(bond(height),atom(height))
              height=height+1
              goto 10
            endif

c-----------------------------------------------------------------------
c  If the cut2 atom is found, there must be a loop to it, which
c  disallows a partitioning in this fashion.
c-----------------------------------------------------------------------

          elseif(partition(i12(bond(height),atom(height))).eq.2) then
            write(6,*) 'Bond to be split is a member of a loop!'
            stop
          endif

c-----------------------------------------------------------------------
c At this point the current attachment on 'atom' and branches up from it
c have all been checked.  Code arrives here if the attachment point
c point atom's partition was already = 1 (which occurs within a loop on
c the 1 side of the structure, or if the attachment point atom had only
c one attachment point (e.g. H))
c
c Either proceed to the next attachment point:
c-----------------------------------------------------------------------

 30       if(bond(height).lt.n12(atom(height))) then
            bond(height)=bond(height)+1
            goto 10

c-----------------------------------------------------------------------
c Or all attachment points this atom have been checked, go down tree:
c-----------------------------------------------------------------------

          elseif(height.ne.1) then
            height=height-1
            goto 30

c-----------------------------------------------------------------------
c If all attachment points of all atoms above height 1 have been checked
c we'll get here.  Proceed to other non-cut2 attachment points on cut1.
c-----------------------------------------------------------------------

          endif
        endif
      enddo

c-----------------------------------------------------------------------
c Repeat for the side of # 2
c-----------------------------------------------------------------------

      do i=1,n12(cut2)
      bond(1)=1
        if(i12(i,cut2).ne.cut1) then
          atom(1)=i12(i,cut2)
 110      partition(atom(height))=2
          if(partition(i12(bond(height),atom(height))).eq.0) then
            if(n12(i12(bond(height),atom(height))).eq.1) then
              partition(i12(bond(height),atom(height)))=2
            else
              bond(height+1)=1
              atom(height+1)=i12(bond(height),atom(height))
              height=height+1
              goto 110
            endif
          elseif(partition(i12(bond(height),atom(height))).eq.1) then
            write(6,*)'ERROR: Bond to be split loops back on itself!'
            stop
          endif
 130      if(bond(height).lt.n12(atom(height))) then
            bond(height)=bond(height)+1
            goto 110
          elseif(height.ne.1) then
            height=height-1
            goto 130
          endif
        endif
      enddo

      return
      end
