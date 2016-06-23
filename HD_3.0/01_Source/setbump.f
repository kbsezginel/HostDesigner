c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'setbump' - sets non-standard bumping lengths to guests
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine setbump

      implicit none
      include 'params.i'
      include 'constants.i'
      include 'hosta.i'
      include 'hostb.i'
      include 'mguest.i'

      integer i,j,k,l,pa,pb,drive
      logical included
      real dx,dy,dz,dist2
      
      pa=na-lguest
      pb=nb-lguest
      
c      write (6,*) ndrivea, ndriveb
      
      do i=1,lguest
        gstn12a(i)=0
        gstn12b(i)=0
      enddo

c-----------------------------------------------------------------------     
c  loop over all bonds from atoms in host b guest to atoms outside it
c-----------------------------------------------------------------------     

      do i=1,lguest
        do j=1,n12b(i+pb)
          if(i12b(j,i+pb).lt.pb) then

c-----------------------------------------------------------------------     
c  dist2 is the scaled square of the distance; this is stored squared 
c  to save time in the bumpcheck routines (which use the majority of the
c  time.)  It is scaled down to allow closer contacts than the bonding
c  distance to the guest atoms.
c-----------------------------------------------------------------------     

c-----------------------------------------------------------------------     
c  loop over all drives
c-----------------------------------------------------------------------     

            do drive=1,ndriveb

              dx = xb(i12b(j,i+pb),drive,1)-xb(i+pb,drive,1)
              dy = yb(i12b(j,i+pb),drive,1)-yb(i+pb,drive,1) 
              dz = zb(i12b(j,i+pb),drive,1)-zb(i+pb,drive,1)
              dist2=scalevdw22*(dx*dx+dy*dy+dz*dz)

c-----------------------------------------------------------------------     
c  Now code loops over all atoms which are symmetry equivalent
c  to the atom (i+pb) in hostb, to set bump thresholds for the atoms
c  equivalently bonded to it.
c
c  Because each atom may have more than one bond to the guest, or
c  be symmetry equivalent to more than one atom in the guest, there is
c  a check to see if that atom is already on the bumping list; if so, 
c  the bump threshhold for that atom is set to the lower threshhold. 
c-----------------------------------------------------------------------     

              do k=1,nequivb
                included=.false.
                do l=1,gstn12a(mapeqb(k,i))
                  if(gsti12a(l,mapeqb(k,i)).eq.i12b(j,i+pb)) then
                    gstbmpa(l,drive,mapeqb(k,i))=
     &                         min(gstbmpa(l,drive,mapeqb(k,i)),dist2)
                    included=.true.
                  endif
                enddo
                if(.not.included) then
                  gstn12a(mapeqb(k,i))=gstn12a(mapeqb(k,i))+1
                  gsti12a(gstn12a(mapeqb(k,i)),mapeqb(k,i))=i12b(j,i+pb)
                  gstbmpa(gstn12a(mapeqb(k,i)),drive,mapeqb(k,i))=dist2
                endif
              enddo
            enddo

c-----------------------------------------------------------------------     
c  end of loop over bonds
c-----------------------------------------------------------------------     

          endif
        enddo
      enddo

c-----------------------------------------------------------------------     
c  Now do the inverse process, looping over the atoms in the guest of 
c  complex a and the symmetry equivalents in complex b, for those bumps:
c-----------------------------------------------------------------------     

      do i=1,lguest
        do j=1,n12a(i+pa)
          if(i12a(j,i+pa).lt.pa) then
            do drive=1,ndrivea
              dx = xa(i12a(j,i+pa),drive,1)-xa(i+pa,drive,1)
              dy = ya(i12a(j,i+pa),drive,1)-ya(i+pa,drive,1) 
              dz = za(i12a(j,i+pa),drive,1)-za(i+pa,drive,1)
              dist2=scalevdw22*(dx*dx+dy*dy+dz*dz)
              do k=1,nequiva
                included=.false.
                do l=1,gstn12b(mapeqa(k,i))
                  if(gsti12b(l,mapeqa(k,i)).eq.i12a(j,i+pa)) then
                    gstbmpb(l,drive,mapeqa(k,i))=
     &                    min(gstbmpb(l,drive,mapeqa(k,i)),dist2)
                    included=.true.
                  endif
                enddo
                if(.not.included) then
                  gstn12b(mapeqa(k,i))=gstn12b(mapeqa(k,i))+1
                  gsti12b(gstn12b(mapeqa(k,i)),mapeqa(k,i))=i12a(j,i+pa)
                  gstbmpb(gstn12b(mapeqa(k,i)),drive,mapeqa(k,i))=dist2
                endif
              enddo
            enddo
          endif
        enddo
      enddo

      return
      end
