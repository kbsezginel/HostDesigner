c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'metalshape' - compares a potential hit's metal coordination
c     sphere with members of an idealized shape specified in the input,
c     rejecting poorly matched cases. tkf
c
c     passed variables:
c     integer    jk         serial number of current drive of hosta
c     integer    mn         serial number of current drive of hostb
c     logical    shapely    set to .false. if molecule fails threshhold
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine metalshape(shapely,jk,mn)

      implicit none
      include 'params.i'
      include 'hosta.i'
      include 'constants.i'
      include 'moleculed.i'
      include 'mguest.i'

      integer  i,j,ncompare,guest1,guest2,nideal,pairing(2,7)
      integer  jk,mn
      
      real testx(7),testy(7),testz(7)
      real idealx(7),idealy(7),idealz(7)
      real rmsd, needrmsd, avelen
      
      logical shapely
      
      rmsd=123456.0e0
      
c-----------------------------------------------------------------------
c Determine the serial numbers of the two host metal placeholders:
c-----------------------------------------------------------------------
      
      guest1=na-1
      guest2=nd
 
c-----------------------------------------------------------------------
c Save the coordinates of the first metal's substituents
c-----------------------------------------------------------------------

      ncompare=0
      do i=1,n12d(guest1)
        ncompare=ncompare+1
        testx(ncompare)=xd(i12d(i,guest1))
        testy(ncompare)=yd(i12d(i,guest1))
        testz(ncompare)=zd(i12d(i,guest1))
      enddo
 
c-----------------------------------------------------------------------
c Move the second set of test coordinates as if the two metal 
c placeholders were at the same point.
c-----------------------------------------------------------------------

      do i=1,n12d(guest2)
        ncompare=ncompare+1
        testx(ncompare)=xd(i12d(i,guest2))+xd(guest1)-xd(guest2)
        testy(ncompare)=yd(i12d(i,guest2))+yd(guest1)-yd(guest2)
        testz(ncompare)=zd(i12d(i,guest2))+zd(guest1)-zd(guest2)
      enddo

c-----------------------------------------------------------------------
c Find the average length of bonds to the guest, and set the needed 
c rmsd to pass (which is scaled by this length of bonds)
c bump distances are used because they are already computed, but must be
c scaled back to bond distances from the square of bump distances
c-----------------------------------------------------------------------

      avelen=0.0e0
      do i=1,gstn12a(1)
        avelen=avelen+sqrt(gstbmpa(i,jk,1))
      enddo
      do i=1,gstn12b(1)
        avelen=avelen+sqrt(gstbmpb(i,mn,1))
      enddo
      avelen=avelen/(ncompare*sqrt(scalevdw22))
      needrmsd=avelen*shapethresh
 
c-----------------------------------------------------------------------
c Define an ideal geometry to fit to:
c-----------------------------------------------------------------------

      idealx(1)=avelen
      idealy(1)=0.0e0
      idealz(1)=0.0e0
      If (metshape(1:4).eq.'TETR') then
         nideal=4
         idealx(2)=-avelen/3.0e0
         idealy(2)=avelen*sqrt(8.0e0)/3.0e0
         idealz(2)=0.0e0         
         idealx(3)=-avelen/3.0e0
         idealy(3)=-avelen*sqrt(2.0e0)/3.0e0
         idealz(3)=avelen*sqrt(6.0e0)/3.0e0        
         idealx(4)=-avelen/3.0e0
         idealy(4)=-avelen*sqrt(2.0e0)/3.0e0
         idealz(4)=-avelen*sqrt(6.0e0)/3.0e0        
      Elseif(metshape(1:4).eq.'SQPL') then
         nideal=4
         idealx(2)=0.0e0
         idealy(2)=avelen
         idealz(2)=0.0e0         
         idealx(3)=-avelen
         idealy(3)=0.0e0
         idealz(3)=0.0e0      
         idealx(4)=0.0e0
         idealy(4)=-avelen
         idealz(4)=0.0e0        
      Elseif(metshape(1:4).eq.'TBPY') then
         nideal=5
         idealx(2)=0.0e0
         idealy(2)=avelen
         idealz(2)=0.0e0
         idealx(3)=0.0e0
         idealy(3)=-avelen/2.0e0
         idealz(3)=avelen*sqrt(3.0e0)/2.0e0
         idealx(4)=0.0e0
         idealy(4)=-avelen/2.0e0
         idealz(4)=-avelen*sqrt(3.0e0)/2.0e0
         idealx(5)=-avelen
         idealy(5)=0.0e0
         idealz(5)=0.0e0
      Elseif(metshape(1:4).eq.'SQPY') then
c-----------------------------------------------------------------------
c   This is not quite a fixed shape, there is a degree of freedom.
c   This angle used is arbitrary, equal to a ~96. angle to the apex.
c-----------------------------------------------------------------------
         nideal=5
         idealx(2)=-avelen/10.0e0
         idealy(2)=avelen*sqrt(99.0e0)/10.0e0
         idealz(2)=0.0e0
         idealx(3)=-avelen/10.0e0
         idealy(3)=0.0e0
         idealz(3)=avelen*sqrt(99.0e0)/10.0e0
         idealx(4)=-avelen/10.0e0
         idealy(4)=-avelen*sqrt(99.0e0)/10.0e0
         idealz(4)=0.0e0
         idealx(5)=-avelen/10.0e0
         idealy(5)=0.0e0
         idealz(5)=-avelen*sqrt(99.0e0)/10.0e0
      Elseif(metshape(1:4).eq.'OCTA') then
         nideal=6
         idealx(2)=0.0e0
         idealy(2)=avelen
         idealz(2)=0.0e0
         idealx(3)=-avelen
         idealy(3)=0.0e0
         idealz(3)=0.0e0
         idealx(4)=0.0e0
         idealy(4)=-avelen
         idealz(4)=0.0e0
         idealx(5)=0.0e0
         idealy(5)=0.0e0
         idealz(5)=avelen
         idealx(6)=0.0e0
         idealy(6)=0.0e0
         idealz(6)=-avelen
      Else
         write(6,*) 'Error: metal geometry flag not recognized.'
         write(6,*) 'Use one of NONE,TETR,SQPL,SQPY,TBPY, or OCTA'
         goto 9999
      Endif

      if(ncompare.gt.nideal) then
         write(6,*) 'Error: too many metal connections for ',metshape
         goto 9999
      endif

c-----------------------------------------------------------------------
c  Prepare the test matrix for the superimposition of ideal and tested,
c  and the initial ideal matrix
c-----------------------------------------------------------------------

      do i=1,nideal
         pairing(1,i)=i
         pairing(2,i)=i
      enddo

c-----------------------------------------------------------------------
c  Now the actual comparison.   Simplest way is to iterate through all
c  pairings of ideal and tested, but symmetry enables us to skip some.
c  This, however, requires specific coding for each shape.
c  We must try every arrangement as if there were 6 different atoms-
c  so, for a full octahedron, of the 6!=720 possible ways to place six
c  atoms, we eliminate a factor of 24 of these for rotations of cube,
c  leaving us only 30 possibilities to iterate through.
c-----------------------------------------------------------------------

      If (metshape(1:4).eq.'TETR') then
        call superimpose(nideal,idealx,idealy,idealz,ncompare
     &      ,testx,testy,testz,ncompare,pairing,rmsd)
        if(rmsd.lt.needrmsd) goto 1000
        if(ncompare.eq.4) then
          pairing(1,3)=4
          pairing(1,4)=3
          call superimpose(nideal,idealx,idealy,idealz,ncompare
     &      ,testx,testy,testz,ncompare,pairing,rmsd)
          if(rmsd.lt.needrmsd) goto 1000
        endif

c-----------------------------------------------------------------------

      Elseif (metshape(1:4).eq.'SQPL') then
        call superimpose(nideal,idealx,idealy,idealz,ncompare
     &      ,testx,testy,testz,ncompare,pairing,rmsd)
        if(rmsd.lt.needrmsd) goto 1000
        if(ncompare.ge.3) then
          do 100 i=1,2
            pairing(1,2)=i+2
            pairing(1,3)=mod((i+1),3)+2
            pairing(1,4)=mod((i+2),3)+2
            call superimpose(nideal,idealx,idealy,idealz,ncompare
     &        ,testx,testy,testz,ncompare,pairing,rmsd)
            if(rmsd.lt.needrmsd) goto 1000
 100      continue
        endif

c-----------------------------------------------------------------------

      Elseif (metshape(1:4).eq.'TBPY') then
        if(ncompare.eq.2) then
           call superimpose(nideal,idealx,idealy,idealz,ncompare
     &      ,testx,testy,testz,ncompare,pairing,rmsd)
          if(rmsd.lt.needrmsd) goto 1000
        else    
          do 200 j=2,0,-1
            pairing(1,1)=j+1
            pairing(1,2)=mod((j+1),3)+1
            do 250 i=2,0,-1
              pairing(1,3+i)=mod((j+2),3)+1
              if(ncompare.ge.4) then
                pairing(1,3+mod((i+1),3))=5
                pairing(1,3+mod((i+2),3))=4
                call superimpose(nideal,idealx,idealy,idealz,ncompare
     &            ,testx,testy,testz,ncompare,pairing,rmsd)
                if(rmsd.lt.needrmsd) goto 1000
              endif
              pairing(1,3+mod((i+1),3))=4
              pairing(1,3+mod((i+2),3))=5
              call superimpose(nideal,idealx,idealy,idealz,ncompare
     &          ,testx,testy,testz,ncompare,pairing,rmsd)
              if(rmsd.lt.needrmsd) goto 1000
 250        continue
 200      continue
        endif

        pairing(1,1)=4
        pairing(1,4)=1
        call superimpose(nideal,idealx,idealy,idealz,ncompare
     &    ,testx,testy,testz,ncompare,pairing,rmsd)
        if(rmsd.lt.needrmsd) goto 1000
        if(ncompare.ge.4) then
          pairing(1,2)=3
          pairing(1,3)=2
          call superimpose(nideal,idealx,idealy,idealz,ncompare
     &      ,testx,testy,testz,ncompare,pairing,rmsd)
          if(rmsd.lt.needrmsd) goto 1000
        endif

c-----------------------------------------------------------------------
        
      Elseif (metshape(1:4).eq.'SQPY') then
        if(ncompare.eq.2) then
          call superimpose(nideal,idealx,idealy,idealz,ncompare
     &      ,testx,testy,testz,ncompare,pairing,rmsd)
          if(rmsd.lt.needrmsd) goto 1000
          pairing(1,1)=3
          call superimpose(nideal,idealx,idealy,idealz,ncompare
     &      ,testx,testy,testz,ncompare,pairing,rmsd)        
          if(rmsd.lt.needrmsd) goto 1000
        else
          do 300 j=0,4
            pairing(1,j+1)=1
            pairing(1,1+mod((j+1),5))=2
            do 350  i=0,2
              pairing(1,1+mod((j+2),5))=3+i
              if(ncompare.ge.4) then
                pairing(1,1+mod((j+3),5))=3+mod(i+2,3)
                pairing(1,1+mod((j+4),5))=3+mod(i+1,3)
                call superimpose(nideal,idealx,idealy,idealz,ncompare
     &            ,testx,testy,testz,ncompare,pairing,rmsd)
                if(rmsd.lt.needrmsd) goto 1000
              endif
              pairing(1,1+mod((j+3),5))=3+mod(i+1,3)
              pairing(1,1+mod((j+4),5))=3+mod(i+2,3)
              call superimpose(nideal,idealx,idealy,idealz,ncompare
     &          ,testx,testy,testz,ncompare,pairing,rmsd)
              if(rmsd.lt.needrmsd) goto 1000
 350        continue
 300      continue
        endif

c-----------------------------------------------------------------------

      Elseif (metshape(1:4).eq.'OCTA') then
        if(ncompare.eq.2) then
          call superimpose(nideal,idealx,idealy,idealz,ncompare,
     &      testx,testy,testz,ncompare,pairing,rmsd)
        if(rmsd.lt.needrmsd) goto 1000
        else
          do 400 j=0,4
            pairing(1,2)=2+j
            pairing(1,3)=2+mod(j+1,5)
            if(ncompare.eq.3) then
              call superimpose(nideal,idealx,idealy,idealz,ncompare
     &          ,testx,testy,testz,ncompare,pairing,rmsd)
              if(rmsd.lt.needrmsd) goto 1000
            else
              do 450 i=0,2
                pairing(1,4+i)=2+mod(j+2,5)
                if(ncompare.ge.5) then
                  pairing(1,4+mod(i+1,3))=2+mod(j+4,5)
                  pairing(1,4+mod(i+2,3))=2+mod(j+3,5)
                  call superimpose(nideal,idealx,idealy,idealz,ncompare
     &              ,testx,testy,testz,ncompare,pairing,rmsd)
                  if(rmsd.lt.needrmsd) goto 1000
                endif
                pairing(1,4+mod(i+1,3))=2+mod(j+3,5)
                pairing(1,4+mod(i+2,3))=2+mod(j+4,5)
                call superimpose(nideal,idealx,idealy,idealz,ncompare
     &            ,testx,testy,testz,ncompare,pairing,rmsd)
                if(rmsd.lt.needrmsd) goto 1000
 450          continue
            endif
 400      continue
        endif
      Endif

c-----------------------------------------------------------------------

      shapely=.false.
      return

c-----------------------------------------------------------------------
      
 1000 shapely=.true.
      return
      
 9999 end
