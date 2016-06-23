c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'readlink' -- reads a link from the link library 
c     
c     passed variables:
c     integer        numlink      number of links read
c     integer        ii           number of links already analyzed
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine readlink (numlink, ii)
     
      implicit none
      include 'params.i'
      include 'linkage.i'
      include 'control.i'
      
      character*80 record
      character*4 init4
      character*1 junk
      
      integer i, j, k, numlink, ii
      logical flagchiral, flagprochi, flagsym
      

      if(useclass) then
  5      read (8, '(a4)', err = 200, end = 100) init4
         if (init4.ne.'CLAS') goto 5
      end if

 25   read (8, '(a80)', err = 200, end = 100) record

      if(record(1:4).eq.'CLAS') then
         read (8, '(a80)', err = 200, end = 100) record
      end if

      if(record(1:4).ne.'LINK') then
        close(8)
        write (6, 10) 
        stop
      endif

      read (record(5:12),  *, err = 200, end = 200) linknum(ii)

      do i=1,37
        linkname(ii)(i:i)= record(i+37:i+37)
      enddo

  10  format (/,'subroutine readlibrary - expected LINK in linklib')
  109 continue
               
c-----------------------------------------------------------------------
c     read attachments
c-----------------------------------------------------------------------

      read (8, '(a80)') record
      read (record, *, err = 200, end = 200) looserga, 
     & ibonderga,atlabga,ivalenga,
     & classga,sclassga,torrefga

      read (8, '(a80)') record
      read (record, *, err = 200, end = 200) loosergb,
     & ibondergb,atlabgb,ivalengb,
     & classgb,sclassgb,torrefgb

      read (8, '(a80)') record
      read (record, *, err = 200, end = 200) conng, symmetry,
     &   mirrorg, flagchiral, flagprochi, flagsym,
     &   confeg(ii), nrot(ii)

      read (8, '(a80)') record
      read (record, *, err = 200, end = 200) distbg, angag, angbg,
     &   twistg

      read (8, '(a80)') record
      read (record, *, err = 200, end =200) ng(ii)

c-----------------------------------------------------------------------
c     exit if current link is longer than maxconn
c     (this assumes that the library is in ascending conng order. 
c      The current library is, and users are instructed to continue this)
c     skip link if it fails nrot test 
c     skip link if it fails minconn test 
c     skip link if it fails maxconfe test                        tkf
c     skip link if it fails chiral, prochiral, or asymmetry tests
c-----------------------------------------------------------------------
      
      if(conng.gt.maxconn) then
         goto 100
      endif

      if((maxconfe.lt.confeg(ii)).or.(conng.lt.minconn).or.
     &   (maxnrot.lt.nrot(ii))) then
          goto 300
      endif

      if(nochiral.and.flagchiral) then
         goto 300
      endif

      if(noprochiral.and.flagprochi) then
         goto 300
      endif
      
      if(noasym) then
        if(.not.flagsym) then
           if(flagprochi) then
              goto 777
           else
              goto 300
           endif
         endif
      endif

  777 continue

c-----------------------------------------------------------------------
c     read atoms
c-----------------------------------------------------------------------

      do i = 1, ng(ii)

        do j=1, MAXVAL
          i12g(j,i) = 0
        end do

        read (8, '(a80)') record
        read (record, *, err = 200, end = 20) k, atlabg(i), 
     &      xg(i), yg(i), zg(i), itypeg(i),
     &      (i12g(j,i),j=1,MAXVAL)

   20 continue
      enddo

c-----------------------------------------------------------------------
c     count and sort attached atoms
c-----------------------------------------------------------------------

      do i = 1, ng(ii)
        n12g(i) = 0
        do j = MAXVAL, 1, -1
          if (i12g(j,i) .ne. 0) then
             n12g(i) = j
             go to 30
          end if
        end do
   30   continue

        call sortinteger (n12g(i), i12g(1,i))

      enddo

      numlink = numlink + 1

      return
      
c-----------------------------------------------------------------------
c     we should arrive here only if 'linklib' finished library
c-----------------------------------------------------------------------      
      
  100 continue

      write (6, *) 'Error in library read - no usable links'

c-----------------------------------------------------------------------      
c  The other possibility is that the looping structure is messed up-
c  the most likely cause of this error is the countlink is using parameters
c  different than here to select usable links.
c-----------------------------------------------------------------------      

      close(8)
      stop

c-----------------------------------------------------------------------
c     An error in a read line in the middle of a link goes here
c-----------------------------------------------------------------------
      
  200 write (6, 210) numlink, linknum(ii)
  210 format (/,'subroutine readlibrary - error reading link ', 2i6) 
      close(8)
      stop
      
c-----------------------------------------------------------------------
c     link is unacceptable, skip the rest of it
c-----------------------------------------------------------------------      
      
  300 continue

      numlink = numlink + 1
      do i = 1, ng(ii)
        read (8, '(a1)') junk
      end do
      goto 25

      end
