c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'countlib' -- counts the number of useable links in frag library
c     
c     passed variables:
c     integer        totlinkused      number of links to be used
c
c     Note this number may or may not be the number of links in the 
c     library - usually it will be lower due to keyword based filters
c-----------------------------------------------------------------------     

      subroutine countlib (totlinkused)
      
      implicit none
      include 'params.i'
      include 'linkage.i'
      include 'control.i'
      
      character*80 record
      character*4 init4
      character*1 junk
      
      integer i, ii
      logical flagchiral, flagprochi, flagsym
      
      real         tconfeg

      integer      tng, tconng
      integer      tnrot, totlink, totlinkused

      logical      tsymmetry, tmirrorg

      totlink = 0

      open (unit = 8, file = linklib, status = 'old')

c-----------------------------------------------------------------------
c     Start loop over all links to be read
c-----------------------------------------------------------------------

      do 800 ii = 1, MAXLINKS

        if(useclass) then
  5        read (8, '(a4)', err = 200, end = 100) init4
           if (init4.ne.'CLAS') goto 5
        end if

 25     read (8, '(a80)', err = 200, end = 100) record

        if(record(1:4).eq.'CLAS') then
           read (8, '(a80)', err = 200, end = 100) record
        end if

        if(record(1:4).ne.'LINK') then
          close(8)
          write (6, 10) 
          stop
        endif
         
  10    format (/,'subroutine countlib - expected LINK in linklib')
  109   continue
               
        totlink = totlink + 1

c-----------------------------------------------------------------------
c     read flags
c-----------------------------------------------------------------------

        read (8, '(a1)') junk     
        read (8, '(a1)') junk

        read (8, '(a80)') record
        read (record, *, err = 200, end = 200) tconng, tsymmetry,
     &     tmirrorg, flagchiral, flagprochi, flagsym,
     &     tconfeg, tnrot

        read (8, '(a1)') junk
      
        read (8, '(a80)') record
        read (record, *, err = 200, end =200) tng

c-----------------------------------------------------------------------
c     exit if current link is longer than maxconn
c     (this assumes that the library is in ascending conng order. 
c      The current library is, and users are instructed to continue this)
c     skip link if it fails nrot test
c     skip link if it fails minconn test
c     skip link if it fails maxconfe test
c     skip link if it fails chiral, prochiral, or asymmetry tests    tkf
c-----------------------------------------------------------------------
      

        if(tconng.gt.maxconn) goto 100

        do i = 1, tng
          read (8, '(a1)') junk
        end do

        if((maxconfe.lt.tconfeg).or.(tconng.lt.minconn).or.
     &     (maxnrot.lt.tnrot)) goto 25

        if(nochiral.and.flagchiral) goto 25

        if(noprochiral.and.flagprochi) goto 25

c-----------------------------------------------------------------------
c   Next is the test for no asymmetric links.  Note that if the
c   noprochiral control flag has not been set, then also allow prochiral
c   links to be accepted as noasym group 
c-----------------------------------------------------------------------

        if(noasym) then
           if(.not.flagsym) then
              if(flagprochi) then
                 goto 800 
              else
                 goto 25
              endif
           endif   
        endif

c       if((.not.flagsym.and.noasym) goto 25    !original code setting

  800 continue

      write(6,1234) MAXLINKS
 1234 format('Can use only first ',i6,' usable links in library')

      close(8)

      return
      
c-----------------------------------------------------------------------
c     we should arrive here only if 'linklib' finished library
c-----------------------------------------------------------------------      
      
  100 continue

      totlinkused = ii - 1

      close(8)

      return

c-----------------------------------------------------------------------
c     An error in a read line in the middle of a link goes here
c-----------------------------------------------------------------------
      
  200 write (6, 210)
  210 format (/,'subroutine countlib - error reading link') 
      close(8)
      stop
      
c-----------------------------------------------------------------------

      end
