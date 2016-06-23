c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'fetchrecord' - reads up to a 120 character line from any type of
c     text file
c
c     passed variables:
c     character  record       120 character string
c     integer    ftype        type of file (1=unix,2=mac,3=win)
c     integer    iunit        io port
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine fetchrecord(iunit,ftype,record)
      implicit none
      integer i, test, fgetc, iunit, numch
      integer ftype
      character ch
      character*120 record

c-----------------------------------------------------------------------
c zero out record
c-----------------------------------------------------------------------

    5 do i = 1, 120
         record(i:i) = ' '
      enddo
      numch=0
            
c-----------------------------------------------------------------------
c load the record
c-----------------------------------------------------------------------

   10 test=fgetc(iunit,ch)
 
      if(test.eq.0) then

         if(ichar(ch).ge.32) then

            numch=numch+1
            record(numch:numch) = ch
            goto 10
     
         elseif(ichar(ch).eq.10) then
            return

         elseif(ichar(ch).eq.13) then

            if(ftype.eq.2) then
               return
            else
               test=fgetc(iunit,ch)
               return
            endif
            
         endif
         
      else
      
         return
         
      endif
      
      end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'filetype' - returns an integer specifying type of text file
c
c     passed variables:
c     character  filename     name of file to test
c     integer    ftype        type of file (1=unix,2=mac,3=win)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine filetype(filename, ftype)

      implicit none
      integer test, fgetc, ftype, nunit
      character ch
      character*20 filename
      
      open(unit=25,file=filename)
      
   10 test = fgetc(25,ch)

      if(ichar(ch).ge.32) then
         goto 10
      elseif(ichar(ch).eq.10) then    
         ftype = 1                    ! we have a unix file
         close(25)
         return
      elseif(ichar(ch).eq.13) then
         test = fgetc(25,ch)
         if(ichar(ch).ne.10.or.test.eq.-1) then
            ftype = 2                 ! we have a mac file
            close(25)
            return
         else
            ftype = 3                 ! we have a windows file
            close(25)
            return
         endif
      elseif(test.eq.-1) then
         write(6,'(a)')'subroutine filetype: ERROR'
         write(6,'(a,a)')'nothing in file named ',filename 
         stop
      endif

      return
      end
