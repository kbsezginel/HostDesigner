c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'summary' -- writes a run summary file 
c
c     passed variables:
c     integer numtot       total number of structures examined
c     integer numlink      number of links read from library
c     integer numlinkused  number of links used in building
c     integer numhits      number of hits in storage array
c     real    setuptime    elapsed time to initialize program
c     real    buildtime    elapsed time of the structure building 
c     real    writetime    elapsed time to write viewing file
c     integer mmtime       elapsed time in MM runs
c     integer conftime     elapsed time doing searches
c     char    date         the date character string
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine summary (numtot, bnumtot, numlink, numlinkused,numhits, 
     & setuptime, buildtime, writetime, mmtime, conftime,
     & date, adummy, bdummy)
      
      implicit none
      include 'params.i'
      include 'control.i'
      include 'hosta.i'
      include 'hostb.i'
      
      integer  numlink, numlinkused, numhits, ndrives
      integer  mmtime, conftime
      integer  length, i, adummy, bdummy
      integer  numtot, bnumtot
      
      real  buildtime, writetime, cpuper, buildmin
      real  setuptime, totaltime
      
      character*1  foo
      character*24 date
      character*80 line

      logical flag1
      
      cpuper = buildtime/(real(numtot)+real(1000000000*bnumtot))
      buildmin = buildtime/60.0e0
      totaltime = buildtime + writetime + real(mmtime) + setuptime 
     & + conftime
      
      open (unit = 9, file = archive(1:index(archive,' ') - 1)//
     &    '.summ', status = 'unknown')

      write(9,10)
      write(9,10)
      write(9,'(a)') 'HOSTDESIGNER RUN SUMMARY'
      write(9,'(a,a)') 'Date of run:  ', date
      write(9,10)
      write(9,10)
      write(9,'(a,a)') 'Name of Host A:  ', titlea 
      write(9,10)
      if(dowhat.eq.'LINK') then
         write(9,'(a,a)') 'Name of Host B:  ', titleb 
         write(9,10)
      endif
      write(9,'(a)') 'Input control file follows:'
      write(9,'(a)') ' '

c-----------------------------------------------------------------------     
c The control file is reread and written verbatum, except that blank
c lines and spaces at the end of lines, if any, are omitted.
c-----------------------------------------------------------------------     

      inquire(file='control', exist=flag1)
      if(flag1) then
         open(unit=7, file='control', status='old')
      else 
         open(unit=7, file='control.txt', status='old')
      endif
  50  read (7,'(a80)', err = 100, end = 100) line
      write (9,'(a)') line
      goto 50
 100  write(9,'(a)') ' '
      close(7)

      if(drivea) then
         write(9,'(a)') 'Drive input for host A follows:'
         open (unit = 25, file = hosta, status = 'old')
         do i=1,na+3+adummy+nattacha*(1+nsymattacha)
            read (25,'(a1)') foo
         enddo
         read (25,'(a80)', err=110, end = 110) line
         if(line(1:1).eq.'D') then
            read (line,*) foo,ndrives
         else
            read (line,*) ndrives
         endif
         do i=1,ndrives
            read (25,'(a80)', err = 110, end = 110) line
            length=len(line)
            do while (line(length:length).eq.' ')
               length=length-1
            enddo
            write (9,'(a)') line(1:length)
          enddo
 110     write(9,'(a)') ' '
         close(25)
      endif

      if((driveb).and.((hosta.ne.hostb).or.(.not.drivea))) then
         write(9,'(a)') 'Drive input for host B follows:'
         open (unit = 26, file = hostb, status = 'old')
         do i=1,nb+3+bdummy+nattachb
            read (26,'(a1)') foo
         enddo
         read (26,'(a80)', err=120, end = 120) line
         if(line(1:1).eq.'D') then
            read (line,*) foo,ndrives
         else
            read (line,*) ndrives
         endif
         do i=1,ndrives
            read (26,'(a80)', err = 120, end = 120) line
            length=len(line)
            do while (line(length:length).eq.' ')
               length=length-1
            enddo
            write (9,'(a)') line(1:length)
         enddo
 120     write(9,'(a)') ' '
         close(26)
      endif


      write(9,10)     
      write(9,'(a,5x,i8)') '               Number of links read: ',
     &   numlink  
      write(9,'(a,5x,i8)') '               Number of links used: ',
     &   numlinkused  
      if(bnumtot.eq.0) then
         write(9,'(a,4x,i9)') 'Total number of structures examined: ',
     &   numtot 
      else
         write(9,'(a,i5,i9.9)') 'Total number of structures examined:',
     &   bnumtot, numtot 
      endif
      write(9,'(a,i13)') '        Number of structures stored: ',
     &   numhits 
      write(9,'(a,5x,F8.4)') '  Time spent (sec) initializing run: ',
     &   setuptime
      write(9,'(a,x,F12.4)') 'Time spent (min) building molecules: ',
     &   buildmin
      write(9,'(a,x,F12.4)') 'Time spent (sec) building molecules: ',
     &   buildtime
      write(9,'(a,4x,F9.7)') '     build time (sec) per structure: ',
     &   cpuper
      write(9,'(a,5x,F8.4)') ' Time spent (sec) writing molecules: ',
     &   writetime
      write(9,'(a,x,F12.4)') '                   Total time (sec): ',
     &   totaltime
      write(9,'(a,x,F12.4)') '                   Total time (min): ',
     &   totaltime/60.0e0
      write(9,10)


      if(mmtime.ne.0)
     &   write(9,'(a,x,F12.4)') '    Time spent (min) doing MM calcs: ',
     &   real(mmtime)/60.0e0
      if(conftime.ne.0)
     &   write(9,'(a,x,F12.4)') '    Time spent (min) doing searches: ',
     &   real(conftime)/60.0e0
   
   10 format('*********************************************',
     &'*********************************')
     
      close(9)
      
      return
      end
