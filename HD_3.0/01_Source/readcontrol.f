c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'readcontrol' -- reads control file 
c-----------------------------------------------------------------------     

      subroutine readcontrol

      implicit none
      include 'params.i'
      include 'constants.i'
      include 'linkage.i'
      include 'hosta.i'
      include 'hostb.i'
      include 'control.i'

      integer ftype
      character*20 filename
      character*120 record, keyword
      logical andflag, flag1
      integer length, length2

c-----------------------------------------------------------------------     
c The file named control is composed of keywords separated by spaces.  
c If the keyword needs to be set to a value, it must be followed 
c immediately by an = and the value to be set.  The order the keywords 
c are given in is irrelevant.

c KEYWORDS:  
c LINK and OVER:  exactly 4 characters, either is required somewhere in input
c AND: if the line needs continuation, it must contain the word AND
c mirrora:  mirrora defaults to on, this keyword turns it off
c mirrorb:  mirrorb defaults to on, this keyword turns it off
c useclass:  useclass defaults to off, this keyword turns it on
c metshape:  allowed values are 'NONE', TETR','SQPL,'TBPY','SQPY','OCTA'.
c            Defaults to 'NONE'.
c linklib:  location of the link library, optionally with path.
c            Defaults to 'LIBRARY'.
c minconn:  minimum number of atoms in the minimal path through the link.
c            Defaults to 0
c maxconn:  maximum number of atoms in the minimal path through the link.
c            Defaults to 9
c maxconfe: maximum energy of link conformer w.r.t. global minimum conformer.
c            Defaults to 100.
c maxnrot: maximum number of rotatable bonds in link conformer
c            Defaults to 9
c numkeep:  number of structures to save in storage array.
c            Defaults to MAXHITS
c numview:  number of structures to write to output
c            Defaults to 100
c hosta:  name of the file describing the first fragment in linker, 
c         or the initial complex in overlay.
c            Defaults to 'hosta'
c hostb:  name of the file describing the second linker fragment.
c            Defaults to 'hostb'
c out:   name of the prefix for the output files.
c            Defaults to 'out'
c drivea:  logical indicating whether or not to drive hosta. 
c            Defaults to .false.
c driveb:  logical indicating whether or not to drive hostb.
c            Defaults to .false.
c testdrive:   writes out drivin input structures instead of proceeding
c            Defaults to .false.
c maxrmsd: maximum amount rmsd may be and still be stored.  
c            Defaults to 100
c nochiral: excludes all chiral links
c
c noprochiral: excludes all chiral and prochiral links
c
c noasym: excludes all asymmetric links
c 
c xyz: turns on simple xyz output file format, default is false
c
c notype: turns off reading of atom types in input files and writing
c           of atom types to output files.  Defaults to .false.
c
c tightness: used in linker to control the rejection on bad torsions
c            Defaults to 1.0, allowing variance of 20 deg
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     Set defaults 
c-----------------------------------------------------------------------

      andflag=.false.
      dowhat='NONE'
      mirrora=.true.
      mirrorb=.true.
      useclass=.false.
      metshape='NONE'
      linklib='LIBRARY'
      minconn=0
      maxconn=20
      maxnrot=20
      maxconfe=100.0e0
      numkeep=MAXHITS
      numview=100
      hosta='hosta'
      hostb='hostb'
      archive='out'
      drivea=.false.
      driveb=.false.
      testdrive=.false.
      maxrmsd=100.0e0
      nochiral=.false.
      noprochiral=.false.
      noasym=.false.
      notype=.false.
      xyz=.false.
      tightness=1.0e0

c-----------------------------------------------------------------------
c     Open and read input 'control' file
c-----------------------------------------------------------------------

      inquire(file='control', exist=flag1)
      if(flag1) then
        filename='control'
        call filetype(filename,ftype)
        open(unit=7, file=filename, status='old')
      else 
        filename='control.txt'
        inquire(file=filename, exist=flag1)
        if(flag1) then
          call filetype(filename,ftype)
          open(unit=7, file=filename, status='old')
        else
          write(6,*) 'Main program - cannot find control file.'
          write(6,*) 'neither control nor control.txt in directory'
          stop
        endif
      endif

          
c-----------------------------------------------------------------------
c     Top of loop over lines of code
c-----------------------------------------------------------------------

 10   call fetchrecord(7,ftype,record)
      length=len(record)
      
c-----------------------------------------------------------------------
c     Top of loop over keywords; the code reads each line *right* to *left*,
c     returning here after each keyword and starting on the next
c-----------------------------------------------------------------------

 20   do while (record(length:length).eq.' ')
        length=length-1
      enddo

c-----------------------------------------------------------------------
c     check if line is finished 
c     (this is only necessary if lines start with spaces)
c-----------------------------------------------------------------------

      if(length.lt.2) then
        length2=1
        goto 40
      endif

c-----------------------------------------------------------------------
c     'length' is now equal to the column number of the right hand side
c     of the keyword string, set by finding the first non-space character
c-----------------------------------------------------------------------

      length2=length
      do while (record(length2:length2).ne.' ')
        length2=length2-1
        if(length2.eq.0) goto 30
      enddo
 30   length2=length2+1

c-----------------------------------------------------------------------
c     the code then looks to the left of position 'length' for the next 
c     space in 'record', and then sets length2 to be the character to the 
c     right of that space.
c-----------------------------------------------------------------------

      keyword=record(length2:length)

c-----------------------------------------------------------------------
c     the variable 'keyword' has been set equal to a string from input,
c     defined to be whatever parts are separated by spaces.
c     the following logic determines which allowed keyword it matches, or
c     gives an error message.
c-----------------------------------------------------------------------

      if(keyword.eq.'LINK') then
        if(dowhat.eq.'OVER') then
          write(6,*) 'Main program - error in control file.'
          write(6,*) 'only one LINK or OVER keyword may be used'
          stop
        endif
        dowhat='LINK'
      elseif(keyword.eq.'OVER') then
        if(dowhat.eq.'LINK') then
          write(6,*) 'Main program - error in control file.'
          write(6,*) 'only one LINK or OVER keyword may be used'
          stop
        endif
        dowhat='OVER'
      elseif(keyword(1:6).eq.'hosta=') then
        hosta=keyword(7:26)
      elseif(keyword(1:6).eq.'hostb=') then
        hostb=keyword(7:26)
      elseif(keyword(1:4).eq.'out=') then
        archive=keyword(5:24)
      elseif(keyword.eq.'AND') then
        andflag=.true.
      elseif(keyword.eq.'mirroraoff') then
        mirrora=.false.
      elseif(keyword.eq.'mirrorboff') then
        mirrorb=.false.
      elseif(keyword.eq.'useclass') then
        useclass=.true.
      elseif(keyword.eq.'drivea') then
        drivea=.true.
      elseif(keyword.eq.'driveb') then
        driveb=.true.
      elseif(keyword.eq.'testdrive') then
        testdrive=.true.
      elseif(keyword(1:9).eq.'metshape=') then
        if(keyword(10:13).eq.'TETR') then
          metshape='TETR'
        elseif(keyword(10:13).eq.'SQPL') then
          metshape='SQPL'
        elseif(keyword(10:13).eq.'SQPY') then
          metshape='SQPY'
        elseif(keyword(10:13).eq.'TBPY') then
          metshape='TBPY'
        elseif(keyword(10:13).eq.'OCTA') then
          metshape='OCTA'
        elseif(keyword(10:13).ne.'NONE') then
          write(6,*) 'Main program - error in control file.'
          write(6,*) 'metshape= set to a keyword not recognized'
          stop
        endif
      elseif(keyword(1:8).eq.'minconn=') then
        keyword=keyword(9:120)
        read (keyword, *, err = 100, end = 100) minconn
      elseif(keyword(1:8).eq.'maxconn=') then
        keyword=keyword(9:120)
        read (keyword, *, err = 100, end = 100) maxconn
      elseif(keyword(1:8).eq.'maxnrot=') then
        keyword=keyword(9:120)
        read (keyword, *, err = 100, end = 100) maxnrot
      elseif(keyword(1:9).eq.'maxconfe=') then
        keyword=keyword(10:120)
        read (keyword, *, err = 100, end = 100) maxconfe
      elseif(keyword(1:8).eq.'linklib=') then
        linklib=keyword(9:120)
      elseif(keyword(1:8).eq.'numkeep=') then
        keyword=keyword(9:120)
        read (keyword, *, err = 100, end = 100) numkeep
      elseif(keyword(1:8).eq.'numview=') then
        keyword=keyword(9:120)
        read (keyword, *, err = 100, end = 100) numview 
      elseif(keyword(1:8).eq.'maxrmsd=') then
        keyword=keyword(9:120)
        read (keyword, *, err = 100, end = 100) maxrmsd
      elseif(keyword.eq.'nochiral') then
        nochiral=.true.
      elseif(keyword.eq.'noprochiral') then
        noprochiral=.true.
      elseif(keyword.eq.'noasym') then
        noasym=.true.
      elseif(keyword.eq.'xyz') then
        xyz=.true.
      elseif(keyword.eq.'notype') then
        notype=.true.
      elseif(keyword(1:10).eq.'tightness=') then
        keyword=keyword(11:120)
        read (keyword, *, err = 100, end = 100) tightness 

      else
        write(6,*)'subroutine readcontrol-error in control file.'
        write(6,*) record(length2:length),' keyword not recognized'
        stop
      endif

c-----------------------------------------------------------------------
c     end of loops over keywords and lines- determines whether further
c     input needs to be read in and where to start looking for it
c-----------------------------------------------------------------------

 40   if(length2.ne.1) then
        length=length2-1
        goto 20
      elseif(andflag) then
        andflag=.false.
        goto 10
      endif
      
      close(7)

      if(dowhat.eq.'NONE') then
        write(6,*)'subroutine readcontrol-error in control file.'
        write(6,*) 'Either the OVER or the LINK keyword must be used'
        stop
      endif

      inquire(file=hosta, exist=flag1)
      if(.not.flag1) then
         write(6,*)'subroutine readcontrol'
         write(6,*) 'Cannot find file ',hosta,' in current directory.'
         stop
      endif

      if(dowhat.eq.'LINK') then
        inquire(file=hostb, exist=flag1)
        if(.not.flag1) then
          write(6,*)'subroutine readcontrol.'
          write(6,*) 'Cannot find file ',hostb,' in current directory.'
          stop
        endif
      endif

      if(numview.gt.numkeep) then
         write(6,*)'subroutine readcontrol-error in control file.'
         write(6,*) 'Problem:  numview > numkeep'
         stop
      endif
      
      if(numkeep.gt.MAXHITS) then
         write(6,*)'subroutine readcontrol-error in control file.'
         write(6,*) 'Problem:  numkeep > MAXHITS'
         stop
      endif
        
      if (linklib.eq.'LIBRARY') then
        length=len(HDCONPATH)
        do while (HDCONPATH(length:length).eq.' ')
          length=length-1
        enddo
        if(HDCONPATH(length:length).eq.'/') then
          linklib= HDCONPATH(1:length) // 'LIBRARY'
        else
          linklib= HDCONPATH(1:length) // '/LIBRARY'
        endif  
        inquire(file=linklib, exist=flag1)
        if(.not.flag1) then
           write(6,*) 'Main program - error in LIBRARY read'
           write(6,*) 'looking for file ',linklib
           write(6,*) 'Make sure environment variable HD_DIR is set'
           write(6,*) 'to a directory containing a LIBRARY file'
           write(6,*) 'or specify a library file with a different name'
           stop
        endif
      endif

      goto 9999

  100 write(6,*) 'Subroutine readcontrol - error in control file.'
      write(6,*) record(length2:length),' variable not recognized'
      stop

 9999 return
      end
