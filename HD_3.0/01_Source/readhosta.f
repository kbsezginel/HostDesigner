c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'readhosta' -- reads host A input file 
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine readhosta(adummy)
      
      implicit none
      include 'params.i'
      include 'constants.i'
      include 'hosta.i'
      include 'linkage.i'
      include 'mguest.i'
      include 'control.i'
      
      character*120 record
      character*1 dtype(10)
      character*1 split(10)
      character*1 char1
      integer i, j, k, m, adummy, ftype, mm, nn
      integer macnum, ndrives, nstep(10),di(10),dj(10),dk(10)
      integer numsym(10)
      logical kill, drivinga
      real xa0(MAXATOMA),ya0(MAXATOMA),za0(MAXATOMA)
      real dstart(10),dend(10),dstep(10),dcurr(10)

      na=0
      adummy=0

c-----------------------------------------------------------------------
c     open the file 'hosta.inp'
c-----------------------------------------------------------------------

      call filetype(hosta,ftype)
      open (unit = 25, file = hosta, status = 'old')
      
      kill = .true.
      drivinga= .true.

c-----------------------------------------------------------------------
c     read title line
c-----------------------------------------------------------------------

      call fetchrecord(25,ftype,record)
      read (record,'(a60)', err = 100, end = 100) titlea

c-----------------------------------------------------------------------
c     read atoms
c-----------------------------------------------------------------------

      call fetchrecord(25,ftype,record)
      read (record, *, err = 100, end = 5) na, lguest
      goto 10
    5 lguest=1
      if(na.le.0) goto 100

   10 if(lguest.gt.MAXGST) then
         write(6,*)'subroutine readhosta-error in ',hosta,' file.'
         write(6,*) 'guest may not be larger than ',MAXGST, 'atoms!'
         stop
      endif

      if(lguest.le.0) then
         write(6,*)'subroutine readcontrol-error in ',hosta,' file.'
         write(6,*)'The guest must have a positive number of atoms!'
         stop
      endif

      if((lguest.ne.1).and.(metshape.ne.'NONE')) then
         write(6,*)'metshape cannot be used with multi-atom guests.'
         stop
      endif

      do 15 i = 1, na

         do j=1, MAXVAL
            i12a(j,i) = 0
         enddo

         call fetchrecord(25,ftype,record)
         if(notype) then
            read (record, *, err = 100, end = 15) atlaba(i), k, xa0(i),    
     &      ya0(i), za0(i), (i12a(j,i),j=1,MAXVAL)
         else
            read (record, *, err = 100, end = 15) atlaba(i), k, xa0(i),    
     &      ya0(i), za0(i), itypea(i), (i12a(j,i),j=1,MAXVAL)
         end if
         
   15 continue

c-----------------------------------------------------------------------
c     read attachments - modified by Slava 2005
c-----------------------------------------------------------------------
c     Added a feature to read symmetrical attachment locations:
c
c     nsymattacha - the number of additional symmetrical attachments
c                   2 - for tripods, 3 for tetrapods and so on
c
c     symloosera(nsymattacha, nattacha)  - the serial numbers of 
c                                          hydrogens that will be lost 
c                                          during symmetrical overlay
c     After each normal attachment specification, user shoud provide 
c     serial numbers of symmetrically related hydrogen atoms - each on 
c     separate line. 
c-----------------------------------------------------------------------

      call fetchrecord(25,ftype,record)
      read (record, *, err = 100, end = 20) nattacha, nsymattacha
      goto 25 
   20 nsymattacha=0
      if(nattacha.le.0) goto 100
 
   25 if(nattacha.gt.MAXATTACHA) then
         write(6,*)'subroutine readhosta - error in host file ',hosta
         write(6,'(A5,I3,A20)')'Only ',MAXATTACHA,' attachments allowed'
         stop
      endif

      if(nsymattacha.gt.MAXSYMATTACHA) then
         write(6,*)'subroutine readhosta - error in host file ',hosta
         write(6,'(A5,I3,A43)')'Only ',MAXSYMATTACHA,
     &   ' additional symmetrical attachments allowed'
         stop
      endif
    
c-----------------------------------------------------------------------
c      Attachment points are read in.  If there are symmetric locations
c      their serial numbers are read and serial numbers of the atoms 
c      connected with these hydrogens are identified
c----------------------------------------------------------------------- 

      do i = 1, nattacha

         nrota(i)=0
         call fetchrecord(25,ftype,record)
         read (record, *, err = 100, end = 30) loosera(i),
     &   needlaba(i), needvala(i), nrota(i)
   30    ibondera(i)=i12a(1,loosera(i))

         if(nsymattacha.gt.0) then
            do j = 1, nsymattacha
               call fetchrecord(25,ftype,record) 
               read (record, *, err = 100, end = 100) symloosera(j,i) 
               symibondera(j,i) = i12a(1,symloosera(j,i))
            enddo
         endif

      enddo

c-----------------------------------------------------------------------
c     read drives
c-----------------------------------------------------------------------

      if(drivea) then

         call fetchrecord(25,ftype,record)
         if(record(1:1).eq.'D') then
            read (record, *, err = 100, end = 100) char1, ndrives
         else
            read (record, *, err = 100, end = 100) ndrives
         endif
      
         if((ndrives.lt.0).or.(ndrives.gt.10)) then
            write (6, 1000) hosta
            stop
         endif
 1000 FORMAT('Error in ',A20,': number of drives to run must be 0-10')

         do 40 i = 1, ndrives

            call fetchrecord(25,ftype,record)
            read (record, *, err = 100, end = 35) dtype(i), split(i),
     &      dstart(i),dend(i), dstep(i), di(i), dj(i), dk(i)
   35       continue

            if(dstart(i).gt.dend(i)) then 
               write(6,*)'ERROR: drive must start before it ends'
               stop
            elseif(abs(dstart(i)-dend(i)).le.0.00001E0) then
               nstep(i)=1
            else
               if(dstep(i).le.0.E0) then
                  write(6,*)'ERROR: drive step size must be positive'
                  stop
               endif
               nstep(i)=nint((dend(i)-dstart(i))/dstep(i))+1
               if(abs(dstart(i)+(nstep(i)-1)*dstep(i)-dend(i))
     &         .ge.(0.01)) write(6,*)
     &'Warning: drive length not an integer number of steps- Rounding.'
               dend(i)=dstart(i)+(nstep(i)-1)*dstep(i)
            endif

            if(dtype(i).ne.'B'.and.dtype(i).ne.'A'.and.dtype(i).ne.
     &         'D'.and.dtype(i).ne.'+') then
                write(6,*) 'ERROR: drive type must be B, A, D, or +'
                stop
            endif

   40    continue

c-----------------------------------------------------------------------
c     Using symmetrical drives
c
c     To use symmetrical drives the user should start a new line 
c     with '+'. The program reads the line and assigns the driving
c     type based on the info of the last non-symmetric drive 
c     To distinguish between symmetrical and nonsymmetrical drives
c     the program uses numsym(ndrives) array:
c     numsym(ndrives) = 0 - normal (nonsymmetrical drive),
c                       1- the first symmetrical drive,
c                       2- the second symmetrical drive
c                       etc., etc.,
c-----------------------------------------------------------------------

         do i=1, ndrives
            numsym(i) = 0
         enddo

         j=1
         do i=1,ndrives
            if(dtype(i).eq.'+') then
               mm = i - 1
   45          if (dtype(i).ne.dtype(mm)) then
                  numsym(i) = numsym(i) + 1
               else 
                  mm = mm - 1
                  numsym(i) = numsym(i) + 1
                  goto 45
               endif
            else
               j=j*nstep(i)
            endif
         enddo
 
         do i=1, ndrives
            if (numsym(i).ne.0) then
               dtype(i) = dtype(i-numsym(i))
            endif
         enddo   

         if(j.gt.MAXDRIVEA) then
            write(6,*)'ERROR: drive positions exceeds maximum.'
            write(6,*)MAXDRIVEA,' are allowed, ',hosta,' requests ',j
            stop
         endif

      endif

c-----------------------------------------------------------------------
c     write testdrive title line
c     This part of the code was modified to account for sym drives 
c-----------------------------------------------------------------------
 
      if(testdrive.and.drivea) then

         do i=1,ndrives
            dcurr(i)=dstart(i)
         enddo
         i=0
         k=ndrives

   50    continue
         i = i + 1 

c-----------------------------------------------------------------------
c      Drive the deepest loop. Write the necessary information out 
c-----------------------------------------------------------------------

         do mm = 1, nstep(k) - 1
            write(driveinfoa(i),55)(dcurr(j),j=1,10)
   55       format('_Drive:_',10F7.2)
            do j=9,9+7*ndrives
               if(driveinfoa(i)(j:j).eq.' ')driveinfoa(i)(j:j)='_'
            enddo
            do j=ndrives,10
               driveinfoa(i)((9+7*j):(15+7*j))='       '
            enddo
            do nn = 0, numsym(k) 
               dcurr(k-nn) = dcurr(k-nn) + dstep(k-nn)
            enddo
            i = i  + 1 
         enddo

c-----------------------------------------------------------------------
c    Write the necessary inormation out for the last deepest drive
c-----------------------------------------------------------------------

         write(driveinfoa(i),55)(dcurr(j),j=1,10)

         do j=9,9+7*ndrives
            if(driveinfoa(i)(j:j).eq.' ')driveinfoa(i)(j:j)='_'
         enddo

         do j=ndrives,10
            driveinfoa(i)((9+7*j):(15+7*j))='       '
         enddo

c---------------------------------------------------------------------
c reset deepest drive to starting position
c---------------------------------------------------------------------        

         do nn = 0, numsym(k)
            dcurr(k-nn) = dstart(k-nn)
         enddo
         j = numsym(k)
         k = k - 1 - j
        
   60    if(k.eq.0) goto 65     

c---------------------------------------------------------------------
c     normal driving (all loops but deepest)
c---------------------------------------------------------------------

         if(dcurr(k).lt.dend(k)) then
          
            do nn = 0, numsym(k)
               dcurr(k-nn) = dcurr(k-nn) + dstep(k-nn)
            enddo
            k=ndrives
            goto 50

         else 

            j = numsym(k)
            do nn = 0, numsym(k)
               dcurr(k-nn) = dstart(k-nn)
            enddo
            k = k -1 -j
            goto 60

         endif 

   65    continue

      endif

c-----------------------------------------------------------------------
c     count and sort attached atoms
c-----------------------------------------------------------------------

      do i = 1, na
         n12a(i) = 0
         do j = MAXVAL, 1, -1
            if (i12a(j,i) .ne. 0) then
               n12a(i) = j
               goto 70
            endif
         enddo
   70    continue
         call sortinteger (n12a(i), i12a(1,i))
      enddo

c-----------------------------------------------------------------------
c     drive hosta 
c     align each hosta with z along its respective attachment point, 
c     and assign connlist numbers according to connectivity proximity
c-----------------------------------------------------------------------

      if(drivea) then

         call driver(ndrives,dtype,dstart,dstep,nstep,split,na,xa0,ya0,
     &   za0,di,dj,dk,n12a,i12a,nattacha,ibondera,loosera,drivinga,
     &   ndrivea,numsym)

      else

         ndrivea=1
         k = 1
         do i=1,nattacha
            call alignzpass(na, xa0, ya0, za0, ibondera(i), loosera(i),
     &      k, xa(1,1,i), ya(1,1,i),za(1,1,i))
         enddo

      endif
      
      if(testdrive) then
         close(25)          
         return
      endif
     
c-----------------------------------------------------------------------
c     Fill iguest array
c-----------------------------------------------------------------------

      do i = 1, na - lguest
         iguesta(i)=0
      enddo

      do i = na-lguest+1, na
         iguesta(i)=1
      enddo

c-----------------------------------------------------------------------
c     Remove dummy atom(s) from host a coordinates.  This will
c     remove any number of dummy atoms, updating the serial numbers
c     other variables.
c-----------------------------------------------------------------------

   75 continue

      do i=1,na

         if((atlaba(i).eq.'DU').or.(atlaba(i).eq.'Du')) then

c-----------------------------------------------------------------------
c    first remove i from the connection arrays: (n12a and i12a)
c-----------------------------------------------------------------------

            do j=1,na
               do k=1,n12a(j)
                  if(i12a(k,j).eq.i) then
                     do m=k,n12a(j)-1
                        i12a(m,j)=i12a(m+1,j)
                     enddo
                     n12a(j)=n12a(j)-1
                  endif
                  if(i12a(k,j).gt.i) i12a(k,j)=i12a(k,j)-1
               enddo
            enddo

c-----------------------------------------------------------------------
c     decrement all the post-i atoms 
c-----------------------------------------------------------------------

            do j=i,na-1
               atlaba(j)=atlaba(j+1)
               itypea(j)=itypea(j+1)
               iguesta(j) = iguesta(j+1)
               do k=1,ndrivea
                  do m=1,nattacha
                     xa(j,k,m)=xa(j+1,k,m)
                     ya(j,k,m)=ya(j+1,k,m)
                     za(j,k,m)=za(j+1,k,m)
                  enddo
               enddo
               do k=1,n12a(j+1)
                  i12a(k,j)=i12a(k,j+1)
               enddo
               n12a(j)=n12a(j+1)
            enddo

c-----------------------------------------------------------------------
c     correct attachment points and symmetrical attachment points
c-----------------------------------------------------------------------

            do j=1,nattacha
               if(ibondera(j).gt.i) ibondera(j)=ibondera(j)-1
               if(loosera(j).gt.i) loosera(j)=loosera(j)-1
            
               do mm=1,nsymattacha
                  if(symibondera(mm,j).gt.i) symibondera(mm,j)=
     &            symibondera(mm,j)-1
                  if(symloosera(mm,j).gt.i) symloosera(mm,j)=
     &            symloosera(mm,j)-1
               enddo
            enddo

            na=na-1
            adummy=adummy+1
            goto 75

         endif

      enddo

c-----------------------------------------------------------------------
c     if an error occurs reading 'hosta.inp' then what?
c-----------------------------------------------------------------------

      kill = .false.

  100 close(25)
      
      if (kill) then
         write (6, 2000) hosta
         stop
      endif
 2000 FORMAT (/,'subroutine readhosta - Error reading file ',A20)

c-----------------------------------------------------------------------
c     assign atomic number and vdw radius based on label
c     This should be done only after removing the dummy atoms
c-----------------------------------------------------------------------

      do i = 1, na

         do j = 1, MAXAT 
            if(atlaba(i).eq.atomname(j)) then
               atnuma(i) = j
               atwta(i) = xmass(j)
               rada(i) = vdwrad(j)
               goto 200
            elseif ((j.eq.MAXAT).and.(atlaba(i).ne.atomname(j))) then
               write (6, '(a)') 'subroutine readhostA - bad atom label'
     &         ,i,atlaba(i)
               stop
            endif
         enddo
  200 continue

      enddo

c-----------------------------------------------------------------------
c     Set connlist and bond type to the each attachment point. 
c-----------------------------------------------------------------------

      do i=1,nattacha
         call connlist(na, n12a, i12a, ibondera(i), lista(1,i), 1)
         call assignclass(na,xa(1,1,i),ya(1,1,i),za(1,1,i),n12a,i12a,
     &   atlaba,iguesta,torrefa(i),ibondera(i),loosera(i),
     &   classa(i),sclassa(i),altrefa(i))
      enddo
 
      return
      end
