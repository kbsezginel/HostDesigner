c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'readhostb' -- reads host B input file 
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine readhostb(bdummy)
      
      implicit none
      include 'params.i'
      include 'constants.i'
      include 'hostb.i'
      include 'linkage.i'
      include 'mguest.i'
      include 'control.i'
      
      character*120 record
      character*1 dtype(10)
      character*1 split(10)
      character*1 char1
      integer i, j, k, m, bdummy, ftype
      integer macnum, ndrives, nstep(10),di(10),dj(10),dk(10)
      integer numsym(10) 
      logical kill, drivinga
      real xb0(MAXATOMB),yb0(MAXATOMB),zb0(MAXATOMB)
      real dstart(10),dend(10),dstep(10),dcurr(10)

      nb=0
      bdummy=0
      
c-----------------------------------------------------------------------
c     open the file 'hostb.inp'
c-----------------------------------------------------------------------

      call filetype(hostb,ftype)
      open (unit = 25, file = hostb, status = 'old')
      
      kill = .true.
      drivinga= .false.

c-----------------------------------------------------------------------
c     read title line
c-----------------------------------------------------------------------

      call fetchrecord(25,ftype,record)
      read (record,'(a60)', err = 100, end = 100) titleb

c-----------------------------------------------------------------------
c     read atoms
c-----------------------------------------------------------------------

      call fetchrecord(25,ftype,record)
      read (record, *, err = 100, end = 5) nb, m

      if(m.ne.lguest) then
         write(6,*)'subroutine readhostb-inconsistent guest length.'
         write(6,*)'hosta and hostb must have the same guest length'
         stop
      endif

    5 continue

      do 10 i = 1, nb

         do j=1, MAXVAL
            i12b(j,i) = 0
         enddo

         call fetchrecord(25,ftype,record)
         if(notype) then
            read (record, *, err = 100, end = 10) atlabb(i), k, xb0(i), 
     &      yb0(i), zb0(i), (i12b(j,i),j=1,MAXVAL)
         else
            read (record, *, err = 100, end = 10) atlabb(i), k, xb0(i), 
     &      yb0(i), zb0(i), itypeb(i), (i12b(j,i),j=1,MAXVAL)
         end if

   10    continue
      
c-----------------------------------------------------------------------
c     read attachments
c-----------------------------------------------------------------------

      call fetchrecord(25,ftype,record)
      read (record, *, err = 100, end = 100) nattachb
      
      if(nattachb.gt.MAXATTACHB) then
         write(6,*)'subroutine readhostb-error in host file ',hostb
         write(6,'(a5,i3,a20)')'Only ',MAXATTACHB,' attachments allowed'
         stop
      endif
      
      do i = 1, nattachb
         nrotb(i)=0
         call fetchrecord(25,ftype,record)
         read (record, *, err = 100, end = 15) looserb(i),
     &   needlabb(i), needvalb(i), nrotb(i)
   15    ibonderb(i)=i12b(1,looserb(i))
      enddo

c-----------------------------------------------------------------------
c     read drives
c-----------------------------------------------------------------------

      if(driveb) then

         call fetchrecord(25,ftype,record)
         if(record(1:1).eq.'D') then
            read (record, *, err = 100, end = 100) char1, ndrives
         else
            read (record, *, err = 100, end = 100) ndrives
         endif
      
         if((ndrives.lt.0).or.(ndrives.gt.10)) then
            write (6, 1000) hostb
            stop
         endif
 1000 FORMAT('Error in ',A20,': number of drives to run must be 0-10')
      
         do 25 i = 1, ndrives

            call fetchrecord(25,ftype,record)

            read (record, *, err = 100, end = 20) dtype(i), split(i),
     &      dstart(i),dend(i), dstep(i), di(i), dj(i), dk(i)
   20       continue

            if(dstart(i).gt.dend(i)) then 
               write(6,*)'ERROR: drive must start before it ends'
               stop
            elseif(abs(dstart(i)-dend(i)).le.0.00001E0) then
               nstep(i)=1
            else
               if(dstep(i).le.0.) then
                  write(6,*)'ERROR: drive step size must be positive'
                  stop
               endif
               nstep(i)=nint((dend(i)-dstart(i))/dstep(i))+1
               if(abs(dstart(i)+(nstep(i)-1)*dstep(i)-dend(i)).ge.
     &         (0.01E0)) write(6,*) 'Warning: drive increment not an
     &         integer number of steps- Rounding.'

               dend(i)=dstart(i)+(nstep(i)-1)*dstep(i)

            endif

   25    continue

         j=1
         do i=1,ndrives
            j=j*nstep(i)
         enddo
         if(j.gt.MAXDRIVEB) then
            write(6,*)'ERROR: maximum drive positions exceeded.'
            write(6,*)MAXDRIVEB,' are allowed, ',hostb,' requests ',j
            stop
         endif

      endif

c-----------------------------------------------------------------------
c     There is no symmetrical drives for hostb
c-----------------------------------------------------------------------

      do i = 1, 10 
         numsym(i) = 0
      enddo  

c-----------------------------------------------------------------------
c     write testdrive title line
c-----------------------------------------------------------------------

      if(testdrive.and.driveb) then

         do i=1,ndrives
            dcurr(i)=dstart(i)
         enddo
         i=0
   30    k=ndrives
         i=i+1
         write(driveinfob(i),2000)(dcurr(j),j=1,10)
 2000 FORMAT('_Drive:_',10F7.2)

         do j=9,9+7*ndrives
            if(driveinfob(i)(j:j).eq.' ')driveinfob(i)(j:j)='_'
         enddo
         do j=ndrives,10
            driveinfob(i)((9+7*j):(15+7*j))='       '
         enddo
   35    continue
         if(dcurr(k).lt.dend(k)) then
            dcurr(k)=dcurr(k)+dstep(k)
            goto 30
         elseif(k.ne.1) then
            dcurr(k)=dstart(k)
            k=k-1
            goto 35
         endif

      endif

c-----------------------------------------------------------------------
c     count and sort attached atoms
c-----------------------------------------------------------------------

      do i = 1, nb
         n12b(i) = 0
         do j = MAXVAL, 1, -1
            if (i12b(j,i) .ne. 0) then
               n12b(i) = j
               goto 40
            endif
         enddo
   40    continue
         call sortinteger (n12b(i), i12b(1,i))
      enddo

c-----------------------------------------------------------------------
c     drive hostb if needed. 
c     align each hostb with z along its respective attachment point, 
c     and assign connlist numbers according to connectivity proximity
c-----------------------------------------------------------------------

      if(driveb) then

         call driver(ndrives,dtype,dstart,dstep,nstep,split,nb,xb0,yb0,
     &   zb0,di,dj,dk,n12b,i12b,nattachb,ibonderb,looserb,drivinga,
     &   ndriveb,numsym)

      else

         ndriveb=1
         k = -1
         do i=1,nattachb
            call alignzpass(nb, xb0, yb0, zb0, ibonderb(i), looserb(i),
     &      k, xb(1,1,i), yb(1,1,i), zb(1,1,i))
         enddo

      endif

      if(testdrive) then
         close(25)
         return
      endif

c-----------------------------------------------------------------------
c     Fill iguest array
c-----------------------------------------------------------------------

      do i = 1, nb - lguest
         iguestb(i)=0
      enddo

      do i = nb-lguest+1, nb
         iguestb(i)=1
      enddo

c-----------------------------------------------------------------------
c     Remove dummy atom(s) from host b coordinates.  This will
c     remove any number of dummies as it reruns itself each time it 
c     finds and removes a dummy atom. tkf
c-----------------------------------------------------------------------

   45 continue

      do i=1,nb

         if((atlabb(i).eq.'DU').or.(atlabb(i).eq.'Du')) then

            do j=1,nb
               do k=1,n12b(j)
                  if(i12b(k,j).eq.i) then
                     do m=k,n12b(j)-1
                        i12b(m,j)=i12b(m+1,j)
                     enddo
                     n12b(j)=n12b(j)-1
                  endif
                  if(i12b(k,j).gt.i) i12b(k,j)=i12b(k,j)-1
               enddo
            enddo

            do j=i,nb-1
               atlabb(j)=atlabb(j+1)
               itypeb(j)=itypeb(j+1)
               iguestb(j)=iguestb(j+1)
               do k=1,ndriveb
                  do m=1,nattachb
                     xb(j,k,m)=xb(j+1,k,m)
                     yb(j,k,m)=yb(j+1,k,m)
                     zb(j,k,m)=zb(j+1,k,m)
                  enddo
               enddo
               do k=1,n12b(j+1)
                  i12b(k,j)=i12b(k,j+1)
               enddo
               n12b(j)=n12b(j+1)
            enddo

            do j=1,nattachb
               if(ibonderb(j).gt.i) ibonderb(j)=ibonderb(j)-1
               if(looserb(j).gt.i) looserb(j)=looserb(j)-1            
            enddo

            nb=nb-1
            bdummy=bdummy+1
            goto 45

         endif

      enddo

c-----------------------------------------------------------------------
c     if an error occurs reading 'hostb.inp' then what?
c-----------------------------------------------------------------------

      kill = .false.

  100 close(25)
      
      if (kill) then
         write (6, 3000)
         stop
      endif
 3000 FORMAT (/,'subroutine readhostb - Error reading hostb input')

c-----------------------------------------------------------------------
c     assign atomic number and vdw radius based on label
c-----------------------------------------------------------------------

      do i = 1, nb
         do j = 1, MAXAT 
            if(atlabb(i).eq.atomname(j)) then
               atnumb(i) = j
               atwtb(i) = xmass(j)
               radb(i) = vdwrad(j)
               goto 50
            else if ((j.eq.MAXAT).and.(atlabb(i).ne.atomname(j))) then
               write (6, '(a)') 'subroutine readhostB - bad atom label'
               stop
            endif
         enddo
   50 continue
      enddo

c-----------------------------------------------------------------------
c     Set connlist and bond type to the each attachment point. 
c-----------------------------------------------------------------------

      do i=1,nattachb
         call connlist(nb,n12b, i12b, ibonderb(i), listb(1,i), 1)
         call assignclass(nb,xb(1,1,i),yb(1,1,i),zb(1,1,i),n12b,i12b,
     &   atlabb,iguestb,torrefb(i),ibonderb(i),looserb(i),
     &   classb(i),sclassb(i),altrefb(i))
      enddo

      return
      end
