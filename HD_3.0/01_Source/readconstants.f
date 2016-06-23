c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'readconstants' -- reads atom properties and dihedral angles from
c     a data file called CONSTANTS
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine readconstants
      
      implicit none
      include 'params.i'
      include 'constants.i'
     
      character*120 record
      character*120 conpath
      character*2   trash
      
      integer length, k
      integer class1,class2,sclass1,sclass2,nrots,i

      real angle, energy

      logical flag1

c-----------------------------------------------------------------------
c     open the file 'CONSTANTS'
c-----------------------------------------------------------------------

      length=len(HDCONPATH)
      do while (HDCONPATH(length:length).eq.' ')
        length=length-1
      enddo
      if(HDCONPATH(length:length).eq.'/') then
        conpath= HDCONPATH(1:length) // 'CONSTANTS'
      else
        conpath= HDCONPATH(1:length) // '/CONSTANTS'
      endif  
      inquire(file=conpath, exist=flag1)
      if(.not.flag1) then
         write(6,*) 'Main program - error in CONSTANTS read'
         write(6,*) 'do not see file ',conpath
         write(6,*) 'Make sure environment variable HD_DIR is set'
         write(6,*) 'to a directory containing a CONSTANTS file'
         stop
      endif

      open (unit = 20, file = conpath , status = 'old')
      
      natoms = 0
      
   10 read (20,'(a120)', err = 100, end = 200) record

      if (record(1:1).ne.'#') then
      
        if (record(1:1).eq.'A') then
        
           natoms = natoms + 1
           
           if(natoms.gt.MAXAT) then
           write(6,'(a)') 'subroutine readconstants - too many atoms'
           write(6,'(a)') 'increase parameter MAXAT'
           stop
           endif
           
           read(record, *, err = 100, end = 100) trash, 
     &       atomname(natoms), k, xmass(natoms), elneg(natoms),
     &       vdwrad(natoms), covrad(natoms)
       
           go to 10
           
        else if (record(1:1).eq.'D') then
        
           read (record, *, err = 100, end = 100) trash, 
     &     class1,sclass1,class2,sclass2,nrots
     
           nrotval(class1,sclass1,class2,sclass2)=nrots
           nrotval(class2,sclass2,class1,sclass1)=nrots
           
           do i=1,nrots
             read (20,*, err = 100, end = 100)angle,energy
             rotval(class1,sclass1,class2,sclass2,i)=angle
             rotval(class2,sclass2,class1,sclass1,i)=angle
             roten(class1,sclass1,class2,sclass2,i)=energy
             roten(class2,sclass2,class1,sclass1,i)=energy
           enddo
           go to 10
           
        else if (record(1:1).eq.'T') then
                
           read (record, *, err = 100, end = 100) trash, scalevdw12, 
     &     HHpass2, scalevdw22, toldist, tolphi, rmstol, moitol,
     &     shapethresh, enrot
     
c-----------------------------------------------------------------------
c     To speed execution, some constants used in bumpcheck are pre-squared
c-----------------------------------------------------------------------

           scalevdw12=scalevdw12*scalevdw12
           scalevdw22=scalevdw22*scalevdw22
           HHpass2=HHpass2*HHpass2
           go to 10
     
        else
        
          write(6,'(a)') 'subroutine readconstants - bad CONSTANTS file'
          write(6,'(a)') 'line contains unreadable initial character'
          stop
                
        end if

      end if
      
      goto 10      

c-----------------------------------------------------------------------
c     if an error occurs reading 'CONSTANTS' then what?
c-----------------------------------------------------------------------
      
  100 close(20)
      write (6, 110)
  110 format (/,'subroutine readconstants - Error reading CONSTANTS')
      stop

  200 close(20)
            
      return
      end
