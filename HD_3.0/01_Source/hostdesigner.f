c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions
c are met:
c
c 1. Redistributions of source code must retain the above copyright
c    notice, this list of conditions and the following disclaimer.
c
c 2. Redistributions in binary form must reproduce the above copyright
c    notice, this list of conditions and the following disclaimer in the
c    documentation and/or other materials provided with the distribution.
c
c THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
c IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
c OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
c IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
c INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
c NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
c DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
c THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
c (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
c THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c This program and accompanying subroutines are the code HostDesigner,
c Version 3.  HostDesigner is a de novo molecular builder for general
c application in supramolecular chemistry (see User's Manual).

c Development of HostDesigner has been supported in part by the Critical
c Materials Institute, an Energy Innovation Hub funded by the U.S
c Department of Energy (DOE), Office of Energy Efficiency and Renewable
c Energy, Advanced Manufacturing Office, in part by the Chemical Sciences,
c Office of Basic Energy Sciences, Office of Science, DOE, and in part by
c the Laboratory Directed Research and Development Program, Fundamental
c Science Division of the Pacific Northwest National Laboratory.
c
c The lead author of this software is Ben Hay, a computational chemist
c with over 25 years experience in the design and evaluation of organic
c molecules for varied applications in supramolecular chemistry.
c
c Code development took place during the period of 2000 - 2006 at
c Pacific Northwest National Laboratory. Contributions were made by the 
c following individuals:
c
c John Nicholas played a key role in the initial code design and the 
c development of the data structures used in the prototype code 
c named HostBuilder
c
c Tim Firman rewrote HostBuilder, produced the first working
c version of HostDesigner, and added a number of new features and 
c modifications that significantly increased the speed of execution
c
c Vyacheslav Bryantsev added the ability to use symmetric sets of 
c bonding vectors allowing the code to build symmetric tripodal and
c macrocyclic architectures
c
c Further code development took place during the period of 2013 - 2015 
c at Oak Ridge National Laboratory and the Supramolecular Design 
c Institute.  Contributions were made by:
c
c Billy Wayne McCann assisted with generalizing the code to treat
c any type of single bond, updating the algorithms for assigning 
c bond lengths and rotor types, as well as greatly expanding the 
c dihedral angle data in the CONSTANTS file
c
c Several subroutines (inertia.f, superimpose.f, center.f, jacobi.f, 
c and quatfit.f) were obtained from the Tinker website 
c dasher.wustl.edu/tinker and are provided with permission from Jay
c Ponder
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     program HostDesigner -- main
c-----------------------------------------------------------------------     

      program HostDesigner
      
      implicit none
      include 'params.i'
      include 'constants.i'
      include 'hosta.i'
      include 'hostb.i'
      include 'mguest.i'
      include 'linkage.i'
      include 'storage.i'
      include 'control.i'
      include 'moleculec.i'
      include 'moleculed.i'
      include 'rotateme.i'

      integer i,j,k
      integer numlink, numlinkused
      integer numhits
      integer adummy, bdummy
      integer mmtime, conftime
      integer numtot,bnumtot
      integer time 

      real elapsed1, elapsed2, elapsed3, elapsed4
      real buildtime, writetime, setuptime
 
      character*24 date
      character*120 cmd

      logical samehost

c-----------------------------------------------------------------------
c     Initialize time Calculation
c-----------------------------------------------------------------------

      mmtime=0
      conftime=0
      elapsed1=0.0e0
      elapsed2=0.0e0
      elapsed3=0.0e0
      elapsed4=0.0e0
      buildtime=0.0e0
      writetime=0.0e0
      setuptime=0.0e0

      call gettime(elapsed1)

c-----------------------------------------------------------------------
c     Initialize counters and variables
c-----------------------------------------------------------------------

      numhits = 0
      numtot  = 0
      bnumtot  = 0
      PI = 4.0e0 * atan(1.0e0)

c-----------------------------------------------------------------------
c     Specify path to the constants file and default library file
c-----------------------------------------------------------------------

      call getenv ('HD_DIR',HDCONPATH)

c-----------------------------------------------------------------------
c     Read the control file 
c-----------------------------------------------------------------------

      call readcontrol

c-----------------------------------------------------------------------
c     Read the CONSTANTS file
c-----------------------------------------------------------------------

      call readconstants

c-----------------------------------------------------------------------
c     Read the host file(s) and call the linker or the overlay routine
c-----------------------------------------------------------------------

      call readhosta(adummy)

      call gettime(elapsed2)
      
      if(dowhat.eq.'LINK') then
         call readhostb(bdummy)

         samehost=(hosta.eq.hostb)
         if(testdrive) then
            call writedrives(samehost)
            stop
         endif
         call guestsym
         call setbump

         call linker (numtot, bnumtot, numhits, numlinkused, numlink, 
     &                samehost) 

      else 
      
         if(testdrive) then
            call writedrives(.true.)
            stop
         endif
         call overlay(numtot, bnumtot, numhits, numlinkused, numlink)

      end if
      
c-----------------------------------------------------------------------
c     SORT THE RMSD ARRAY IF NOT FULL         
c-----------------------------------------------------------------------
      
      if(numhits.lt.numkeep) call sortreal(numhits, stormsd, pointer)

c-----------------------------------------------------------------------
c     GET TIME SPENT BUILDING MOLECULES
c-----------------------------------------------------------------------

      call gettime(elapsed3)
      
      setuptime = elapsed2 - elapsed1
      buildtime = elapsed3 - elapsed2

c-----------------------------------------------------------------------
c     MAKE THE VIEW  FILES - will write an ordered list of a 
c     maximum number of hits specified in the control file in ascending
c     order from the best hit to the worst hit kept.           
c-----------------------------------------------------------------------

      if(numhits.lt.numview) numview=numhits
      
      call makeview
      call makeview2(numhits)
      
      call gettime(elapsed4)

      writetime = elapsed4 - elapsed3

c-----------------------------------------------------------------------
c     WRITE THE SUMMARY  FILE and TIDY DIRECTORY IF NEEDED
c-----------------------------------------------------------------------

      call fdate(date) 

      call summary(numtot, bnumtot, numlink, numlinkused, numhits,
     &   setuptime, buildtime, writetime, mmtime, conftime,
     &   date, adummy, bdummy)

      write(6,*)'------------------------------------------------------'
      write(6,*)'    HostDesigner runs successfully! '
      write(6,*)'------------------------------------------------------'

c-----------------------------------------------------------------------
c     Done 
c-----------------------------------------------------------------------

      end
