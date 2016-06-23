c----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'overlay' -- reads one host-guest structure and checks linkage 
c     structures from a library to superimpose two host bond vectors 
c     with two link bond vectors
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine overlay(numtot, bnumtot, numhits, totlinkused, numlink)

      implicit none
      include 'params.i'
      include 'constants.i'
      include 'control.i'
      include 'hosta.i'
      include 'linkage.i'
      include 'moleculec.i'
      include 'rotateme.i'

      integer i, j, ii, jj, kk, mm
      integer ia, ja, ka, la, ig, jg, kg, lg, ilost,jlost,klost,llost
      integer nn, nrotab
      integer ivalga,ivalgb,clga,clgb,sclga,sclgb
      integer totlinkused, numhits, numlink
      integer nfit, ifit(2,4)
      integer numtot, bnumtot
      integer jlostc, klostc, ic, lc
      integer cla1,cla2
      integer equiv_class, equiv_sclass

      real xguse(MAXATOMG),yguse(MAXATOMG),zguse(MAXATOMG)
      real xause(MAXATOMA),yause(MAXATOMA),zause(MAXATOMA)
      real dx, dy, dz, scaler
      real rmsd, delta1, delta2
      real dist, dist1a, dist2a, dist1g, dist2g
      real torqa, torqg, deltorq
      real torqc1, torqc2, tthreshh
      real rrmsd,rx(MAXATOMC),ry(MAXATOMC),rz(MAXATOMC)
      real rxg(MAXATOMG),ryg(MAXATOMG),rzg(MAXATOMG)
      real energy, elapsed, updatetime,energysave
      real diffless, diffcurr, rotendebug
 
      character*2 alabga,alabgb
           
      logical passo, failtorsion, lastpass


      tthreshh = 20.0e0/tightness
      
c-----------------------------------------------------------------------
c     CHECK FOR AN EVEN NUMBER OF ATTACHMENT PAIRS
c-----------------------------------------------------------------------

      if (nattacha.eq.0) then
      
         write(6,*)  'Error in Overlay Subroutine:'
         write(6,*)  'Host A input has no connections specified'
         stop
         
      elseif (mod(nattacha,2).ne.0) then
      
         write(6,*)  'Error in Overlay Subroutine:'
         write(6,*)  'Host A input has an odd number of connections'
         stop
         
      end if

c-----------------------------------------------------------------------
c     INITIALIZE VARIABLES 
c-----------------------------------------------------------------------
      
      rrmsd=1000.0e0
      energy=0.0e0
      energysave=energy
      write(6,*) 'Starting OVERLAY run'
      call gettime(elapsed)
      updatetime=elapsed+3.0e0

      call countlib(totlinkused)
      
      open (unit = 8, file = linklib, status = 'old')
      numlink=0    

c-----------------------------------------------------------------------
c     ENTER THE DO LOOP OVER LINKS
c-----------------------------------------------------------------------

      do 800 ii = 1, totlinkused

         call readlink(numlink,ii)

c-----------------------------------------------------------------------
c     This section provides runtime status updates.
c-----------------------------------------------------------------------

         call gettime(elapsed)
         if(elapsed.ge.updatetime) then
            write(6,'(i8,a,i8)') ii,
     &           ' links analyzed so far out of ',totlinkused
            updatetime=elapsed+5.0e0
         endif

c-----------------------------------------------------------------------
c     assign atomic number and vdw radius to link based on label
c-----------------------------------------------------------------------

         do i = 1, ng(ii)

            do j = 1, MAXAT
         
               if((i.eq.looserga).or.(i.eq.loosergb)) then
                  radg(i) = 0.0e0
                  goto 40
               endif
         
               if(atlabg(i).eq.atomname(j)) then
                  atnumg(i) = j
                  atwtg(i) = xmass(j)
                  radg(i) = vdwrad(j)
                  go to 40
               elseif ((j.eq.MAXAT).and.(atlabg(i).ne.atomname(j)))then
               write (6,*)'atom label used in library not in CONSTANTS!'
               stop
               endif

            enddo

   40    continue
         enddo

c-----------------------------------------------------------------------
c     ENTER THE DO LOOP OVER ATTACHMENT PAIRS
c-----------------------------------------------------------------------

      do 700 jj = 1, nattacha, 2

         if((nrot(ii)+nrota(jj)+nrota(jj+1)).gt.maxnrot) goto 700

         lastpass=symmetry
         nrotab=nrota(jj)+nrota(jj+1)
         jlost = loosera(jj)
         klost = loosera(jj+1)
         ia = ibondera(jj)
         ja = torrefa(jj)
         ka = torrefa(jj+1)
         la = ibondera(jj+1)        
        
c-----------------------------------------------------------------------
c     If the link does not have symmetry, we will return here to swap
c     the ends of the link.  The above variables do not vary with 
c     this switch.
c-----------------------------------------------------------------------

  900    continue

c-----------------------------------------------------------------------
c     set link-based variables (toggle as needed for switching ends of link)
c-----------------------------------------------------------------------

         if(symmetry.eqv.lastpass) then
            ivalga=ivalenga
            ivalgb=ivalengb
            alabga=atlabga
            alabgb=atlabgb
            clga=  classga
            clgb=  classgb
            sclga= sclassga
            sclgb= sclassgb
            ilost= looserga
            llost= loosergb
            ig =   torrefga
            jg =   ibonderga
            kg =   ibondergb
            lg =   torrefgb
         else
            ivalga=ivalengb
            ivalgb=ivalenga
            alabga=atlabgb
            alabgb=atlabga
            clga=  classgb
            clgb=  classga
            sclga= sclassgb
            sclgb= sclassga
            ilost= loosergb
            llost= looserga
            ig =   torrefgb
            jg =   ibondergb
            kg =   ibonderga
            lg =   torrefga
         endif

c-----------------------------------------------------------------------
c     VERIFY THE LINK (and this orientation) HAS APPROPRIATE ATOM TYPES
c-----------------------------------------------------------------------

         if(ng(ii).gt.0) then

            if(((needlaba(jj).ne.'XX').and.(needlaba(jj).ne.alabga))
     &        .or. 
     &      ((needlaba(jj+1).ne.'XX').and.(needlaba(jj+1).ne.alabgb)))
     &         goto 600

            if (((needvala(jj).ne.0).and.(needvala(jj).ne.ivalga))
     &       .or.
     &      ((needvala(jj+1).ne.0).and.(needvala(jj+1).ne.ivalgb)))
     &         goto 600


c-----------------------------------------------------------------------
c     GET THE CONNECTION LISTS FOR EACH BOND
c-----------------------------------------------------------------------

            call connlist(ng(ii), n12g(1), i12g(1,1), jg, listg, 1)
            call connlist(ng(ii), n12g(1), i12g(1,1), kg, listg0,-1)
      
         endif

c-----------------------------------------------------------------------
c     BUILD THE OVERLAY MATRIX
c-----------------------------------------------------------------------

         ifit(1,1) = ia
         ifit(1,2) = jlost
         ifit(1,3) = klost
         ifit(1,4) = la
         ifit(2,1) = ilost
         ifit(2,2) = jg
         ifit(2,3) = kg
         ifit(2,4) = llost
         
         nfit = 4
    
c-----------------------------------------------------------------------
c     loop over drives : this would be more efficient outside of the 
c     loop over links, but it's harder to do it that way and speed isn't
c     crucial in overlay
c-----------------------------------------------------------------------

      do 300 kk=1,ndrivea
 
c-----------------------------------------------------------------------
c     Load coordinates
c-----------------------------------------------------------------------
         do i = 1, na
            xause(i) = xa(i,kk,jj)
            yause(i) = ya(i,kk,jj)
            zause(i) = za(i,kk,jj)
         end do

c-----------------------------------------------------------------------
c     This block adds treatment of the null link
c-----------------------------------------------------------------------

         if(ng(ii).eq.0) then

            equiv_class = classa(jj+1)
            equiv_sclass = sclassa(jj+1)
      
            call getbondl(ia,la,n12a,i12a,atlaba,iguestc,scaler)

c-----------------------------------------------------------------------
c     move the first hydrogen on the host
c-----------------------------------------------------------------------
         
            dx = xause(jlost) - xause(ia)
            dy = yause(jlost) - yause(ia)
            dz = zause(jlost) - zause(ia)
            dist = sqrt(dx*dx+dy*dy+dz*dz)
            scaler = scaler/dist
            xause(jlost) = scaler*dx + xause(ia)
            yause(jlost) = scaler*dy + yause(ia)
            zause(jlost) = scaler*dz + zause(ia)
         
c-----------------------------------------------------------------------
c     move the second hydrogen on the host
c-----------------------------------------------------------------------

            dx = xause(klost) - xause(la)
            dy = yause(klost) - yause(la)
            dz = zause(klost) - zause(la)
            dist = sqrt(dx*dx+dy*dy+dz*dz)
            scaler = scaler/dist
            xause(klost) = scaler*dx + xause(la)
            yause(klost) = scaler*dy + yause(la)
            zause(klost) = scaler*dz + zause(la)

c-----------------------------------------------------------------------
c     get distances on the ends of bonds
c-----------------------------------------------------------------------

            dx = xause(jlost) - xause(la)
            dy = yause(jlost) - yause(la)
            dz = zause(jlost) - zause(la)
            dist1a = sqrt(dx*dx+dy*dy+dz*dz)
            dx = xause(klost) - xause(ia)
            dy = yause(klost) - yause(ia)
            dz = zause(klost) - zause(ia)
            dist2a = sqrt(dx*dx+dy*dy+dz*dz)

c-----------------------------------------------------------------------
c     make the molecule if possible
c-----------------------------------------------------------------------

            if((dist1a+dist2a).lt.toldist) then
               numtot = numtot + 1
               rmsd = sqrt((dist1a+dist2a)/2)
               call makeover(ia,jlost,klost,la,ilost,jg,kg,llost,
     &         ja,ka,ig,lg,xguse,yguse,zguse,xause,yause,zause,
     &         ii,jj,0)

               failtorsion = .false.
               energy=(1+nrota(jj)+nrota(jj+1))*enrot

               if ((classa(jj).gt.3).and.(classa(jj+1).gt.3)) then
                  call equivrotor(classa(jj+1),sclassa(jj+1))
               endif

               nrotvala=nrotval(classa(jj),sclassa(jj),classa(jj+1),
     &               sclassa(jj+1))
       
               if (nrotvala.le.4) then
                  failtorsion = .true.

                  call torsionval(xc,yc,zc,torrefaca,ibonderaca,
     &            ibonderacb,torrefacb,torqc1)
                  torqc1 = torqc1*180.0e0/PI
                  diffless = 360.0e0
                  mm = 0
                  do nn = 1, nrotvala
                     diffcurr = abs(torqc1-rotval(classa(jj),
     &               sclassa(jj),classa(jj+1),sclassa(jj+1),nn))
                     if (diffcurr.le.diffless) then
                        diffless = diffcurr 
                        mm = nn
                     endif
                  enddo

                  if (diffless.le.tthreshh) then
                     failtorsion=.false.
                     energy=energy+roten(classa(jj),sclassa(jj),
     &               classa(jj+1),sclassa(jj+1),mm)
                  endif
               endif

            endif

            classa(jj+1) = equiv_class
            sclassa(jj+1) = equiv_sclass

            if(failtorsion) goto 300
 
            goto 250

         endif

c-----------------------------------------------------------------------
c    We reload xg every time; the subsequent logic requires the xg bond
c    lengths to be 1.0 angstroms and writes over them a longer distance.
c    It's inefficient, but it works to simply reload every time.
c-----------------------------------------------------------------------

         do i = 1, ng(ii)
            xguse(i) = xg(i)
            yguse(i) = yg(i)
            zguse(i) = zg(i)
         end do
   
c-----------------------------------------------------------------------
c     OBTAIN THE LINK-INVARIANT HOST DISTANCE
c-----------------------------------------------------------------------

         dx = xause(ia) - xause(la)
         dy = yause(ia) - yause(la)
         dz = zause(ia) - zause(la)
 
         dist2a = sqrt(dx*dx+dy*dy+dz*dz)

c-----------------------------------------------------------------------
c     SET THE "BONDER-LOOSER" BOND LENGTHS. THESE VALUES ARE SET
c     AT 1.000 IN THE LINK LIBRARY AND ARE AT C-H DISTANCES IN THE 
c     HOST MOLECULES.
c-----------------------------------------------------------------------
         
c-----------------------------------------------------------------------
c     move the first dummy atom on the link
c-----------------------------------------------------------------------

         call getbondo(ia,jlost,n12a,i12a,atlaba,iguesta,
     &                 jg,ilost,n12g,i12g,atlabg,scaler)

         dx = xguse(ilost) - xguse(jg)
         dy = yguse(ilost) - yguse(jg)
         dz = zguse(ilost) - zguse(jg)
         
         xguse(ilost) = scaler*dx + xguse(jg)
         yguse(ilost) = scaler*dy + yguse(jg)
         zguse(ilost) = scaler*dz + zguse(jg)
         
c-----------------------------------------------------------------------
c     then move the first hydrogen on the host
c-----------------------------------------------------------------------
         
         dx = xause(jlost) - xause(ia)
         dy = yause(jlost) - yause(ia)
         dz = zause(jlost) - zause(ia)
         
         dist = sqrt(dx*dx+dy*dy+dz*dz)
         
         scaler = scaler/dist
         
         xause(jlost) = scaler*dx + xause(ia)
         yause(jlost) = scaler*dy + yause(ia)
         zause(jlost) = scaler*dz + zause(ia)
         
c-----------------------------------------------------------------------
c     we move the second dummy atom on the link
c-----------------------------------------------------------------------

         call getbondo(la,klost,n12a,i12a,atlaba,iguesta,
     &                 kg,llost,n12g,i12g,atlabg,scaler)

         dx = xguse(llost) - xguse(kg)
         dy = yguse(llost) - yguse(kg)
         dz = zguse(llost) - zguse(kg)
         
         xguse(llost) = scaler*dx + xguse(kg)
         yguse(llost) = scaler*dy + yguse(kg)
         zguse(llost) = scaler*dz + zguse(kg)
                  
c-----------------------------------------------------------------------
c     then we move the second hydrogen on the host
c-----------------------------------------------------------------------

         dx = xause(klost) - xause(la)
         dy = yause(klost) - yause(la)
         dz = zause(klost) - zause(la)
         
         dist = sqrt(dx*dx+dy*dy+dz*dz)
         
         scaler = scaler/dist
         
         xause(klost) = scaler*dx + xause(la)
         yause(klost) = scaler*dy + yause(la)
         zause(klost) = scaler*dz + zause(la)
        
c-----------------------------------------------------------------------
c     OBTAIN THE LINK DISTANCES
c-----------------------------------------------------------------------

         dist1g = distbg
   
         dx = xguse(ilost) - xguse(llost)
         dy = yguse(ilost) - yguse(llost)
         dz = zguse(ilost) - zguse(llost)

         dist2g = sqrt(dx*dx+dy*dy+dz*dz)

c-----------------------------------------------------------------------
c     COMPUTE THE LINK-DEPENDENT HOST DISTANCE
c-----------------------------------------------------------------------

         dx = xause(jlost) - xause(klost)
         dy = yause(jlost) - yause(klost)
         dz = zause(jlost) - zause(klost)

         dist1a = sqrt(dx*dx+dy*dy+dz*dz)

c-----------------------------------------------------------------------
c     COMPARE THE TWO DISTANCES
c-----------------------------------------------------------------------

         delta1 = abs(dist1a - dist1g)
         delta2 = abs(dist2a - dist2g)
         
         if((delta1.gt.toldist).or.(delta2.gt.toldist)) then
           go to 300
         end if
                                               
c-----------------------------------------------------------------------
c     COMPARE THE DIHEDRAL ANGLES, if necessary, change the chirality on
c     the link to have the correct sign.  Note that the torsion angles
c     torqa and torqg range from 0 to 360 deg.
c-----------------------------------------------------------------------

         call torsionval(xause, yause, zause, ia,jlost,klost,la,torqa)

         torqa = torqa*180.0e0/PI
         torqg = twistg
      
         deltorq = min(abs(torqa-torqg),abs(abs(torqa-torqg)-360))
         
         if (deltorq.gt.tolphi) then
      
            torqg = -torqg
            deltorq = min(abs(torqa-torqg),abs(abs(torqa-torqg)-360))
         
            if (deltorq.gt.tolphi) then
               goto 300
            endif
         
            do i = 1, ng(ii)
               zguse(i) = -zguse(i)
            enddo
                  
         end if
         
c-----------------------------------------------------------------------
c     OVERLAY
c-----------------------------------------------------------------------

         call superimpose(na,xause,yause,zause,ng(ii),xguse,yguse,
     &      zguse,nfit,ifit,rmsd)

c-----------------------------------------------------------------------
c     MAKE THE NEW MOLECULE AND BUMPCHECK IT.
c-----------------------------------------------------------------------
          
         numtot = numtot + 1

         call makeover(ia, jlost, klost, la, ilost, jg, kg, llost,
     &     ja,ka,ig,lg, xguse,yguse,zguse, xause,yause,zause,ii,jj,0)

         call bumpcheckover(passo,ii)
                    
         if(passo) then

c-----------------------------------------------------------------------
c     CHECK TORSION ABOUT THE FIRST CREATED BOND
c-----------------------------------------------------------------------

            failtorsion = .false.
            energy=confeg(ii)+((nrot(ii)+nrota(jj)+nrota(jj+1))*enrot)
            nrotvala=nrotval(classa(jj),sclassa(jj),clga,sclga)
            nrotvalb=nrotval(clgb,sclgb,classa(jj+1),sclassa(jj+1))
       
            if (classa(jj).le.14) then

               failtorsion = .true.

               call torsionval(xc, yc, zc, torrefaca, ibonderaca, 
     &         ibondergca, torrefgca, torqc1)
               torqc1 = torqc1*180.0e0/PI

               diffless = 360.0e0
               mm = 0

               do nn = 1, nrotvala

                  diffcurr = abs(torqc1-rotval(classa(jj),sclassa(jj),
     &            clga,sclga,nn))

                  if (diffcurr.le.diffless) then
                     diffless = diffcurr 
                     mm = nn
                  endif

               enddo

               if (diffless.le.tthreshh) then
                  failtorsion=.false.
                  energy=energy+roten(classa(jj),sclassa(jj),
     &            clga,sclga,mm)
               endif

            endif

            if (failtorsion) goto 300

c-----------------------------------------------------------------------
c     CHECK TORSION ABOUT THE SECOND CREATED BOND
c-----------------------------------------------------------------------
    
            if (classa(jj+1).le.14) then

               failtorsion = .true.

               call torsionval(xc, yc, zc, torrefgcb, ibondergcb, 
     &         ibonderacb, torrefacb, torqc2)
               torqc2 = torqc2*180.0e0/PI

               diffless = 360.0
               mm = 0

               do nn = 1, nrotvalb

                  diffcurr = abs(torqc2-rotval(clgb,sclgb,classa(jj+1),
     &            sclassa(jj+1),nn))

                  if (diffcurr.le.diffless) then
                     diffless = diffcurr
                     mm = nn
                  endif

               enddo

               if (diffless.le.tthreshh) then
                  failtorsion=.false.
                  energy=energy +roten(clgb,sclgb,classa(jj+1),
     &            sclassa(jj+1),mm)
               endif

            endif

            if (failtorsion) goto 300

c-----------------------------------------------------------------------
c     TEMPORARILY STORE THE HIT RMSD AND COORDINATES
c-----------------------------------------------------------------------

            if(rmsd.lt.rrmsd) then
               rrmsd=rmsd
               energysave=energy
               do i=1,nc
                  rx(i)=xc(i)
                  ry(i)=yc(i)
                  rz(i)=zc(i)
               enddo

c-----------------------------------------------------------------------
c     need to save the current coordinates of the link because they
c     can be inverted from their original state during the drives
c-----------------------------------------------------------------------

               if(nsymattacha.gt.0) then
                  do i = 1, ng(ii)
                     rxg(i) = xguse(i)
                     ryg(i) = yguse(i)
                     rzg(i) = zguse(i)
                  enddo
               endif

            endif

         endif

c-----------------------------------------------------------------------
c     handling the null link case here
c-----------------------------------------------------------------------

  250 continue
      if(ng(ii).eq.0) then
         if(rmsd.lt.rrmsd) then
            rrmsd=rmsd
            energysave=energy
            do i=1,nc
               rx(i)=xc(i)
               ry(i)=yc(i)
               rz(i)=zc(i)
            enddo
         endif
      endif

c-----------------------------------------------------------------------
c     BOTTOM OF THE DO LOOP OVER DRIVES
c-----------------------------------------------------------------------

  300 continue

      if(rrmsd.le.999.9e0) then

c-----------------------------------------------------------------------
c     IF THERE ARE SYMMETRIC ATTACHMENTS FOR LINKS, ADD THEM NOW
c     (feature added by Slava 2006)
c-----------------------------------------------------------------------

         if(nsymattacha.gt.0) then

c-----------------------------------------------------------------------
c     PUT THE BEST COORDINATES INTO MOLECULE C    
c-----------------------------------------------------------------------
             
            do i = 1, nc
               xc(i) = rx(i) 
               yc(i) = ry(i)
               zc(i) = rz(i)  
            enddo

            do mm = 1, nsymattacha 

c-----------------------------------------------------------------------
c     UPDATE ATTACHMENT INFORMATION     
c-----------------------------------------------------------------------

               jlostc = symlooserc(mm,jj)
               klostc = symlooserc(mm,jj+1)
               ic = symibonderc(mm,jj)
               lc = symibonderc(mm,jj+1)

c-----------------------------------------------------------------------
c     UPDATE ATTACHMENT INFORMATION     
c-----------------------------------------------------------------------

               if(ng(ii).eq.0) then

                  call makeover(ic,jlostc,klostc,lc,ilost,jg,kg,llost,
     &            0,0,0,0, rxg,ryg,rzg, xc,yc,zc, ii,jj,mm)

c-----------------------------------------------------------------------
c      ADJUST C-H BOND LENGTHS FOR ATTACHMENTS TO MOLECULE C       
c-----------------------------------------------------------------------

               else

                  call getbondo(ic,jlostc,n12c,i12c,atlabc,iguestc,
     &                          jg,ilost,n12g,i12g,atlabg,scaler)

                  dx = xc(jlostc) - xc(ic)
                  dy = yc(jlostc) - yc(ic)
                  dz = zc(jlostc) - zc(ic) 
                  dist = sqrt(dx*dx + dy*dy + dz*dz)
                  scaler = scaler/dist
                  xc(jlostc) = scaler*dx + xc(ic)
                  yc(jlostc) = scaler*dy + yc(ic)
                  zc(jlostc) = scaler*dz + zc(ic)      


                  call getbondo(lc,klostc,n12c,i12c,atlabc,iguestc,
     &                          kg,llost,n12g,i12g,atlabg,scaler)

                  dx = xc(klostc) - xc(lc)
                  dy = yc(klostc) - yc(lc)
                  dz = zc(klostc) - zc(lc) 
                  dist = sqrt(dx*dx + dy*dy + dz*dz)
                  scaler = scaler/dist
                  xc(klostc) = scaler*dx + xc(lc)
                  yc(klostc) = scaler*dy + yc(lc)
                  zc(klostc) = scaler*dz + zc(lc)      

c-----------------------------------------------------------------------
c     OVERLAY 
c-----------------------------------------------------------------------

                  ifit(1,1) = ic
                  ifit(1,2) = jlostc
                  ifit(1,3) = klostc
                  ifit(1,4) = lc
                  ifit(2,1) = ilost
                  ifit(2,2) = jg
                  ifit(2,3) = kg
                  ifit(2,4) = llost
                  nfit = 4

                  call superimpose(nc,xc,yc,zc,ng(ii),rxg,
     &            ryg,rzg,nfit,ifit,rmsd) 

                  call connlist(nc, n12c, i12c, ic, listc, 1)
                  call connlist(nc, n12c, i12c, lc, listc0, -1)

                  call makeover(ic,jlostc,klostc,lc,ilost,jg,kg,llost,
     &            0,0,0,0, rxg,ryg,rzg, xc,yc,zc, ii,jj,mm)

c-----------------------------------------------------------------------
c     BUMPCHECK AT END OF NSYMATTACHA LOOP
c-----------------------------------------------------------------------

                  if(mm.eq.nsymattacha) then 
                     call bumpcheckover(passo,ii) 
                     if(.not.passo) goto 500
                  endif

               endif

               do i = 1, nc
                  rx(i) = xc(i)
                  ry(i) = yc(i)
                  rz(i) = zc(i)  
               enddo

            enddo

         endif

c-----------------------------------------------------------------------
c     STORE BEST HIT FROM DRIVES
c-----------------------------------------------------------------------

         if(rmsd.lt.maxrmsd) then
            call hitcheck(nc,rx,ry,rz,atwtc,itypec,iguestc,atlabc,ii,
     &        n12c,i12c,rrmsd,numhits,0,energysave,nrotab)
         endif

  500    rrmsd=1000.0e0

      endif

c-----------------------------------------------------------------------
c     IF LINK NOT SYMMETRIC, THEN TRY OTHER LINK ATTACHMENT
c-----------------------------------------------------------------------

  600 continue

      if(.not.lastpass) then
         lastpass = .true.
         go to 900 
      endif

c-----------------------------------------------------------------------
c     BOTTOM OF THE DO LOOP OVER ATTACHMENT PAIRS
c-----------------------------------------------------------------------

  700 continue  

c-----------------------------------------------------------------------
c     BOTTOM OF THE DO LOOP OVER LINKS
c-----------------------------------------------------------------------
                 
      if(numtot.ge.1000000000)then
         numtot=numtot-1000000000
         bnumtot=bnumtot+1
      endif

  800 continue
  
      close(8)

c-----------------------------------------------------------------------
c     Back to main
c-----------------------------------------------------------------------
     
      return
      end
