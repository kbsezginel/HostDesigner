c----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'linker' -- reads two host-guest structures and connects them 
c     together using linkage structures from a library to create
c     new host structures.  The new structures are evaluated by the
c     degree of overlay between the two guests.
c-----------------------------------------------------------------------     

      subroutine linker(numtot, bnumtot, numhits, totlinkused, numlink,
     &                  samehost)

      implicit none
      include 'params.i'
      include 'constants.i'
      include 'control.i'
      include 'hosta.i'
      include 'hostb.i'
      include 'linkage.i'
      include 'moleculec.i'
      include 'moleculed.i'
      include 'mguest.i'

      integer i, j, ii, jj, jk, kk, mm, mn, nn
      integer iswap, jswap, mswap,clswap,sclswap,trefswap
      integer totlinkused, numlink
      integer nrots, numhits
      integer rcp
      integer rcptype(MAXATOMD,MAXRCP)
      integer rcpguest(MAXATOMD,MAXRCP)
      integer rcpn12(MAXATOMD,MAXRCP)
      integer rcpi12(MAXVAL,MAXATOMD,MAXRCP)
      integer numtot,bnumtot
      integer nrotab, orig_sclassa, orig_sclassb
      integer orig_torrefc, orig_idihedc, orig_ldihedd
      integer equiv_classb, equiv_sclassb

      real a,b
      real rmsd, energy, energyc
      real cdihed, ddihed
      real elapsed, updatetime
      real bondlc, bondld
      real rcprmsd(MAXRCP)
      real rcpen(MAXRCP)
      real rcpx(MAXATOMD,MAXRCP)
      real rcpy(MAXATOMD,MAXRCP)
      real rcpz(MAXATOMD,MAXRCP)
      real rcpatwt(MAXATOMD,MAXRCP)

      character*2 swapchar
      character*2 rcpatlab(MAXATOMD,MAXRCP)

      logical passc, passd, passe
      logical domirrora, domirrorb, donesymm
      logical samehost, shapely

c-----------------------------------------------------------------------
c     Initialize variables and counters 
c-----------------------------------------------------------------------

      write(6,*) ' Starting LINKER run'

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
     &      ' links analyzed so far out of',totlinkused
           updatetime=elapsed+5.0e0
        endif
     
c-----------------------------------------------------------------------
c     assign atomic number and vdw radius to link based on label
c-----------------------------------------------------------------------

        do i = 1, ng(ii)
          do j = 1, MAXAT
         
            if((i.eq.looserga).or.(i.eq.loosergb)) then
              radg(i) = 0.0e0
              go to 40
            endif
         
            if(atlabg(i).eq.atomname(j)) then
              atnumg(i) = j
              atwtg(i) = xmass(j)
              radg(i) = vdwrad(j)
              go to 40
            elseif ((j.eq.MAXAT).and.(atlabg(i).ne.atomname(j)))then
              write (6,*)'atom label',atlabg(i),
     &          'in library is not in CONSTANTS!'
              stop
            endif
          end do
   40   continue
        end do

c-----------------------------------------------------------------------
c     Top of link symmetry loop: If the link does not have symmetry, and
c     the hosts or hosts' attachment points are inequivalent, we will 
c     return to 750 to check the alternative connectivity.
c-----------------------------------------------------------------------

        donesymm=symmetry

  750   if(ng(ii).gt.0) call connlistg(ng(ii),n12g(1),i12g(1,1), 
     &        ibonderga, listg, 1)
        
c-----------------------------------------------------------------------
c     Orient link's first attachment point along the negative Z-axis
c
c     Note: this does edit the coordinates of the link (though it's only
c     translation and rotation, and it's only done a maximum of twice per
c     link)
c-----------------------------------------------------------------------

        if(ng(ii).gt.0)call alignz(ng(ii), xg(1),yg(1),zg(1),
     &              ibonderga, looserga, -1)

c-----------------------------------------------------------------------
c     ENTER THE DO LOOP FOR ATTACHMENT SITES ON HOST A 
c-----------------------------------------------------------------------

        do 700 jj = 1, nattacha

c-----------------------------------------------------------------------
c     Check if the link provides the correct atom bonding types.  Host A
c     always attaches to first attachment site on the link bonder.
c  
c     (a) the atom label must match unless host spec is XX (any atom) 
c     (b) the valence must match unless host spec is 0 (any valence).
c-----------------------------------------------------------------------

         if(((needlaba(jj).ne.'XX').and.(needlaba(jj).ne.atlabga))
     &     .or.((needvala(jj).ne.0).and.(needvala(jj).ne.ivalenga)))
     &      goto 700
         if((nrot(ii)+nrota(jj)).gt.maxnrot) goto 700
        
c-----------------------------------------------------------------------
c     initialize rcprmsd array
c-----------------------------------------------------------------------

         do i=1,8*MAXTORVAL*MAXTORVAL*MAXATTACHB
           rcprmsd(i)=1000.0e0
         enddo

c-----------------------------------------------------------------------
c     create molecule C 
c
c     remove terminal atoms
c     renumber link atoms and combine the bond lists
c     assign all molecule C variables
c-----------------------------------------------------------------------
          
         call makec(ii, jj)   

         if(ng(ii).ne.0) then
            call getbondl(jdihedc,kdihedc,n12c,i12c,atlabc,iguestc,
     &                    bondlc)
         end if

         orig_sclassa = sclassa(jj)
         orig_idihedc = idihedc
         orig_torrefc = torrefc

         call connlist(nc, n12c, i12c, ibonderc, listc, 1)

         do 675 jk=1,ndrivea
           if(ng(ii).eq.0) then
             do i=1,na
               xc(i)=xa(i,jk,jj)
               yc(i)=ya(i,jk,jj)
               zc(i)=za(i,jk,jj)
             enddo
           else 
             do i=1,loosera(jj)-1
               xc(i)=xa(i,jk,jj)
               yc(i)=ya(i,jk,jj)
               zc(i)=za(i,jk,jj)
             enddo
             do i=loosera(jj),na-1
               xc(i)=xa(i+1,jk,jj)
               yc(i)=ya(i+1,jk,jj)
               zc(i)=za(i+1,jk,jj)
             enddo
             do i=1,looserga-1
               xc(i+na-1)=xg(i)
               yc(i+na-1)=yg(i)
               zc(i+na-1)=zg(i)+bondlc
             enddo
             do i=looserga,ng(ii)-1
               xc(i+na-1)=xg(i+1)
               yc(i+na-1)=yg(i+1)
               zc(i+na-1)=zg(i+1)+bondlc
             enddo
           endif
          
c-----------------------------------------------------------------------
c     TOP OF MIRRORA LOOP
c-----------------------------------------------------------------------

          domirrora = .false.

 650      continue

c-----------------------------------------------------------------------
c     Swap the subclass for the rotor when classes are pyramidal
c     (amine, phosphine, sulfoxide) and subclass is 2 or 3
c-----------------------------------------------------------------------

          if ((classa(jj).eq.6).or.(classa(jj).eq.8).or.
     &        (classa(jj).eq.9)) then
             if (orig_sclassa.eq.2) then
                if(domirrora) then
                   sclassa(jj) = 3
                else
                   sclassa(jj) = 2
                endif
             else if (orig_sclassa.eq.3) then
                if(domirrora) then
                   sclassa(jj) = 2
                else
                   sclassa(jj) = 3
                endif
             endif
          endif
                   
c-----------------------------------------------------------------------
c     Swap the torsion reference for the rotor when classes are 
c     pyramidal (amine, phosphine) with subclass 1 or 4
c-----------------------------------------------------------------------

          if ((classa(jj).eq.6).or.(classa(jj).eq.9)) then
             if ((sclassa(jj).eq.1).or.(sclassa(jj).eq.4)) then
                if(domirrora) then
                   idihedc = alt_idihedc
                   if (nc.eq.na) then
                      torrefc = alt_torrefc
                   else
                      torrefc = orig_torrefc
                   endif
                else
                   idihedc = orig_idihedc
                endif
             endif
          endif

c-----------------------------------------------------------------------
c     Assign the number of rotamers
c-----------------------------------------------------------------------

          if(nc.eq.na) then
            nrots=1
          else
            nrots=nrotval(classa(jj),sclassa(jj),classga,sclassga)
          endif

c-----------------------------------------------------------------------
c     ENTER MOLECULE C DIHEDRAL DO LOOP
c-----------------------------------------------------------------------

         do 600 kk = 1,nrots

c-----------------------------------------------------------------------
c     read in coordinates to fill xc array
c-----------------------------------------------------------------------

            if (nc.eq.na) then
              energyc=confeg(ii)
              if(domirrora) then
                do i=1,nc
                  xc(i)=-xa(i,jk,jj)
                  yc(i)=ya(i,jk,jj)
                  zc(i)=za(i,jk,jj)
                enddo
              endif

            else

c-----------------------------------------------------------------------
c     find energy; skip loop if total energy with new bond is too high
c-----------------------------------------------------------------------

              energyc=confeg(ii)+
     &        roten(classa(jj),sclassa(jj),classga,sclassga,kk)
              if(energyc.gt.maxconfe) goto 600
              if(domirrora) then
                do i=1,loosera(jj)-1
                  xc(i)=-xa(i,jk,jj)
                  yc(i)=ya(i,jk,jj)
                  zc(i)=za(i,jk,jj)
                enddo
                do i=loosera(jj),na-1
                  xc(i)=-xa(i+1,jk,jj)
                  yc(i)=ya(i+1,jk,jj)
                  zc(i)=za(i+1,jk,jj)
                enddo
              else
                do i=1,loosera(jj)-1
                  xc(i)=xa(i,jk,jj)
                  yc(i)=ya(i,jk,jj)
                  zc(i)=za(i,jk,jj)
                enddo
                do i=loosera(jj),na-1
                  xc(i)=xa(i+1,jk,jj)
                  yc(i)=ya(i+1,jk,jj)
                  zc(i)=za(i+1,jk,jj)
                enddo
              endif
              do i=1,looserga-1
                xc(i+na-1)=xg(i)
                yc(i+na-1)=yg(i)
                zc(i+na-1)=zg(i)+bondlc
              enddo
              do i=looserga,ng(ii)-1
                xc(i+na-1)=xg(i+1)
                yc(i+na-1)=yg(i+1)
                zc(i+na-1)=zg(i+1)+bondlc
              enddo
      
              a =
     &        rotval(classa(jj),sclassa(jj),classga,sclassga,kk)  

              if(classa(jj).ne.16) then
                call torsionval(xc, yc, zc, idihedc, jdihedc,
     &                 kdihedc, ldihedc, cdihed)
                call setdihedral(na, nc, xc, yc, a, cdihed)
              endif
          
              call bumpcheckc(passc)
                    
              if(.not.passc) then
          
                 numtot = numtot + 1

                 goto 600
            
              endif

              call alignz (nc, xc, yc, zc, ibonderc, looserc, 1)

            endif

c-----------------------------------------------------------------------
c     ENTER THE DO LOOP FOR ATTACHMENT SITES ON HOST B
c
c     If Host A and Host B are different, necessary to do whole range.
c
c     If host A and Host B are identical, then:
c
c        if the link has no attachment symmetry, the loop that 
c        would otherwise toggle this is skipped, the price of which
c        is that the asymmetric columns here are different. It's faster
c        to do the full n^2 here once than do an upper diagonal twice.
c
c        if the link has attachment symmetry, there is no need to do
c        every combination of jj and mm; the upper and lower diagonals
c        are identical even when skipping the attachment symmetry.
c-----------------------------------------------------------------------

            if(samehost.and.symmetry) then
              i=jj
            else
              i=1
            endif

            do 500 mm = i, nattachb

c-----------------------------------------------------------------------
c     Check if the link provides the correct atom bonding types.  Host B
c     always attaches to second attachment site on the link bonder.
c     (a) the atom label must match unless host spec is XX (any atom) 
c     (b) the valence must match unless host spec is 0 (any valence).
c     (c) the attachment does not match the torsion of C being used
c         (used only if the label of the attachment matters to C)
c-----------------------------------------------------------------------

         if(((needlabb(mm).ne.'XX').and.(needlabb(mm).ne.atlabgb))
     &     .or.((needvalb(mm).ne.0).and.(needvalb(mm).ne.ivalengb)))
     &     goto 500

         if((nrot(ii)+nrota(jj)+nrotb(mm)).gt.maxnrot) goto 500

c-----------------------------------------------------------------------
c     create molecule D 
c
c     remove terminal atoms
c     renumber link atoms and combine the bond lists
c-----------------------------------------------------------------------

              call maked(looserb(mm), mm)
       
              orig_sclassb = sclassb(mm)
              orig_ldihedd = ldihedd 

c-----------------------------------------------------------------------
c     Set correct bond length by moving atoms of host B with respect to 
c     a fixed molecule C.  This is done after maked to avoid changing zb.
c     the xd array should someday be eliminated for speed reasons.
c-----------------------------------------------------------------------

              do i=1,looserc-1
                xd(i)=xc(i)
                yd(i)=yc(i)
                zd(i)=zc(i)
              enddo
              do i=looserc,nc-1
                xd(i)=xc(i+1)
                yd(i)=yc(i+1)
                zd(i)=zc(i+1)
              enddo

            do 475 mn=1,ndriveb
              do i=1,looserb(mm)-1
                xd(i+nc-1)=xb(i,mn,mm)
                yd(i+nc-1)=yb(i,mn,mm)
                zd(i+nc-1)=zb(i,mn,mm)
              enddo
              do i=looserb(mm),nb-1
                xd(i+nc-1)=xb(i+1,mn,mm)
                yd(i+nc-1)=yb(i+1,mn,mm)
                zd(i+nc-1)=zb(i+1,mn,mm)
              enddo

c-----------------------------------------------------------------------
c     get bond length and translate hostb along z-axis in positive dir
c       * the molecule c portion runs from 1 to nc-1
c       * the hostb portion runs from nc to nd
c-----------------------------------------------------------------------

              call getbondl(jdihedd,kdihedd,n12d,i12d,atlabd,iguestd,
     &             bondld)

              do i = nc, nd
                 zd(i) = zd(i) + bondld
              enddo

c-----------------------------------------------------------------------
c     TOP OF MIRRORB LOOP
c     All of the above is independent of the coords of host b.
c     Since the molecule is aligned so that the formed bond is on the
c     z-axis a simple y -> -y transformation mirrors the host B part
c     of D.  This is analogous to the mirror A loop. tkf
c-----------------------------------------------------------------------

              domirrorb = .false.

 450          continue

c-----------------------------------------------------------------------
c     Swap the subclass for the rotor when classes are pyramidal
c     (amine, phosphine, sulfoxide) and subclass is 2 or 3
c-----------------------------------------------------------------------

              if ((classb(mm).eq.6).or.(classb(mm).eq.8).or.
     &            (classb(mm).eq.9)) then
                 if (orig_sclassb.eq.2) then
                    if(domirrorb) then
                       sclassb(mm) = 3
                    else
                       sclassb(mm) = 2
                    endif
                 else if (orig_sclassb.eq.3) then
                    if(domirrorb) then
                       sclassb(mm) = 2
                    else
                       sclassb(mm) = 3
                    endif
                 endif
              endif

c-----------------------------------------------------------------------
c     Swap the torsion reference for the rotor when classes are 
c     pyramidal (amine, phosphine) with subclass 1 or 4
c-----------------------------------------------------------------------

              if ((classb(mm).eq.6).or.(classb(mm).eq.9)) then
                 if ((sclassb(mm).eq.1).or.(sclassb(mm).eq.4)) then
                    if(domirrorb) then
                       ldihedd = alt_ldihedd
                    else
                       ldihedd = orig_ldihedd
                    endif
                 endif
              endif

c-----------------------------------------------------------------------
c     Rotor data in the CONSTANTS file is available for classes 1 thru 3
c     with all other classes.  However, in the event that we have the
c     NULL link (attach hosta directly to hostb) there is a possibility
c     that we will have cases where there is no data.  To get around 
c     this problem, dihedral data is temporarily assigned by 
c     equivalencing the hostb rotor to be one of the rotors in class 1
c     thru 3 that it most resembles
c-----------------------------------------------------------------------

              if (nc.eq.na) then
                 if ((classgb.gt.3).and.(classb(mm).gt.3)) then
                    equiv_classb = classb(mm)
                    equiv_sclassb = sclassb(mm)
                    call equivrotor(classb(mm),sclassb(mm))
                 endif
              end if

c-----------------------------------------------------------------------
c     ENTER MOLECULE D DIHEDRAL DO LOOP
c-----------------------------------------------------------------------

           do 400 nn=1,nrotval(classgb,sclassgb,classb(mm),sclassb(mm))

               numtot = numtot + 1

c-----------------------------------------------------------------------
c     skip loop if total energy with new bonds is too high
c     energy determined here is also stored if d is saved
c-----------------------------------------------------------------------

               energy = energyc +
     &         roten(classgb,sclassgb,classb(mm),sclassb(mm),nn)+
     &         (nrot(ii)+nrota(jj)+nrotb(mm))*enrot
               if(energy.gt.maxconfe) goto 400
                         
               b = rotval(classgb,sclassgb,classb(mm),sclassb(mm),nn)

               if(classb(mm).ne.16) then
                 call torsionval(xd, yd, zd, idihedd, jdihedd,
     &               kdihedd, ldihedd, ddihed)
                 call setdihedral(nc, nd, xd, yd, b, ddihed)
               endif

               i=na-lguest-1

               call getrmsd(i,nd,xd,yd,zd,rmsd)

               if(rmsd.lt.maxrmsd) then

                  call bumpcheckd(passd,jk,mn)

                  if(passd) then  

                    if(metshape.ne.'NONE') then
                      call metalshape(shapely,jk,mn)
                      if(.not.shapely) goto 400
                    endif

c-----------------------------------------------------------------------
c  rcp, or the recipe number, is unique for a given set of values of 
c  kk, mm, nn, domirrora, domirrorb.
c
c  rcp doesn't distinguish between different jk and mn values, therefore
c  within a single jj loop if two hits have the same rcp number they
c  were created with the same procedure with the exception of drives.
c  hitcheck will only save one instance of a given recipe number.
c
c  Different values of ii, jj, and the mirroring of the link are outside
c  this loop, and are always unique for the purposes of this check.
c
c  in this way the near duplicates differing only in the drive rotation 
c  will not be saved.  This is intended to choose the 'best' hit from 
c  all the drives with a given recipe.  tkf
c
c  if MAXTORVAL is ever changed from 6, two numbers below will need to
c  be altered (from 24 and 144 to 4*MAXTORVAL and 4*MAXTORVAL*MAXTORVAL)
c-----------------------------------------------------------------------

                    rcp=4*(kk-1)+24*(nn-1)+144*(mm-1)+1
                    if(domirrora)rcp=rcp+1
                    if(domirrorb)rcp=rcp+2

                    if(rcprmsd(rcp).gt.rmsd) then
                      rcprmsd(rcp)=rmsd
                      rcpen(rcp)=energy
                      do i = 1, nd
                        rcpx(i,rcp) = xd(i)
                        rcpy(i,rcp) = yd(i)
                        rcpz(i,rcp) = zd(i)
                        rcpatwt(i,rcp) = atwtd(i)
                        rcptype(i,rcp) = ityped(i)
                        rcpguest(i,rcp) = iguestd(i)
                        rcpatlab(i,rcp) = atlabd(i)
                        rcpn12(i,rcp) = n12d(i)
                        do j = 1, n12d(i)
                          rcpi12(j,i,rcp) = i12d(j,i)
                        end do
                      end do
                    endif             
                  endif           
                endif
                
c-----------------------------------------------------------------------
c     BOTTOM OF MOLECULE D DIHEDRAL DO LOOP (nn)
c-----------------------------------------------------------------------
      
  400         continue

c-----------------------------------------------------------------------
c     If class and subclass were changed as a result of the NULL link
c     and missing rotor values, reset the class and subclass for hostb
c-----------------------------------------------------------------------

              if (nc.eq.na) then
                 if ((classgb.gt.3).and.(equiv_classb.gt.3)) then
                    classb(mm) = equiv_classb
                    sclassb(mm) = equiv_sclassb
                 endif
              end if
  
c-----------------------------------------------------------------------
c     BOTTOM OF MIRROR HOST B LOOP
c-----------------------------------------------------------------------

              if(mirrorb.neqv.domirrorb) then
                 do i = nc, nd
                    yd(i) = -yd(i)
                 enddo
                 domirrorb=.true.
                 goto 450
              endif

c Here is what was there, seems like too much work since all we need to
c do is invert the sign on the hostb portion of the molecule
c               do i=1,looserb(mm)-1
c                 xd(i+nc-1)=xb(i,mn,mm)
c                 yd(i+nc-1)=-yb(i,mn,mm)
c                 zd(i+nc-1)=zb(i,mn,mm)
c               enddo
c               do i=looserb(mm),nb-1
c                 xd(i+nc-1)=xb(i+1,mn,mm)
c                 yd(i+nc-1)=-yb(i+1,mn,mm)
c                 zd(i+nc-1)=zb(i+1,mn,mm)
c               enddo
c               call setbondlength(itypeb(ibonderb(mm)),
c    &            itypec(ibonderc),nb-1,zd(nc))

c-----------------------------------------------------------------------
c     BOTTOM OF THE DRIVE HOST B LOOP  (mn)
c-----------------------------------------------------------------------
 
  475       continue

c-----------------------------------------------------------------------
c     BOTTOM OF ATTACHMENTS ON HOST B DO LOOP (mm)
c-----------------------------------------------------------------------
      
  500       continue

c-----------------------------------------------------------------------
c     BOTTOM OF MOLECULE C DIHEDRAL DO LOOP (kk)
c-----------------------------------------------------------------------
      
  600     continue
 
c-----------------------------------------------------------------------
c     BOTTOM OF MIRROR HOST A LOOP
c-----------------------------------------------------------------------

          if((mirrora.neqv.domirrora).and.(ng(ii).ne.0)) then
            domirrora=.true.
            goto 650
          endif

c-----------------------------------------------------------------------
c     BOTTOM OF THE DRIVE HOST A LOOP (jk)
c-----------------------------------------------------------------------
 
  675    continue

c-----------------------------------------------------------------------
c     BOTTOM OF ATTACHMENTS ON HOST A DO LOOP (jj)
c
c     store (if apropos) the hits generated during this loop. 
c     this holds the link and the attachment on host A constant,
c     but includes all other variations.
c
c     if MAXTORVAL is ever changed from 6, the 144s below will need
c     to be altered (to 4*MAXTORVAL*MAXTORVAL)
c-----------------------------------------------------------------------

         do i=1,144*nattachb
           if(rcprmsd(i).le.maxrmsd) then

c-----------------------------------------------------------------------
c     check is needed again as new, lower rmsd hits can be saved during
c     this loop that would knock the current hit out of top <numkeep>.
c     information such as mm about the drive can be extracted from rcp:
c             foo=(rcpx(nd,i)-rcpx(na-1,i))*(rcpx(nd,i)-rcpx(na-1,i))+
c     &           (rcpy(nd,i)-rcpy(na-1,i))*(rcpy(nd,i)-rcpy(na-1,i))+
c     &           (rcpz(nd,i)-rcpz(na-1,i))*(rcpz(nd,i)-rcpz(na-1,i))
c             foo=sqrt(foo)-.0001
c             if(foo.gt.rcprmsd(i)) then
c                write(6,*)foo,rcprmsd(i)
c             endif
c-----------------------------------------------------------------------

             mm=(i-mod(i,144))/144+1
             nrotab=nrota(jj)+nrotb(mm)
             call hitcheck(nd,rcpx(1,i),rcpy(1,i),rcpz(1,i),
     &          rcpatwt(1,i),rcptype(1,i),rcpguest(1,i),
     &          rcpatlab(1,i),ii,rcpn12(1,i),rcpi12(1,1,i),
     &          rcprmsd(i),numhits,0,rcpen(i),nrotab)
           endif
           rcprmsd(i)=1000.0e0
         enddo

        if(numtot.ge.1000000000) then
          numtot=numtot-1000000000
          bnumtot=bnumtot+1
        endif

  700   continue
  
c-----------------------------------------------------------------------
c     if the link does not have symmetry, then we try the alternate
c     attachment of host A and host B to the link. However, don't do if
c     host A and host B are identical.
c-----------------------------------------------------------------------

        if(.not.donesymm) then
      
           if(.not.samehost) then
      
                iswap = looserga
                jswap = ibonderga
                mswap = ivalenga
             swapchar = atlabga
               clswap = classga
              sclswap = sclassga
             trefswap = torrefga
             
             looserga = loosergb
            ibonderga = ibondergb
             ivalenga = ivalengb
              atlabga = atlabgb
              classga = classgb
             sclassga = sclassgb
             torrefga = torrefgb
         
             loosergb = iswap
            ibondergb = jswap
             ivalengb = mswap
              atlabgb = swapchar
              classgb = clswap
             sclassgb = sclswap
             torrefgb = trefswap
         
             donesymm = .true.
          
             goto 750 
          
           endif

        endif

c-----------------------------------------------------------------------
c     BOTTOM OF THE DO LOOP OVER LINKS
c-----------------------------------------------------------------------
              
  800 continue
  
      close(8)

c-----------------------------------------------------------------------
c     Back to main
c-----------------------------------------------------------------------

      return
      end
