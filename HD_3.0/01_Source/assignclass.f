c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'assignclass' -- defines class and subclass of a rotating group
c     and identifies the terminal atom that will be used to define the
c     torsion angles.
c
c                 LIST OF CLASSES AND SUBCLASSES
c
c     Note that unless otherwise specified Y is any non-H atom and X is
c     any atom. Refer to user manual for further information.
c
c     Class 1 - Tetrahedral
c       subclass 1:  -CH2Z
c       subclass 2:  -CHZ2
c       subclass 3:  -CX3
c     Class 2 - Trigonal Planar (alkene-like)
c       subclass 1:  -C(H)=C(H)(X) where H is cis
c       subclass 2:  -C(H)=C(Z)(X) where Z is cis
c       subclass 3:  -C(Z)=C(H)(X) where H is cis
c       subclass 4:  -C(Z)=C(Z)(X) where Z is cis
c     Class 3 - Trigonal Planar - arene
c       subclass 1:  -(arene) two ortho H
c       subclass 2:  -(arene) one ortho H and one ortho Z
c       subclass 3:  -(arene) two ortho Z
c     Class 4 - Amide Nitrogen
c       subclass 1:  -N(H)C(=O)X attaching cis to carbonyl
c       subclass 2:  -N(H)C(=O)X attaching trans to carbonyl
c       subclass 3:  -N(Z)C(=O)X attaching cis to carbonyl
c       subclass 4:  -N(Z)C(=O)X attaching trans to carbonyl
c     Class 5 - Ether-like
c       subclass 1:  -OZ
c       subclass 2:  -OH
c       subclass 3:  -SZ 
c       subclass 4:  -SH 
c     Class 6 - Amine
c       subclass 1:  -NH2
c       subclass 2:  -NHZ
c       subclass 3:  -NHZ, enantiomer
c       subclass 4:  -NZ2
c     Class 7 - Amide Carbon
c       subclass 1:  -C(=O)N(H)(X), where H is trans to carbonyl
c       subclass 2:  -C(=O)N(Z)(X), where Z is trans to carbonyl
c     Class 8 - Oxides
c       subclass 1:  Phosphine Oxide : -P(=O)Z2
c       subclass 2:  Sulfoxide : -S(=O)X
c       subclass 3:  Sulfoxide : -S(=O)X, enantiomer
c       subclass 4:  Sulfone   : -S(=O)2Z
c     Class 9 - Phosphines
c       subclass 1:  -PH2
c       subclass 2:  -PHZ
c       subclass 3:  -PHZ, enantiomer
c       subclass 4:  -PZ2 
c     Class 10 - Aldehyde and Ketone
c       subclass 1:  -C(=O)H
c       subclass 2:  -C(=O)Z
c       subclass 1:  -C(=S)H
c       subclass 2:  -C(=S)Z
c     Class 11 - Imine Nitrogen
c       subclass 1:  -N=C(H)(X), where H is cis
c       subclass 2:  -N=C(Z)(X), where Z is cis
c       subclass 3:  -N=N-X, where X is trans
c     Class 12 - Imine Carbon
c       subclass 1:  -C(Z)=N-X, where X is trans
c       subclass 2:  -C(H)=N-X, where X is trans
c       subclass 3:  -C(Z)=N-X, where X is cis
c       subclass 4:  -C(H)=N-X, where X is cis
c     Class 13 - Planar heterocycle - unsubstituted ortho position
c       subclass 1:  1 ortho unsubstituted + 1 ortho H
c       subclass 2:  2 ortho unsubstituted
c       subclass 3:  1 ortho unsubstituted + 1 ortho Z
c     Class 14 - Acid, ester, phenoxy
c       subclass 1:  -CO2 anion
c       subclass 2:  -C(=O)O-X
c       subclass 3:  -OC(=O)X
c       subclass 4:  –OPh
c     Class 15 - No torsional preference
c       subclass 1:  Linear, Hypervalent
c     Class 16 - Terminal case
c       subclass 1:  bonding atom attached to only a loser atom
c
c-----------------------------------------------------------------------
c
c     CAUTION:  Three arrays in constants.i must be updated if you 
c     add more classes.  These are as follows:
c
c     nrotval(#classes,#subclasses,#classes,#subclasses)  
c     rotval(#classes,#subclasses,#classes,#subclasses,#minima)  
c     roten(#classes,#subclasses,#classes,#subclasses,#minima)  
c
c     At this time #classes=16, #subclasses=4, and #minima=6
c
c     The number of minima is also defined in params.i as MAXTORVAL
c     Currently set at 6.
c     
c-----------------------------------------------------------------------
c
c                           META
c
c     This subroutine is called by readhosta.F (line 563) and readhostb.F
c     (line 396) as 
c
c     call assignclass(na,xa(1,1,i),ya(1,1,i),za(1,1,i),n12a,i12a,
c    &                atlaba,iguesta,torrefa(i),ibondera(i),loosera(i),
c    &                classa(i),sclassa(i),altrefa(i))
c
c     from readhosta.F. From readhostb.F, change the 'a' suffixes to b.
c     Compare this to the dummy variable which are used in this subroutine
c
c     subroutine assignclass(natom, x, y, z, n12, i12, atlab, guest, 
c    &                       ii, jj, kk, class, sub, alt_ii)
c
c     Therefore the actual --> dummy variable mappings are:
c       na          --> natom
c       xa(1,1,i)   --> xr
c       ya(1,1,i)   --> yr
c       za(1,1,i)   --> zr
c       n12a        --> n12r
c       i12a        --> i12r
c       atlaba      --> atlabr
c       iguesta     --> guest
c       torrefa(i)  --> ii
c       ibondera(j) --> jj
c       loosera(i)  --> kk
c       classa(i)   --> class
c       sclassa(i)  --> sub
c       altrefa(i)  --> alt_ii
c
c-----------------------------------------------------------------------
c      
c     P A S S E D   V A R I A B L E S:
c
c     integer   natom       number of atoms 
c     real      xr, yr, rz  coordinates of the atoms
c     integer   n12r        array of number attached atoms
c     integer   i12r        array of attached atoms
c     char      atlabr      atom labels
c     integer   guest       guest flags
c     integer   ii          serial number of driver atom
c     integer   jj          serial number of attach atom
c     integer   kk          serial number of the dummy atom to be lost
c     integer   class       the class of the rotor
c     integer   sub         the subclass of the rotor
c     integer   alt_ii      alternative driver atom to be used when
c                           a host is mirrored -  used only for 
c                           asymmetric trigonal pyramidal hosts.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine assignclass(natom, xr, yr, zr, n12r, i12r, atlabr, 
     &                       guest, ii, jj, kk, class, sub, alt_ii)

      implicit none
      include 'params.i'

c  passed variables      
      integer natom                 
      real xr(MAXATOMD), yr(MAXATOMD), zr(MAXATOMD)
      integer n12r(MAXATOMD)         
      integer i12r(MAXVAL,MAXATOMD)
      character*2 atlabr(MAXATOMD) 
      integer guest(MAXATOMD)      
      integer ii                   
      integer jj                    
      integer kk                 
      integer class, sub            
      integer alt_ii               

c  local variables used to hold the molecule
      integer n                   
      integer n12(MAXATOMD) 
      integer i12(MAXVAL,MAXATOMD)
      real x(MAXATOMD), y(MAXATOMD), z(MAXATOMD)
      character*2 atlab(MAXATOMD)

c  other LOCAL variables
      integer i, j           
      integer ntoss
      integer valency 
      integer nsmall, nmedium, nlarge
      integer idsmall(MAXVAL), idmedium(MAXVAL), idlarge(MAXVAL)
      integer val1, val2, vala1, vala2
      integer id1, id2, ida1, ida2
      integer idcis1, idcis2, idcisa1, idcisa2
       
      real ang, dist2p, yloser

      logical copl1, copl2, copla1, copla2
      logical cis1, cis2, cisa1, cisa2
      logical flat1, flat2, flata1, flata2
      logical sqpl
      
c-----------------------------------------------------------------------
c                END OF VARIABLE DECLARATIONS
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     Initialize variables.
c-----------------------------------------------------------------------
      copl1  = .false.
      copl2  = .false.
      copla1 = .false.
      copla2 = .false.
      cis1   = .false.
      cis2   = .false.
      cisa1  = .false.
      cisa2  = .false.
      flat1  = .false.
      flat2  = .false.
      flata1 = .false.
      flata2 = .false.
      sqpl   = .false.

      do i = 1, MAXATOMD
         x(i) = 0.0e0         
         y(i) = 0.0e0         
         z(i) = 0.0e0         
         n12(i) = 0        
         do j = 1, MAXVAL
            i12(j,i) = 0
         enddo
      enddo

      nsmall = 0
      nmedium = 0
      nlarge = 0
      do i = 1, MAXVAL
         idsmall(i) = 0
         idmedium(i) = 0
         idlarge(i) = 0
      end do
       
      ii = 0
      alt_ii = 0

c-----------------------------------------------------------------------
c     Copy molecule into local arrays used only in assignclass.
c     While doing this, remove any guest atoms.  The guest is
c     always at the end of the input fragment. 
c-----------------------------------------------------------------------

c  count the number of guest atoms
      ntoss = 0
      do i = 1, natom
         if(guest(i).eq.1) then
            ntoss = ntoss+1 
         endif
      enddo

      n = natom - ntoss

c  copy the passed molecule into the local molecule
      do i = 1, n
         x(i) = xr(i)
         y(i) = yr(i)
         z(i) = zr(i)
         atlab(i) = atlabr(i)
         n12(i) = n12r(i)
         do j = 1, n12(i)
            i12(j,i) = i12r(j,i)
         enddo 
      end do

c  remove any connections to guest atoms
      do i = 1, n
         do j = 1, n12(i)
            if(i12(j,i).gt.n) then
               i12(j,i) = 0
               n12(i) = n12(i) - 1
            endif
         enddo
      enddo

c-----------------------------------------------------------------------
c                          BEGIN MAIN BODY
c     The primary test is the valency and geometry of the link 
c     atom (jj).  Subsequent tests are based on the nature of the atoms
c     alpha to jj, beta to jj, and sometimes the identity of the element
c     jj.  The cases are terminal, linear, bent, trigonal planar, 
c     trigonal pyramidal, tetrahedral, square planar, and hypervalent. 
c-----------------------------------------------------------------------

      valency = n12(jj)

c-----------------------------------------------------------------------
c     Valency of 0: Not possible
c-----------------------------------------------------------------------

      if (valency.eq.0) then
          write(*,*) "*** ABORTING FROM ASSIGNCLASS.F. ***"
          write(*,*) '*** Impossible input - host fragment is a ',
     &     'single atom with no attachment point ***'
          stop
      end if

c-----------------------------------------------------------------------
c     Valency of 1: Terminal case
c
c     User has provided an input fragment in which the host consists of
c     a single atom with one bonding vector.  In this case, it is not
c     possible to define a torsion angle to the bond made between this
c     host atom and a linkage.
c-----------------------------------------------------------------------

      if (valency.eq.1) then
          class = 16
          sub   = 1
          ii    = 1
          return
      end if

c-----------------------------------------------------------------------
c     In what follows, sizes are specified as follows:
c      - Small = H atom
c      - Medium = any other atom with a valency of 2 or less
c      - Large = any other atom with a valency of 3 or more
c-----------------------------------------------------------------------
 
      do j = 1, valency
         if (i12(j,jj).eq.kk) cycle
         if (atlab(i12(j,jj)).eq.'H ') then
            nsmall = nsmall + 1
            idsmall(nsmall) = i12(j,jj)
         else if (n12(i12(j,jj)).le.2) then
            nmedium = nmedium + 1
            idmedium(nmedium) = i12(j,jj)
         else
            nlarge = nlarge + 1
            idlarge(nlarge) = i12(j,jj)
         endif
      end do

c-----------------------------------------------------------------------
c     Valency of 2
c-----------------------------------------------------------------------

      if (valency.eq.2) then
         do j = 1, 2
            if (i12(j,jj).eq.kk) cycle
            ii = i12(j,jj)
         end do
         call angleval(x,y,z,kk,jj,ii,ang)
         if (ang.ge.150.0e0) then                                !linear
            if (ang.eq.180.0e0) then
               xr(kk) = x(kk) + 0.02e0
               yr(kk) = y(kk) + 0.02e0
            end if
            class = 15
            sub = 1
            return
         else                                                      !bent
            if(nsmall.eq.1) then
               if(atlab(jj).eq.'O ') then
                  class = 5
                  sub = 2
               else
                  class = 5
                  sub = 4
               end if
               return
            else if ( (nmedium.eq.1).and.
     &              (atlab(ii).eq.'N ').or.
     &              (atlab(ii).eq.'P ').or.
     &              (atlab(ii).eq.'As') ) then
               call coplanar(ii,jj,kk,n12,i12,x,y,z,
     &              copla1,cisa1,idcisa1)
               if (cisa1) then
                  if (atlab(idcisa1).eq.'H ') then
                     class = 11
                     sub = 1
                  else
                     class = 11
                     sub = 2
                  end if
                  return
               else
                  class = 11
                  sub = 3
                  return
               end if
            else if (nmedium.eq.1) then
               class = 11
               sub = 3
               return
            else                               !attachment must be large
               if(n12(ii).eq.3) then
                  call dist2plane(ii,i12(1,ii),i12(2,ii),
     &                 i12(3,ii),x,y,z,dist2p)
                  if (dist2p.le.0.15e0) then 
                     call coplanar(ii,jj,kk,n12,i12,x,y,z,
     &                    copla1,cisa1,idcisa1)
                  end if
               else
                  dist2p = 1.0e0                 
               end if
               if (dist2p.le.0.015e0) then          !trigonal and planar
                  if (atlab(jj).eq.'O ') then
                     nlarge = 0
                     do j = 1, 3
                        if (i12(j,ii).eq.jj) cycle
                        if (n12(i12(j,ii)).ge.3) then
                           nlarge = nlarge + 1
                        end if
                     end do
                     if (nlarge.eq.2) then
                        class = 14
                        sub = 4
                     else 
                        class = 14
                        sub = 3
                     end if
                     return
                  end if 
                  if (atlab(idcisa1).eq.'H ') then
                     class = 11
                     sub = 1
                  else
                     class = 11
                     sub = 2
                  end if
                  return
               else                 !≥ 3 beta attachments and not planar
                  if (atlab(jj).eq.'O ') then
                     class = 5
                     sub = 1
                  else
                     class = 5
                     sub = 3
                  end if
                  return
               end if
            end if
         end if      
      end if

c-----------------------------------------------------------------------
c     Valency of 3
c
c     Trivalent rotors are first classed as pyramidal or planar
c
c     Subsequent assignments based on properties of attached groups
c        size - small, medium, or large
c     and when rotor is planar
c        planarity of attached groups
c        coplanarity of attached groups with bonding group
c        size of cis attachments
c-----------------------------------------------------------------------

      if (valency.eq.3) then

         call dist2plane(jj,i12(1,jj),i12(2,jj),i12(3,jj),x,y,z,dist2p)

c-----------------------------------------------------------------------
c  Pyramidal cases - distance from jj to plane defined by attached
c  atoms is greater than 0.15 angstroms 
c-----------------------------------------------------------------------

      if (dist2p.ge.0.15e0) then                        !pyramidal rotor

         if (nsmall.eq.2) then
            if (atlab(jj).eq.'N ') then
               class = 6
            else
               class = 9
            end if
            sub = 1
            ida1 = idsmall(1)
            ida2 = idsmall(2)
            call pyramid(jj,kk,ida1,ida2,x,y,z,yloser)
            if (yloser.gt.0.0e0) then
               ii = ida1
               alt_ii = ida2
            else
               ii = ida2
               alt_ii = ida1
            end if
            return
         else if (nsmall.eq.1) then
            if (atlab(jj).eq.'N ') then
               class = 6
            else
               class = 9
            end if 
            if (nmedium.eq.1) then
               ida2 = idmedium(1)
            else
               ida2 = idlarge(1)
            end if
            ida1 = idsmall(1)
            call pyramid(jj,kk,ida1,ida2,x,y,z,yloser)
            if (yloser.gt.0.0e0) then
               sub = 2
            else
               sub = 3
            end if
            ii = ida2
            return
         else if (nmedium.eq.1.and.nlarge.eq.1) then
            class = 8
            ida1 = idlarge(1)
            ida2 = idmedium(1)
            call pyramid(jj,kk,ida1,ida2,x,y,z,yloser)
            if (yloser.gt.0.0e0) then
               sub = 2
            else
               sub = 3
            end if
            ii = ida2
            return
         else                                !Either 2 medium or 2 large
            if (atlab(jj).eq.'N ') then
               class = 6
            else
               class = 9
            end if
            sub = 4
            if (nmedium.eq.2) then
               ida1 = idmedium(1)
               ida2 = idmedium(2)
            else
               ida1 = idlarge(1)
               ida2 = idlarge(2)
            end if
            if (yloser.gt.0.0e0) then
               ii = ida2
               alt_ii = ida1
            else
               ii = ida1
               alt_ii = ida2
            end if
               return
          end if

      else                                                 !planar rotor

c-----------------------------------------------------------------------
c     Characterization of planar trivalent rotors
c
c     highest priority alpha attachment is assigned to a1
c     priority is in this order:
c        first by valence
c        next coplanar has higher priority
c        next a cis substituent is higher priority than no cis 
c-----------------------------------------------------------------------

      if (nsmall.eq.2) then
         ida1 = idsmall(1)
         vala1 = 1
         ida2 = idsmall(2)
         vala2 = 1

      else if (nsmall.eq.1.AND.nmedium.eq.1) then
         ida1 = idmedium(1)
         vala1 = n12(ida1)
         if(vala1.eq.2) then
            call coplanar(ida1,jj,kk,n12,i12,x,y,z,copla1,cisa1,idcisa1)
         end if
         ida2 = idsmall(1) 
         vala2 = 1

      else if (nsmall.eq.1.AND.nlarge.eq.1) then
         ida1 = idlarge(1)
         vala1 = n12(ida1)
         if (vala1.eq.3) then
            call dist2plane(ida1, i12(1,ida1), i12(2,ida1),
     &           i12(3,ida1), x, y, z, dist2p)
            if (dist2p.le.0.15e0) then
               flata1 = .true.
               call coplanar(ida1, jj, kk, n12, i12, x, y, z, 
     &              copla1,cisa1,idcisa1)
            end if
         end if
         ida2 = idsmall(1)
         vala2 = 1

      else if (nmedium.eq.2) then
         id1 = idmedium(1)
         val1 = n12(id1)
         if (val1.eq.2) then
            call coplanar(id1,jj,kk,n12,i12,x,y,z,copl1,cis1,idcis1)
         end if
         id2 = idmedium(2)
         val2 = n12(id2)
         if (val2.eq.2) then
            call coplanar(id2,jj,kk,n12,i12,x,y,z,copl2,cis2,idcis2)
         end if

      else if (nmedium.eq.1.AND.nlarge.eq.1) then
         ida1 = idlarge(1)
         vala1 = n12(ida1)
         if (vala1.eq.3) then
            call dist2plane(ida1, i12(1,ida1), i12(2,ida1),
     &           i12(3,ida1), x, y, z, dist2p)
            if (dist2p.le.0.15e0) then
               flata1 = .true.
               call coplanar(ida1, jj, kk, n12, i12, x, y, z, 
     &              copla1,cisa1,idcisa1)
            end if
         end if
         ida2 =  idmedium(1)
         vala2 = n12(ida2)
         if (vala2.eq.2) then
            call coplanar(ida2,jj,kk,n12,i12,x,y,z,copla2,cisa2,idcisa2)
         end if

      else
         id1 = idlarge(1)
         val1 = n12(id1)
         if (val1.eq.3) then
            call dist2plane(id1, i12(1,id1), i12(2,id1),
     &           i12(3,id1), x, y, z, dist2p)
            if (dist2p.le.0.15e0) then
               flat1 = .true.
               call coplanar(id1, jj, kk, n12, i12, x, y, z, 
     &              copl1,cis1,idcis1)
            end if
         end if
         id2 = idlarge(2)
         val2 = n12(id2)
         if (val2.eq.3) then
            call dist2plane(id2, i12(1,id2), i12(2,id2),
     &           i12(3,id2), x, y, z, dist2p)
            if (dist2p.le.0.15e0) then
               flat2 = .true.
               call coplanar(id2, jj, kk, n12, i12, x, y, z, 
     &              copl2,cis2,idcis2)
            end if
         end if
      end if

c  if the groups are both medium or both large, assign priority order
      if (nmedium.eq.2.or.nlarge.eq.2) then
         if ( (flat1.and.cis1) .or.
     &        (flat1.and.(.not.cis2)) .or.
     &        (.not.flat2) ) then
            ida1 = id1
            vala1 = val1
            flata1 = flat1
            copla1 = copl1
            cisa1 = cis1
            idcisa1 = idcis1
            ida2 = id2
            vala2 = val2
            flata2 = flat2
            copla2 = copl2
            cisa2 = cis2
            idcisa2 = idcis2
         else 
            ida1 = id2
            vala1 = val2
            flata1 = flat2
            copla1 = copl2
            cisa1 = cis2
            idcisa1 = idcis2
            ida2 = id1
            vala2 = val1
            flata2 = flat1
            copla2 = copl1
            cisa2 = cis1
            idcisa2 = idcis1
         endif
      endif

c-----------------------------------------------------------------------
c     Assignments for planar trivalent rotors
c-----------------------------------------------------------------------

c  Case 1: 2 H
      if (nsmall.eq.2) then
         class = 14
         sub = 1
         ii = idsmall(1)
         return
      end if

c  Case 2: H + terminal non-H atom
      if (nsmall.eq.1.and.vala1.eq.1) then
         if (atlab(ida1).eq.'S ') then
            class = 10
            sub = 3
         else
            class = 10
            sub = 1
         end if
         ii = ida1
         return
      end if

c  Case 3: H + divalent attachment
      if (nsmall.eq.1.and.vala1.eq.2) then
         if ( (atlab(ida1).eq.'N ').or.
     &        (atlab(ida1).eq.'P ').or.
     &        (atlab(ida1).eq.'As') ) then
            if (cisa1) then
               class = 12
               sub = 4
            else
               class = 12
               sub = 2
            end if
         else
            class = 14
            sub = 2
         end if
            ii = ida1
            return
      end if

c  Case 4: H + ≥ trivalent attachment
      if (nsmall.eq.1.and.vala1.ge.3) then
         if (.not.flata1) then                         !large not planar
            class = 14
            sub = 2
            ii = ida1
            return
         else if (cisa1) then                        !planar, cis present
            if (atlab(idcisa1).eq.'H ') then                !cis is small
               class = 2                           
               sub = 1
               ii = ida1
               return
            else
               if (atlab(jj).eq.'N ') then                 !amide N rotor
                  if ( (n12(idcisa1).le.2).and.
     &            (atlab(idcisa1).ne.'H ') ) then          !cis is medium 
                     class = 4
                     sub = 1
                  else if (n12(idcisa1).ge.3) then         !cis is large 
                     class = 4
                     sub = 2
                  end if
                  ii = ida1
                  return
               else                    
                  class = 2
                  sub = 2
                  ii = ida1
                  return
               end if  
            end if
         end if
      end if

c  Case 5: 2 non-H groups, both monovalent or 1 monovalent and 1 divalent
      if (nmedium.eq.2.and.vala1.eq.1.and.vala2.eq.1 )then
         class = 14
         sub = 1
         ii = ida1
         return
      else if (nmedium.eq.2.and.vala1.eq.1.and.vala2.eq.2) then
         class = 14
         sub = 2
         ii = ida1
         return
      else if (nmedium.eq.2.and.vala1.eq.2.and.vala2.eq.1) then
         class = 14
         sub = 2
         ii = ida2
         return
      end if

c  Case 6: 2 divalent groups
      if (nmedium.eq.2.and.vala1.eq.2.and.vala2.eq.2) then
         if ( (copla1.and.copla2) ) then                  !both coplanar
            if ( (.not.cisa1).and.(.not.cisa2) ) then        !both trans 
               class = 13
               sub = 2
            else
               class = 14
               sub = 1
            end if
            ii = ida1
            return
         else
            class = 14
            sub = 1
            ii = ida1
            return
         end if
      end if

c  Case 7: 1 ≥ trivalent and 1 non-H monovalent
      if (nmedium.eq.1.and.vala2.eq.1.and.nlarge.eq.1) then
         if (.not.flata1) then                         !large not planar
            if(atlab(ida2).eq.'S ') then
               class = 10
               sub = 4
            else
               class = 10
               sub = 2
            end if
            ii = ida2
            return
         else if (atlab(idcisa1).eq.'H ') then            !planar, cis H
            class = 7
            sub = 1
            ii = ida2
            return
         else                                         !planar, cis not-H
            class = 7
            sub = 2
            ii = ida2
            return
         end if
      end if

c  Case 8: 1 ≥ trivalent and 1 divalent
      if (nmedium.eq.1.and.vala2.eq.2.and.nlarge.eq.1) then
         if ( (atlab(ida2).eq.'N ').or.
     &        (atlab(ida2).eq.'P ').or.
     &        (atlab(ida2).eq.'As') ) then
              flata2 = .true.
         end if
         if ((.not.flata1).and.(.not.flata2)) then
            class = 14
            sub = 1
            ii = ida1
            return
         else if ((flata1).and.(.not.flata2)) then
            if (atlab(idcisa1).eq.'H ') then
               class = 2
               sub = 3
               ii = ida1
               return
            else
               class = 2
               sub = 4
               ii = ida1
               return
            end if
         else if ((.not.flata1).and.(flata2)) then
            if(cisa2) then
               class = 12
               sub = 3
               ii = ida2
               return
            else
               class = 12
               sub = 1
               ii = ida2
               return
            endif
         else                                               !both planar
            if (cisa2) then
               if (atlab(idcisa1).eq.'H ') then
                  class = 2
                  sub = 3
                  ii = ida1
                  return
               else
                  class = 2
                  sub = 4
                  ii = ida1
                  return
               end if
            else
               if (atlab(idcisa1).eq.'H ') then
                  class = 13
                  sub = 1
                  ii = ida2
                  return
               else
                  class = 13
                  sub = 3
                  ii = ida2
                  return
               end if
            end if
         end if
      end if

c  Case 9:  2 groups ≥ trivalent
      if (nlarge.eq.2) then
         if ((.not.flata1).and.(.not.flata2)) then
            class = 14
            sub = 1
            ii = ida1
            return
         else if ((flata1).and.(.not.flata2)) then
            if (atlab(idcisa1).eq.'H ') then
               class = 2
               sub = 3
               ii = ida1
               return
            else
               if (atlab(jj).eq.'N ') then                !amide N rotor
                  if ( (n12(idcisa1).le.2).and.
     &            (atlab(idcisa1).ne.'H ') ) then         !cis is medium 
                     class = 4
                     sub = 3
                     ii = ida1
                     return
                  else if (n12(idcisa1).ge.3) then         !cis is large 
                     class = 4
                     sub = 4
                     ii = ida1
                     return
                  end if
               else                    
                  class = 2
                  sub = 4
                  ii = ida1
                  return
               end if
            end if
         else                                               !both planar
  
            if ( (atlab(idcisa1).eq.'H ').and.
     &           (atlab(idcisa2).eq.'H ') ) then
               class = 3
               sub = 1
               ii = ida1
               return
            else if ( (atlab(idcisa1).eq.'H ').and.
     &                (atlab(idcisa2).ne.'H ') ) then
               class = 3
               sub = 2
               ii = ida2
               return
            else if ( (atlab(idcisa1).ne.'H ').and.
     &                (atlab(idcisa2).eq.'H ') ) then
               class = 3
               sub = 2
               ii = ida1
               return
            else
               class = 3
               sub = 3
               ii = ida1
               return
            end if
         end if
      end if

c-----------------------------------------------------------------------
      end if                                ! End of planar cases
      end if                                ! End of Valency = 3 section
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     Valency of 4
c-----------------------------------------------------------------------

      if (valency.eq.4) then

         do j = 1, valency
            if (i12(j,jj).eq.kk) cycle
            ii = i12(j,jj)
            call angleval(x,y,z,kk,jj,ii,ang)
            if (ang.ge.155.0e0) then
               sqpl = .true.
            else
               id1 = ii
            end if
         end do

         if (sqpl) then                                   !square planar
            class = 3
            sub = 3
            ii = id1
            return

         else                                               !tetrahedral

c-----------------------------------------------------------------------
c     Tetrahedral rotors are characterized by the relative size of their 
c     substituents, for example:
c      - Three same size groups, e.g., t-butyl rotor
c      - Two large groups one small group, e.g., i-propyl rotor
c      - Two small groups, one large group, e.g., ethyl rotor
c-----------------------------------------------------------------------

c Case 1 - all same size
            if ((nsmall.eq.3).or.(nmedium.eq.3).or.(nlarge.eq.3)) then
               class = 1
               sub = 3
               if(nsmall.ne.0) then
                  ii = idsmall(1)
               else if (nmedium.ne.0) then
                  ii = idmedium(1)
               else
                  ii = idlarge(1)
               end if 
               return

c Case 2 - two small 
            else if (nsmall.eq.2) then
               class = 1
               sub = 1
               if (nmedium.eq.1) then
                  ii = idmedium(1)
               else
                  ii = idlarge(1)
               endif
               return

c Case 3 - one small
            else if ((nsmall.eq.1)) then
               class = 1
               sub = 2
               ii = idsmall(1)
               return

c Case 4 - two medium, one large
            else if (nmedium.eq.2) then
               if ( (atlab(jj).eq.'S ').or.
     &              (atlab(jj).eq.'Se').or.
     &              (atlab(jj).eq.'Te').or.
     &              (atlab(jj).eq.'P ').or.
     &              (atlab(jj).eq.'As')) then
                  class = 8
                  sub = 4
                  ii = idlarge(1)
               else
                  class = 1
                  sub = 1
                  ii = idlarge(1)
               end if
               return
          
c Case 5 - one medium, two large
            else if (nmedium.eq.1) then
               if ( (atlab(jj).eq.'S ').or.
     &              (atlab(jj).eq.'Se').or.
     &              (atlab(jj).eq.'Te').or.
     &              (atlab(jj).eq.'P ').or.
     &              (atlab(jj).eq.'As')) then
                  class = 8
                  sub = 1
                  ii = idmedium(1)
               else
                  class = 1
                  sub = 2
                  ii = idmedium(1)
               end if
               return

            end if
         end if                              !End of tetrahedral section
      end if                                 !End of Valency = 4 section

c-----------------------------------------------------------------------
c     Valency of ≥ 5: Hypervalent cases 
c
c     When there are 5 or more attachments, then it is assumed that 
c     there are no dihedral preferences about the bond made to this 
c     atom.
c-----------------------------------------------------------------------

      do j = 1, valency
         if (i12(j,jj).eq.kk) cycle  !skip loser atom kk
         call angleval(x,y,z,kk,jj,i12(j,jj),ang)
         if (ang.lt.160.0e0) then
            ii = i12(j,jj)
            exit
         end if
      end do
      class = 15
      sub = 1
      return

c-----------------------------------------------------------------------
c     No more cases possible 
c-----------------------------------------------------------------------
      end



c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'coplanar' -- determines if an attachment bound to planar group 
c     through bond ii-jj is coplanar.  The attachment is either divalent
c     or planar trivalent.  Also determines if an attachment
c     on ii is cis to jj-kk bond.  
c
c     passed variables:
c     integer ii           serial number of attached atom
c     integer jj           serial number of attached planar atom
c     integer kk           serial number of bonding vector to jj
c     integer n12, i12     attachment arrays
c     logical copl         .true. if coplanar
c     logical cis          .true. if cis group present
c     integer idcis        serial number of cis atom
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine coplanar(ii,jj,kk,n12,i12,x,y,z,copl,cis,idcis)

      implicit none
      include 'params.i'
      include 'constants.i'

      integer ii, jj, kk, j
      integer n12(MAXATOMD), i12(MAXVAL,MAXATOMD)
      integer idcis
      integer numii, numjj, idii(2), idjj(2)

      real x(MAXATOMD), y(MAXATOMD), z(MAXATOMD)
      real ang

      logical copl, cis

c  Count attachments on ii and jj and remember their serial numbers
      numii = 0
      do j = 1, n12(ii)
         if(i12(j,ii).eq.jj) cycle
         numii = numii + 1
         idii(numii) = i12(j,ii)
      enddo
      numjj = 0
      do j = 1, n12(jj)
         if(i12(j,jj).eq.ii) cycle
         numjj = numjj + 1
         idjj(numjj) = i12(j,jj)
      enddo

c  Stop for unallowed attachments
      if(numjj.lt.1.or.numjj.gt.2) then
         print*,'Expected numjj = 1 or 2'
         print*,'Obtained numjj = ',numjj
         stop
      end if
      if(numii.lt.1.or.numii.gt.2) then
         print*,'Expected numjj = 1 or 2'
         print*,'Obtained numjj = ',numjj
         stop
      end if

c  Coplanar means that the dihedral angle between an attachment
c  on ii and jj-kk is within ± 10 deg of 0 or 180 

      call torsionval(x,y,z,idii(1),ii,jj,kk,ang)
      ang = ang*180.0e0/pi

      if (ang.gt.180.0e0) then
          ang = 360.0e0 - ang
      end if

      copl = .false.
      if(ang.ge.0.0e0.and.ang.le.10.0e0) then
         copl = .true.
      else if (ang.gt.170.0e0.and.ang.le.180) then
         copl = .true.
      end if

      if(ang.ge.0.0e0.and.ang.le.90.0e0) then
         cis = .true.
         idcis = idii(1)
      else if (ang.gt.90.0e0.and.ang.le.180.0e0) then
         if (numii.eq.2) then
            cis = .true.
            idcis = idii(2)
         else
            cis = .false.
            idcis = 0
         end if
      else
         cis = .false.
         idcis = 0
      end if

      return
      end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'dist2plane' -- returns distance of atom ii from a plane
c                     defined by atoms jj,kk,ll
c
c     passed variables:
c     integer ii           serial number of test atom
c     integer jj           serial number of first plane atom
c     integer kk           serial number of second plane atom
c     integer ll           serial number of third plane atom
c     real    x,y,z        atomic coordinates
c     real    dist2pl      distance ii from plane jj,kk,ll
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine dist2plane(ii,jj,kk,ll,x,y,z,dist2pl)

      implicit none
      include 'params.i'

      integer ii, jj, kk, ll

      real x(MAXATOMD), y(MAXATOMD), z(MAXATOMD)
      real v1x, v1y, v1z, v2x, v2y, v2z, nvx, nvy, nvz
      real magnormvect, dist2pl

c     first vector in plane
      v1x = x(ll)-x(jj)
      v1y = y(ll)-y(jj)
      v1z = z(ll)-z(jj)

c     second vector in plane
      v2x = x(kk)-x(jj)
      v2y = y(kk)-y(jj)
      v2z = z(kk)-z(jj)

c     normal vector to plane
      nvx = v1y*v2z - v1z*v2y
      nvy = v1z*v2x - v1x*v2z
      nvz = v1x*v2y - v1y*v2x

      magnormvect = sqrt(nvx*nvx + nvy*nvy + nvz*nvz)

      dist2pl = (nvx*(x(ii)-x(jj)) + nvy*(y(ii)-y(jj)) +
     &           nvz*(z(ii)-z(jj)))/magnormvect
     
      dist2pl = abs(dist2pl)

      return
      end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'pyramid' -- orients a pyramidal rotating group such that the 
c     first attached atom is on the origin, the second attached atom is
c     on the +z-axis, and the central bonding atom is on the x-axis.  
c     Returns the y value for the atom to be lost.
c
c     passed variables:
c     integer jj           serial number of bonding atom
c     integer kk           serial number of atom to be lost
c     integer ll           serial number of first attached atom
c     integer mm           serial number of second attached atom
c     real    x,y,z        coordinates
c     real    yloser       y coordinate for atom to be lost
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine pyramid(jj,kk,ll,mm,x,y,z,yloser)

      implicit none
      include 'params.i'

      integer jj, kk, ll, mm

      real x(MAXATOMD), y(MAXATOMD), z(MAXATOMD)
      real dx(5), dy(5), dz(5), yloser
      real ang, phi
     
c     Atom numbering used in orientation is as follows:
c
c     1 = dummy
c     2 = first attached atom
c     3 = second attached atom
c     4 = central bonding atom
c     5 = atom to be lost

c initially, we don't know where dummy goes, doesn't matter as will
c reasign it after aligning to the z-axis
      dx(1) = 0.0e0
      dy(1) = 0.0e0
      dz(1) = 0.0e0

c this first attached atom will end up on the origin
      dx(2) = x(ll)
      dy(2) = y(ll)
      dz(2) = z(ll)

c this second attached atom will end up on the +z axis
      dx(3) = x(mm)
      dy(3) = y(mm)
      dz(3) = z(mm)

c this central bonding atom will end up in the xy-plane   
      dx(4) = x(jj)
      dy(4) = y(jj)
      dz(4) = z(jj)

      dx(5) = x(kk)
      dy(5) = y(kk)
      dz(5) = z(kk)

      call alignz(5,dx,dy,dz,2,3,1)

c reposition the dummy
      dx(1) = 1.0e0
      dy(1) = 0.0e0
      dz(1) = 0.0e0

      call torsionval(dx,dy,dz,1,2,3,4,phi)
      ang = 0.0e0
      call setdihedral(1,5,dx,dy,ang,phi)

c     The sign of the y coordinate of the loser atom, the fifth atom in
c     the dx,dy,dz arrays, provides the criterion for assignments

c     When the y value is positive the lost atom (bonding vector) is    
c     pointing down when viewing the x-z plane from the +z direction

c     When the y value is negative the lost atom (bonding vector) is
c     pointing up when viewing the x-z plane from the +z direction

      yloser = dy(5)

      return
      end
