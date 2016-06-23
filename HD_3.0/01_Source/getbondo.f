c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'getbondo' -- returns length of bond that will be formed between
c                   atom ii in the input fragment and atom jj in the 
c                   linking fragment
c
c     passed variables:
c     integer ii           bonding atom in input
c     integer ilost        lost atom in input
c     integer n12i, i12i   attachment arrays in input
c     character atlabi     atom labels in input
c     integer iguest       flags for guest atoms
c     integer jj           bonding atom in link
c     integer jlost        lost atom in link
c     integer n12j, i12j   attachment arrays in link
c     character atlabj     atom labels in link
c     real bondl           bond length
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine getbondo(ii,ilost,n12i,i12i,atlabi,iguest,
     &                    jj,jlost,n12j,i12j,atlabj,bondl)

      implicit none
      include 'params.i'
      include 'constants.i'

      integer i, j, k, ii, jj, kk
      integer ilost, jlost
      integer elemi, elemj
      integer n12i(MAXATOMD), i12i(MAXVAL,MAXATOMD)
      integer iguest(MAXATOMD)
      integer n12j(MAXATOMD), i12j(MAXVAL,MAXATOMD)
      integer nc1, nc2, nnit1, nnit2, noxo1, noxo2
      integer nattc, ncc, natto, nco
      integer nc_on_att, boost
      integer val1, val2

      real elnegi, elnegj, radi, radj
      real bondl

      character*2  atlabi(MAXATOMD)
      character*2  atlabj(MAXATOMD)
      character*2  atlabo

c-----------------------------------------------------------------------
c     Get constants for bonded atoms ii and jj, here natoms refers to 
c     the number of elements stored in the CONSTANTS file
c-----------------------------------------------------------------------

      elemi = 0
      elemj = 0
      radi = 0.0e0
      radj = 0.0e0
      do i = 1, natoms
         if (atlabi(ii).eq.atomname(i)) then
            elemi = i
            elnegi = elneg(i)
            radi = covrad(i)
         endif
         if (atlabj(jj).eq.atomname(i)) then
            elemj = i
            elnegj = elneg(i)
            radj = covrad(i) 
         end if
      end do

c-----------------------------------------------------------------------
c     If do not find an element, then stop
c-----------------------------------------------------------------------

      if (elemi.eq.0) then
         print*, 'Stopped in subroutine GETBONDLENGTH:'
         print*, 'Do not find element ',atlabi(ii),' in CONSTANTS file'
         stop
      end if
      if (elemj.eq.0) then
         print*, 'Stopped in subroutine GETBONDLENGTH:'
         print*, 'Do not find element ',atlabj(jj),' in CONSTANTS file'
         stop
      end if

c-----------------------------------------------------------------------
c     Compute default single bond length using the Blom-Haaland method:
c     J. Mol. Struct. 1985, 128, 21
c-----------------------------------------------------------------------
      
      if ((elnegi.ne.0.0e0).and.(elnegj.ne.0.0e0)) then
         bondl = radi + radj - 0.085e0*(abs(elnegi - elnegj))**1.40e0
      else
         bondl = radi + radj
      end if

c-----------------------------------------------------------------------
c     If neither atom in the bond is a carbon, then return
c-----------------------------------------------------------------------

      if ((atlabi(ii).ne.'C ').and.(atlabj(jj).ne.'C ')) then
         return
      end if

c-----------------------------------------------------------------------
c     If this is a carbon-carbon bond . . . 
c-----------------------------------------------------------------------

      if ((atlabi(ii).eq.'C ').and.(atlabj(jj).eq.'C ')) then

c-----------------------------------------------------------------------
c     Get the valence on each of the bonded atoms, ignore any attached
c     guest atoms on the input fragment
c-----------------------------------------------------------------------

         val1 = 0
         do j = 1, n12i(ii)
            if (iguest(i12i(j,ii)).eq.1) cycle
            val1 = val1 + 1
         end do

         val2 = n12j(jj)

c-----------------------------------------------------------------------
c     Count number of carbon atoms attached to each bonded atom,
c     excluding the one that defines the bond in question
c     Count number of carbon atoms on carbon atoms attached to each
c     bonded atom.  Boost = boost + 1 when the following is true:
c        -bonded atom is sp2
c        -the attached carbon has ≥ 2 carbon attached to it
c-----------------------------------------------------------------------

         boost = 0

         nc1 = 0
         noxo1 = 0
         nnit1 = 0
         do j = 1, n12i(ii)
            if((i12i(j,ii).eq.ilost).or.(iguest(i12i(j,ii)).eq.1)) cycle
            if (atlabi(i12i(j,ii)).eq.'N ') then
               nnit1 = nnit1 + 1
            else if (atlabi(i12i(j,ii)).eq.'O ') then
               noxo1 = noxo1 + 1
            else if (atlabi(i12i(j,ii)).eq.'X ') then
               nc1 = nc1+1
            else if (atlabi(i12i(j,ii)).eq.'C ') then
               nc1 = nc1+1
               k = i12i(j,ii)
               nc_on_att = 0
               do kk = 1, n12i(k)
                  if (i12i(kk,k).eq.ii) cycle
                  if (atlabi(i12i(kk,k)).eq.'C ') then
                     nc_on_att = nc_on_att + 1
                  end if
               end do
               if ((val1.eq.3).and.(nc_on_att.eq.2)) then
                  boost = boost + 1
               end if
            end if
         end do

         nc2 = 0
         noxo2 = 0
         nnit2 = 0
         do j = 1, n12j(jj)
            if (i12j(j,jj).eq.jlost) cycle
            if (atlabj(i12j(j,jj)).eq.'N ') then
               nnit2 = nnit2 + 1
            else if (atlabj(i12j(j,jj)).eq.'O ') then
               noxo2 = noxo2 + 1
            else if (atlabj(i12j(j,jj)).eq.'X ') then
               nc2 = nc2+1
            else if (atlabj(i12j(j,jj)).eq.'C ') then
               nc2 = nc2+1
               k = i12j(j,jj)
               nc_on_att = 0
               do kk = 1, n12j(k)
                  if (i12j(kk,k).eq.jj) cycle
                  if (atlabj(i12j(kk,k)).eq.'C ') then
                     nc_on_att = nc_on_att + 1 
                  end if
               end do
               if ((val2.eq.3).and.(nc_on_att.eq.2)) then
                  boost = boost + 1
               end if
            end if
         end do

c-----------------------------------------------------------------------
c     Adjust the bond length due to hybridization followed by sterics
c-----------------------------------------------------------------------

c Default Blom-Haaland = 1.580 for a C-C single bond

c Case 1: both sp3
         if ((val1.eq.4).and.(val2.eq.4)) then        
c Correction to give the MM3 = 1.536 for least hindered case 
            bondl = bondl - 0.044                          ! gives 1.536
c Steric corrections based on MM3 calculations
            if ((nc1.eq.3).and.(nc2.eq.3)) then
               bondl = bondl + 0.040                       ! gives 1.576
            else if ( ((nc1.eq.3).and.(nc2.eq.2)).or.
     &                ((nc1.eq.2).and.(nc2.eq.3)) ) then
               bondl = bondl + 0.025                       ! gives 1.561
            else if ( ((nc1.eq.3).and.(nc2.eq.1)).or.
     &                ((nc1.eq.1).and.(nc2.eq.3)) ) then
               bondl = bondl + 0.012                       ! gives 1.548
            else if ((nc1.eq.2).and.(nc2.eq.2)) then
               bondl = bondl + 0.016                       ! gives 1.552
            else if ( ((nc1.eq.2).and.(nc2.eq.1)).or.
     &              ((nc1.eq.1).and.(nc2.eq.2)) ) then
               bondl = bondl + 0.007                       ! gives 1.543
            end if

c Case 2: sp3 + sp2
         else if ( ((val1.eq.4).and.(val2.eq.3)).or.
     &             ((val1.eq.3).and.(val2.eq.4)) ) then
c MM3 = 1.506, hybridization correction for least hindered case
            bondl = bondl - 0.074                          ! gives 1.506 
c Steric corrections based on MM3 calculations
            if ( ((nc1.eq.3).and.(nc2.eq.2)).or.
     &           ((nc2.eq.3).and.(nc1.eq.2)) ) then
               bondl = bondl + 0.018                       ! gives 1.524
               if(boost.ne.0) then
                  bondl = bondl + 0.0075*boost
               end if
            else if ( ((nc1.eq.3).and.(nc2.eq.1)).or.
     &                ((nc1.eq.1).and.(nc2.eq.3)) ) then
               bondl = bondl + 0.009                       ! gives 1.515
            else if ( ((nc1.eq.3).and.(nnit2.eq.1).and.(noxo2.eq.1)).or.
     &           ((nc2.eq.3).and.(nnit1.eq.1).and.(noxo1.eq.1)) ) then
               bondl = bondl + 0.009                       ! gives 1.523
            else if ((nc1.eq.2).and.(nc2.eq.2)) then
               bondl = bondl + 0.010                       ! gives 1.516
            else if ( ((nc1.eq.2).and.(nc2.eq.1)).or.
     &                ((nc1.eq.1).and.(nc2.eq.2)) ) then
               bondl = bondl + 0.004                       ! gives 1.510
            else if ( ((nc1.le.1).and.(nc2.eq.2)).or.
     &                ((nc1.eq.2).and.(nc2.le.1)) ) then
               bondl = bondl + 0.005                       ! gives 1.511
            end if
c Amide correction
            if ( ((nnit1.eq.1).and.(noxo1.eq.1)).or.
     &           ((nnit2.eq.1).and.(noxo2.eq.1)) ) then
               bondl = bondl + 0.008                       ! gives 1.514
            end if

c Case 3: sp3 + sp
         else if ( ((val1.eq.4).and.(val2.eq.2)).or.
     &             ((val1.eq.2).and.(val2.eq.4)) ) then
c MM3 = 1.474, hybridization correction for least hindered case
            bondl = bondl - 0.106                          ! gives 1.474
c Steric corrections based on MM3 calculations
            if ((nc1.eq.3).or.(nc2.eq.3)) then
               bondl = bondl + 0.005                       ! gives 1.479
            else if ((nc1.eq.2).or.(nc2.eq.2)) then
               bondl = bondl + 0.003                       ! gives 1.477
            end if
 
c Case 4: both sp2
         else if ((val1.eq.3).and.(val2.eq.3)) then
c MM3 = 1.464, hybridization correction for least hindered case
            bondl = bondl - 0.116                          ! gives 1.464
c Steric corrections based on MM3 calculations
            if ((nc1.eq.2).and.(nc2.eq.2)) then
               bondl = bondl + 0.022                       ! gives 1.486
            else if (((nc1.eq.2).and.(nc2.eq.1)).or.
     &               ((nc1.eq.1).and.(nc2.eq.2))) then
               bondl = bondl + 0.008                       ! gives 1.472
            else if (((nnit1.eq.1).and.(noxo1.eq.1)).or.
     &               ((nnit2.eq.1).and.(noxo2.eq.1))) then
               bondl = bondl + 0.012                       ! gives 1.476
               if ((nc1.eq.2).or.(nc2.eq.2)) then
                  bondl = bondl + 0.025                    ! gives 1.501
               end if
            else if ((noxo1.eq.2).or.(noxo2.eq.2)) then
               bondl = bondl + 0.022                       ! gives 1.486
            end if
            if (boost.ne.0) then
                bondl = bondl + 0.0075*boost
                if ((boost.eq.2).and.
     &             (((noxo1.eq.1).and.(nnit1.eq.0)).or.
     &              ((noxo2.eq.1).and.(nnit2.eq.0)))) then
                   bondl = bondl + 0.026                   ! gives 1.513
                end if
            end if

c Case 5: sp2 + sp
         else if ( ((val1.eq.3).and.(val2.eq.2)).or.
     &             ((val1.eq.2).and.(val2.eq.3)) ) then
c MM3 = 1.434, hybridization correction for least hindered case
            bondl = bondl - 0.146                          ! gives 1.434 

c Case 6: both sp
         else if ((val1.eq.2).and.(val2.eq.2)) then
c MM3 = 1.394, hybridization correction for least hindered case
            bondl = bondl - 0.186                          ! gives 1.394 

         end if

c-----------------------------------------------------------------------
c     One of the bonded atoms is a carbon atom
c-----------------------------------------------------------------------

      else 

c-----------------------------------------------------------------------
c     Characterize the bonded atoms
c-----------------------------------------------------------------------

         nattc = 0
         ncc = 0
         natto = 0
         nco = 0
         boost = 0

         if (atlabi(ii).eq.'C ') then
            do j = 1, n12i(ii)
               if (iguest(i12i(j,ii)).eq.1) cycle
               nattc = nattc + 1
            end do
            if ( (elnegj.eq.0.0e0).and.(nattc.eq.3) ) then
               bondl = bondl - 0.065e0
               return
            end if
            natto = n12j(jj)
            atlabo = atlabj(jj)
            do j = 1, n12i(ii)
               if ( (iguest(i12i(j,ii)).eq.1).or.
     &              (i12i(j,ii).eq.ilost) ) cycle
               if (atlabi(i12i(j,ii)).eq.'C ') then
                  ncc=ncc+1
                  k = i12i(j,ii)
                  nc_on_att = 0
                  do kk = 1, n12i(k)
                     if (i12i(kk,k).eq.ii) cycle
                     if (atlabi(i12i(kk,k)).eq.'C ') then
                        nc_on_att = nc_on_att + 1
                     end if
                  end do
                  if ((nattc.eq.3).and.(nc_on_att.eq.2)) then
                     boost = boost + 1
                  end if
               end if
            end do
            do j = 1, natto
               if(i12j(j,jj).eq.jlost) cycle
               if (atlabj(i12j(j,jj)).eq.'C ') then
                  nco=nco+1
               end if
            end do
         else if (atlabj(jj).eq.'C ') then
            nattc = n12j(jj)
            if ( (elnegi.eq.0.0e0).and.(nattc.eq.3) ) then
               bondl = bondl - 0.065e0
               return
            end if
            do j = 1, n12i(ii)
               if (iguest(i12i(j,ii)).eq.1) cycle
               natto = natto + 1
            end do
            atlabo = atlabi(ii)
            do j = 1, nattc
               if(i12j(j,jj).eq.jlost) cycle
               if (atlabj(i12j(j,jj)).eq.'C ') then
                  ncc=ncc+1
                  k = i12j(j,jj)
                  nc_on_att = 0
                  do kk = 1, n12j(k)
                     if (i12j(kk,k).eq.jj) cycle
                     if (atlabj(i12j(kk,k)).eq.'C ') then
                        nc_on_att = nc_on_att + 1
                     end if
                  end do
                  if ((nattc.eq.3).and.(nc_on_att.eq.2)) then
                     boost = boost + 1
                  end if
               end if
            end do
            do j = 1, natto
               if ( (iguest(i12i(j,ii)).eq.1).or.
     &              (i12i(j,ii).eq.ilost) ) cycle
               if (atlabi(i12i(j,ii)).eq.'C ') then
                  nco=nco+1
               end if
            end do
         end if

c-----------------------------------------------------------------------
c     Adjust the bond length due to hybridization followed by sterics
c-----------------------------------------------------------------------

         if (nattc.eq.4) then                              ! sp3 carbon

c Case 1: Bound to N
            if (atlabo.eq.'N ') then

c Default Blom-Haaland = 1.481 for a C-N single bond

c    Case A: N has 4 attached atoms
               if (natto.eq.4) then
c    Correction to give the MM3 = 1.518 for least hindered case
                  bondl = bondl + 0.037                    ! gives 1.518
c    Steric corrections based on MM3 calculations
                  if ((ncc.eq.3).and.(nco.eq.3)) then
                     bondl = bondl + 0.045                 ! gives 1.563
                  else if (((ncc.eq.3).and.(nco.eq.2)).or.
     &                    ((ncc.eq.2).and.(nco.eq.3))) then
                     bondl = bondl + 0.028                 ! gives 1.546
                  else if (((ncc.eq.3).and.(nco.eq.1)).or.
     &                    ((ncc.eq.1).and.(nco.eq.3))) then
                     bondl = bondl + 0.013                 ! gives 1.531
                  else if ((ncc.eq.2).and.(nco.eq.2)) then
                     bondl = bondl + 0.016                 ! gives 1.537
                  else if (((ncc.eq.2).and.(nco.eq.1)).or.
     &                    ((ncc.eq.1).and.(nco.eq.2))) then
                     bondl = bondl + 0.008                 ! gives 1.526
                  end if

c    Case B: N has 3 attached atoms
               else if (natto.eq.3) then
c    Correction to give the MM3 = 1.462 for least hindered case
                  bondl = bondl - 0.019                    ! gives 1.462
c    Steric corrections based on MM3 calculations
                  if ((ncc.eq.3).and.(nco.eq.2)) then
                     bondl = bondl + 0.022                 ! gives 1.484
                  else if ((ncc.eq.3).and.(nco.le.1)) then
                     bondl = bondl + 0.011                 ! gives 1.473
                  else if (ncc.eq.2) then
                     bondl = bondl + 0.007                 ! gives 1.469
                  end if

c    Case C: N has 2 attached atoms
               else if (natto.eq.2) then
c    Correction to give the MM3 = 1.447 for least hindered case
                  bondl = bondl - 0.034                    ! gives 1.447
c    Steric corrections based on MM3 calculations
                  if(ncc.eq.3) then
                     bondl = bondl + 0.009                 ! gives 1.456
                  else if (ncc.eq.2) then
                     bondl = bondl + 0.005                 ! gives 1.452
                  end if
               end if

c Case 2: Bound to O
            else if (atlabo.eq.'O ') then
c Blom-Haaland = 1.425 for a C-O single bond
c Correction to give the MM3 = 1.421 for least hindered case
               bondl = bondl - 0.004                       ! gives 1.421
c Steric corrections based on MM3 calculations
               if (ncc.eq.3) then
                  bondl = bondl + 0.008                    ! gives 1.429
               else if (ncc.eq.2) then
                  bondl = bondl + 0.005                    ! gives 1.426
               end if

c Case 3: Bound to P
            else if (atlabo.eq.'P ') then
c Blom-Haaland = 1.853 for a C-P single bond

c    Case A: P has 4 attached atoms
               if (natto.eq.4) then
c    Correction to give the MM3 = 1.812 for least hindered case
                  bondl = bondl - 0.041                    ! gives 1.812
c    Steric corrections based on MM3 calculations
                  if (ncc.eq.3) then
                     bondl = bondl + 0.028                 ! gives 1.840
                  else if (ncc.eq.2) then
                     bondl = bondl + 0.017                 ! gives 1.829
                  end if

c    Case B: P has 3 attached atoms
               else if (natto.eq.3) then
c    Correction to give the MM3 = 1.855 for least hindered case
                  bondl = bondl + 0.002                    ! gives 1.855
c    Steric corrections based on MM3 calculations
                  if (ncc.eq.3) then
                     bondl = bondl + 0.024                 ! gives 1.879
                  else if (ncc.eq.2) then
                     bondl = bondl + 0.015                 ! gives 1.870
                  end if
               end if

c Case 4: Bound to S       
            else if (atlabo.eq.'S ') then
c Blom-Haaland = 1.818 for a C-S single bond

c    Case A: S has 4 attached atoms
               if (natto.eq.4) then
c    Correction to give the MM3 = 1.787 for least hindered case
                  bondl = bondl - 0.031                    ! gives 1.787
c    Steric corrections based on MM3 calculations
                  if (ncc.eq.3) then
                     bondl = bondl + 0.030                 ! gives 1.817
                  else if (ncc.eq.2) then
                     bondl = bondl + 0.015                 ! gives 1.802
                  end if

c    Case B: S has 3 attached atoms
               else if (natto.eq.3) then
c    Correction to give the MM3 = 1.812 for least hindered case
                  bondl = bondl - 0.006                    ! gives 1.812
c    Steric corrections based on MM3 calculations
                  if (ncc.eq.3) then
                     bondl = bondl + 0.034                 ! gives 1.846
                  else if (ncc.eq.2) then
                     bondl = bondl + 0.015                 ! gives 1.827
                  end if

c    Case C: S has 2 attached atoms
               else if (natto.eq.2) then
c    Correction to give the MM3 = 1.816 for least hindered case
                  bondl = bondl - 0.002                    ! gives 1.816
                  if (ncc.eq.3) then
                     bondl = bondl + 0.024                 ! gives 1.840
                  else if (ncc.eq.2) then
                     bondl = bondl + 0.014                 ! gives 1.830
                  end if
               end if
            end if

         else if (nattc.eq.3) then                         ! sp2 carbon

c Case 1: Bound to N
            if (atlabo.eq.'N ') then

c    Case A: N has 4 attached atoms
               if (natto.eq.4) then
c    Correction to give the MMFF94 = 1.459 for least hindered case
                  bondl = bondl - 0.022                    ! gives 1.459
c    Steric corrections based on MMFF94 calculations
                  if (nco.eq.3) then
                     bondl = bondl + 0.065                 ! gives 1.524
                  else if (nco.eq.2) then
                     bondl = bondl + 0.021                 ! gives 1.479
                  end if

c    Case B: N has 3 attached atoms
               else if (natto.eq.3) then
c    Correction to give the MM3 = 1.405 for least hindered case
                  bondl = bondl - 0.076                    ! gives 1.405
c    Steric corrections based on MM3 calculations
                  bondl = bondl + 0.007*boost
                  if (nco.eq.2) then
                     bondl = bondl + 0.005                 ! gives 1.410
                  end if

c    Case C: N has 2 attached atoms
               else if (natto.eq.2) then
c    Correction to give the MM3 = 1.399 for least hindered case
                  bondl = bondl - 0.082                    ! gives 1.399
c    Steric corrections based on MM3 calculations
                  if (ncc.eq.2) then
                    bondl = bondl + 0.012
                  end if
                  bondl = bondl + 0.009*boost
               end if

c Case 2: Bound to O
            else if (atlabo.eq.'O ') then
c Correction to give the MM3 = 1.362 for least hindered case
               bondl = bondl - 0.063                       ! gives 1.362
c Steric corrections based on MM3 calculations
               if (ncc.eq.2) then
                  bondl = bondl + 0.005                    ! gives 1.367
                  if (boost.eq.2) then
                     bondl = bondl + 0.036                 ! gives 1.403
                 end if
               end if
          
c Case 3: Bound to P
            else if (atlabo.eq.'P ') then

c    Case A: P has 4 attached atoms
               if (natto.eq.4) then
c    Correction to give the MM3 = 1.799 for least hindered case
                  bondl = bondl - 0.054                    ! gives 1.799
c    Steric corrections based on MM3 calculations
                  if (ncc.eq.2) then
                     bondl = bondl + 0.010                 ! gives 1.809
                     bondl = bondl + 0.015*boost
                  end if
                
c    Case B: P has 3 attached atoms
               else if (natto.eq.3) then
c    Correction to give the MM3 = 1.817 for least hindered case
                  bondl = bondl - 0.036                    ! gives 1.817
                  if (ncc.eq.2) then
                     bondl = bondl + 0.015                 ! gives 1.832
                     bondl = bondl + 0.005*boost  
                  end if
               end if

c Case 4: Bound to S       
            else if (atlabo.eq.'S ') then

c    Case A: S has 4 attached atoms
               if (natto.eq.4) then
c    Correction to give the MM3 = 1.769 for least hindered case
                  bondl = bondl - 0.049                    ! gives 1.769
                  if (ncc.eq.2) then
                     bondl = bondl + 0.015                 ! gives 1.784
                     bondl = bondl + 0.015*boost
                  end if

c    Case B: S has 3 attached atoms
               else if (natto.eq.3) then
c    Correction to give the MM3 = 1.752 for least hindered case
                  bondl = bondl - 0.066                    ! gives 1.752
                  if (ncc.eq.2) then
                     bondl = bondl + 0.054                 ! gives 1.806
                     bondl = bondl + 0.015*boost
                  end if

c    Case C: S has 2 attached atoms
               else if (natto.eq.2) then
c    Correction to give the MM3 = 1.729 for least hindered case
                  bondl = bondl - 0.089                    ! gives 1.729
                  if (ncc.eq.2) then
                     bondl = bondl + 0.031                 ! gives 1.760
                     bondl = bondl + 0.020*boost           
                  end if

               end if
            end if

         else if (nattc.eq.2) then                         ! sp carbon
c Due to lack of parameters, assignments here are based on x-ray data

            if (atlabo.eq.'N ') then
               bondl = bondl - 0.151                       ! gives 1.330

            else if (atlabo.eq.'O ') then
               bondl = bondl - 0.142                       ! gives 1.283

            else if (atlabo.eq.'P ') then
               bondl = bondl - 0.083                       ! gives 1.770 

            else if (atlabo.eq.'S ') then
               bondl = bondl - 0.114                       ! gives 1.704

            end if

         end if

c-----------------------------------------------------------------------
c     bottom of decision tree
c-----------------------------------------------------------------------

      end if

      return
      end 
