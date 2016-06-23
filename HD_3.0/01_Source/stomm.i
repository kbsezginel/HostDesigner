c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     stomm.i -- molecular information pertaining to mmruns
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c     mmnatom      total number of atoms
c     mmatlab      atom label
c     mmguest      tag to identify guest atoms
c     mmx          x-coordinates 
c     mmy          y-coordinates
c     mmz          z-coordinates
c     mmtype       force field atom type
c     mmn12        total number of atoms directly bonded to each atom
c     mmi12        serial numbers of atoms directly bonded to each atom
c     lignatom     number of atoms in ligand
c     ligatlab     atom lables in ligand
c     ligx         x-coordinates of the ligand
c     ligy         y-coordinates of the ligand
c     ligz         z-coordinates of the ligand
c     ligtype      force field atom type 
c     lign12       total number of atoms directly bonded to each atom
c     ligi12       serial numbers of atoms directly bonded to each atom
c     mmnsb        number of searchable bonds
c     mmsbi        ith atom of searchable bond
c     mmsbj        jth atom of searchable bond
c     lignsb        number of searchable bonds
c     ligsbi        ith atom of searchable bond
c     ligsbj        ith atom of searchable bond

c     key          pointer giving the order of taking hits
c     rotbonds     number of rotatable bonds
c     linkid       serial number of the link
c     miss1        =0 if OK, =1 if missing params for host-guest
c     miss2        =0 if OK, =1 if missing params for host
c     fail1        =0 if OK, =1 if calc on host-guest not complete
c     fail2        =0 if OK, =1 if calc on host not complete
c     deltaE1      change in energy on optimization of host-guest
c     deltaE2      change in energy of host after removing guest
c     enlig        estimated conformer energy of the link
c     enbond       estimated energy of the bonds formed
c     e2lig        steric energy of host after optimization
c     h2lig        heat of formation of host after optimization
c     identity     string giving the link name
c     istart       start of the link
c     iend         end of the link

c     entot        ???
c     rmsd0
c     rmsd1
c     rmsd2
c     rmsd3
c     crmsd
c     startwith
c     endwith

c     confkey
c     mmrank
c     hdrank
c     rotme
c     serlink
c     missin1
c     missin2
c     failin1
c     failin2  
c     Ecomp
c     Econf
c     Eest
c     Etot
c     srchfail
c     icount

c     nprint       number of conformers to be printed
c     connatom     number of atoms
c     conx         coordinates
c     cony         coordinates
c     conz         coordinates
c     contype      force field atom type 
c     conn12       total number of atoms directly bonded to each atom
c     coni12       serial numbers of atoms directly bonded to each atom
c     conelist     steric energies of conformers
c     conhlist     heats of formation of conformers

      character*2  mmatlab, ligatlab, conatlab
      character*37 identity
      
      logical srchfail      

      real  mmx, mmy, mmz
      real  ligx, ligy, ligz 
      real  conx, cony, conz
      real  deltaE1, deltaE2, e2lig, h2lig
      real  enlig, enbond, conelist, conhlist
      real rmsd0, rmsd1, rmsd2, rmsd3, crmsd, entot
      real Ecomp, Econf, Eest, Etot
      
      integer  mmnatom, mmtype, mmguest, mmn12, mmi12
      integer  lignatom, ligtype, lign12, ligi12
      integer  connatom, contype, conn12, coni12
      integer  key, rotbonds, miss1, miss2, fail1, fail2
      integer  linkid, nprint, istart, iend
      integer  startwith, endwith
      integer confkey, mmrank, hdrank, rotme, serlink
      integer missin1, missin2, failin1, failin2, icount
      integer mmnsb, mmsbi, mmsbj, lignsb, ligsbi, ligsbj


      common /mmmol/ mmx(MAXMM,MAXATOMD), mmy(MAXMM,MAXATOMD), 
     &        mmz(MAXMM,MAXATOMD), ligx(MAXMM,MAXATOMD), 
     &        ligy(MAXMM,MAXATOMD), ligz(MAXMM,MAXATOMD),
     &        conx(MAXCONF,MAXSRCH,MAXATOMD),
     &        cony(MAXCONF,MAXSRCH,MAXATOMD),
     &        conz(MAXCONF,MAXSRCH,MAXATOMD),
     &        conelist(MAXCONF,MAXSRCH), conhlist(MAXCONF,MAXSRCH)
    
      common /mmmol2/ mmnatom(MAXMM), mmtype(MAXMM,MAXATOMD),
     &        mmguest(MAXMM,MAXATOMD), lignatom(MAXMM),
     &        ligtype(MAXMM,MAXATOMD), deltaE1(MAXMM),
     &        deltaE2(MAXMM), e2lig(MAXMM), h2lig(MAXMM),
     &        enlig(MAXMM), enbond(MAXMM),
     &        connatom(MAXSRCH), contype(MAXSRCH,MAXATOMD),
     &        rmsd0(MAXMM), rmsd1(MAXMM), rmsd2(MAXMM),
     &        rmsd3(MAXMM), entot(MAXMM), crmsd(MAXSRCH),
     &        Ecomp(MAXSRCH), Econf(MAXSRCH), Eest(MAXSRCH),
     &        Etot(MAXSRCH) 
 
      common /mmconn/ mmn12(MAXMM,MAXATOMD),lign12(MAXMM,MAXATOMD),
     &        mmi12(MAXMM,MAXVAL,MAXATOMD), 
     &        ligi12(MAXMM,MAXVAL,MAXATOMD),
     &        conn12(MAXSRCH,MAXATOMD),
     &        coni12(MAXSRCH,MAXVAL,MAXATOMD),
     &        key(MAXHITS), rotbonds(MAXMM), miss1(MAXMM),
     &        miss2(MAXMM), fail1(MAXMM), fail2(MAXMM),
     &        linkid(MAXMM), nprint(MAXSRCH),
     &        istart(MAXMM), iend(MAXMM),
     &        startwith, endwith,
     &        confkey(MAXSRCH), mmrank(MAXSRCH), hdrank(MAXSRCH), 
     &        rotme(MAXSRCH), serlink(MAXSRCH), missin1(MAXSRCH),
     &        missin2(MAXSRCH), failin1(MAXSRCH),
     &        failin2(MAXSRCH), icount, mmnsb, mmsbi(MAXSBD), 
     &        mmsbj(MAXSBD), lignsb(MAXMM), ligsbi(MAXMM,MAXATOMD),
     &        ligsbj(MAXMM,MAXATOMD)

      common /mmchar/ mmatlab(MAXMM,MAXATOMD),
     &        ligatlab(MAXMM,MAXATOMD),
     &        conatlab(MAXSRCH,MAXATOMD),
     &        identity(MAXMM)

      common /conflog/ srchfail(MAXSRCH)

