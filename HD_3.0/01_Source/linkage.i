c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     linkage.i -- input parameters for current link
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c     titleg       descriptive header of link
c     ng           total number of atoms
c     atnumg       atomic number
c     atwtg        atomic weight
c     atlabg       atom label
c     xg           x-coordinates 
c     yg           y-coordinates
c     zg           z-coordinates
c     radg         van der Waals radii
c     itypeg       MM3 atom type
c     n12g         total number of atoms directly bonded to each atom
c     i12g         serial numbers of atoms bonded to each atom
c     looserga     serial number of first atom to be replaced by host
c     loosergb     serial number of second atom to be replaced by host
c     ibonderga    serial number of atom that forms first bond with host
c     ibondergb    serial number of atom that forms second bond with host
c     atlabga      atom label of first connecting atom
c     atlabgb      atom label of second connecting atom
c     ivalenga     hybridization of first connecting atom
c     ivalengb     hybridization of second connecting atom
c     conng        connectivity in linker
c     symmetry     true if swapping connections yields same link
c     mirrorg      true if link lacks a superimposable mirror image
c     confeg       conformational energy of link
c     distbg       distance between the two connecting atoms
c     angag        angle for loosera-ibondera-distbg
c     angbg        angle for looserb-ibonderb-distbg
c     twistg       dihedral angle of loosera-ibondera-ibonderb-looserb
c     listg        connectivity list
c     listg0       connectivity list used by overlay routine
c     linknum      number as input in library of link ii

c     useclass     logical to use class or not
c     minconn      minimum connectivity allowed
c     maxconn      maximum connectivity allowed
c     maxconfe     maximum conformational energy allowed 
c     maxnrot      maximum rotatable bonds allowed in link

      character*37 linkname
      character*2  atlabg, atlabga, atlabgb

      real       xg, yg, zg, radg, atwtg
      real       distbg, angag, angbg, twistg
      real       maxconfe, confeg

      integer      ng, itypeg, atnumg
      integer      n12g,i12g, conng, listg, listg0
      integer      looserga, loosergb, ibonderga, ibondergb
      integer      ivalenga, ivalengb
      integer      minconn, maxconn, maxnrot
      integer      classga, classgb, sclassga, sclassgb
      integer      nrot, torrefga, torrefgb
      integer      linknum

      logical      symmetry, mirrorg, useclass

      common /linkcha1/ atlabg(MAXATOMG), 
     &         atlabga, atlabgb
      common /linkcha2/ linkname(MAXLINKS)
           
      common /comlink1/ xg(MAXATOMG), yg(MAXATOMG),
     &         zg(MAXATOMG),atwtg(MAXATOMG), radg(MAXATOMG)
     
      common /comlink2/ ng(MAXLINKS), itypeg(MAXATOMG), 
     &         atnumg(MAXATOMG),looserga, 
     &         loosergb, ibonderga, 
     &         ibondergb,
     &         minconn, maxconn, maxnrot, conng,
     &         ivalenga, ivalengb,
     &         nrot(MAXLINKS), linknum(MAXLINKS)
    
      common /comlink3/ maxconfe, confeg(MAXLINKS)

      common /overlayg/ distbg, twistg, 
     &         angag, angbg

      common /connn12g/ n12g(MAXATOMG), 
     &         i12g(MAXVAL,MAXATOMG), 
     &         listg(MAXATOMG), listg0(MAXATOMG),
     &         classga, sclassga,
     &         classgb, sclassgb,
     &         torrefga, torrefgb

      common /linkcont/ symmetry, mirrorg ,useclass

