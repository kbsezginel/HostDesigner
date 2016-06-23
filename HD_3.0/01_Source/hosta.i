c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     hosta.i -- input parameters for host A
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c     titlea       descriptive header of host A file
c     na           total number of atoms
c     atnuma       atomic number
c     atwta        atomic weight
c     atlaba       atom label
c     iguesta      tag to identify guest atoms
c     xa           x-coordinates 
c     ya           y-coordinates
c     za           z-coordinates
c     rada         van der Waals radii
c     itypea       MM3 atom type
c     n12a         total number of atoms directly bonded to each atom
c     i12a         serial numbers of atoms bonded to each atom
c     nattacha     total number of attachment sites
c     ndrivea      total number of drives per attachment site
c     loosera      serial number of atom to be replaced by link
c     ibondera     serial number of atom that forms bond with link
c     needlaba     linker must attach with this atom label
c     needvala     linker must attach with this atom valence
c     mirrora      .true. if chiral
c     lista        connectivity list
c     drivea       .true. if host a is to be driven
c     nrota        number rotatable bonds between guest and attachment
c     nsymattacha  number of additional symmetrical attachment sites
c     symloosera   serial number of symm atom to be replaced by link
c     symibondera  serial number of symm atom that forms bond with link
c     classa       class of rotor
c     subclassa    subclass of rotor
c     torrefa      atom used to define rotor dihedral angle
c     altrefa      alternate atom used to define rotor dihedral angle


      character*90 driveinfoa
      character*60 titlea     
      character*2  atlaba, needlaba
      
      real  xa, ya, za, rada, atwta
      
      integer      na, itypea, iguesta, atnuma
      integer      n12a,i12a, lista
      integer      nattacha, ndrivea, loosera, ibondera, needvala
      integer      classa, sclassa, torrefa, altrefa, nrota
      integer      nsymattacha, symloosera, symibondera
      
      logical      mirrora, drivea


      common /hostaa1/ xa(MAXATOMA,MAXDRIVEA,MAXATTACHA), 
     &                 ya(MAXATOMA,MAXDRIVEA,MAXATTACHA), 
     &                 za(MAXATOMA,MAXDRIVEA,MAXATTACHA), 
     &                 rada(MAXATOMA), atwta(MAXATOMA)
     
      common /hostachar/ atlaba(MAXATOMA), titlea, needlaba(MAXATTACHA),
     &        driveinfoa(MAXDRIVEA)
     
      common /hostaa2/ na, atnuma(MAXATOMA), iguesta(MAXATOMA), 
     &        itypea(MAXATOMA)
    
      common /conn12a/ n12a(MAXATOMA), i12a(MAXVAL,MAXATOMA), 
     &        lista(MAXATOMA,MAXATTACHA)

      common /attacha/ nattacha, ndrivea, loosera(MAXATTACHA), 
     &        ibondera(MAXATTACHA), needvala(MAXATTACHA), 
     &        classa(MAXATTACHA),sclassa(MAXATTACHA),
     &        torrefa(MAXATTACHA),altrefa(MAXATTACHA),nrota(MAXATTACHA),
     &        nsymattacha, symloosera(MAXSYMATTACHA, MAXATTACHA),
     &        symibondera(MAXSYMATTACHA, MAXATTACHA)

      common /loga/ mirrora, drivea
