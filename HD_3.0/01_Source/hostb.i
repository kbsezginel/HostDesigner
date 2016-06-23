c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     hostb.i -- input parameters for host B
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c     titleb       descriptive header of host B file
c     nb           total number of atoms
c     atnumb       atomic number
c     atwtb        atomic weight
c     atlabb       atom label
c     iguestb      tag to identify guest atoms
c     xb           x-coordinates 
c     yb           y-coordinates
c     zb           z-coordinates
c     radb         van der Waals radii
c     itypeb       MM3 atom type
c     n12b         total number of atoms directly bonded to each atom
c     i12b         serial numbers of atoms bonded to each atom
c     nattachb     total number of attachment sites
c     ndriveb      total number of drives per attachment site
c     looserb      serial number of atom to be replaced by link
c     ibonderb     serial number of atom that forms bond with link
c     needlabb     linker must attach with this atom label
c     needvalb     linker must attach with this atom valence
c     mirrorb      .true. if chiral
c     listb        connectivity values
c     driveb       .true. if host b is to be driven
c     nrotb        number rotatable bonds between guest and attachment
c     classb       class of rotor
c     subclassb    subclass of rotor
c     torrefb      atom used to define rotor dihedral angle
c     altrefb      alternate atom used to define rotor dihedral angle


      character*90 driveinfob
      character*60 titleb     
      character*2  atlabb, needlabb
      
      real   xb, yb, zb, radb, atwtb

      integer      nb, itypeb, iguestb, atnumb
      integer      n12b, i12b, listb
      integer      nattachb, ndriveb, looserb, ibonderb, needvalb
      integer      classb, sclassb, torrefb, altrefb, nrotb
      
      logical      mirrorb, driveb


      common /hostbb1/ xb(MAXATOMB,MAXDRIVEB,MAXATTACHB),
     &                 yb(MAXATOMB,MAXDRIVEB,MAXATTACHB),
     &                 zb(MAXATOMB,MAXDRIVEB,MAXATTACHB), 
     &                 radb(MAXATOMB), atwtb(MAXATOMB)
     
      common /hostbchar/ titleb, atlabb(MAXATOMB), needlabb(MAXATTACHB),
     &        driveinfob(MAXDRIVEB)
     
      common /hostbb2/ nb, atnumb(MAXATOMB), iguestb(MAXATOMB), 
     &        itypeb(MAXATOMB)

      common /conn12b/ n12b(MAXATOMB), i12b(MAXVAL,MAXATOMB), 
     &        listb(MAXATOMB,MAXATTACHB)

      common /attachb/ nattachb, ndriveb, looserb(MAXATTACHB), 
     &        ibonderb(MAXATTACHB), needvalb(MAXATTACHB), 
     &        classb(MAXATTACHB),sclassb(MAXATTACHB),
     &        torrefb(MAXATTACHB),altrefb(MAXATTACHB),nrotb(MAXATTACHB)
    
      common /logb/ mirrorb, driveb
