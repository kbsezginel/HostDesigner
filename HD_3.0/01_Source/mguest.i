c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     mguest.i -- parameters multi-atom guests
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c     lguest     number of atoms in guest
c     nequiva    number of symmetry equivalent placements of b onto a
c     mapeqa     serial number lists of all equivalent maps of b onto a
c     nequivb    number of symmetry equivalent placements of a onto b
c     mapeqb     serial number lists of all equivalent maps of a onto b
c     gstn12a    total number of atoms of b directly or symmetry equivalent 
c                    to those bonded to each guest atom on a 
c     gsti12a    serial numbers of atoms of b directly or symmetry equivalent 
c                    to those bonded to each guest atom on a 
c     gstbmpa    maximum length allowed from guest a atoms to bound atoms
c                    or atoms symmetry equivalent on the complex fragment b
c     gstn12b    total number of atoms of a directly or symmetry equivalent 
c                    to those bonded to each guest atom on b 
c     gsti12b    serial numbers of atoms of a directly or symmetry equivalent 
c                    to those bonded to each guest atom on b 
c     gstbmpb    maximum length allowed from guest b atoms to bound atoms
c                    or atoms symmetry equivalent on the complex fragment a

      real  gstbmpa,gstbmpb

      integer           lguest,gstn12a,gsti12a,mapeqa,nequiva
      integer           gstn12b,gsti12b,mapeqb,nequivb
      
      common /guesti/ lguest,nequiva,nequivb,gstn12a(MAXGST),
     &                gstn12b(MAXGST),gsti12a(MAXVAL,MAXGST),
     &                gsti12b(MAXVAL,MAXGST),
     &                mapeqa(1000,MAXGST),mapeqb(1000,MAXGST)

      common /guestd/ gstbmpa(MAXATOMA,MAXDRIVEA,MAXGST)
     &               ,gstbmpb(MAXATOMB,MAXDRIVEB,MAXGST)

