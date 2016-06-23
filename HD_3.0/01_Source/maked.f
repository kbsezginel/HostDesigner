c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'maked' -- creates MOLECULE D from HOST B + MOLECULE C
c
c     remove terminal atoms
c     renumber link atoms and combine the bond lists
c     assign all molecule c variables
c     generate an exclusion list for molecule c
c
c     passed variables:
c     integer   lostb - current value of looserb(mm) 
c     integer   mm
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     
      
      subroutine maked(lostb, mm)
      
      implicit none
      
      include 'params.i'
      include 'hosta.i'
      include 'hostb.i'
      include 'moleculec.i'
      include 'moleculed.i'
      
      integer i, j
      integer lostb, mm

c-----------------------------------------------------------------------
c     move molecule C into molecule D
c-----------------------------------------------------------------------

      do i = 1, looserc - 1
      
         radd(i) = radc(i)
         atnumd(i) = atnumc(i)
         atwtd(i) = atwtc(i)
         ityped(i) = itypec(i)
         atlabd(i) = atlabc(i)
         iguestd(i) = iguestc(i)
         
         listd(i) = listc(i)
         
         n12d(i) = n12c(i)
         do j = 1, n12d(i)
            i12d(j,i) = i12c(j,i)
         end do
         
      end do

c-----------------------------------------------------------------------
c     assign serial numbers of the molecule C half of the rotatable bond 
c     that will be used in the dihedral driver 
c-----------------------------------------------------------------------

      idihedd = torrefc
      jdihedd = ibonderc
      
      if(idihedd.gt.looserc) idihedd = idihedd - 1
      if(jdihedd.gt.looserc) jdihedd = jdihedd - 1
                   
c-----------------------------------------------------------------------
c     remove terminal atom from molecule C list
c-----------------------------------------------------------------------

      do i = looserc, nc - 1
      
         radd(i) = radc(i+1)
         atnumd(i) = atnumc(i+1)
         atwtd(i) = atwtc(i+1)
         ityped(i) = itypec(i+1)
         atlabd(i) = atlabc(i+1)
         iguestd(i) = iguestc(i+1)
         
         listd(i) = listc(i+1)
         
         n12d(i) = n12c(i+1)
         do j = 1, n12d(i)
            i12d(j,i) = i12c(j,i+1)
         end do
         
      end do      
         
      nd = nc - 1

      do i = 1, nd

         do j = 1, n12d(i)

            if(i12d(j,i).gt.looserc) then
               i12d(j,i) = i12d(j,i) - 1
            end if
            
         end do
         
      end do
      
c-----------------------------------------------------------------------
c     move host B into molecule D
c-----------------------------------------------------------------------

      do i = 1, lostb - 1
      
         radd(i + nd) = radb(i)
         atnumd(i + nd) = atnumb(i)
         atwtd(i + nd) = atwtb(i)
         ityped(i + nd) = itypeb(i)
         atlabd(i + nd) = atlabb(i)
         iguestd(i + nd) = iguestb(i)
         
         listd(i + nd) = listb(i,mm)
         
         n12d(i + nd) = n12b(i)
         do j = 1, n12d(i + nd)
            i12d(j, i + nd) = i12b(j, i) + nd
         end do
         
      end do

c-----------------------------------------------------------------------
c     remove lost hydrogen atom from host B part
c-----------------------------------------------------------------------

      do i = lostb, nb - 1
      
         radd(i+nd) = radb(i + 1)
         atnumd(i+nd) = atnumb(i + 1)
         atwtd(i+nd) = atwtb(i + 1)
         ityped(i+nd) = itypeb(i + 1)
         atlabd(i+nd) = atlabb(i + 1)
         iguestd(i+nd) = iguestb(i + 1)
         
         listd(i+nd) = listb(i + 1, mm)
         
         n12d(i+nd) = n12b(i + 1)
         do j = 1, n12d(i+nd)
            i12d(j,i+nd) = i12b(j, i + 1) + nd
         end do
         
      end do        
      
      do i = nd + 1, nd + nb - 1
      
         do j = 1, n12d(i)

            if(i12d(j,i).gt.nd+lostb) then
               i12d(j,i) = i12d(j,i) - 1
            end if
            
         end do
            
      end do
            
      nd = nd + nb - 1 
      
c-----------------------------------------------------------------------
c     assign serial numbers of the host B half of the rotatable bond 
c     that will be used in the dihedral driver 
c-----------------------------------------------------------------------

      kdihedd = ibonderb(mm) + nc - 1
      ldihedd = torrefb(mm) + nc - 1
      alt_ldihedd = altrefb(mm) + nc -1
      
      if(ibonderb(mm).gt.lostb) kdihedd = kdihedd - 1
      if(torrefb(mm).gt.lostb) ldihedd = ldihedd - 1
      if(altrefb(mm).gt.lostb) alt_ldihedd = alt_ldihedd - 1
                  
c-----------------------------------------------------------------------
c     fix connectivity of bonding atoms between molecule C and host B
c-----------------------------------------------------------------------

      do i = 1, n12d(jdihedd)
      
         if(i12d(i,jdihedd).eq.looserc) then
         
            i12d(i,jdihedd) = kdihedd
            
            go to 30
            
         end if
         
      end do
      
   30 continue
      
      do i = 1, n12d(kdihedd)
      
         if(i12d(i,kdihedd).eq.lostb + nc - 1) then
         
            i12d(i,kdihedd) = jdihedd
            
            go to 40

         end if
         
      end do
      
   40 continue
      
      return
      end
