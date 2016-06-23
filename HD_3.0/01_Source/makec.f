c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'makec' -- creates MOLECULE C from HOST A + LINK
c
c     remove terminal atoms
c     renumber link atoms and combine the bond lists
c     assign all molecule c variables
c     generate an exclusion list for molecule c
c
c     passed variables:
c     integer   losta     current value of loosera(jj) 
c     integer   ii        number of link to be used
c     integer   jj        number of attachment to link ii to be used
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine makec(ii,jj)
      
      implicit none
      
      include 'params.i'
      include 'hosta.i'
      include 'linkage.i'
      include 'moleculec.i'
      
      integer i, j, ii, jj
      integer losta

      losta = loosera(jj)

c-----------------------------------------------------------------------
c     move host A into molecule C
c-----------------------------------------------------------------------

      do i = 1, na
      
         radc(i) = rada(i)
         atnumc(i) = atnuma(i)
         atwtc(i) = atwta(i)
         itypec(i) = itypea(i)
         atlabc(i) = atlaba(i)
         iguestc(i) = iguesta(i)
         
         listc0(i) = lista(i,jj)
         
         n12c(i) = n12a(i)
         do j = 1, n12c(i)
            i12c(j,i) = i12a(j,i)
         enddo

      enddo

c-----------------------------------------------------------------------
c     if we get the NULL link, then ...
c-----------------------------------------------------------------------
      
      if(ng(ii).eq.0) then
      
         idihedc = 0
         jdihedc = 0
         kdihedc = 0
         ldihedc = 0
         
         looserc = losta
         ibonderc = ibondera(jj)
         torrefc= torrefa(jj)
         alt_torrefc = altrefa(jj)
                  
         atlabgb = atlabc(ibonderc)
         ivalengb = n12a(ibonderc) - 1
         classgb = classa(jj)
         sclassgb = sclassa(jj)

         nc = na
         
         return
         
      end if
            
c-----------------------------------------------------------------------
c     assign serial numbers of the host A half of the rotatable bond 
c     that will be used in the dihedral driver 
c-----------------------------------------------------------------------

      idihedc = torrefa(jj)      
      alt_idihedc = altrefa(jj)
      jdihedc = ibondera(jj)

      if(idihedc.gt.losta) idihedc = idihedc - 1
      if(alt_idihedc.gt.losta) alt_idihedc = alt_idihedc -1
      if(jdihedc.gt.losta) jdihedc = jdihedc - 1
             
c-----------------------------------------------------------------------
c     remove the lost hydrogen atom from host A list 
c-----------------------------------------------------------------------

      do i = losta, na
      
         radc(i) = radc(i+1)
         atnumc(i) = atnumc(i+1)
         atwtc(i) = atwtc(i+1)
         itypec(i) = itypec(i+1)
         atlabc(i) = atlabc(i+1)
         iguestc(i) = iguestc(i+1)
         
         listc0(i) = listc0(i+1)
         
         n12c(i) = n12c(i+1)
         do j = 1, n12c(i)
            i12c(j,i) = i12c(j,i+1)
         end do
         
      end do   
      
      nc = na - 1
      
      do i = 1, nc
         do j = 1, n12c(i)
            if(i12c(j,i).gt.losta) then
               i12c(j,i) = i12c(j,i) - 1
            end if
         end do
      end do
                  
c-----------------------------------------------------------------------
c     move link into molecule C
c-----------------------------------------------------------------------

      do i = 1, ng (ii)
      
         radc(i + nc) = radg(i)
         atnumc(i + nc) = atnumg(i)
         atwtc(i + nc) = atwtg(i)
         itypec(i + nc) = itypeg(i)
         atlabc(i + nc) = atlabg(i)
         iguestc(i + nc) = 0
         
         listc0(i + nc) = listg(i)
         
         n12c(i + nc) = n12g(i)
         do j = 1, n12c(i + nc)
            i12c(j,i + nc) = i12g(j,i) + nc
         end do
         
      end do
      
c-----------------------------------------------------------------------
c     remove terminal atom from link list 
c-----------------------------------------------------------------------

      do i = nc + looserga, nc + ng(ii)
      
         radc(i) = radc(i + 1)
         atnumc(i) = atnumc(i + 1)
         atwtc(i) = atwtc(i + 1)
         itypec(i) = itypec(i + 1)
         atlabc(i) = atlabc(i + 1)
         iguestc(i) = iguestc(i + 1)
         
         listc0(i) = listc0(i+1)
         
         n12c(i) = n12c(i+1)
         do j = 1, n12c(i)
            i12c(j,i) = i12c(j,i+1)
         end do
         
      end do         
      
      do i = nc + 1, nc + ng(ii) - 1
      
         do j = 1, n12c(i)
         
            if(i12c(j,i).gt.nc+looserga) then
               i12c(j,i) = i12c(j,i) - 1
            end if
            
         end do
         
      end do
               
      nc = nc + ng(ii) - 1 
            
c-----------------------------------------------------------------------
c     assign serial numbers of the linker half of the rotatable bond 
c     that will be used in the dihedral driver 
c-----------------------------------------------------------------------

      kdihedc = ibonderga + na - 1
      ldihedc = torrefga  + na - 1

      if(ibonderga.gt. looserga) kdihedc = kdihedc - 1
      if(torrefga .gt. looserga) ldihedc = ldihedc - 1
                   
c-----------------------------------------------------------------------
c     fix connectivity lists of bonding atoms between host A and link
c-----------------------------------------------------------------------

      do i = 1, n12c(jdihedc)
      
         if(i12c(i,jdihedc).eq.losta) then
         
            i12c(i,jdihedc) = kdihedc
            
            go to 30
            
         end if
         
      end do
      
   30 continue
      
      do i = 1, n12c(kdihedc)
      
         if(i12c(i,kdihedc).eq.looserga + na - 1) then
         
            i12c(i,kdihedc) = jdihedc
            
            go to 40

         end if
         
      end do
      
   40 continue
         
c-----------------------------------------------------------------------
c     reassign serial numbers of the attachment site for B  
c-----------------------------------------------------------------------
      
      looserc  = loosergb + na - 1
      ibonderc = ibondergb + na - 1
 
      if (loosergb .gt. looserga) looserc  = looserc  - 1
      if (ibondergb.gt. looserga) ibonderc = ibonderc - 1

      if(torrefgb.ne.looserga) then
        torrefc  = torrefgb + na - 1
        if (torrefgb .gt. looserga) torrefc  = torrefc  - 1
      else
        torrefc = ibondera(jj)
      endif
      
      return

      end
