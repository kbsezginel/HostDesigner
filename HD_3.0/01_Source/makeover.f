c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'makeover' -- creates molecule formed by overlaying
c     a link with the host, stores this in arrays for molecule C.
c
c     if the null link, form one bond between b1a to b2a
c     else two bonds will be formed - b1a to b1g and b2a to b2g
c
c     passed variables:
c     integer   b1a    serial number of atom attached to first hydrogen
c     integer   h1     serial number of first hydrogen lost by host
c     integer   h2     serial number of 2nd hydrogen lost by host
c     integer   b2a    serial number of atom attached to 2nd hydrogen
c     integer   d1     serial number of first dummy atom lost by link
c     integer   b1g    serial number of atom attached to first dummy
c     integer   b2g    serial number of atom attached to 2nd dummy
c     integer   d2     serial number of 2nd dummy atom lost by link
c     integer   t1a    serial number of torsion reference point for b1a
c     integer   t2a    serial number of torsion reference point for b2a
c     integer   t1g    serial number of torsion reference point for b1g
c     integer   t2g    serial number of torsion reference point for b2g
c     real      xgu,ygu,zgu   coordinates for the link
c     real      xu,yu,zu      coordinates for the fragment
c     integer   ii     number of link to be used
c     integer   jj     number of current set of attachment points
c     integer   mm     number of symmetrical linkage attachment points
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------     

      subroutine makeover(bb1a,h1,h2,bb2a,d1,bb1g,bb2g,d2,tt1a,tt2a,
     &           tt1g,tt2g,xgu,ygu,zgu,xu,yu,zu,ii,jj,mm)
      
      implicit none
      
      include 'params.i'
      include 'hosta.i'
      include 'linkage.i'
      include 'moleculec.i'
      
      integer i, j, spot1, spot2, ii, jj, mm
      integer h1, b1a, h2, b2a, d1, b1g, d2, b2g, t1a, t2a, t1g, t2g
      integer bb1a, bb2a, bb1g, bb2g, tt1a, tt2a, tt1g, tt2g
      integer ncc, ngg
      
      real xgu(MAXATOMG),ygu(MAXATOMG),zgu(MAXATOMG)
      real xu(MAXATOMC),yu(MAXATOMC),zu(MAXATOMC)

      ngg = ng(ii)

      if(mm.eq.0) then
         ncc = na
      else
         ncc = nc
      endif
      
      b1a= bb1a
      b2a= bb2a
      b1g= bb1g
      b2g= bb2g
      t1a= tt1a
      t2a= tt2a
      t1g= tt1g
      t2g= tt2g
      
c-----------------------------------------------------------------------
c     move host A into molecule C, done first time this routine is used
c-----------------------------------------------------------------------

      if(mm.eq.0) then

         do i = 1, ncc

            xc(i) = xu(i)
            yc(i) = yu(i)
            zc(i) = zu(i)
            radc(i) = rada(i)
            atnumc(i) = atnuma(i)
            atwtc(i) = atwta(i)
            itypec(i) = itypea(i)
            atlabc(i) = atlaba(i)
            iguestc(i) = iguesta(i)
         
            n12c(i) = n12a(i)
            do j = 1, n12c(i)
               i12c(j,i) = i12a(j,i)
            end do

            listc(i) = lista(i,jj)
            listc0(i) = -lista(i,jj+1)
         
         end do

      endif 

c-----------------------------------------------------------------------
c     remove first lost hydrogen
c-----------------------------------------------------------------------

      spot1 = max(h1,h2)

      do i = spot1, ncc
         xc(i) = xc(i+1)
         yc(i) = yc(i+1)
         zc(i) = zc(i+1)
         radc(i) = radc(i+1)
         atnumc(i) = atnumc(i+1)
         atwtc(i) = atwtc(i+1)
         itypec(i) = itypec(i+1)
         atlabc(i) = atlabc(i+1)
         iguestc(i) = iguestc(i+1)
         
         listc(i) = listc(i+1)
         listc0(i) = listc0(i+1)
         
         n12c(i) = n12c(i+1)
         do j = 1, n12c(i)
            i12c(j,i) = i12c(j,i+1)
         end do
         
      end do       
      
      do i = 1, ncc - 1
      
         do j = 1, n12c(i)

            if(i12c(j,i).gt.spot1) then
               i12c(j,i) = i12c(j,i) - 1
            else if (i12c(j,i).eq.spot1) then
               i12c(j,i) = -1
            end if
            
         end do
         
      end do
      
c-----------------------------------------------------------------------
c     remove second lost hydrogen
c-----------------------------------------------------------------------

      spot2 = min(h1,h2)

      do i = spot2, ncc - 1

         xc(i) = xc(i+1)
         yc(i) = yc(i+1)
         zc(i) = zc(i+1)
         radc(i) = radc(i+1)
         atnumc(i) = atnumc(i+1)
         atwtc(i) = atwtc(i+1)
         itypec(i) = itypec(i+1)
         atlabc(i) = atlabc(i+1)
         iguestc(i) = iguestc(i+1)
         
         listc(i) = listc(i+1)
         listc0(i) = listc0(i+1)
         
         n12c(i) = n12c(i+1)
         do j = 1, n12c(i)
            i12c(j,i) = i12c(j,i+1)
         end do
         
      end do         
      
      do i = 1, ncc - 2
      
         do j = 1, n12c(i)
         
            if(i12c(j,i).gt.spot2) then
               i12c(j,i) = i12c(j,i) - 1
            else if (i12c(j,i).eq.spot2) then
               i12c(j,i) = -1
            end if
            
         end do
         
      end do
    
      ncc = ncc - 2
      
c-----------------------------------------------------------------------
c     update position of the bonding and reference atoms in frag1
c-----------------------------------------------------------------------

      if(b1a.gt.spot1) b1a = b1a - 1
      if(b1a.gt.spot2) b1a = b1a - 1
      if(b2a.gt.spot1) b2a = b2a - 1
      if(b2a.gt.spot2) b2a = b2a - 1
      ibonderaca = b1a
      ibonderacb = b2a

      if(mm.eq.0) then
         if(t1a.gt.spot1) t1a = t1a - 1
         if(t1a.gt.spot2) t1a = t1a - 1
         if(t2a.gt.spot1) t2a = t2a - 1
         if(t2a.gt.spot2) t2a = t2a - 1
         torrefaca = t1a
         torrefacb = t2a 
      endif

c-----------------------------------------------------------------------
c     if we will be doing symmetrical links, then we must transfer the
c     information from hosta into moleculec, and then update the 
c     positions if needed
c-----------------------------------------------------------------------

      if(nsymattacha.gt.0) then 

         if(mm.eq.0) then

            do i = 1, nsymattacha
               symlooserc(i,jj)=symloosera(i,jj)
               symibonderc(i,jj)=symibondera(i,jj)
               symlooserc(i,jj+1)=symloosera(i,jj+1)
               symibonderc(i,jj+1)=symibondera(i,jj+1)
            enddo

         endif

c-----------------------------------------------------------------------
c     update position of symmetrical attachment points in frag1
c-----------------------------------------------------------------------

         do i=1, nsymattacha

            if(symlooserc(i,jj).gt.spot1) 
     &      symlooserc(i,jj) = symlooserc(i,jj) - 1
            if(symlooserc(i,jj).gt.spot2)
     &      symlooserc(i,jj) = symlooserc(i,jj) - 1
        
            if(symibonderc(i,jj).gt.spot1)
     &      symibonderc(i,jj) = symibonderc(i,jj) - 1
            if(symibonderc(i,jj).gt.spot2)
     &      symibonderc(i,jj) = symibonderc(i,jj) - 1

            if(symlooserc(i,jj+1).gt.spot1)
     &      symlooserc(i,jj+1) = symlooserc(i,jj+1) - 1
            if(symlooserc(i,jj+1).gt.spot2)
     &      symlooserc(i,jj+1) = symlooserc(i,jj+1) - 1

            if(symibonderc(i,jj+1).gt.spot1)
     &      symibonderc(i,jj+1) = symibonderc(i,jj+1) - 1
            if(symibonderc(i,jj+1).gt.spot2)
     &      symibonderc(i,jj+1) = symibonderc(i,jj+1) - 1

         enddo 

      endif
                        
c-----------------------------------------------------------------------
c     if it is the null link, make the bond between b1a and b2a
c     and return
c-----------------------------------------------------------------------

      if(ngg.eq.0) then

         do i = 1, n12c(b1a)
            if(i12c(i,b1a).eq.-1) then
               i12c(i,b1a) = b2a
               go to 30
            end if
         end do
   30    continue
         do i = 1, n12c(b2a)
            if(i12c(i,b2a).eq.-1) then
               i12c(i,b2a) = b1a
               go to 40
            end if
         end do
   40    continue

         nc = ncc

         return

      endif

c-----------------------------------------------------------------------
c     move link into molecule C, excluding the last two atoms which are
c     the dummy atoms that will be lost
c-----------------------------------------------------------------------

      do i = 1, ngg - 2
      
         xc(i + ncc) = xgu(i)
         yc(i + ncc) = ygu(i)
         zc(i + ncc) = zgu(i)
         radc(i + ncc) = radg(i)
         atnumc(i + ncc) = atnumg(i)
         atwtc(i + ncc) = atwtg(i)
         itypec(i + ncc) = itypeg(i)
         atlabc(i + ncc) = atlabg(i)
         iguestc(i + ncc) = 0
         
         listc(i + ncc) = listg(i)
         listc0(i + ncc) = listg0(i)
         
         n12c(i + ncc) = n12g(i)
         do j = 1, n12c(i + ncc)
            i12c(j,i + ncc) = i12g(j,i) + ncc
         end do
         
      end do
      
      nc = ncc + ngg - 2 
            
c-----------------------------------------------------------------------
c     update the position of the bonding and reference atoms in link
c-----------------------------------------------------------------------

      b1g = b1g + ncc     
      b2g = b2g + ncc
      ibondergca = b1g
      ibondergcb = b2g
     
      if(mm.eq.0) then
         if((t1g.eq.d2).and.(t2g.eq.d1)) then
            t1g = b2a
            t2g = b1a
            torrefgca = t1g
            torrefgcb = t2g
         else
            t1g = t1g + ncc
            t2g = t2g + ncc
            torrefgca = t1g
            torrefgcb = t2g
         endif
      endif

c-----------------------------------------------------------------------
c     correct the connectivity lists on bonding atoms
c-----------------------------------------------------------------------
      
      do i = 1, n12c(b1a)
      
         if(i12c(i,b1a).eq.-1) then
            i12c(i,b1a) = b1g
            go to 50
         end if
         
      end do
      
   50 continue
   
      do i = 1, n12c(b2a)
      
         if(i12c(i,b2a).eq.-1) then
            i12c(i,b2a) = b2g
            go to 60
         end if
         
      end do
     
   60 continue
      
      do i = 1, n12c(b1g)
      
         if(i12c(i,b1g).eq.d1+ncc) then
            i12c(i,b1g) = b1a
            go to 70
         end if
         
      end do
      
   70 continue
      
      do i = 1, n12c(b2g)
      
         if(i12c(i,b2g).eq.d2+ncc) then
            i12c(i,b2g) = b2a
            go to 80
         end if
         
      end do
     
   80 continue

      return      
      end
