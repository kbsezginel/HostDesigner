c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'equivrotor' -- based on similar geometric features, this routine
c     assigns a class and subclass for any rotor to be the same as one
c     of the hydrocarbon classes
c
c     passed variables:
c     integer class        class of the rotor
c     integer sclass       subclass of the rotor
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine equivrotor(class,sclass)

      implicit none

      integer class, sclass

      select case(class)
      case(4)
          select case(sclass)
          case(1,2)
              class  = 2
              sclass = 2
          case(3,4)
              class  = 2
              sclass = 4
          end select
      case(5)
          class  = 1
          sclass = 1
      case(6)
          select case(sclass)
          case(1,4)
              class  = 1
              sclass = 3
          case(2,3)
              class  = 1
              sclass = 1
          end select
      case(7)
          select case(sclass)
          case(1)
              class  = 2
              sclass = 2
              return
          case(2)
              class  = 2
              sclass = 4
          end select
      case(8)
          select case(sclass)
          case(1,2,3)
              class  = 1
              sclass = 2
          case(4)
              class  = 1
              sclass = 1
          end select
      case(9)
          class  = 1
          sclass = 3
      case(10)
          select case(sclass)
          case(1,3)
             class  = 2
             sclass = 1
          case(2,4)
             class  = 2
             sclass = 2
          end select   
      case(11)
          select case(sclass)
          case(1,3)
             class  = 2
             sclass = 1
          case(2)
             class  = 2
             sclass = 2
          end select   
      case(12)
          select case(sclass)
          case(1,3)
             class  = 2
             sclass = 2
          case(2,4)
             class  = 2
             sclass = 1
          end select   
      case(13)
          select case(sclass)
          case(1,2)
              class  = 3
              sclass = 1
          case(3)
              class  = 3
              sclass = 2
          end select
      case(14)
          select case(sclass)
          case(1)
              class  = 2
              sclass = 2
          case(2)
              class  = 2
              sclass = 2
          case(3)
              class  = 1
              sclass = 1
          case(4)
              class  = 1
              sclass = 1
          end select
      case(15)
          class  = 1
          sclass = 3
      case(16)
          class  = 1
          sclass = 3
      case default
          class  = 1
          sclass = 3
      end select 

      return
      end
