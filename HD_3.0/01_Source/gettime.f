c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                    Copyright © 2014 by Ben Hay
c                  Supramolecular Design Institute
c
c             This code is part of HostDesigner, Version 3
c                (see LICENSE.txt for conditions of use)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subroutine gettime - reports back the elapsed time in seconds
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine gettime(time)

      implicit none
      
      real time, times(2), etime 

      time = etime(times)

      return
      end
