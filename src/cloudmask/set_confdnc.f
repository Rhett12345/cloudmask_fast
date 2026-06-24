      subroutine set_confdnc(confdnc,testbits)

      implicit none
      save

c--------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for setting output "bit" flags according
c     to final confidence of clear sky.
c
c!Input Parameters:
c confdnc       Unobstructed fov confidence
c
c!Output Parameters:
c testbits      6 cell byte array containing cloud mask bit results
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c--------------------------------------------------------------------

c ... scalar arguments
      real confdnc

c ... array arguments
      byte testbits(6)

c ... local scalars
      ! (none)

      if(confdnc .gt. 0.99) then
        call set_bit(testbits,1)
        call set_bit(testbits,2)
      else if(confdnc .gt. 0.95) then
        call set_bit(testbits,2)
      else if(confdnc .gt. 0.66) then
        call set_bit(testbits,1)
      end if

      return
      end
