      subroutine set_unused_bits(testbits)

      implicit none
      save

c-------------------------------------------------------------------
C!F77 
c
c!Description:
c     Set (as yet) unused bits. 
c     These bits are:
c          Spare bit (bit 23)
c          The cloud adjancency bit (bit 12)
c          The temporal consistency test (bit 24)
c          Spare bits (30-31)
c
c!Input parameters:
c testbits      6 cell byte array containing cloud mask bit results
c
c!Output Parameters:
c None.
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c-------------------------------------------------------------------

c ... arguments ..
      byte testbits(6)

c ... external subroutines
      external set_bit

c     Set temporal consistency test bit
      call set_bit(testbits,24)

c     Set cloud adjancency bit (a post-launch product)
      call set_bit(testbits,12) 

c     Set spare bits
      call set_bit(testbits,31)

      return
      end
