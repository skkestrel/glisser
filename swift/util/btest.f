      logical function btest (flag, bit)

      integer*4 flag, bit, tmp, remainder

      tmp = flag/2**bit
      remainder = mod(tmp, 2)
      btest = remainder .eq. 1

      return
      end
