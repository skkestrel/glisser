      integer*4 function ibset (flag, bit)

      integer*4 flag, bit, new_flag
      logical btest
      external btest

      if (bit .ge. 32) then
         new_flag = flag
      else
         if (btest(flag,bit)) then
            new_flag = flag
         else
            new_flag = flag + 2**bit
         end if
      end if
      ibset = new_flag

      return
      end
