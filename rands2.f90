!****************************************************************
! **********************************************************************
        double precision FUNCTION rands (SEED)
! **********************************************************************

!-----  this is a special function for random number generation
!        on 32-bit machines that do not support long integer
!        multiplication and truncation.  the technique used is to do
!        the multiplication and addition in parts, by splitting all
!       integers in a 'high' and a 'low' part.  the algorithm is
!-----        exact, and should give machine-independent results.

!-----        the algorithm implemented is (following d.e. knuth):
!        seed = seed*1592653589 + 453816691
!        if (seed.lt.0) seed = seed + 1 + 2147483647
!-----        note that 1592653589 = 48603*2**15 + 30485

! 32768 = 2^15, 65536 = 2^16, 4.65661287308E-10 = 2^(-31)

        INTEGER SEED, I1, I2

        I2 = MOD (SEED, 32768) * 30485
        I1 = MOD (SEED / 32768 * 30485, 65536) + MOD (MOD (SEED, 32768)   &
         * 48603, 65536) + I2 / 32768 + MOD (SEED / 32768, 2) * 32768 +  &
          13849
        I2 = MOD (I2, 32768) + 12659
        I1 = I1 + I2 / 32768
        I2 = MOD (I2, 32768)
        I1 = MOD (I1, 65536)
        SEED = I1 * 32768 + I2
        rands = REAL(I1 * 32768 + I2) * 4.65661287308E-10

        RETURN
        END

