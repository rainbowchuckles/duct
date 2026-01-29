# ifndef _IUPAC_
# define _IUPAC_

#     define          LNEL 9

      character*NEL   lcel
      dimension       lcel(LNEL)
      data            lcel(1) / 'O ' /,
     1                lcel(2) / 'H ' /,
     2                lcel(3) / 'C ' /,
     3                lcel(4) / 'N ' /,
     4                lcel(5) / 'AR' /,
     5                lcel(6) / 'CL' /,
     6                lcel(7) / 'NE' /,
     7                lcel(8) / 'HE' /,
     8                lcel(9) / 'E ' /

      real*8           lwel
      dimension       lwel(LNEL)
      data            lwel(1) /  0.01600d0 /,
     1                lwel(2) /  0.00101d0 /,
     2                lwel(3) /  0.01201d0 /,
     3                lwel(4) /  0.01401d0 /,
     4                lwel(5) /  0.03995d0 /,
     5                lwel(6) /  0.03545d0 /,
     6                lwel(7) /  0.02018d0 /,
     7                lwel(8) /  0.0040026d0 /,
     8                lwel(9) /  5.4858d-7 /

# endif
