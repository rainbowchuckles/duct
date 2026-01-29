#  ifndef _US76ATM_
#  define _US76ATM_

#  include <cstdlib>
#  include <cstdio>
#  include <cmath>
#  include <cassert>

// include <typd.h>

/* US76 atmosphere up to 86km
   z                     geometric altitude km 
   d                     Density/Sea Level density
   t                     Temperature
   p                     Pressure/Sea Level pressure
*/

   void atm( double z,double &d,double &t,double &p );

#  endif
