#ifndef _TWINLAWS_DEF_H
#define _TWINLAWS_DEF_H

   void yyy_cell2tg(
      double cell[6], double &sc_tol,
      int &ng, int uu_g[][3][3], int u_g[][3],
      int &lc, int &nc, int &nc2, int uu_c[][3][3], double sc_c[],
      int &ivb, int &ierr
   );

#endif
