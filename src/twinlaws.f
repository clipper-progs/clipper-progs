
C
C     TWINLAWS
C     Copyright (C) 2004-2008 Andrey Lebedev
C
C     This code is distributed under the terms and conditions of the
C     CCP4 Program Suite Licence Agreement as a CCP4 Application.
C     A copy of the CCP4 licence can be obtained by writing to the
C     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
C
C
c     Code for calculating twinning operators from first principles
c     extracted from Sfcheck
c
c ====================================================================
      subroutine yyy_cell2tg(
     + cell, sc_tol, ng, uu_g, u_g, lc, nc, nc2, uu_c, sc_c, ivb,ierr)
      implicit   none
c
c     cell = cell
c     sc_tol = limit acceptable value of: (.05)
c    ( coordinate error on the surface of a sphere )/ ( radius of the sphere )
c          to estimate error, compared are:
c          1) operator that fits translations, but is not orthogonal
c          2) its best fit to an orthogonal operator
c     ng = No of symops
c     uu_g = rotation matrixes ( h(j)* uu_g(j,i)* x(i) !!! )
C     x'(j) = uu_g(j,1)* x(1) + uu_g(j,2)* x(2) + uu_g(j,3)* x(3)
c     u_g = symops' translations (*12)
c     lc = dimension of array for storage of c ( = 48)
c     nc = total number of non-equivalent twin operators
c     nc - 1       
c     uu_c = matrices of the above operators ( h(j)* uu_c(j,i)* x(i) !!! ) (*12)
C     trace =-1 --> 2; =0 --> 3; =1 --> 4 ; =2 --> 6 ; =3 --> [1] 
c     sc_c = scores for operators e ( in the above sense of errors on the
c     surface of the sphere), i.e. sc_c < sc_tol except those operators
c     from the factor group, that are generated from the low score
c     generators to close the factor group.
c
c     ivb = 0 : no output,
c     ivb = 1 : a table of traces of cell-group operators (vs SG-ops and classses)
c     ivb = 2 : same and operators of cell-group etc.
c     ivb = 3 : test output
c

      real*8     cell(6), muu(3,3), mvv(3,3)
      real*8     celt(6)
      integer*4  uSv(3,3), vSu(3,3)

      real*8     sc_tol
      real*8     sc_eps 
      save       sc_eps

      integer*4  ng
      integer*4  u_g(3,ng)
      integer*4  uu_g(3,3,ng)

      integer*4  nh
      integer*4  vv_h(3,3,24)
      integer*4  uu_h(3,3,24)
      integer*4  h_hh(24,24)
      integer*4  h_h(24)
      integer*4  tr_h(24)
      integer*4  p_h(24)
      real*8     sc_h(24)

      integer*4  np
      integer*4  lc, nc, nc2, ic
      integer*4  h_c(24)
      integer*4  h_pc(24,24)
      integer*4  uu_c(3,3,lc)
      real*8     sc_c(lc)

      integer*4  i, j, ivb, ierr
      logical*4  ok

c --------------------------------------------------------------------
      sc_eps = 1.0D-04
      nc     = 0
      nc2    = 0
      
c      write(88,*) 'ng = ', ng
c      write(88,*) ' '
      do i=1,ng
c        write(88,'(3(i6,2x))') u_g(1,i), u_g(2,i), u_g(3,i)
      enddo
c      write(88,*) ' '
      do i=1,ng
         do j=1,3
c            write(88,'(3(i6,2x))') uu_g(j,1,i), uu_g(j,2,i), uu_g(j,3,i)
         enddo
c         write(88,*) ' '
      enddo

c      write(88,*) cell(1), cell(2), cell(3), cell(4), cell(5), cell(6)
      call yyy_cell2met( cell, muu, ierr )
      if(ierr.ne.0) return
      !write(88,*) 'pt 1'
      call yyy_met2cell( celt, muu, ierr )
      !write(88,*) 'pt 2'
      if(ierr.ne.0) return
      !write(88,*) 'pt 3'
      call yyy_find_base( ng, uu_g, u_g, uSv, ierr )
      !write(88,*) 'pt 4'
      if(ierr.ne.0) return
      call yyy_short_base( muu, uSv, mvv, vSu, ierr )
      !write(88,*) 'pt 5'
      if(ierr.ne.0) return

      call yyy_cell_group( ivb, sc_tol, mvv, nh, vv_h, sc_h, ierr )
      !write(88,*) 'pt 6'
      if(ierr.ne.0) return
      call yyy_test_group( nh, vv_h, h_hh, h_h, tr_h, ok )

      if( .not. ok ) then
c        stop'yyy_cell2tg:a'
        ierr = 1
        return
      endif
      !write(88,*) 'pt 7'

      if( ivb .ge. 2 )then
         call yyy_write_group( 'vv_h', nh, vv_h, h_hh, h_h, tr_h )
      endif
            !write(88,*) 'pt 8'

      call yyy_transform_u_to_v( vSu, uSv, nh, vv_h, uu_h )
      call yyy_test_group12( nh, uu_h, h_hh, h_h, tr_h, ok )
                  !write(88,*) 'pt 9'

      if( .not. ok ) then
c        stop'yyy_cell2tg:b'
        ierr = 1
        return
      endif
                  !write(88,*) 'pt 10 ',ok


      if( ivb .ge. 2 )then
         call yyy_write_group( 'uu_h', nh, uu_h, h_hh, h_h, tr_h )
      endif
                  !write(88,*) 'pt 11'

      call yyy_map_subgroup( ng, uu_g, nh, uu_h, p_h, np, h_pc, ok )
                  !write(88,*) 'pt 12 ',ok

      if( .not. ok )                                            return

      call yyy_right_classes(
     +   nh, h_hh, p_h, tr_h, sc_h, sc_eps,
     +   nc, h_c,
     +   np, h_pc, nc2, ok , ierr)
      !write(88,*) 'pt 13 ',ierr
      if(ierr.ne.0) return
      !write(88,*) 'pt W'

      if( .not. ok )                                            return
      if( ivb .ge. 1 )then
         call yyy_write_classes( np, nc, h_c, h_pc, tr_h, sc_h, ivb )
      endif

      if( nc .gt. lc ) then
c        stop'yyy_cell2tg:c'
        ierr = 1
        return
      endif
      !write(88,*) 'pt X nc = ',nc

      do ic = 1,nc
         !write(88,*) 'ic', ic
         sc_c(ic) = sc_h(h_c(ic))
         do i = 1,3
            do j = 1,3
               !write(88,*) i,j
               uu_c(j,i,ic) = uu_h(j,i,h_c(ic))
            enddo
         enddo
      enddo
c     write(88,*) ' '
c     write(88,*) 'lc = ',lc,' nc = ',nc,'sc_tol = ',sc_tol
      do i=1,nc
c        write(88,*) 'score = ', sc_c(i)
         do j=1,3
c           write(88,'(3(i6,2x))') uu_c(j,1,i), uu_c(j,2,i), uu_c(j,3,i)
         enddo
c        write(88,*) ' '
      enddo
      !write(88,*), 'exiting twinlaws, nc = ',nc
      end

c ====================================================================
      subroutine yyy_write_classes(
     +   np, nc, h_c, h_pc, tr_h, sc_h, ivb )
      implicit   none

      integer*4  np, ip
      integer*4  nc, ic
      integer*4  h_c(24)
      integer*4  h_pc(24,24)

      integer*4  tr_h(24)
      real*8     sc_h(24)

      integer*4  ivb

c --------------------------------------------------------------------
      if( ivb .ge. 2 )then

c         write(*,'(78(''-''))')
c         write(*,'(a22,a44)')
c     +      'CG-op vs PG-op\\class,',
c     +      'and selected CG-op and score*1000 vs class:'
c        write(*,'(a)')
c         write(*,'((6x,24i3))') ( ic, ic = 1,nc )
c         write(*,'(a)')
         do ip = 1,np
c            write(*,'((i3,3x,24i3))') ip, ( h_pc(ip,ic), ic = 1,nc )
         enddo
c         write(*,'(a)')
c         write(*,'((a3,3x,24i3))') 'sel', ( h_c(ic), ic = 1,nc )
c         write(*,'((a3,3x,24i3))') 'sco',
c     +      ( int( 0.5D+00 + 1000* sc_h(h_c(ic)) ), ic = 1,nc )

      endif
      if( ivb .ge. 1 )then

c         write(*,'(a)')
c         write(*,'(78(''-''))')
c         write(*,'(a26,a48)')
c     +      'tr(CG-op) vs PG-op\\class,',
c     +      'and tr(selected CG-op) and score*1000 vs class:'
c         write(*,'(a)')
c         write(*,'((6x,24i3))') ( ic, ic = 1,nc )
c         write(*,'(a)')
         do ip = 1,np
c            write(*,'((i3,3x,24i3))')
c     +         ip, ( tr_h(h_pc(ip,ic)), ic = 1,nc )
         enddo
c         write(*,'(a)')
c         write(*,'((a3,3x,24i3))') 'sel', ( tr_h(h_c(ic)), ic = 1,nc )
c         write(*,'((a3,3x,24i3))') 'sco',
c     +      ( int( 0.5D+00 + 1000* sc_h(h_c(ic)) ), ic = 1,nc )
c         write(*,'(a)')

      endif
      if( ivb .ge. 2 )then

c         write(*,'(78(''-''))')

      endif
      end

c ====================================================================
      subroutine yyy_right_classes(
     +   nh, h_hh, p_h, tr_h, sc_h, sc_eps,
     +   nc, h_c,
     +   np, h_pc, nc2, ok, ierr )
      implicit   none

c input:   nh, h_hh, p_h, tr_h
c          np, h_pc(1:np,1)

c output:  nc, h_c
c          h_pc(1:np,2:nc), nc2, ok

      integer*4  nh, ih, jh
      integer*4  h_hh(24,24)
      integer*4  p_h(24)
      integer*4  tr_h(24)
      real*8     sc_h(24)

      integer*4  np, ip
      integer*4  nc, ic, nc2
      integer*4  h_c(24)
      integer*4  h_pc(24,24)

      real*8     sc_min, sc_max, sc_pc, sc_eps
      integer*4  tr_min, tr_max, tr_pc
      integer*4  ip_min, ip_max, ierr
      logical*4  ok

c --------------------------------------------------------------------
      ok = .false.

      do ip = 1,np
         ih = h_pc(ip,1)
         p_h(ih) = ip
      enddo
      ic = 1
      do ih = 1,nh
         if( p_h(ih) .eq. 0 )then
            ic = ic + 1
            do ip = 1,np
               jh = h_hh(h_pc(ip,1),ih)
               if( p_h(jh) .ne. 0 )                             return
               p_h(jh) = ip
               h_pc(ip,ic) = jh
            enddo
         endif
      enddo
      nc = ic
      if( np* nc .ne. nh )                                      return
      if( h_pc(1,1) .ne. 1 )                                    return
      h_c(1) = 1
      nc2 = 0
      do ic = 2,nc
         tr_min = + 4
         tr_max = - 2
         ip_min = 0
         ip_max = 0
         do ip = 1,np
            tr_pc = tr_h(h_pc(ip,ic))
            if( tr_max .le. tr_pc )then
               tr_max = tr_pc
               ip_max = ip
            endif
            if( tr_min .ge. tr_pc )then
               tr_min = tr_pc
               ip_min = ip
            endif
         enddo
         if( tr_min .lt. -1 .or. tr_max .gt. 3 )                return
         h_c(ic) = h_pc(ip_min,ic)
         if( tr_min .eq. -1 ) nc2 = nc2 + 1
      enddo
      do ic = 1,nc
         sc_max = 0
         do ip = 1,np
            sc_pc = sc_h(h_pc(ip,ic))
            sc_max = max( sc_max, sc_pc )
         enddo
         sc_min = sc_max
         do ip = 1,np
            sc_pc = sc_h(h_pc(ip,ic))
            sc_min = min( sc_min, sc_pc )
         enddo

         if( sc_max - sc_min .gt. sc_eps ) then
c           stop'yyy_right_classes:a'
           ierr = 1
           return
         endif

      enddo

      ok = .true.
      end

c ====================================================================
      subroutine yyy_map_subgroup(
     +   ng, uu_g,
     +   nh, uu_h, p_h,
     +   np, h_p,
     +   ok )
      implicit   none

      integer*4  ng, ig
      integer*4  uu_g(3,3,ng)

      integer*4  nh, ih
      integer*4  uu_h(3,3,24)
      integer*4  p_h(24)

      integer*4  np, ip
      integer*4  h_p(24)

      integer*4  i, j
      integer*4  uu_x(3,3), uu_y(3,3), d_x
      logical*4  ok

c --------------------------------------------------------------------
      do ih = 1,nh
         p_h(ih) = 0
      enddo
      ip = 0
      do ig = 1,ng
         do i = 1,3
            do j = 1,3
               uu_x(j,i) = uu_g(j,i,ig)* 12
            enddo
         enddo
         call yyy_invert_int( uu_x, uu_y, d_x )
         if( d_x .lt. 0 )then
            do i = 1,3
               do j = 1,3
                  uu_x(j,i) = - uu_x(j,i)
               enddo
            enddo
         endif
         ih = 0
         ok = .false.
         do while( .not. ok .and. ih .lt. nh )
            ih = ih + 1
            call yyy_the_same_ops( uu_x, uu_h(1,1,ih), ok )
         enddo
         if( .not. ok )                                         return
         if( p_h(ih) .eq. 0 )then
            ip = ip + 1
            h_p(ip) = ih
            p_h(ih) = ip
         endif
      enddo
      np = ip
      end

c ====================================================================
      subroutine yyy_test_group12(
     +   nh, uu_h, h_hh, h_h, tr_h, ok )
      implicit   none

c Collects tables of multiplication and reciprocal elements
c thus checking integrity of a group vv_h(,,)

      integer*4  nh, ih, jh, kh, i, j
      integer*4  uu_h(3,3,24)
      integer*4  uu_x(3,3)
      integer*4  uu_y(3,3)
      integer*4  h_hh(24,24)
      integer*4  h_h(24)
      integer*4  tr_h(24)
      integer*4  tr_x
      logical*4  ok

c --------------------------------------------------------------------
      do ih = 1,nh
         do jh = 1,nh
            call yyy_multiply_ops( uu_h(1,1,jh), uu_h(1,1,ih), uu_x )
            kh = h_hh(jh,ih)
            do i = 1,3
               do j = 1,3
                  uu_y(j,i) = uu_h(j,i,kh)* 12
               enddo
            enddo
            call yyy_the_same_ops( uu_x, uu_y, ok )
            if( .not. ok )                                      return
         enddo
      enddo
      ok = .false.
      do ih = 1,nh
         tr_x = 0
         do i = 1,3
            tr_x = tr_x + uu_h(i,i,ih)
         enddo
         if( tr_x .ne. tr_h(ih)* 12 )                           return
      enddo
      ok = .true.
      end

c ====================================================================
      subroutine yyy_find_base( ng, uu_g, u_g, u_b, ierr )
      implicit   none

      integer*4  ng, ig
      integer*4  uu_g(3,3,ng)
      integer*4  u_g(3,ng)

      integer*4  u(3), v(3), u_b(3,3), r(3,3), d, tr
      integer*4  i, j, ib, ierr
      logical*4  found

c --------------------------------------------------------------------
      do i = 1,3
         do j = 1,3
            u_b(j,i) = 0
            r(j,i) = 0
         enddo
         u_b(i,i) = 12
         r(i,i) = 1
      enddo
      d = 12
      ig = 0
      found = .true.
      do while( found .and. ig .lt. ng )
         ig = ig + 1
         tr = 0
         do i = 1,3
            tr = tr + uu_g(i,i,ig)
         enddo
         if( tr .eq. 3 )then
            do i = 1,3
               v(i) = 0
               do j = 1,3
                  v(i) = v(i) + r(i,j)* u_g(j,ig)
               enddo
            enddo
            ib = 0
            call yyy_divide_int( 3, v, d, found )
            do while( .not. found .and. ib .lt. 3 )
               ib = ib + 1
               do i = 1,3
                  u(i) = u_b(i,ib)
                  u_b(i,ib) = u_g(i,ig)
               enddo
               call yyy_invert_int( u_b, r, d )
               do i = 1,3
                  v(i) = 0
                  do j = 1,3
                     v(i) = v(i) + r(i,j)* u(j)
                  enddo
               enddo
               call yyy_divide_int( 3, v, d, found )
               if( .not. found )then
                  do i = 1,3
                     u_b(i,ib) = u(i)
                  enddo
               endif
            enddo
         endif
      enddo
      if( .not. found ) then
c        stop'find_base:a'
        ierr = 1
        return
      endif

      end

c ====================================================================
      subroutine yyy_short_base( muu, uSv, mvv, vSu, ierr )
      implicit   none

      real*8     muu(3,3), mvv(3,3), sc_t(13), sc_min, sc_max
      logical*4  found, ok
      integer*4  v_b(3,3), uSv(3,3), vSu(3,3), b_t(13), detb, xx
      integer*4  i, j, icy, it, jt, ib, ierr
      integer*4  v_t(3,26) 
      data  v_t /
     +                       1, 0, 0,
     +                       0, 1, 0,
     +                       0, 0, 1,
     +                       0, 1, 1,
     +                       0, 1,-1,
     +                       1, 0, 1,
     +                      -1, 0, 1,
     +                       1, 1, 0,
     +                       1,-1, 0,
     +                       1, 1, 1,
     +                      -1, 1, 1,
     +                       1,-1, 1,
     +                       1, 1,-1,
     +                         39* 0 /

c --------------------------------------------------------------------
      do i = 1,3
         do j = 1,3
            mvv(j,i) = muu(j,i)/ 144
         enddo
      enddo
      do it = 14,26
         do i = 1,3
            v_t(i,it) = - v_t(i,it-13)
         enddo
      enddo
      call yyy_transform_m( uSv, mvv )
      icy = 0
      found = .false.
      do while( .not. found )
         icy = icy + 1

         if( icy .ge. 1000 ) then
c          stop'short_base:a'
           ierr = 1
           return
         endif

         sc_max = 0
         do it = 1,13
            b_t(it) = 0
            sc_t(it) = 0
            do i = 1,3
               do j = 1,3
                  sc_t(it) = sc_t(it) + v_t(j,it)* mvv(j,i)* v_t(i,it)
               enddo
            enddo
            sc_max = max( sc_max, sc_t(it) )
         enddo

         found = .true.
         do ib = 1,3
            sc_min = sc_max* 2 + 1
            jt = 0
            do it = 1,13
               if( b_t(it) .eq. 0 .and. sc_t(it) .lt. sc_min )then
                  do i = 1,3
                     v_b(i,ib) = v_t(i,it)
                  enddo
                  if( ib .eq. 3 )then
                     call yyy_invert_int( v_b, vSu, detb )
                     call yyy_divide_int( 9, vSu, detb, ok )
                  else
                     ok = .true.
                  endif
                  if( ok )then
                     sc_min = sc_t(it)
                     jt = it
                  endif
               endif
            enddo

            if( jt .le. 0 .or. jt .gt. 13 ) then
c             stop'short_base:b'
              ierr = 1
              return
            endif

            b_t(jt) = ib
            do i = 1,3
               v_b(i,ib) = v_t(i,jt)
            enddo
            found = found .and. jt .le. 3
         enddo
         call yyy_transform_b( v_b, uSv )
         call yyy_transform_m( v_b, mvv )
      enddo
      call yyy_invert_int( uSv, vSu, xx )
      if( .not. found .or. xx .eq. 0 ) then 
c       stop'short_base:c'
        ierr = 1
        return
      endif
      detb = xx/ 12
      if( xx .ne. detb* 12 ) then
c       stop'short_base:d'
        ierr = 1
        return
      endif
      do i = 1,3
         do j = 1,3
            xx = vSu(j,i)/ detb
            if( vSu(j,i) .ne. xx* detb ) then
c             stop'short_base:e'
              ierr = 1
              return
            endif
            vSu(j,i) = xx
         enddo
      enddo
      end

c ====================================================================
      subroutine yyy_score_ops( mvv, no, vv_o, sc_o, ierr )
      implicit   none

      real*8     mvv(3,3), a(3,3), r(3), b(3,3)
      real*8     of(3,3), fo(3,3)

      integer*4  no, io, ierr
      real*8     sc_o(no)
      integer*4  vv_o(3,3,no)
      real*8     qo(3,3), oqo(3,3)
      integer*4  i, j, k

c --------------------------------------------------------------------
      call yyy_eigv3( mvv, a, r, b, ierr )
      if(ierr.ne.0) return
      do i = 1,3
         r(i) = sqrt( r(i) )
      enddo
      do i = 1,3
         do j = 1,3
            of(j,i) = 0
            do k = 1,3
               of(j,i) = of(j,i) + a(j,k)* r(k)* b(k,i)
            enddo
         enddo
      enddo
      do i = 1,3
         r(i) = 1/ r(i)
      enddo
      do i = 1,3
         do j = 1,3
            fo(j,i) = 0
            do k = 1,3
               fo(j,i) = fo(j,i) + a(j,k)* r(k)* b(k,i)
            enddo
         enddo
      enddo
      do io = 1,no
         do i = 1,3
            do j = 1,3
               qo(j,i) = 0
               do k = 1,3
                  qo(j,i) = qo(j,i) + vv_o(j,k,io)* fo(k,i)
               enddo
            enddo
         enddo
         do i = 1,3
            do j = 1,3
               oqo(j,i) = 0
               do k = 1,3
                  oqo(j,i) = oqo(j,i) + of(j,k)* qo(k,i)
               enddo
            enddo
         enddo
         call yyy_eigv3( oqo, a, r, b, ierr )
         if(ierr.ne.0) return
         sc_o(io) = 0
         do i = 1,3
            sc_o(io) = sc_o(io) + ( 1 - r(i) )**2
         enddo
         sc_o(io) = sqrt( sc_o(io)/ 3 )
      enddo
      end

c ====================================================================
      subroutine yyy_transform_u_to_v( uSv, vSu, np, uu_p, vv_p )
      implicit   none

      integer*4  uSv(3,3), vSu(3,3)
      integer*4  np, ip
      integer*4  uu_p(3,3,24)
      integer*4  vv_p(3,3,24)
      integer*4  vu(3,3)

c --------------------------------------------------------------------
      do ip = 1,np
         call yyy_multiply_ops( vSu, uu_p(1,1,ip), vu )
         call yyy_multiply_ops( vu, uSv, vv_p(1,1,ip) )
      enddo
      end

c ====================================================================
      subroutine yyy_divide_int( n, v, d, ok )
      implicit   none

      integer*4  n, v(n), d, i, j
      logical*4  ok

c --------------------------------------------------------------------
      ok = .false.
      if( d .eq. 0 )                                            return
      do i = 1,n
         j = v(i)/ d
         if( v(i) .ne. d* j )                                   return
         v(i) = j
      enddo
      ok = .true.
      end

c ====================================================================
      subroutine yyy_transform_b( db, b )
      implicit   none

      integer*4  b(3,3), db(3,3), bdb(3,3)
      integer*4  i, j, k

c --------------------------------------------------------------------
      do i = 1,3
         do j = 1,3
            bdb(j,i) = 0
            do k = 1,3
               bdb(j,i) = bdb(j,i) + b(j,k)* db(k,i)
            enddo
         enddo
      enddo
      do i = 1,3
         do j = 1,3
            b(j,i) = bdb(j,i)
         enddo
      enddo
      end

c ====================================================================
      subroutine yyy_transform_m( db, m )
      implicit   none

      real*8     m(3,3), mdb(3,3)
      integer*4  db(3,3)
      integer*4  i, j, k

c --------------------------------------------------------------------
      do i = 1,3
         do j = 1,3
            mdb(j,i) = 0
            do k = 1,3
               mdb(j,i) = mdb(j,i) + m(j,k)* db(k,i)
            enddo
         enddo
      enddo
      do i = 1,3
         do j = 1,3
            m(j,i) = 0
            do k = 1,3
               m(j,i) = m(j,i) + db(k,j)* mdb(k,i)
            enddo
         enddo
      enddo
      end

c ====================================================================
      subroutine yyy_cell2met( cell, met, ierr )
      implicit   none

      real*8     cell(6), met(3,3), d2r
      logical*4  rad, deg
      integer*4  i, j, i1, i2, ierr
      integer*4  m3(6)
      data  m3 / 1, 2, 3, 1, 2, 3 /

c --------------------------------------------------------------------
      rad = .true.
      deg = .true.
      do i = 4,6
         rad = rad .and. cell(i) .lt. 3.2D+00 .and. cell(i) .gt. 0
         deg = deg .and. cell(i) .gt. 3.2D+00 .and. cell(i) .lt. 180
      enddo
      if( rad )then
         d2r = 1.0D+00
      else if( deg )then
         d2r = atan( 1.0D+00 )/ 45
      else
c        stop'yyy_cell2met:a'
        ierr = 1
        return
      endif
      do i = 1,3
         i1 = m3(i+1)
         i2 = m3(i+2)
         met(i,i) = 1
         met(i1,i2) = cos( cell(i+3)* d2r )
         met(i2,i1) = met(i1,i2)
      enddo
      do i = 1,3
         do j = 1,3
            met(j,i) = cell(j)* met(j,i)* cell(i)
         enddo
      enddo
      end

c ====================================================================
      subroutine yyy_met2cell( cell, met, ierr )
      implicit   none

      real*8     test
      real*8     cell(6), met(3,3), r2d
      integer*4  i, i1, i2, tol , ierr
      integer*4  m3(6) 

      data  m3/ 1, 2, 3, 1, 2, 3 /
      data  tol / 10000 /
c --------------------------------------------------------------------
      do i = 1,3
         if( met(i,i) .le. 0 ) then
c           stop'yyy_met2cell:a'
           ierr = 1
           return
         endif
         cell(i) = sqrt( met(i,i) )
      enddo
      do i = 1,3
         i1 = m3(i+1)
         i2 = m3(i+2)
         cell(i+3) = met(i1,i2)/ ( cell(i1)* cell(i2) )
         test = cell(i+3) - met(i2,i1)/ ( cell(i1)* cell(i2) ) + tol
         if( test .ne. tol* 1.0D+00 ) then
c          stop'yyy_met2cell:b'
           ierr = 1
           return
         endif
      enddo
      r2d = 45/ atan( 1.0D+00 )
      do i = 4,6
         cell(i) = r2d* acos( cell(i) )
      enddo
      end

c ====================================================================
      subroutine yyy_eigv3( w, a, r, b, ierr )
      implicit   none

      real*8     w(3,3)
      real*8     m(3,3), r(3), a(3,3), b(3,3), ab(2,2), abc, maxab
      integer*4  ia, ja, ierr
      integer*4  ib, jb

      integer*4  c3(5) 
      integer*4  i1, i2, i3, j1, icy
c     integer*4  k1
      real*8     dp, op, tp, cp, sp
      real*8     dm, om, tm, cm, sm
      real*8     cosa, sina, testa
      real*8     cosb, sinb, testb
      real*8     cosx, sinx, tanx, m2, m3
      real*8     tmax, tplus
      logical*4  more(3)
      data c3 / 1, 2, 3, 1, 2 /
c --------------------------------------------------------------------
      do i1 = 1,3
         do j1 = 1,3
            m(j1,i1) = w(j1,i1)
            a(j1,i1) = 0
            b(j1,i1) = 0
         enddo
         a(i1,i1) = 1
         b(i1,i1) = 1
      enddo
      i1 = 1
      icy = 0
      more(1) = .true.
      more(2) = .true.
      more(3) = .true.
      do while( more(1) .or. more(2) .or. more(3) )
         icy = icy + 1
         if( icy .gt. 40 ) then
c          stop'yyy_eigv3:a'
           ierr = 1
           return
         endif
         i2 = c3(i1+1)
         i3 = c3(i1+2)
         dp = m(i2,i2) + m(i3,i3)
         dm = m(i2,i2) - m(i3,i3)
         op = m(i3,i2) + m(i2,i3)
         om = m(i3,i2) - m(i2,i3)
         tp = sqrt( om**2 + dp**2 )
         tm = sqrt( op**2 + dm**2 )

c handling zeroes:

         tmax = max( tp, tm )

         tplus = tmax + dp
         if( tplus .eq. tmax ) dp = 0

         tplus = tmax + dm
         if( tplus .eq. tmax ) dm = 0

         tplus = tmax + op
         if( tplus .eq. tmax ) op = 0

         tplus = tmax + om
         if( tplus .eq. tmax ) om = 0

         tplus = tmax + tp
         if( tplus .eq. tmax ) tp = 0

         tplus = tmax + tm
         if( tplus .eq. tmax ) tm = 0

c ordering:

         if( i1 .ne. 2 ) tm = - tm

c sines, cosines:

         if( tm .eq. 0 )then
            cp = 1
            sp = 0
         else
            cp = dm/ tm
            sp = op/ tm
         endif
         if( tp .eq. 0 )then
            cm = 1
            sm = 0
         else
            cm = dp/ tp
            sm = om/ tp
         endif
         ab(1,1) = ( cm + cp )/ 2.0D+00
         ab(2,2) = ( cm - cp )/ 2.0D+00
         ab(2,1) = ( sp + sm )/ 2.0D+00
         ab(1,2) = ( sp - sm )/ 2.0D+00
         m(i2,i2) = ( tp + tm )/ 2
         m(i3,i3) = ( tp - tm )/ 2
         m(i3,i2) = 0
         m(i2,i3) = 0

         maxab = 0
         do ib = 1,2
            do ia = 1,2
               abc = abs( ab(ia,ib) )
               if( maxab .lt. abc )then
                  ja = ia
                  jb = ib
                  maxab = abc
               endif
            enddo
         enddo

         tanx = ab(3-ja,jb)/ ab(ja,jb)
         cosx = 1/ sqrt( 1 + tanx**2 )
         sinx = tanx* cosx
         if( ab(ja,jb) .lt. 0 )then
            cosx = - cosx
            sinx = - sinx
         endif

         testa = 1 + tanx
         if( ja .eq. 1 )then
            cosa = cosx
            sina = sinx
         else
            cosa = sinx
            sina = cosx
         endif

         tanx = ab(ja,3-jb)/ ab(ja,jb)
         cosx = 1/ sqrt( 1 + tanx**2 )
         sinx = tanx* cosx

         testb = 1 + tanx
         if( jb .eq. 1 )then
            cosb = cosx
            sinb = sinx
         else
            cosb = sinx
            sinb = cosx
         endif

         m2 = m(i2,i1)
         m3 = m(i3,i1)
         m(i2,i1) = cosa* m2 + sina* m3
         m(i3,i1) = cosa* m3 - sina* m2
         m2 = m(i1,i2)
         m3 = m(i1,i3)
         m(i1,i2) = cosb* m2 + sinb* m3
         m(i1,i3) = cosb* m3 - sinb* m2

         do j1 = 1,3
            m2 = a(j1,i2)
            m3 = a(j1,i3)
            a(j1,i2) = cosa* m2 + sina* m3
            a(j1,i3) = cosa* m3 - sina* m2
            m2 = b(i2,j1)
            m3 = b(i3,j1)
            b(i2,j1) = cosb* m2 + sinb* m3
            b(i3,j1) = cosb* m3 - sinb* m2
         enddo
         more(i1) = testa .ne. 1 .or. testb .ne. 1
         i1 = c3(i1+1)
      enddo
      do i1 = 1,3
         r(i1) = m(i1,i1)
      enddo
c     do i1 = 1,3
c        do j1 = 1,3
c           m(j1,i1) = 0
c           do k1 = 1,3
c              m(j1,i1) = m(j1,i1) + a(j1,k1)* r(k1)* b(k1,i1)
c           enddo
c        enddo
c     enddo
      end

c ====================================================================
cA:
c ====================================================================
      subroutine yyy_write_group( ch, nh, vv_h, h_hh, h_h, tr_h )
      implicit   none

      character  ch*(*)
      integer*4  nh, ih, jh
      integer*4  vv_h(3,3,24)
      integer*4  h_hh(24,24)
      integer*4  h_h(24)
      integer*4  tr_h(24)
      integer*4  i, j

c --------------------------------------------------------------------
c      write(*,'(a)') '--------------------------------------'
c      write(*,'(a)') ch
c      write(*,'(a)')
c      write(*,'(2a6)') 'no', 'tr'
      do ih = 1,nh
c         write(*,'(a)')
c         write(*,'((2i6,3i4))')
c     +      ih, tr_h(ih), ( vv_h(1,i,ih), i = 1,3 )
c         write(*,'((12x,3i4))')
c     +      ( ( vv_h(j,i,ih), i = 1,3 ), j = 2,3 )
      enddo
c      write(*,'(a)')
c      write(*,'((6x,24i3))') ( ih, ih = 1,nh )
c      write(*,'(a)')
      do jh = 1,nh
c         write(*,'((i3,3x,24i3))') jh, ( h_hh(jh,ih), ih = 1,nh )
      enddo
c      write(*,'(a)')
c      write(*,'((6x,24i3))') ( ih, ih = 1,nh )
c      write(*,'(a)')
c      write(*,'((6x,24i3))') ( h_h(ih), ih = 1,nh )
c      write(*,'(a)')
      end

c ====================================================================
      subroutine yyy_test_group(
     +   nh, vv_h, h_hh, h_h, tr_h, found )
      implicit   none

c Collects tables of multiplication and reciprocal elements
c thus checking integrity of a group vv_h(,,)

      integer*4  nh, ih, jh, kh, i
      integer*4  vv_h(3,3,24)
      integer*4  vv_x(3,3)
      integer*4  h_hh(24,24)
      integer*4  h_h(24)
      integer*4  tr_h(24)
      logical*4  found

c --------------------------------------------------------------------
      do ih = 1,nh
         do jh = 1,nh
            call yyy_multiply_ops( vv_h(1,1,jh), vv_h(1,1,ih), vv_x )
            kh = 0
            found = .false.
            do while( .not. found .and. kh .lt. nh )
               kh = kh + 1
               call yyy_the_same_ops( vv_x, vv_h(1,1,kh), found )
            enddo
            if( .not. found )                                   return
            h_hh(jh,ih) = kh
         enddo
      enddo
      found = .false.
      if( h_hh(1,1) .ne. 1 )                                    return
      do ih = 1,nh
         jh = 0
         do while( .not. found .and. jh .lt. nh )
            jh = jh + 1
            found = h_hh(jh,ih) .eq. 1
         enddo
         if( .not. found )                                      return
         h_h(ih) = jh
         found = .false.
         do while( jh .lt. nh )
            jh = jh + 1
            if( h_hh(jh,ih) .eq. 1 )                            return
         enddo
      enddo
      do ih = 1,nh
         if( h_h(h_h(ih)) .ne. ih )                             return
      enddo
      do ih = 1,nh
         tr_h(ih) = 0
         do i = 1,3
            tr_h(ih) = tr_h(ih) + vv_h(i,i,ih)
         enddo
      enddo
      found = .true.
      end

c ====================================================================
      subroutine yyy_cell_group(
     +   ivb, sc_tol, mvv, nh, vv_h, sc_h, ierr )
      implicit   none

c Produces a group with generators vv_o(io) such as sc_o(io) < sc_tol.
c An operator vv_o(jo) with sc_o(jo) < sc_tol is not present in the group
c if it is not consistent with an operator vv_o(io) and vv_o(io) is present
c in the group and sc_o(io) < sc_o(jo).

      real*8     sc_tol, mvv(3,3)

      integer*4  mo, no, ierr
      parameter  ( mo = 504 )
      integer*4  vv_o(3,3,mo)
      real*8     sc_o(mo)
      integer*4  o_s(mo*2)
      integer*4  h_s(mo)

      integer*4  nh
      integer*4  s_h(mo)
      integer*4  vv_h(3,3,24)
      real*8     sc_h(24)

      integer*4  ivb

c --------------------------------------------------------------------
      call yyy_generate_ops( mo, no, vv_o, ierr )
      if(ierr.ne.0) return

      call yyy_score_ops( mvv, no, vv_o, sc_o, ierr )
      if(ierr.ne.0) return

      sc_o(1) = - 1
      call yyy_sort_key( no, o_s, sc_o )
      sc_o(1) = 0

      call yyy_collect_group(
     +   ivb, sc_tol,
     +   no, h_s, o_s, vv_o, sc_o,
     +   nh, s_h, vv_h, sc_h, ierr )
      if(ierr.ne.0) return

      end

c ====================================================================
      subroutine yyy_collect_group(
     +   ivb, sc_tol,
     +   no, h_s, o_s, vv_o, sc_o,
     +   nh, s_h, vv_h, sc_h, ierr )
      implicit   none

c Multiplication table is not collected as elements of group
c are collected in order different from sorting order required

      real*8     sc_tol

      integer*4  no, ns, is, ts, ierr
      integer*4  o_s(no)
      integer*4  h_s(no)

      integer*4  vv_o(3,3,no)
      real*8     sc_o(no)

      integer*4  nh, ih, jh, kh, th
      integer*4  s_h(no)
      integer*4  vv_h(3,3,24)
      real*8     sc_h(24)

      integer*4  vv_x(3,3)
      logical*4  found
      integer*4  ivb

      integer*4  i, j

c --------------------------------------------------------------------
      if( ivb .ge. 3 )then
c         write(*,'(32(''-''))')
      endif

      ns = no
      do is = 1,no
         h_s(is) = 0
      enddo
      nh = 0
      is = 1
      found = .true.
      do while( found )
         ih = nh + 1
         jh = nh + 1
         kh = 1
         h_s(is) = ih
         s_h(ih) = is
         do while( found .and. jh .le. ih )
            call yyy_multiply_ops(
     +         vv_o(1,1,o_s(s_h(kh))),
     +         vv_o(1,1,o_s(s_h(jh))),
     +         vv_x )
            ts = 0
            found = .false.
            do while( .not. found .and. ts .lt. no )
               ts = ts + 1
               call yyy_the_same_ops( vv_x, vv_o(1,1,o_s(ts)), found )
            enddo
            if( found .and. h_s(ts) .eq. 0 )then
               ih = ih + 1
               s_h(ih) = ts
               h_s(ts) = ih
            endif

            if( ivb .ge. 3 )then
c               write(*,'(i4,a2,i3,a1)') ts, '(', h_s(ts), ')'
            endif

            th = kh
            kh = jh
            jh = th
            if( kh .lt. jh )then
               kh = kh + 1
            else if( kh .eq. jh )then
               kh = 1
               jh = jh + 1
            endif
         enddo
         if( found )then
            nh = ih
         else
            do jh = nh+1,ih
               h_s(s_h(jh)) = 0
            enddo
         endif

         if( ivb .ge. 3 )then
c            write(*,'(32(''-''))')
c            write(*,'((256i4))') ( s_h(jh), jh = 1,ih )
c            write(*,'((256i4))') ( s_h(jh), jh = 1,nh )
c            write(*,'(32(''-''))')
         endif

         found = .false.
         do while( .not. found .and. is .lt. ns )
            is = is + 1
            found = h_s(is) .eq. 0
         enddo
         found = found .and. sc_o(o_s(is)) .le. sc_tol
      enddo
      if( nh .gt. 24 ) then
c       stop'yyy_collect_group:a'
        ierr = 1
        return
      endif
      ih = 0
      do is = 1,no
         if( h_s(is) .ne. 0 )then
            ih = ih + 1
            if( ih .gt. nh ) then
c             stop'yyy_collect_group:b'
              ierr = 1
              return
            endif
            s_h(ih) = is
            sc_h(ih) = sc_o(o_s(is))
            do i = 1,3
               do j = 1,3
                  vv_h(j,i,ih) = vv_o(j,i,o_s(is))
               enddo
            enddo
         endif
      enddo

      if( ivb .ge. 3 )then
c         write(*,'((256i4))') ( s_h(jh), jh = 1,nh )
c         write(*,'(32(''-''))')
      endif
      end

c ====================================================================
      subroutine yyy_generate_ops( mo, no, vv_o, ierr )
      implicit   none

      integer*4  it, jt, kt, ierr
      integer*4  v_t(3,26) 

      integer*4  mo, no
      integer*4  vv_o(3,3,mo)

      integer*4  i, j, k, d, tr, trr, m
      integer*4  w(3,3), ww(3,3), www(3,3), wr(3,3)
      logical*4  found

      data   v_t /
     +                       1, 0, 0,
     +                       0, 1, 0,
     +                       0, 0, 1,
     +                       0, 1, 1,
     +                       0, 1,-1,
     +                       1, 0, 1,
     +                      -1, 0, 1,
     +                       1, 1, 0,
     +                       1,-1, 0,
     +                       1, 1, 1,
     +                      -1, 1, 1,
     +                       1,-1, 1,
     +                       1, 1,-1,
     +                         39* 0 /
c --------------------------------------------------------------------
      do it = 1,13
         do i = 1,3
            v_t(i,it+13) = - v_t(i,it)
         enddo
      enddo
      no = 0
      do it = 1,26
         do jt = 1,26
            do kt = 1,26
               do i = 1,3
                  w(i,1) = v_t(i,it)
                  w(i,2) = v_t(i,jt)
                  w(i,3) = v_t(i,kt)
               enddo
               do i = 1,3
                  do j = 1,3
                     ww(j,i) = w(j,i)
                  enddo
               enddo
               call yyy_invert_int( w, wr, d )
               tr = 0
               trr = 0
               do i = 1,3
                  tr = tr + w(i,i)
                  trr = trr + wr(i,i)
               enddo
               found = tr .eq. trr .and. d .eq. 1
               m = 0
               do while( found .and. m .lt. tr + 3 )
                  m = m + 1
                  do i = 1,3
                     do j = 1,3
                        www(j,i) = 0
                        do k = 1,3
                           www(j,i) = www(j,i) + ww(j,k)* w(k,i)
                        enddo
                     enddo
                  enddo
                  do i = 1,3
                     do j = 1,3
                        ww(j,i) = www(j,i)
                     enddo
                  enddo
                  do i = 1,3
                     do j = 1,3
                        found = found .and. abs( www(j,i) ) .le. 1
                     enddo
                  enddo
               enddo
               if( found )then
                  no = no + 1

                  if( no .gt. mo ) then
c                   stop'yyy_generate_ops:a'
                    ierr = 1
                    return
                  endif

                  do i = 1,3
                     do j = 1,3
                        vv_o(j,i,no) = w(j,i)
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo
      end

c ====================================================================
      subroutine yyy_multiply_ops( a, b, c )
      implicit   none

      integer*4  a(3,3), b(3,3), c(3,3)
      integer*4  i, j, k

c --------------------------------------------------------------------
      do i = 1,3
         do j = 1,3
            c(j,i) = 0
            do k = 1,3
               c(j,i) = c(j,i) + a(j,k)* b(k,i)
            enddo
         enddo
      enddo
      end

c ====================================================================
      subroutine yyy_the_same_ops( a, b, found )
      implicit   none

      integer*4  a(3,3), b(3,3)
      integer*4  i, j
      logical*4  found

c --------------------------------------------------------------------
      found = .false.
      do i = 1,3
         do j = 1,3
            if( a(j,i) .ne. b(j,i) )                            return
         enddo
      enddo
      found = .true.
      end

c ====================================================================
      subroutine yyy_sort_key( nc, c, rc )
      implicit   none

      integer*4  ia, ja
      integer*4  ib, jb
      integer*4  ic, jc, kc, nc
      integer*4  c(nc*2)
      real*8     rc(nc)

c --------------------------------------------------------------------
      do ic = 1,nc
         c(ic) = ic
         c(ic+nc) = ic                              ! for nc = 1
      enddo
      jc = nc* 2
      kc = 1
      do while ( kc .lt. nc )
         ia = 1
         ib = 1 + kc
         ic = 1 + nc
         do while( ic .le. jc )
            ja = min( nc, ia + kc - 1 )
            jb = min( nc, ib + kc - 1 )
            do while( ia .le. ja .and. ib .le. jb )
               if( rc(c(ib)) - rc(c(ia)) .ge. 0 )then
                  c(ic) = c(ia)
                  ia = ia + 1
               else
                  c(ic) = c(ib)
                  ib = ib + 1
               endif
               ic = ic + 1
            enddo
            do while( ia .le. ja )
               c(ic) = c(ia)
               ia = ia + 1
               ic = ic + 1
            enddo
            do while( ib .le. jb )
               c(ic) = c(ib)
               ib = ib + 1
               ic = ic + 1
            enddo
            ia = ia + kc
            ib = ib + kc
         enddo
         do ic = 1,nc
            c(ic) = c(ic+nc)
         enddo
         kc = kc* 2
      enddo
      end

c ====================================================================
      subroutine yyy_invert_int( b, r, d )
      implicit   none

      integer*4  b(3,3), r(3,3), d
      integer*4  c3(5) 
      save       c3
      integer*4  i1, i2, i3
      integer*4  j1, j2, j3
      data   c3/ 1, 2, 3, 1, 2 /
c --------------------------------------------------------------------
      do i1 = 1,3
         i2 = c3(i1+1)
         i3 = c3(i1+2)
         do j1 = 1,3
            j2 = c3(j1+1)
            j3 = c3(j1+2)
            r(i1,j1) = b(j2,i2)* b(j3,i3) - b(j3,i2)* b(j2,i3)
         enddo
      enddo
      d = 0
      do i1 = 1,3
         d = d + b(1,i1)* r(i1,1)
      enddo
      end

c ====================================================================
