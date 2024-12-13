!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM AF
!
subroutine SERIAL_HERMITIAN_diagonalization(n,M,E)
 !
 use pars,           ONLY:SP
 use drivers,        ONLY:l_nl_optics,l_real_time
 use linear_algebra, ONLY:LINEAR_ALGEBRA_error,&
&                         LINEAR_ALGEBRA_WS_reset,LALGEBRA_WS
 use timing_m,       ONLY:timing
#include<y_memory.h>
 !
 integer     :: n
 complex(SP) :: M(n,n)
 real(SP)    :: E(n)
 !
 character(64)     :: subname="SERIAL_HERMITIAN_diagonalization"
 integer           :: lwork
 type(LALGEBRA_WS) :: WS
 !
 lwork=-1
 allocate(WS%v_cmplx(1))
 !
 if(.not.(l_nl_optics.or.l_real_time)) then
   call timing('SERIAL_HERMITIAN_diagonalization',OPR='start')
   YAMBO_ALLOC(WS%v_real,(max(1,3*n-2)))
 else
   allocate(WS%v_real(max(1,3*n-2)))
 endif
 !
 print *, "M = ",M(1:3,1:3)
 print *, "E = ",E(:)
 print *, "WS%v_real = ",WS%v_real(1:10)
#if defined _DOUBLE
 call ZHEEV('V','U',n,M,size(M,1),E,WS%v_cmplx,lwork,WS%v_real,WS%i_fail)
#else
 call CHEEV('V','U',n,M,size(M,1),E,WS%v_cmplx,lwork,WS%v_real,WS%i_fail)
#endif
 !
 lwork=nint(real(WS%v_cmplx(1)))
 deallocate(WS%v_cmplx)
 !
 if(.not.(l_nl_optics.or.l_real_time)) then
   YAMBO_ALLOC(WS%v_cmplx,(lwork))
 else
   allocate(WS%v_cmplx(lwork))
 endif
 !
#if defined _DOUBLE
 call ZHEEV('V','U',n,M,size(M,1),E,WS%v_cmplx,lwork,WS%v_real,WS%i_fail)
 if(WS%i_fail.ne.0) call LINEAR_ALGEBRA_error(subname,'performing ZHEEV')
#else
 call CHEEV('V','U',n,M,size(M,1),E,WS%v_cmplx,lwork,WS%v_real,WS%i_fail)
 if(WS%i_fail.ne.0) call LINEAR_ALGEBRA_error(subname,'performing CHEEV')
#endif
 !
 call LINEAR_ALGEBRA_WS_reset(WS)
 !
 if(.not.(l_nl_optics.or.l_real_time)) then
   call timing('SERIAL_HERMITIAN_diagonalization',OPR='stop')
 endif
 !
end subroutine

subroutine SERIAL_HERMITIAN_chase(n,M,E,neig)
 !
 use pars,           ONLY:SP
 use drivers,        ONLY:l_nl_optics,l_real_time
 use linear_algebra, ONLY:LINEAR_ALGEBRA_error,&
&                         LINEAR_ALGEBRA_WS_reset,LALGEBRA_WS
 use timing_m,       ONLY:timing
#include<y_memory.h>
 !
 integer     :: n,neig,ii
 complex(SP) :: M(n,n)
 real(SP)    :: E(neig),abstol

 integer :: nev, nex, init, deg, flag
 real(SP) :: tol
 character :: mode, opt, qr
 !
 character(64)     :: subname="SERIAL_HERMITIAN_chase"
 complex(SP), allocatable :: Vtemp(:,:)
 real(SP), allocatable :: Etemp(:)

 if(.not.(l_nl_optics.or.l_real_time)) then
   call timing('SERIAL_HERMITIAN_chase',OPR='start')
 endif
 allocate(Vtemp(n,n))
 Vtemp = 0.0
!
 nev = min(neig,n)
 nex = max(0,min(10,n-nev))
 deg = 20
 tol = 1e-10
 mode = 'R'
 opt = 'S'
 qr = 'C'
 print *, M(1:3,1:3)
 print *, Vtemp(1:3,1:3)
 allocate(Etemp(nev+nex))
 Etemp = 0.0
 print *, Etemp(1:10)
 print *, "nev = ", nev, "+ nex = ", nex, "n = ", n
 
#if defined _DOUBLE
 print *, "double go!"
 call zchase_init(n, nev, nex, M, n, Vtemp, Etemp, init)
#else
 print *, "simple go!"
 call cchase_init(n, nev, nex, M, n, Vtemp, Etemp, init)
#endif

 print *, "Initialized"

#if defined _DOUBLE
 call zchase(deg,tol,mode,opt,qr)
#else
 call cchase(deg,tol,mode,opt,qr)
#endif
 
print *, "Solved"

#if defined _DOUBLE
 call zchase_finalize(flag)
#else
 call cchase_finalize(flag)
#endif

 if(flag/=0) print *, "[WARNING] ChASE did not finalize properly"
 !
 M=Vtemp
 do ii = 1,nev
        E(ii) = Etemp(ii)
 end do
 print *, E(:)
 deallocate(Etemp)
 if(.not.(l_nl_optics.or.l_real_time)) then
   call timing('SERIAL_HERMITIAN_chase',OPR='stop')
 endif
 !
end subroutine
