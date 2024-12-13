!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM AF
!
! #CHASE_UPDATE_STAMP
subroutine SERIAL_HERMITIAN_feweigenvalues(n,M,E,neig)
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
 !
 character(64)     :: subname="SERIAL_HERMITIAN_feweigenvalues"
 integer           :: lwork,nfound,info
 type(LALGEBRA_WS) :: WS
 integer, allocatable :: iwork(:),ifail(:)
 real(SP), external :: DLAMCH,SLAMCH
 complex(SP), allocatable :: Vtemp(:,:)
 real(SP), allocatable :: Etemp(:)

 allocate(Etemp(n))
 allocate(Vtemp(n,n))
 Vtemp = 0
 !
 lwork=-1
 allocate(WS%v_cmplx(1))
 !
 if(.not.(l_nl_optics.or.l_real_time)) then
   call timing('SERIAL_HERMITIAN_feweigenvalues',OPR='start')
   YAMBO_ALLOC(WS%v_real,(max(1,7*n)))
 else
   allocate(WS%v_real(max(1,7*n)))
 endif
 !
 abstol = 2*DLAMCH('S')
 allocate(iwork(max(1,5*n)))
        !print *, M(:,:)
 allocate(ifail(n))
 info = 0
 ifail = 0
 E = 0.0
 Vtemp = 0.0
 iwork = 0
 nfound = 0

#if defined _DOUBLE
 abstol = 2*DLAMCH('S')
 call ZHEEVX('V','I','U',n,M,n,1,1,1,neig,abstol,nfound,Etemp,Vtemp&
&       ,n,WS%v_cmplx,lwork,WS%v_real,iwork,ifail,info)
#else
 abstol = 2*SLAMCH('S')
 call CHEEVX('V','I','U',n,M,n,1,1,1,neig,abstol,nfound,Etemp,Vtemp&
&       ,n,WS%v_cmplx,lwork,WS%v_real,iwork,ifail,info)
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
 call ZHEEVX('V','I','U',n,M,n,1,1,1,neig,abstol,nfound,Etemp,Vtemp&
&       ,n,WS%v_cmplx,lwork,WS%v_real,iwork,ifail,info)
 if(info.ne.0) call LINEAR_ALGEBRA_error(subname,'performing ZHEEV')
#else
 call CHEEVX('V','I','U',n,M,n,1,1,1,neig,abstol,nfound,Etemp,Vtemp&
&       ,n,WS%v_cmplx,lwork,WS%v_real,iwork,ifail,info)
 if(info.ne.0) call LINEAR_ALGEBRA_error(subname,'performing CHEEV')
#endif
 !
 call LINEAR_ALGEBRA_WS_reset(WS)
 !
 deallocate(iwork)
 M=Vtemp
 do ii = 1,neig
        E(ii) = Etemp(ii)
 end do
 deallocate(Etemp)
 
 if(.not.(l_nl_optics.or.l_real_time)) then
   call timing('SERIAL_HERMITIAN_feweigenvalues',OPR='stop')
 endif
 !
end subroutine
