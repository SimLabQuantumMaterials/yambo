!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM AF
!
! #CHASE STAMP      

subroutine PARALLEL_HERMITIAN_chase(Mat, V, N, E, neig)
 !
 use pars,           ONLY:SP
 use parallel_int,   ONLY:PP_redux_wait
 use SLK_m,          ONLY:SLK_POOL,SLK_ORTHO
 use linear_algebra, ONLY:LINEAR_ALGEBRA_error,LINEAR_ALGEBRA_WS_reset,LALGEBRA_WS
 use parallel_m,     ONLY:MPI_COMM_WORLD
 use matrix,         ONLY:PAR_matrix
 use timing_m,       ONLY:timing
#include<y_memory.h>
 !
 integer           :: N,neig,i,j,k,nr,nb,nprow,npcol,m,myrow,mycol
 type(PAR_matrix)  :: Mat,V
 real(SP)          :: E(neig)

 integer :: nev, nex, init, deg, flag
 real(SP) :: tol
 character :: mode, opt, qr, grid_major
 !
 character(64)     :: subname="PARALLEL_HERMITIAN_chase"
 complex(SP), allocatable :: Vtemp(:,:)
 real(SP), allocatable :: Etemp(:)

 call timing('SERIAL_HERMITIAN_chase',OPR='start')
!
 nev = min(neig,N)
 nex = max(0,min(10,N-nev))

 allocate(Vtemp(m,nev+nex))
 allocate(Etemp(nev+nex))

 Vtemp = 0.0
 Etemp = 0.0

 deg = 20
 tol = 1e-10
 mode = 'R'
 opt = 'S'
 qr = 'C'
 grid_major = 'R'

 nprow = SLK_ORTHO%grid(1)
 npcol = SLK_ORTHO%grid(2)
 myrow = SLK_ORTHO%coordinate(1)
 mycol = SLK_ORTHO%coordinate(2)
 nb    = Mat%BLCnrows
 m     = Mat%nrows

#if defined _DOUBLE
 print *, "parallel double go!"
 call pzchase_init_blockcyclic(N, nev, nex, nb, nb, Mat%blc(:,:,Mat%I), m, Vtemp, Etemp,&
&                              nprow, npcol, 0, 0, grid_major, MPI_COMM_WORLD, init)
#else
 print *, "parallel simple go!"
 call czchase_init_blockcyclic(N, nev, nex, nb, nb, atM%blc(:,:,Mat%I), m, Vtemp, Etemp,&
&                              nprow, npcol, 0, 0, grid_major, MPI_COMM_WORLD, init)
#endif

 print *, "Initialized"

#if defined _DOUBLE
 call pzchase(deg,tol,mode,opt,qr)
#else
 call pcchase(deg,tol,mode,opt,qr)
#endif
 
print *, "Solved"

#if defined _DOUBLE
 call pzchase_finalize(flag)
#else
 call pcchase_finalize(flag)
#endif

 if(flag/=0) print *, "[WARNING] ChASE did not finalize properly"
 !
 do i = 1,nev
        E(i) = Etemp(i)
 end do

 j = 0
 k = 1
 do i = 1, nev, nb
        nr = nb
        if(nev - i < nb) then
                nr = nev - i  + 1
        end if
        if(mycol == j) then
                V%blc(:,k:k+nr-1,V%I) = Vtemp(:,i:i+nr-1)
                k = k + nr
        end if
        j = MOD(j + 1, nprow)
 end do

 print *, E(:)
 deallocate(Etemp)
 call timing('PARALLEL_HERMITIAN_chase',OPR='stop')
 !
end subroutine
