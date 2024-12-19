!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM AF
!
! #CHASE STAMP      

subroutine PARALLEL_HERMITIAN_chase(Mat, V, N, E, neig)

 use chase_diag
 !
 use pars,           ONLY:SP
 use parallel_int,   ONLY:PP_redux_wait
 use SLK_m,          ONLY:SLK_POOL,SLK_ORTHO
 use linear_algebra, ONLY:LINEAR_ALGEBRA_error,LINEAR_ALGEBRA_WS_reset,LALGEBRA_WS
 use parallel_m,     ONLY:PAR_COM_WORLD
 use matrix,         ONLY:PAR_matrix
 use timing_m,       ONLY:timing
#include<y_memory.h>
 !
 integer           :: N,neig,i,j,k,nr,nb,nprow,npcol,m,myrow,mycol,nl
 type(PAR_matrix)  :: Mat,V
 real(SP)          :: E(neig)

 integer  :: nev, nex, init, deg, flag
 integer  :: comm
 real(SP) :: tol
 character :: mode, opt, qr, grid_major
 !
 character(64)     :: subname="PARALLEL_HERMITIAN_chase"
 complex(SP), allocatable :: Vtemp(:,:)
 real(SP), allocatable :: Etemp(:)

 integer :: myrank_intra, myrank_inter, intra_size, inter_size,&
&           myrank_slk_ortho_comm, slk_ortho_comm_size,ierr

 call timing('PARALLEL_HERMITIAN_chase',OPR='start')
!
 nev = min(neig,N)
 nex = max(0,min(10,N-nev))

 deg = 20
 tol = 1e-10
 mode = 'R'
 opt = 'S'
 qr = 'C'
 grid_major = 'R' !----- Check grid major /!\.

 nprow = SLK_ORTHO%grid(1)
 npcol = SLK_ORTHO%grid(2)
 myrow = SLK_ORTHO%coordinate(1)
 mycol = SLK_ORTHO%coordinate(2)
 nb    = Mat%BLCnrows
 m     = Mat%nrows
 nl    = Mat%ncols
 comm  = SLK_ORTHO%COMM

 allocate(Vtemp(m,nev+nex))
 allocate(Etemp(nev+nex))

 Vtemp = 0.0
 Etemp = 0.0
 
 call mpi_comm_rank(comm,myrank_slk_ortho_comm,ierr)
 call mpi_comm_size(comm,slk_ortho_comm_size  ,ierr)

! print *, "-------- Grid Info : ---------"
! print *, "My rank is ", myrank_slk_ortho_comm, " / ",&
!&                               slk_ortho_comm_size," in SLK_ORTHO_COMM " 
! print *, "Mat% I = ",Mat%I, " nb = ", nb, " m = ", m
! print *, "Rows = ", Mat%rows


#if defined _DOUBLE
! call pzchase_init(N,nev,nex,m,nl,Mat%blc(:,:,Mat%I),m,Vtemp,Etemp,&
!&                  nprow,npcol,grid_major,comm,init)

 call timing('ChASE_INIT',OPR='start')
 call pzchase_init_blockcyclic(N, nev, nex, nb, nb, Mat%blc(:,:,Mat%I), m, Vtemp, Etemp,&
&                              nprow, npcol, grid_major, 0, 0, SLK_ORTHO%COMM, init)
#else
 call pcchase_init_blockcyclic(N, nev, nex, nb, nb, Mat%blc(:,:,Mat%I), m, Vtemp, Etemp,&
&                              nprow, npcol, grid_major, 0, 0, SLK_ORTHO%COMM, init)
#endif
 call timing('ChASE_INIT',OPR='stop')

#if defined _DOUBLE
 call pzchase(deg,tol,mode,opt,qr)
#else
 call pcchase(deg,tol,mode,opt,qr)
#endif
 
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

 call timing('ChASE_Reorganize',OPR='start')
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
 
call timing('ChASE_Reorganize',OPR='stop')

 call PP_redux_wait(E,COMM=SLK_POOL%INTRA_comm)

 deallocate(Etemp)
 call timing('PARALLEL_HERMITIAN_chase',OPR='stop')
 !
end subroutine PARALLEL_HERMITIAN_chase
