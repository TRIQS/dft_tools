
c ******************************************************************************
c  
c   TRIQS: a Toolbox for Research in Interacting Quantum Systems
c  
c   Copyright (C) 2011 by L. Pourovskii, V. Vildosola, C. Martins, M. Aichhorn
c  
c   TRIQS is free software: you can redistribute it and/or modify it under the
c   terms of the GNU General Public License as published by the Free Software
c   Foundation, either version 3 of the License, or (at your option) any later
c   version.
c  
c   TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
c   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
c   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
c   details.
c  
c   You should have received a copy of the GNU General Public License along with
c   TRIQS. If not, see <http://www.gnu.org/licenses/>.
c  
c *****************************************************************************/

        SUBROUTINE orthogonal_wannier
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine orthonormalizes the Wannier-like functions      %%
C %% obtained with the projectors P(icrorb,ik,is), in order to       %%
C %% get a set of "true" Wannier orbitals.                           %%
C %%                                                                 %%
C %% Only the correlated orbitals are treated here.                  %%
C %%                                                                 %%
C %%          THIS VERSION CAN NOT BE USED WITH SPIN-ORBIT           %%
C %% (since the calculation is made independently for up/dn states)  %%
C %%   THIS VERSION CAN BE USED WITH SPIN-POLARIZED INPUT FILES.     %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE almblm_data
        USE common_data
        USE prnt
        USE projections
        USE reps
        IMPLICIT NONE
        COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Dmat, D_orth, D
        INTEGER :: is, ik, l, nbnd, ndim, isrt, nbbot, nbtop
        INTEGER :: icrorb, ind1, ind2, ib, iatom
        INTEGER :: m1, m2, irep
C
        WRITE(buf,'(a)')'Orthonormalization of the projectors...'
        CALL printout(0)
        CALL printout(0)
C
        IF(ncrorb==0) RETURN
C 
C =====================================
C Creation of the overlap matrix Dmat :
C =====================================
C
C -----------------------------------------------------------
C Determination of the dimension ndim of the overlap matrix :
C -----------------------------------------------------------
        ndim=0
C Loop on the correlated orbitals
        DO icrorb=1,ncrorb
          isrt=crorb(icrorb)%sort
          l=crorb(icrorb)%l
C Since this subroutine is used only in the case without SO,
C the correlated ireps can be considered if there are any. (ifsplit=.TRUE.)
          IF(crorb(icrorb)%ifsplit) THEN
C the value of l can not be 0 here, because ifsplit is necessary .FALSE. 
C for s-orbital (restriction in dmftproj.f)
            DO irep=1,reptrans(l,isrt)%nreps
              IF(crorb(icrorb)%correp(irep))
     &         ndim=ndim+reptrans(l,isrt)%dreps(irep)
C The dimension of the irep is added to ndim. 
            ENDDO
          ELSE
C If no particular irep is considered (ifsplit=.FALSE.),
C The whole matrix of the representation is considered. 
           ndim=ndim+2*l+1
          ENDIF
        ENDDO 
C ------------------
C Creation of Dmat :
C ------------------
        ALLOCATE(Dmat(1:ndim,1:ndim))
C 
C =====================================================================
C Computation of the orthonormalized Wannier functions and projectors :  
C =====================================================================
C The computation is performed for each k_point and each spin-value independently 
C because they are good quantum numbers.
        DO ik=1,nk 
          DO is=1,ns
C Only the k-points with inlcuded bands are considered for the projectors.
            IF(.NOT.kp(ik,is)%included) CYCLE
            nbnd=kp(ik,is)%nb_top-kp(ik,is)%nb_bot+1
            nbbot=kp(ik,is)%nb_bot
            nbtop=kp(ik,is)%nb_top
            ALLOCATE(D(1:ndim,1:nbnd))
C
C --------------------------------
C Initialization of the D matrix :
C --------------------------------
C This D matrix of size ndim*nbnd is the complete "projector matrix" 
C which enables to go from the Wannier-like basis |u_orb> to the Bloch states |ik,ib>. 
            ind1=0
            DO icrorb=1,ncrorb
              isrt=crorb(icrorb)%sort
              l=crorb(icrorb)%l
C If l=0, there only possible irep is the whole matrix itself.
              IF (l==0) THEN
               D(ind1+1,1:nbnd)=pr_crorb(icrorb,ik,is)%
     &              mat_rep(1,nbbot:nbtop)
               ind1=ind1+1
              ELSE
C the projectors of the correlated ireps are considered if there are any. (ifsplit=.TRUE.)
               IF(crorb(icrorb)%ifsplit) THEN
C the value of l can not be 0 here, because ifsplit is necessary .FALSE. 
C for s-orbital (restriction in dmftproj.f)
                m1=-l-1
                DO irep=1,reptrans(l,isrt)%nreps
                  IF(crorb(icrorb)%correp(irep)) THEN
                   m2=m1+reptrans(l,isrt)%dreps(irep)
                   ind2=ind1+reptrans(l,isrt)%dreps(irep)
C Since there is no SO, prcrorb%matrep is of size 2*l+1, from -l to l
C (the basis which mix up/dn states are not possible here.)
C The states range from m1+1 to m2 in the irep.   
C The corresponding projector is stored from the line (ind1+1) to the line ind2, in the D matrix.
                   D(ind1+1:ind2,1:nbnd)=pr_crorb(icrorb,ik,is)%
     &               mat_rep(m1+1:m2,nbbot:nbtop)
                   ind1=ind2
                  ENDIF
                  m1=m1+reptrans(l,isrt)%dreps(irep)
                ENDDO
               ELSE
C The projectors of the whole correlated representation is considered. (ifsplit=.FALSE.)
                ind2=ind1+2*l+1
C Since there is no SO, prcrorb%matrep is of size 2*l+1, from -l to l.
C (the basis which mix up/dn states are not possible here.)
C The corresponding projection matrix is stored from the line (ind1+1) to the line ind2, in the D matrix.
                D(ind1+1:ind2,1:nbnd)=pr_crorb(icrorb,ik,is)%
     &            mat_rep(-l:l,nbbot:nbtop)
                ind1=ind2
               ENDIF   ! End of the ifsplit if-then-else
              ENDIF    ! End of the l=0 if-then-else  
            ENDDO      ! End of the icrorb loop
C
C ----------------------------------------
C Computation of the overlap matrix Dmat : 
C ----------------------------------------
C The overlap matrix is stored in Dmat = D*transpose(conjugate(D))
            CALL ZGEMM('N','C',ndim,ndim,nbnd,DCMPLX(1.D0,0.D0),
     &        D,ndim,D,ndim,DCMPLX(0.D0,0.D0),Dmat,ndim)
C
C -------------------------------------------
C Computation of the matrix S = Dmat^{-1/2} : 
C -------------------------------------------
            CALL orthogonal_h(Dmat,ndim,.TRUE.)
C This matrix is stored in Dmat.
C
C -----------------------------------------------
C Computation of the orthonormalized projectors : 
C -----------------------------------------------
C The calculation performed is the following : P=O^(-1/2)*P_tilde.
C Its value is stored in the matrix D_orth (of size ndim*nbnd)

            ALLOCATE(D_orth(1:ndim,1:nbnd))
            CALL ZGEMM('N','N',ndim,nbnd,ndim,DCMPLX(1.D0,0.D0),
     &        Dmat,ndim,D,ndim,DCMPLX(0.D0,0.D0),D_orth,ndim)
            DEALLOCATE(D)
C
C --------------------------------------------------------------------------------
C Storing the value of the orthonormalized projectors in the pr_crorb structures :
C --------------------------------------------------------------------------------
            ind1=0
            DO icrorb=1,ncrorb
              isrt=crorb(icrorb)%sort
              l=crorb(icrorb)%l
C If l=0, there only possible irep is the whole matrix itself.
              IF (l==0) THEN
               pr_crorb(icrorb,ik,is)%mat_rep
     &            (1,nbbot:nbtop)=D_orth(ind1+1,1:nbnd)
               ind1=ind1+1
              ELSE
C the projectors of the correlated ireps are considered if there are any. (ifsplit=.TRUE.)
               IF(crorb(icrorb)%ifsplit) THEN
C the value of l can not be 0 here, because ifsplit is necessary .FALSE. 
C for s-orbital (restriction in dmftproj.f)
                m1=-l-1
                DO irep=1,reptrans(l,isrt)%nreps
                  IF(crorb(icrorb)%correp(irep)) THEN
                   m2=m1+reptrans(l,isrt)%dreps(irep)
                   ind2=ind1+reptrans(l,isrt)%dreps(irep)
C prcrorb%matrep is of size 2*l+1, from -l to l (the basis which mix up/dn states are not possible here.)
C In the D_orth matrix, the corresponding part of the projection matrix ranges from the line (ind1+1) to the line ind2.
C The projector associated to the ireps is stored in the prcrorb%matrep from m1+1 to m2.
                   pr_crorb(icrorb,ik,is)%
     &               mat_rep(m1+1:m2,nbbot:nbtop)=
     &               D_orth(ind1+1:ind2,1:nbnd)
                   ind1=ind2
                  ENDIF
                  m1=m1+reptrans(l,isrt)%dreps(irep)
                ENDDO
               ELSE
C The projectors of the whole correlated representation is considered. (ifsplit=.FALSE.)
                ind2=ind1+2*l+1
C Since there is no SO, prcrorb%matrep is of size 2*l+1, from -l to l.
C (the basis which mix up/dn states are not possible here.)
C In the D_orth matrix, the projection matrix ranges from the line (ind1+1) to the line ind2.
C The projector is stored in the pr_crorb%matrep (from -l to l).
                pr_crorb(icrorb,ik,is)%mat_rep
     &            (-l:l,nbbot:nbtop)=D_orth(ind1+1:ind2,1:nbnd)
                ind1=ind2
               ENDIF ! End of the ifsplit if-then-else
              ENDIF  ! End of the l=0 if-then-else
            ENDDO    ! End of the icrorb loop
C prcrorb%matrep contains now the orthonormalized projectors.
            DEALLOCATE(D_orth)
          ENDDO     ! End of the loop on is
        ENDDO       ! End of the loop on ik
        DEALLOCATE(Dmat)
C
C =============================================================================
C Printing the projectors with k-points 1 and nk in the file fort.18 for test :
C =============================================================================
        DO icrorb=1,ncrorb
          iatom=crorb(icrorb)%atom
          isrt=crorb(icrorb)%sort
          l=crorb(icrorb)%l
          WRITE(18,'()')
          WRITE(18,'(a)') 'apres othonormalizsation'
          WRITE(18,'(a,i4)') 'icrorb = ', icrorb
          WRITE(18,'(a,i4,a,i4)') 'isrt = ', isrt, ' l = ', l
          IF (l==0) THEN
           WRITE(18,'(a,i4)') 'ik = ', 1
           DO ib = kp(1,1)%nb_bot,kp(1,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,1,1)%mat_rep(:,ib)
             IF (ifSP)
     &       WRITE(18,*) pr_crorb(icrorb,1,2)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
           WRITE(18,'(a,i4)') 'ik = ', nk
           DO ib = kp(nk,1)%nb_bot,kp(nk,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,nk,1)%mat_rep(:,ib)
             IF (ifSP)
     &       WRITE(18,*) pr_crorb(icrorb,nk,2)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
           WRITE(18,'(a,i4)') 'ik = ', 1
           DO ib = kp(1,1)%nb_bot,kp(1,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,1,1)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
           WRITE(18,'(a,i4)') 'ik = ', nk
           DO ib = kp(nk,1)%nb_bot,kp(nk,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,nk,1)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
          ELSE
           WRITE(18,'(a,i4)') 'ik = ', 1
           DO ib = kp(1,1)%nb_bot,kp(1,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,1,1)%mat_rep(:,ib)
             IF (ifSP)
     &       WRITE(18,*) pr_crorb(icrorb,1,2)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
           WRITE(18,'(a,i4)') 'ik = ', nk
           DO ib = kp(nk,1)%nb_bot,kp(nk,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,nk,1)%mat_rep(:,ib)
             IF (ifSP)
     &       WRITE(18,*) pr_crorb(icrorb,nk,2)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
          ENDIF
        ENDDO
C
      RETURN
      END
       
 

        SUBROUTINE orthogonal_wannier_SO
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine orthonormalizes the Wannier-like functions      %%
C %% obtained with the projectors P(icrorb,ik,is), in order to       %%
C %% get a set of "true" Wannier orbitals.                           %%
C %%                                                                 %%
C %% Only the correlated orbitals are treated here.                  %%
C %%                                                                 %%
C %%          THIS VERSION MUST BE USED WITH SPIN-ORBIT              %%
C %% (since the calculation for up/dn states is made simultaneously) %%
C %% THIS VERSION CAN NOT BE USED WITHOUT SPIN-POLARIZED INPUT FILES.%%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE almblm_data
        USE common_data
        USE prnt
        USE projections
        USE reps
        IMPLICIT NONE
        COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Dmat, D_orth, D
        INTEGER :: is, ik, l, nbnd, ndim, isrt, nbbot, nbtop
        INTEGER :: icrorb, ind1, ind2, iatom, ib
        INTEGER :: m1, m2, irep
C
        WRITE(buf,'(a)')'Orthonormalization of the projectors...'
        CALL printout(0)
        CALL printout(0)
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if there is no dn part of pr_crorb.
C -------------------------
C
        IF(.not.ifSP) THEN
         WRITE(buf,'(a,a,i2,a)')'The projectors on ',
     &     'the dn states are required for isrt = ',isrt,
     &     ' but there is no spin-polarized input files.'
         CALL printout(0)
         WRITE(buf,'(a)')'END OF THE PRGM'
         CALL printout(0)
         STOP
        ENDIF
C ---------------------------------------------------------------------------------------
C 
C =====================================
C Creation of the overlap matrix Dmat :
C =====================================
C
C -----------------------------------------------------------
C Determination of the dimension ndim of the overlap matrix :
C -----------------------------------------------------------
        ndim=0
C Loop on the correlated orbitals
        DO icrorb=1,ncrorb
          isrt=crorb(icrorb)%sort
          l=crorb(icrorb)%l
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
          IF (l==0) THEN
C Since this subroutine is used only in the case with SO,
C the only irep possible for s-orbital is the matrix itself.
           ndim=ndim+2
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C the projectors of the correlated ireps are considered if there are any. (ifsplit=.TRUE.)
           IF(crorb(icrorb)%ifsplit) THEN
            DO irep=1,reptrans(l,isrt)%nreps
              IF(crorb(icrorb)%correp(irep)) THEN
               ndim=ndim+reptrans(l,isrt)%dreps(irep)
              ENDIF
C The dimension of the irep is added to ndim. 
            ENDDO
           ELSE
C If no particular irep is considered (ifsplit=.FALSE.),
C The whole matrix of the representation is considered. 
            ndim=ndim+2*(2*l+1)
           ENDIF
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
          ELSE
C Since this subroutine is used only in the case with SO,
C the only irep possible for this orbital is the matrix itself.
           ndim=ndim+2*(2*l+1)
          ENDIF
        ENDDO 
C ------------------
C Creation of Dmat :
C ------------------
        ALLOCATE(Dmat(1:ndim,1:ndim))
C 
C =====================================================================
C Computation of the orthonormalized Wannier functions and projectors :  
C =====================================================================
C The computation is performed for each k_point independently 
C because they are still good quantum numbers.
        DO ik=1,nk
C Only the k-points with inlcuded bands are considered for the projectors.
          IF(.NOT.kp(ik,1)%included) CYCLE
          nbnd=kp(ik,1)%nb_top-kp(ik,1)%nb_bot+1
          nbbot=kp(ik,1)%nb_bot
          nbtop=kp(ik,1)%nb_top
C it was checked that nbtop(up)=nbtop(dn) and nbbot(up)=nbbot(dn) 
C for a computation with SO [in set_projections.f]
          ALLOCATE(D(1:ndim,1:nbnd))
C
C --------------------------------
C Initialization of the D matrix :
C --------------------------------
C This D matrix of size ndim*nbnd is the complete "projector matrix" 
C which enables to go from the Wannier-like basis |u_orb> to the Bloch states |ik,ib>. 
          ind1=0
          DO icrorb=1,ncrorb
            isrt=crorb(icrorb)%sort
            l=crorb(icrorb)%l
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
            IF (l==0) THEN
C the only irep possible for s-orbital is the matrix itself.
             DO is=1,ns
C               D(ind1,1:nbnd)=
C Bug correction 8.11.2012
               D(ind1+1,1:nbnd)=
     &           pr_crorb(icrorb,ik,is)%mat_rep(1,nbbot:nbtop)
               ind1=ind1+1
             ENDDO
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
            ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the projection matrix is stored in prcrorb%matrep with is=1.
C the projectors of the correlated ireps are considered if there are any. (ifsplit=.TRUE.)
             IF (crorb(icrorb)%ifsplit) THEN
              m1=0
              DO irep=1,reptrans(l,isrt)%nreps
                IF (crorb(icrorb)%correp(irep)) THEN
                 m2=m1+reptrans(l,isrt)%dreps(irep)
                 ind2=ind1+reptrans(l,isrt)%dreps(irep)
C The states range from m1+1 to m2 in the irep.   
C The corresponding projector is stored from the line (ind1+1) to the line ind2, in the D matrix.
                 D(ind1+1:ind2,1:nbnd)=pr_crorb(icrorb,ik,1)%
     &             mat_rep(m1+1:m2,nbbot:nbtop)
                 ind1=ind2
                ENDIF
                m1=m1+reptrans(l,isrt)%dreps(irep)
              ENDDO
             ELSE
C The projectors of the whole correlated representation is considered. (ifsplit=.FALSE.)
              ind2=ind1+2*(2*l+1)
C The corresponding projection matrix is stored from the line (ind1+1) to the line ind2, in the D matrix.
              D(ind1+1:ind2,1:nbnd)=pr_crorb(icrorb,ik,1)%
     &          mat_rep(1:2*(2*l+1),nbbot:nbtop)
              ind1=ind2
             ENDIF    ! End of the ifsplit if-then-else
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
            ELSE
C the only irep possible for such an orbital is the matrix itself.
             DO is=1,ns
               ind2=ind1+2*l+1
               D(ind1+1:ind2,1:nbnd)=
     &           pr_crorb(icrorb,ik,is)%mat_rep(-l:l,nbbot:nbtop)
               ind1=ind2
             ENDDO
            ENDIF  ! End of the ifmixing if-then-else
          ENDDO    ! End of the icrorb loop
C
C ----------------------------------------
C Computation of the overlap matrix Dmat : 
C ----------------------------------------
C The overlap matrix is stored in Dmat = D*transpose(conjugate(D))
          CALL ZGEMM('N','C',ndim,ndim,nbnd,DCMPLX(1.D0,0.D0),
     &      D,ndim,D,ndim,DCMPLX(0.D0,0.D0),Dmat,ndim)
C
C -------------------------------------------
C Computation of the matrix S = Dmat^{-1/2} : 
C -------------------------------------------
          CALL orthogonal_h(Dmat,ndim,.TRUE.)
C This matrix is stored in Dmat.
C
C -----------------------------------------------
C Computation of the orthonormalized projectors : 
C -----------------------------------------------
C The calculation performed is the following : P=O^(-1/2)*P_tilde.
C Its value is stored in the matrix D_orth (of size ndim*nbnd)
          ALLOCATE(D_orth(1:ndim,1:nbnd))
          CALL ZGEMM('N','N',ndim,nbnd,ndim,DCMPLX(1.D0,0.D0),
     &      Dmat,ndim,D,ndim,DCMPLX(0.D0,0.D0),D_orth,ndim)
          DEALLOCATE(D)
C
C --------------------------------------------------------------------------------
C Storing the value of the orthonormalized projectors in the pr_crorb structures :
C --------------------------------------------------------------------------------
          ind1=0
          DO icrorb=1,ncrorb
            isrt=crorb(icrorb)%sort
            l=crorb(icrorb)%l
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
            IF (l==0) THEN
C the only irep possible for s-orbital is the matrix itself.
             DO is=1,ns
               pr_crorb(icrorb,ik,is)%mat_rep(1,nbbot:nbtop)=
     &           D_orth(ind1+1,1:nbnd)
               ind1=ind1+1
             ENDDO
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
            ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C the projectors of the correlated ireps are considered if there are any. (ifsplit=.TRUE.)
             IF(crorb(icrorb)%ifsplit) THEN
              m1=0
              DO irep=1,reptrans(l,isrt)%nreps
                IF (crorb(icrorb)%correp(irep)) THEN
                 m2=m1+reptrans(l,isrt)%dreps(irep)
                 ind2=ind1+reptrans(l,isrt)%dreps(irep)
C In the D_orth matrix, the corresponding part of the projection matrix ranges from the line (ind1+1) to the line ind2.
C The projector associated to the ireps is stored in the prcrorb%matrep from m1+1 to m2.
                 pr_crorb(icrorb,ik,1)%mat_rep(m1+1:m2,nbbot:nbtop)
     &             =D_orth(ind1+1:ind2,1:nbnd)
                 ind1=ind2
                ENDIF
                m1=m1+reptrans(l,isrt)%dreps(irep)
              ENDDO
             ELSE
C The projectors of the whole correlated representation is considered. (ifsplit=.FALSE.)
              ind2=ind1+2*(2*l+1)
C The corresponding projection matrix is stored from the line (ind1+1) to the line ind2, in the D matrix.
              pr_crorb(icrorb,ik,1)%mat_rep(1:2*(2*l+1),nbbot:nbtop)
     &          =D_orth(ind1+1:ind2,1:nbnd)
              ind1=ind2
             ENDIF    ! End of the ifsplit if-then-else
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
            ELSE
C the only irep possible for this orbital is the matrix itself.
             DO is=1,ns
               ind2=ind1+2*l+1
               pr_crorb(icrorb,ik,is)%mat_rep(-l:l,nbbot:nbtop)
     &           =D_orth(ind1+1:ind2,1:nbnd)
               ind1=ind2
             ENDDO
            ENDIF  ! End of the ifmixing if-then-else
          ENDDO    ! End of the icrorb loop
          DEALLOCATE(D_orth)
        ENDDO  ! End of the loop on ik
        DEALLOCATE(Dmat)
C
C =============================================================================
C Printing the projectors with k-points 1 and nk in the file fort.18 for test :
C =============================================================================
        DO icrorb=1,ncrorb
          iatom=crorb(icrorb)%atom
          isrt=crorb(icrorb)%sort
          l=crorb(icrorb)%l
          WRITE(18,'()')
          WRITE(18,'(a)') 'apres othonormalizsation'
          WRITE(18,'(a,i4)') 'icrorb = ', icrorb
          WRITE(18,'(a,i4,a,i4)') 'isrt = ', isrt, ' l = ', l
          IF (l==0) THEN
           WRITE(18,'(a,i4)') 'ik = ', 1
           DO ib = kp(1,1)%nb_bot,kp(1,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,1,1)%mat_rep(:,ib)
             IF (ifSP)
     &       WRITE(18,*) pr_crorb(icrorb,1,2)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
           WRITE(18,'(a,i4)') 'ik = ', nk
           DO ib = kp(nk,1)%nb_bot,kp(nk,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,nk,1)%mat_rep(:,ib)
             IF (ifSP)
     &       WRITE(18,*) pr_crorb(icrorb,nk,2)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
           WRITE(18,'(a,i4)') 'ik = ', 1
           DO ib = kp(1,1)%nb_bot,kp(1,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,1,1)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
           WRITE(18,'(a,i4)') 'ik = ', nk
           DO ib = kp(nk,1)%nb_bot,kp(nk,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,nk,1)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
          ELSE
           WRITE(18,'(a,i4)') 'ik = ', 1
           DO ib = kp(1,1)%nb_bot,kp(1,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,1,1)%mat_rep(:,ib)
             IF (ifSP)
     &       WRITE(18,*) pr_crorb(icrorb,1,2)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
           WRITE(18,'(a,i4)') 'ik = ', nk
           DO ib = kp(nk,1)%nb_bot,kp(nk,1)%nb_top
             WRITE(18,*) pr_crorb(icrorb,nk,1)%mat_rep(:,ib)
             IF (ifSP)
     &       WRITE(18,*) pr_crorb(icrorb,nk,2)%mat_rep(:,ib)
             WRITE(18,'()')
           ENDDO
          ENDIF
        ENDDO
C
        RETURN
        END
       
 
