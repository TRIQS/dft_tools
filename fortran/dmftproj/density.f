
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

      SUBROUTINE density(tetr,only_corr,qtot,ifprnt)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine computes the density matrices for all the       %%
C %% inlcuded orbitals and the corresponding total charge, within    %%
C %% the energy window for which the projectors were previously      %%
C %% defined.                                                        %%
C %%                                                                 %%
C %% If tetr = .TRUE., the tetrahedron weights are used.             %%
C %%         = .FALSE., a simple point integration is used.          %%
C %%                                                                 %%
C %% If only_corr = .TRUE., the calculation is performed for the     %%
C %%                        correlated orbitals only.                %%
C %%                                                                 %%
C %% If ifprnt = .TRUE., the density matrices are written in the     %%
C %%                     file case.outdmftpr.                        %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definition of the variables :
C ----------------------------
      USE almblm_data
      USE common_data
      USE prnt
      USE projections
      USE reps
      USE symm
      IMPLICIT NONE
      LOGICAL :: tetr, only_corr, ifprnt
      COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: D
      COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: conj_mat, mat
      INTEGER :: is, is1, ik, iorb, l, nbnd, iss 
      INTEGER :: lm, i, icrorb, m1, m2, m, ib 
      INTEGER :: iatom, jatom, jorb, imu, isrt, isym
      INTEGER :: istart, istart1
      INTEGER :: nbbot, nbtop
      REAL(KIND=8) :: q, qtot
C
        WRITE(buf,'(a)')'Evaluation of the density matrices...'
        CALL printout(0)
        CALL printout(0)
C
C Initialization of the variable nsp and of the table crdensmat
C nsp describes the number of block necessary to describe the complete density matrix.
        nsp=ns
        IF(ifSO) nsp=4
C The only possible cases are therefore :
C nsp=1 for a computation without spin-polarized input files, without SO   [block (up/up+dn/dn)]
C nsp=2 for a computation with spin-polarized input files (but without SO) [blocks up/up and dn/dn]
C nsp=4 for a computation (with spin-polarized input files) with SO        [all the blocks]
C
C
C =============================================================
C Computation of the density matrices for correlated orbitals :
C =============================================================
C
C These computations are performed only if there are correlated orbitals in the system.
        IF(ncrorb.NE.0) THEN
C crdensmat is a table of size nsp*ncrorb. 
C For each correlated orbital icrorb, crdensmat(:,icrorb) is the corresponding density matrix.
        IF(.NOT.ALLOCATED(crdensmat)) THEN
         ALLOCATE(crdensmat(nsp,ncrorb))
        ENDIF
        DO icrorb=1,ncrorb
          DO is=1,nsp
            IF(ALLOCATED(crdensmat(is,icrorb)%mat)) 
     &       DEALLOCATE(crdensmat(is,icrorb)%mat)
            ALLOCATE(crdensmat(is,icrorb)%mat(1,1))
            crdensmat(is,icrorb)%mat(1,1)=0.d0
          ENDDO
        ENDDO
C 
C Loop on the correlated orbitals icrorb
C
        DO icrorb=1,ncrorb
          isrt=crorb(icrorb)%sort
          l=crorb(icrorb)%l
C
C -----------------------------------------------------------------------------------
C The s-orbitals are a particular case of a "non-mixing" basis and are treated here :
C -----------------------------------------------------------------------------------
          IF (l==0) THEN
C The field mat of crdensmat has already the good size (1 scalar element). 
C There's no need of a basis change since the representation of an s-orbital is Identity whatever the basis is.
           DO ik=1,nk
             DO iss=1,nsp
C
C Determination of the block indices :
C ------------------------------------
               IF(iss.LE.2) THEN
                is=iss
                is1=iss 
C If iss=1 (up), is=1 and is1=1 -> Description of the up/up block
C If iss=2 (down), is=2 and is1=2 -> Description of the dn/dn block
               ELSE 
                is=iss-2 
                is1=3-is
C If iss=3, is=1 (up) and is1=2 (down) -> Description of the up/dn block
C If iss=4, is=2 (down) and is1=1 (up) -> Description of the dn/up block
               ENDIF
C Only the k-points with included bands are considered for the projectors.
               IF(.NOT.kp(ik,is)%included) CYCLE
               IF(.NOT.kp(ik,is1)%included) CYCLE
               nbnd=kp(ik,is)%nb_top-kp(ik,is)%nb_bot+1
               nbbot=kp(ik,is)%nb_bot
               nbtop=kp(ik,is)%nb_top
C for the diagonal blocks, is=is1 ; thus nbtop and nbtop are the same for is and is1. 
C for the off-diagonal blocks (calculated only when SO is considered), 
C it was checked that nbtop(up)=nbtop(dn) and nbbot(up)=nbbot(dn) [in set_projections.f]
C thus nbtop and nbtop fit again for is and is1.
C
C Computation of the density matrix using the tetrahedron weights for the integration :
C -------------------------------------------------------------------------------------
               IF(tetr) THEN
C The field pr_crorb%mat_rep is used to perform the computation (well defined for s-orbitals)  
                ALLOCATE(mat(1,nbbot:nbtop),conj_mat(1,nbbot:nbtop))
                mat(1,nbbot:nbtop)=
     &            pr_crorb(icrorb,ik,is)%mat_rep(1,nbbot:nbtop)*
     &            SQRT(kp(ik,is)%tetrweight(nbbot:nbtop))
                conj_mat(1,nbbot:nbtop)=CONJG( 
     &            pr_crorb(icrorb,ik,is1)%mat_rep(1,nbbot:nbtop))*
     &            SQRT(kp(ik,is1)%tetrweight(nbbot:nbtop))
C mat = P(icrorb,ik,is)*sqrt(tetrahedron-weight(ik,is))
C conj_mat = conjugate[ P(icrorb,ik,is1)) ]*sqrt(tetrahedron-weight(ik,is1))
                DO ib = nbbot,nbtop
                  crdensmat(iss,icrorb)%mat(1,1)=
     =              crdensmat(iss,icrorb)%mat(1,1)+
     +              mat(1,ib)*conj_mat(1,ib)
                ENDDO
C crdensmat = mat*transpose(conj_mat) which is a matrix of size 1
C The summation over the k-points is done with the do loop ; 
C crdensmat(iss,icrorb) is therefore the block "iss" of the density matrix for the orbital icrorb.
                DEALLOCATE(conj_mat,mat)
C
C Computation of the density matrix using a simple point integration :
C --------------------------------------------------------------------
               ELSE
                DO ib = nbbot, nbtop
                  crdensmat(iss,icrorb)%mat(1,1)=
     =              crdensmat(iss,icrorb)%mat(1,1)+
     +              pr_crorb(icrorb,ik,is)%mat_rep(1,ib)*
     *              CONJG(pr_crorb(icrorb,ik,is1)%mat_rep(1,ib))*
     *              kp(ik,is)%weight
                ENDDO
C crdensmat = P(icrorb,ik,is)*transpose(conjugate(P(icrorb,ik,is1)))*weight(ik,is) which is a matrix of size 1.
C The weight used is a geometric factor associated to k and does not depend on the variable is. 
C That's why we merely multiply by the "weight" each term while summing over the k-points.
               ENDIF
             ENDDO    ! End of the iss loop  
           ENDDO      ! End of the ik loop
C
C -----------------------------------------------------------------------------------------------------
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) ) :
C -----------------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C If this option is used, then ifSO=.TRUE. (because of the restriction in set_ang_trans.f)
C Moreover ifSP=.TRUE. (since ifSO => ifSP, in this version) 
C As a result, we know that nb_bot(up)=nb_bot(dn) and nb_top(up)=nb_top(dn)
C
C The field mat of crdensmat must be resized.
C As the complete spinor rotation approach is used, only one matrix is necessary (with is=1).
           DEALLOCATE(crdensmat(1,icrorb)%mat)
           ALLOCATE(crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1)))
           crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=0.d0
           ALLOCATE(D(1:2*(2*l+1),1:2*(2*l+1)))
C
C Computation of the density matrix using the tetrahedron weights for the integration :
C -------------------------------------------------------------------------------------
           IF(tetr) THEN
            DO ik=1,nk
C Only the k-points with inlcuded bands are considered for the projectors.
              IF(.NOT.kp(ik,1)%included) CYCLE
              nbnd=kp(ik,1)%nb_top-kp(ik,1)%nb_bot+1
              nbbot=kp(ik,1)%nb_bot
              nbtop=kp(ik,1)%nb_top
C A distinction between up and dn states is necessary in order to use the tetrahedron weight.
C As a result, we will transform the projectors back in spherical harmonics basis 
C to perform the calculation and then put the resulting density matrix in the desired basis.
C
              ALLOCATE(mat(1:2*(2*l+1),nbbot:nbtop))
              ALLOCATE(conj_mat(1:2*(2*l+1),nbbot:nbtop))
C The representation of the projectors is put back in the spherical harmonics basis
              mat(1:2*(2*l+1),nbbot:nbtop)=MATMUL(TRANSPOSE(CONJG( 
     =          reptrans(l,isrt)%transmat
     &          (1:2*(2*l+1),1:2*(2*l+1)))),pr_crorb(icrorb,ik,1)%
     &          mat_rep(1:2*(2*l+1),nbbot:nbtop))
C mat = inverse(reptrans)*pr_crorb%mat_rep = <lm|new_i>*proj_{new_i} = proj_{lm} [temporarily]
              conj_mat(1:2*(2*l+1),nbbot:nbtop)=MATMUL( 
     &          TRANSPOSE(CONJG(reptrans(l,isrt)%
     &          transmat(1:2*(2*l+1),1:2*(2*l+1)) )),
     &          pr_crorb(icrorb,ik,1)%mat_rep
     &          (1:2*(2*l+1),nbbot:nbtop))
C conj_mat = inverse(reptrans)*pr_crorb%mat_rep = <lm|new_i>*proj_{new_i} = proj_{lm} [temporarily]
              DO m=1,2*l+1
                mat(m,nbbot:nbtop)=mat(m,nbbot:nbtop)*
     &            SQRT(kp(ik,1)%tetrweight(nbbot:nbtop))
C The first (2*l+1) lines are associated to up states. Hence a multiplication by tetrweight(up)
                mat((2*l+1)+m,nbbot:nbtop)=
     =            mat((2*l+1)+m,nbbot:nbtop)*
     &            SQRT(kp(ik,2)%tetrweight(nbbot:nbtop))
C The last (2*l+1) lines are associated to dn states. Hence a multiplication by tetrweight(dn)
C mat = P(icrorb,ik,is)*sqrt(tetrahedron-weight(ik,is))
                conj_mat(m,nbbot:nbtop)=
     =            CONJG(conj_mat(m,nbbot:nbtop))*
     &            SQRT(kp(ik,1)%tetrweight(nbbot:nbtop))
C The first (2*l+1) lines are associated to up states. Hence a multiplication by tetrweight(up)
                conj_mat((2*l+1)+m,nbbot:nbtop)=
     &            CONJG(conj_mat((2*l+1)+m,nbbot:nbtop))*
     &            SQRT(kp(ik,2)%tetrweight(nbbot:nbtop))
C The last (2*l+1) lines are associated to dn states. Hence a multiplication by tetrweight(dn)
C conj_mat = conjugate[ P(icrorb,ik,is1))*sqrt(tetrahedron-weight(ik,is1)) ]
              ENDDO
              CALL zgemm('N','T',2*(2*l+1),2*(2*l+1),nbnd,
     &          DCMPLX(1.D0,0.D0),mat,2*(2*l+1),conj_mat,2*(2*l+1),
     &          DCMPLX(0.D0,0.D0),D,2*(2*l+1))
              DEALLOCATE(conj_mat,mat)
C D = mat*transpose(conj_mat) is a matrix of size 2*(2*l+1)* 2*(2*l+1)
C
              crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     &          crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))
     &          +D(1:2*(2*l+1),1:2*(2*l+1))
            END DO       ! End of the ik loop 
C The summation over the k-points is done.
C crdensmat(icrorb) is therefore the complete density matrix in spherical harmonic basis.
C
C The density matrix is then put into the desired basis, using reptrans(l,isrt)%transmat 
            crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =        MATMUL(reptrans(l,isrt)%transmat
     &        (1:2*(2*l+1),1:2*(2*l+1)),crdensmat(1,icrorb)%mat
     &        (1:2*(2*l+1),1:2*(2*l+1)) )
            crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =        MATMUL(crdensmat(1,icrorb)%mat
     &        (1:2*(2*l+1),1:2*(2*l+1)),TRANSPOSE( CONJG(
     &        reptrans(l,isrt)%transmat(1:2*(2*l+1),1:2*(2*l+1)))))
C crdensmat = (reptrans)*crdensmat_l*inverse(reptrans) 
C or crdensmat = <new_i|lm> crdensmat_{lm} <lm|new_i>
C crdensmat(icrorb) is now the complete density matrix in the desired basis.
C
C Computation of the density matrix using a simple point integration :
C --------------------------------------------------------------------
           ELSE
C No distinction between up and dn states is necessary because we use a
C geometric factor which depends only of the k-point.
C We can then use directly the projectors in the desired basis.
            DO ik=1,nk
C Only the k-points with inlcuded bands are considered for the projectors.
              IF(.NOT.kp(ik,1)%included) CYCLE
              nbnd=kp(ik,1)%nb_top-kp(ik,1)%nb_bot+1
              nbbot=kp(ik,1)%nb_bot
              nbtop=kp(ik,1)%nb_top
              ALLOCATE(mat(1:2*(2*l+1),nbbot:nbtop))
              mat(1:2*(2*l+1),nbbot:nbtop)=
     =          pr_crorb(icrorb,ik,1)%mat_rep(1:2*(2*l+1),
     &          nbbot:nbtop)
              CALL zgemm('N','C',2*(2*l+1),2*(2*l+1),nbnd,
     &          DCMPLX(1.D0,0.D0),mat,2*(2*l+1),mat,2*(2*l+1),
     &          DCMPLX(0.D0,0.D0),D,2*(2*l+1))
              DEALLOCATE(mat)
C D = P(icrorb,ik,is)*transpose(conjugate(P(icrorb,ik,is))) is a matrix of size 2*(2*l+1) * 2*(2*l+1)
              crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))= 
     =          crdensmat(1,icrorb)%mat(1:2*(2*l+1),1:2*(2*l+1))
     +          +D(1:2*(2*l+1),1:2*(2*l+1))*kp(ik,1)%weight
            ENDDO      ! End of the ik loop
C The summation over the k-points is done in the do loop. 
C The weight used is a geometric factor associated to k and does not depend on the variable is. 
C That's why we merely multiply by the "weight" each term while summing over the k-points.
           ENDIF
           DEALLOCATE(D)
C
C ----------------------------------------------------------------------------------------------
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only) :
C ----------------------------------------------------------------------------------------------
          ELSE
C The field mat of crdensmat must be resized.
C nsp matrices of size (2*l+1) are necessary to represent the whole density matrix.
           DO is=1,nsp
             DEALLOCATE(crdensmat(is,icrorb)%mat)
             ALLOCATE(crdensmat(is,icrorb)%mat(-l:l,-l:l))
             crdensmat(is,icrorb)%mat(-l:l,-l:l)=0.d0
           ENDDO
           ALLOCATE(D(-l:l,-l:l))
C All the computations can be performed in the new basis (using the field pr_crorb%mat_rep)
           DO ik=1,nk
             DO iss=1,nsp
C
C Determination of the block indices :
C ------------------------------------
               IF(iss.LE.2) THEN
                is=iss
                is1=iss 
C If iss=1 (up), is=1 and is1=1 -> Description of the up/up block
C If iss=2 (down), is=2 and is1=2 -> Description of the dn/dn block
               ELSE 
                is=iss-2 
                is1=3-is
C If iss=3, is=1 (up) and is1=2 (down) -> Description of the up/dn block
C If iss=4, is=2 (down) and is1=1 (up) -> Description of the dn/up block
               ENDIF
C Only the k-points with inlcuded bands are considered for the projectors.
               IF(.NOT.kp(ik,is)%included) CYCLE
               IF(.NOT.kp(ik,is1)%included) CYCLE
               nbnd=kp(ik,is)%nb_top-kp(ik,is)%nb_bot+1
               nbbot=kp(ik,is)%nb_bot
               nbtop=kp(ik,is)%nb_top
C for the diagonal blocks, is=is1 ; thus nbtop and nbtop are the same for is and is1. 
C for the off-diagonal blocks (calculated only when SO is considered), 
C it was checked that nbtop(up)=nbtop(dn) and nbbot(up)=nbbot(dn) [in set_projections.f]
C thus nbtop and nbtop fit again for is and is1.
C
C Computation of the density matrix using the tetrahedron weights for the integration :
C -------------------------------------------------------------------------------------
               IF(tetr) THEN
C The representation of the projectors in the desired basis is used (field pr_crorb%mat_rep)
                ALLOCATE(mat(-l:l,nbbot:nbtop))
                ALLOCATE(conj_mat(-l:l,nbbot:nbtop))
                DO m=-l,l
                  mat(m,nbbot:nbtop)=
     =              pr_crorb(icrorb,ik,is)%mat_rep(m,nbbot:nbtop)*
     &              SQRT(kp(ik,is)%tetrweight(nbbot:nbtop))
                  conj_mat(m,nbbot:nbtop)=CONJG(
     &              pr_crorb(icrorb,ik,is1)%mat_rep(m,nbbot:nbtop))
     &              *SQRT(kp(ik,is1)%tetrweight(nbbot:nbtop))
C mat = P(icrorb,ik,is)*sqrt(tetrahedron-weight(ik,is))
C conj_mat = conjugate[ P(icrorb,ik,is1))*sqrt(tetrahedron-weight(ik,is1)) ]
                ENDDO
                CALL zgemm('N','T',(2*l+1),(2*l+1),nbnd,
     &            DCMPLX(1.D0,0.D0),mat,(2*l+1),conj_mat,(2*l+1),
     &            DCMPLX(0.D0,0.D0),D(-l:l,-l:l),(2*l+1))
                DEALLOCATE(conj_mat,mat)
C D = mat*transpose(conj_mat) is a matrix of size (2*l+1)*(2*l+1)
                crdensmat(iss,icrorb)%mat(-l:l,-l:l)=
     =            crdensmat(iss,icrorb)%mat(-l:l,-l:l)+D(-l:l,-l:l)
C The summation over the k-points is done.
C crdensmat(icrorb) is therefore the complete density matrix of the orbital icrorb.
C
C Computation of the density matrix using a simple point integration :
C --------------------------------------------------------------------
               ELSE
                ALLOCATE(mat(-l:l,nbbot:nbtop))
                ALLOCATE(conj_mat(-l:l,nbbot:nbtop))
                mat(-l:l,nbbot:nbtop)=pr_crorb(icrorb,ik,is)
     &            %mat_rep(-l:l,nbbot:nbtop)
                conj_mat(-l:l,nbbot:nbtop)=pr_crorb(icrorb,ik,is1)
     &            %mat_rep(-l:l,nbbot:nbtop)
                CALL zgemm('N','C',(2*l+1),(2*l+1),nbnd,
     &            DCMPLX(1.D0,0.D0),mat,(2*l+1),conj_mat,(2*l+1),
     &            DCMPLX(0.D0,0.D0),D(-l:l,-l:l),(2*l+1))
                DEALLOCATE(mat,conj_mat)
C D = P(icrorb,ik,is)*transpose(conjugate(P(icrorb,ik,is))) is a matrix of size (2*l+1)*(2*l+1)
                crdensmat(iss,icrorb)%mat(-l:l,-l:l)=
     =           crdensmat(iss,icrorb)%mat(-l:l,-l:l)
     +           +D(-l:l,-l:l)*kp(ik,is)%weight
C The weight used is a geometric factor asoociated to k and does not depend on the variable is. 
C That's why we merely multiply by the "weight" each term while summing over the k-points.
               ENDIF
             ENDDO    ! End of the iss loop  
           ENDDO      ! End of the ik loop
           DEALLOCATE(D)
          ENDIF       ! End of the basis if-then-else
C
        ENDDO         ! End of the icrorb loop
C
C
C ===============================
C Symmetrization to the full BZ : 
C ===============================
         CALL symmetrize_mat(crdensmat,crorb,ncrorb)
C
C
C ============================================================================
C Application of the Rloc transformation to go back to the local coordinates :
C ============================================================================
         CALL rotdens_mat(crdensmat,crorb,ncrorb)
C
C
C =================================================================
C Printing the density matrices and the charge in the output file :
C =================================================================
        IF(ifprnt) THEN
         CALL printout(0)
         WRITE(buf,'(a)') '-------------------------------------'
         CALL printout(0)
         WRITE(buf,'(a)') 
     &     'Density Matrices for the Correlated States : '
         CALL printout(0)
C The density matrices and charge are printed for all the corrrelated orbitals
         DO icrorb=1,ncrorb
           CALL printout(0)
           isrt=crorb(icrorb)%sort
           l=crorb(icrorb)%l
C Description of the correlated orbital
           WRITE(buf,'(3(a,i3))')
     &       '  Sort = ',isrt,'  Atom = ',crorb(icrorb)%atom,
     &       '  and Orbital l = ',l
           CALL printout(0)
           q=0d0
C
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
           IF (l==0) THEN
C For a calculation spin-polarized with SO :
            IF (nsp==4) THEN
             WRITE(buf,'(2(2F12.6,2x))')
     &         crdensmat(1,icrorb)%mat(1,1),
     &         crdensmat(3,icrorb)%mat(1,1)
             CALL printout(0)
             WRITE(buf,'(2(2F12.6,2x))')
     &         crdensmat(4,icrorb)%mat(1,1),
     &         crdensmat(2,icrorb)%mat(1,1)
             CALL printout(0)
             q=q+crdensmat(1,icrorb)%mat(1,1)+
     +         crdensmat(2,icrorb)%mat(1,1)
C For a calculation spin-polarized without SO :
            ELSE IF (nsp==2) THEN
             WRITE(buf,'(2(2F12.6,2x))')
     &         crdensmat(1,icrorb)%mat(1,1),DCMPLX(0.D0,0.D0)
             CALL printout(0)
             WRITE(buf,'(2(2F12.6,2x))')
     &         DCMPLX(0.D0,0.D0),crdensmat(2,icrorb)%mat(1,1)
             CALL printout(0)
             q=q+crdensmat(1,icrorb)%mat(1,1)+
     +         crdensmat(2,icrorb)%mat(1,1)
C For a paramagnetic calculation without SO :
            ELSE 
             WRITE(buf,'(2F12.6,2x)') crdensmat(1,icrorb)%mat(1,1)
             CALL printout(0)
             q=q+crdensmat(1,icrorb)%mat(1,1)
            ENDIF
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
           ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the calculation is necessary spin-polarized with SO :
            WRITE(buf,'(a,a)') 'Writing the matrix as : ',
     &        '[ block 1 | block 2 ] with'
            CALL printout(0)
            WRITE(buf,'(a,a)') '                        ',
     &        '[ block 3 | block 4 ]'
            CALL printout(0)
C Printing the different blocks
            ALLOCATE(D(1,1:2*l+1))
            WRITE(buf,'(a,i2,a)') '   # For the block 1 :'
            CALL printout(0)
            DO m=1,2*l+1
              D(1,1:2*l+1)=crdensmat(1,icrorb)%mat(m,1:(2*l+1))
              WRITE(buf,'(7(2F12.6),x)') D(:,:)          
              CALL printout(0)
            ENDDO
            WRITE(buf,'(a,i2,a)') '   # For the block 2 :'
            CALL printout(0)
            DO m=1,2*l+1
              D(1,1:2*l+1)=
     &          crdensmat(1,icrorb)%mat(m,2*l+2:2*(2*l+1))
              WRITE(buf,'(7(2F12.6),x)') D(:,:) 
              CALL printout(0)
            ENDDO
            WRITE(buf,'(a,i2,a)') '   # For the block 3 :'
            CALL printout(0)
            DO m=2*l+2,2*(2*l+1)
              D(1,1:2*l+1)=
     &          crdensmat(1,icrorb)%mat(m,1:(2*l+1))
              WRITE(buf,'(7(2F12.6),x)') D(:,:) 
              CALL printout(0)
            ENDDO
            WRITE(buf,'(a,i2,a)') '   # For the block 4 :'
            CALL printout(0)
            DO m=2*l+2,2*(2*l+1)
              D(1,1:2*l+1)=
     &          crdensmat(1,icrorb)%mat(m,2*l+2:2*(2*l+1))
              WRITE(buf,'(7(2F12.6),x)') D(:,:)
              CALL printout(0)
            ENDDO
            DEALLOCATE(D)
            DO m1=1,2*(2*l+1)
              q=q+crdensmat(1,icrorb)%mat(m1,m1)
            ENDDO
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
           ELSE
            IF(nsp==4) THEN 
             WRITE(buf,'(a,a)') 'Writing the matrix as : ',
     &         '[ block up/up | block up/dn ] with'
             CALL printout(0)
             WRITE(buf,'(a,a)') '                        ',
     &         '[ block dn/up | block dn/dn ]'
             CALL printout(0)
            ELSEIF(nsp==2) THEN
             WRITE(buf,'(a,a)') 'Writing the matrix as : ',
     &         '[ block up/up |      0      ] with'
             CALL printout(0)
             WRITE(buf,'(a,a)') '                        ',
     &         '[      0      | block dn/dn ]'
             CALL printout(0)
            ENDIF
            DO iss=1,nsp
              IF(iss.LE.2) THEN
               is=iss
               is1=iss 
C If iss=1 (up), is=1 and is1=1 -> Description of the up/up block
C If iss=2 (down), is=2 and is1=2 -> Description of the down/down block
               IF (is==1) THEN
                IF ( ifSP.or.ifSO ) THEN
                 WRITE(buf,'(a)') '   # For the Up/Up block     :'  
                ENDIF
                CALL printout(0)
               ELSE
                WRITE(buf,'(a)') '   # For the Down/Down block :' 
                CALL printout(0)
               ENDIF
              ELSE 
               is=iss-2
               is1=3-is
C If iss=3, is=1 (up) and is1=2 (down) -> Description of the up/down block
C If iss=4, is=2 (down) and is1=1 (up) -> Description of the down/up block
               IF (is==1) THEN
                WRITE(buf,'(a)') '   # For the Up/Down block   :'  
                CALL printout(0)
               ELSE
                WRITE(buf,'(a)') '   # For the Down/Up block   :'  
                CALL printout(0)
               ENDIF
              ENDIF
              ALLOCATE(D(1,-l:l))
              DO m=-l,l
                D(1,-l:l)=crdensmat(iss,icrorb)%mat(m,-l:l)
                WRITE(buf,'(7(2F12.6),x)') D(:,:)          
                CALL printout(0)
                IF(is==is1) q=q+crdensmat(iss,icrorb)%mat(m,m)
              ENDDO
              DEALLOCATE(D)
            ENDDO
           ENDIF
C Displaying the charge q of the orbital
           CALL printout(0)
           WRITE(buf,'(a,f10.5)')'The charge of the orbital is : ',q
           CALL printout(0)
         ENDDO
        ENDIF
C
        ENDIF  ! End of the if ncrorb=0 if-then-else
C The calculation is stopped here if the flag only_corr is .TRUE.
        IF (.not.only_corr) THEN
C
C =========================================
C Computation of the total charge density :
C =========================================
C The charge is stored in the variable qtot given in argument.  
        qtot=0d0
        do is=1,ns
          DO ik=1,nk
            if (kp(ik,is)%included) then 
              DO ib=kp(ik,is)%nb_bot,kp(ik,is)%nb_top
                IF(tetr) THEN
                  qtot=qtot+kp(ik,is)%tetrweight(ib)
                ELSE
                  qtot=qtot+kp(ik,is)%weight
                ENDIF
              ENDDO
            endif
          ENDDO
          if (ifSO) exit
        enddo
C
C =================================================================
C Computation of the density matrix for all the included orbitals :
C =================================================================
C The computations are performed with the Theta projectors (pr_orb)
C
C densmat is a table of size nsp*norb. 
C For each included orbital iorb, densmat(:,iorb) is the corresponding density matrix.
        IF(.NOT.ALLOCATED(densmat)) THEN
         ALLOCATE(densmat(nsp,norb))
        ENDIF
        DO iorb=1,norb
          DO is=1,nsp
            IF(ALLOCATED(densmat(is,iorb)%mat)) 
     &       DEALLOCATE(densmat(is,iorb)%mat)
            ALLOCATE(densmat(is,iorb)%mat(1,1))
            densmat(is,iorb)%mat(1,1)=0.d0
          ENDDO
        ENDDO
C
C Loop on the included orbitals iorb
C
        DO iorb=1,norb
          isrt=orb(iorb)%sort
          l=orb(iorb)%l
C
C -----------------------------------------------------------------------------------
C The s-orbitals are a particular case of a "non-mixing" basis and are treated here :
C -----------------------------------------------------------------------------------
          IF (l==0) THEN
C The field mat of densmat has already the good size (1 scalar element). 
C There's no need of a basis change since the representation of an s-orbital is Identity whatever the basis is.
           DO ik=1,nk
             DO iss=1,nsp
C
C Determination of the block indices :
C ------------------------------------
               IF(iss.LE.2) THEN
                is=iss
                is1=iss 
C If iss=1 (up), is=1 and is1=1 -> Description of the up/up block
C If iss=2 (down), is=2 and is1=2 -> Description of the dn/dn block
               ELSE 
                is=iss-2 
                is1=3-is
C If iss=3, is=1 (up) and is1=2 (down) -> Description of the up/dn block
C If iss=4, is=2 (down) and is1=1 (up) -> Description of the dn/up block
               ENDIF
C Only the k-points with inlcuded bands are considered for the projectors.
               IF(.NOT.kp(ik,is)%included) CYCLE
               IF(.NOT.kp(ik,is1)%included) CYCLE
               nbnd=kp(ik,is)%nb_top-kp(ik,is)%nb_bot+1
               nbbot=kp(ik,is)%nb_bot
               nbtop=kp(ik,is)%nb_top
C for the diagonal blocks, is=is1 ; thus nbtop and nbtop are the same for is and is1. 
C for the off-diagonal blocks (calculated only when SO is considered), 
C it was checked that nbtop(up)=nbtop(dn) and nbbot(up)=nbbot(dn) [in set_projections.f]
C thus nbtop and nbtop fit again for is and is1.
C
C Computation of the density matrix using the tetrahedron weights for the integration :
C -------------------------------------------------------------------------------------
               IF(tetr) THEN
C The field pr_orb%matn_rep is used to perform the computation (well defined for s-orbitals)  
                DO i=1,norm_radf(iorb)%n
                  ALLOCATE(mat(1,nbbot:nbtop))
                  ALLOCATE(conj_mat(1,nbbot:nbtop))
                  mat(1,nbbot:nbtop)=
     &              pr_orb(iorb,ik,is)%matn_rep(1,nbbot:nbtop,i)*
     &              SQRT(kp(ik,is)%tetrweight(nbbot:nbtop))
                  conj_mat(1,nbbot:nbtop)=CONJG(
     &              pr_orb(iorb,ik,is1)%matn_rep(1,nbbot:nbtop,i))*
     &              SQRT(kp(ik,is1)%tetrweight(nbbot:nbtop))
C mat = Theta(iorb,ik,is)*sqrt(tetrahedron-weight(ik,is))
C conj_mat = conjugate[ Theta(iorb,ik,is1))*sqrt(tetrahedron-weight(ik,is1)) ]
                  DO ib = nbbot,nbtop
                    densmat(iss,iorb)%mat(1,1)=
     =                densmat(iss,iorb)%mat(1,1)
     &                +mat(1,ib)*conj_mat(1,ib)
                  ENDDO
C densmat = mat*transpose(conj_mat) which is a matrix of size 1
                  DEALLOCATE(conj_mat,mat)
                ENDDO
C The summation over the k-points is done with the do loop ; 
C The summation over the |phi_j> basis is done with the do loop ; 
C densmat(iss,iorb) is therefore the block "iss" of the density matrix for the orbital iorb.
C
C Computation of the density matrix using a simple point integration :
C --------------------------------------------------------------------
               ELSE
                DO i=1,norm_radf(iorb)%n
                  DO ib = nbbot,nbtop
                    densmat(iss,iorb)%mat(1,1)=
     =                densmat(iss,iorb)%mat(1,1)+
     +                pr_orb(iorb,ik,is)%matn_rep(1,ib,i)*
     *                CONJG(pr_orb(iorb,ik,is1)%matn_rep(1,ib,i))*
     *                kp(ik,is)%weight
                  ENDDO
                ENDDO
C densmat = Theta(iorb,ik,is)*transpose(conjugate(Theta(iorb,ik,is1))) which is a matrix of size 1
C The weight used is a geometric factor associated to k and does not depend on the variable is. 
C That's why we merely multiply by the "weight" each term while summing over the k-points.
               ENDIF
             ENDDO    ! End of the iss loop  
           ENDDO      ! End of the ik loop
C
C -----------------------------------------------------------------------------------------------------
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) ) :
C -----------------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C If this option is used, then ifSO=.TRUE. (because of the restriction in set_ang_trans.f)
C Moreover ifSP=.TRUE. (since ifSO => ifSP, in this version) 
C As a result, we know that nb_bot(up)=nb_bot(dn) and nb_top(up)=nb_top(dn)
C
C The field mat of densmat must be resized.
C As the complete spinor rotation approach is used, only one matrix is necessary (with is=1).
           DEALLOCATE(densmat(1,iorb)%mat)
           ALLOCATE(densmat(1,iorb)%mat(1:2*(2*l+1),1:2*(2*l+1)))
           densmat(1,iorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=0.d0
           ALLOCATE(D(1:2*(2*l+1),1:2*(2*l+1)))
C
C Computation of the density matrix using the tetrahedron weights for the integration :
C -------------------------------------------------------------------------------------
           IF(tetr) THEN
            DO ik=1,nk
C Only the k-points with inlcuded bands are considered for the projectors.
              IF(.NOT.kp(ik,1)%included) CYCLE
              nbnd=kp(ik,1)%nb_top-kp(ik,1)%nb_bot+1
              nbbot=kp(ik,1)%nb_bot
              nbtop=kp(ik,1)%nb_top
C A distinction between up and dn states is necessary in order to use the tetrahedron weight.
C As a result, we will transform the projectors back in spherical harmonics basis 
C to perform the calculation and then put the resulting density matrix in the desired basis.
C
C Loop on the |phi_j> basis
              DO i=1,norm_radf(iorb)%n
                ALLOCATE(mat(1:2*(2*l+1),nbbot:nbtop))
                ALLOCATE(conj_mat(1:2*(2*l+1),nbbot:nbtop))
C The representation of the projectors is put back in the spherical harmonics basis
                mat(1:2*(2*l+1),nbbot:nbtop)=MATMUL( 
     &            TRANSPOSE(CONJG(reptrans(l,isrt)%
     &            transmat(1:2*(2*l+1),1:2*(2*l+1)) )),
     &            pr_orb(iorb,ik,1)%matn_rep
     &            (1:2*(2*l+1),nbbot:nbtop,i) )
C mat = inverse(reptrans)*pr_orb%mat_rep = <lm|new_i>*theta_{new_i} = theta_{lm} [temporarily]
                conj_mat(1:2*(2*l+1),nbbot:nbtop)=MATMUL( 
     &            TRANSPOSE(CONJG(reptrans(l,isrt)%
     &            transmat(1:2*(2*l+1),1:2*(2*l+1)) )),
     &            pr_orb(iorb,ik,1)%matn_rep
     &            (1:2*(2*l+1),nbbot:nbtop,i) )
C conj_mat = inverse(reptrans)*pr_orb%mat_rep = <lm|new_i>*theta_{new_i} = theta_{lm} [temporarily]
                DO m=1,2*l+1
                  mat(m,nbbot:nbtop)=mat(m,nbbot:nbtop)*
     &              SQRT(kp(ik,1)%tetrweight(nbbot:nbtop))
C The first (2*l+1) lines are associated to up states. Hence a multiplication by tetrweight(up)
                  mat((2*l+1)+m,nbbot:nbtop)=
     =              mat((2*l+1)+m,nbbot:nbtop)*
     &              SQRT(kp(ik,2)%tetrweight(nbbot:nbtop))
C The last (2*l+1) lines are associated to dn states. Hence a multiplication by tetrweight(dn)
C mat = Theta(icrorb,ik,is)*sqrt(tetrahedron-weight(ik,is))
                  conj_mat(m,nbbot:nbtop)=CONJG( 
     &              conj_mat(m,nbbot:nbtop) )*
     &              SQRT(kp(ik,1)%tetrweight(nbbot:nbtop))
C The first (2*l+1) lines are associated to up states. Hence a multiplication by tetrweight(up)
                  conj_mat((2*l+1)+m,nbbot:nbtop)=
     &              CONJG(conj_mat((2*l+1)+m,nbbot:nbtop))*
     &              SQRT(kp(ik,2)%tetrweight(nbbot:nbtop))
C The last (2*l+1) lines are associated to dn states. Hence a multiplication by tetrweight(dn)
C conj_mat = conjugate[ Theta(icrorb,ik,is1))*sqrt(tetrahedron-weight(ik,is1)) ]
                ENDDO
                CALL zgemm('N','T',2*(2*l+1),2*(2*l+1),nbnd,
     &            DCMPLX(1.D0,0.D0),mat,2*(2*l+1),conj_mat,2*(2*l+1),
     &            DCMPLX(0.D0,0.D0),D,2*(2*l+1))
                DEALLOCATE(conj_mat,mat)
C D = mat*transpose(conj_mat) is a matrix of size 2*(2*l+1)* 2*(2*l+1)
C
                densmat(1,iorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     &            densmat(1,iorb)%mat(1:2*(2*l+1),1:2*(2*l+1))
     &            +D(1:2*(2*l+1),1:2*(2*l+1))
              END DO     ! End of the |phi_j> basis loop
            END DO       ! End of the ik loop 
C The summation over the k-points is done ;
C The summation over the |phi_j> basis is done ; 
C crdensmat(icrorb) is therefore the complete density matrix in spherical harmonic basis.
C
C The density matrix is then put into the desired basis, using reptrans(l,isrt)%transmat 
            densmat(1,iorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =        MATMUL(reptrans(l,isrt)%transmat
     &        (1:2*(2*l+1),1:2*(2*l+1)),densmat(1,iorb)%mat
     &        (1:2*(2*l+1),1:2*(2*l+1)) )
            densmat(1,iorb)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =        MATMUL(densmat(1,iorb)%mat
     &        (1:2*(2*l+1),1:2*(2*l+1)),TRANSPOSE( CONJG(
     &        reptrans(l,isrt)%transmat(1:2*(2*l+1),1:2*(2*l+1)))))
C densmat = (reptrans)*densmat_l*inverse(reptrans) 
C or densmat = <new_i|lm> densmat_{lm} <lm|new_i>
C densmat(iorb) is now the complete density matrix in the desired basis.
C
C Computation of the density matrix using a simple point integration :
C --------------------------------------------------------------------
               ELSE
C No distinction between up and dn states is necessary because we use a
C geometric factor which depends only of the k-point.
C We can then use directly the projectors in the desired basis.
            DO ik=1,nk
C Only the k-points with inlcuded bands are considered for the projectors.
              IF(.NOT.kp(ik,1)%included) CYCLE
C Loop on the |phi_j> basis
              DO i=1,norm_radf(iorb)%n
                nbnd=kp(ik,1)%nb_top-kp(ik,1)%nb_bot+1
                nbbot=kp(ik,1)%nb_bot
                nbtop=kp(ik,1)%nb_top
                ALLOCATE(mat(1:2*(2*l+1),nbbot:nbtop))
                mat(1:2*(2*l+1),nbbot:nbtop)=
     =            pr_orb(iorb,ik,1)%matn_rep
     &            (1:2*(2*l+1),nbbot:nbtop,i)
                CALL zgemm('N','C',2*(2*l+1),2*(2*l+1),nbnd,
     &            DCMPLX(1.D0,0.D0),mat,2*(2*l+1),mat,
     &            2*(2*l+1),DCMPLX(0.D0,0.D0),D,2*(2*l+1))
                DEALLOCATE(mat)
C D = Theta(iorb,ik,is)*transpose(conjugate(Theta(iorb,ik,is))) is a matrix of size 2*(2*l+1) * 2*(2*l+1)
                densmat(1,iorb)%mat(1:2*(2*l+1),1:2*(2*l+1))= 
     =            densmat(1,iorb)%mat(1:2*(2*l+1),1:2*(2*l+1))
     +            +D(1:2*(2*l+1),1:2*(2*l+1))*kp(ik,1)%weight
              END DO   ! End of the |phi_j > basis loop
            ENDDO      ! End of the ik loop
C The summation over the k-points is done ;
C The summation over the |phi_j> basis is done ; 
C The weight used is a geometric factor associated to k and does not depend on the variable is. 
C That's why we merely multiply by the "weight" each term while summing over the k-points.
           ENDIF
           DEALLOCATE(D)
C
C ----------------------------------------------------------------------------------------------
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only) :
C ----------------------------------------------------------------------------------------------
          ELSE
C The field mat of densmat must be resized.
C nsp matrices of size (2*l+1) are necessary to represent the whole density matrix.
           DO is=1,nsp
             DEALLOCATE(densmat(is,iorb)%mat)
             ALLOCATE(densmat(is,iorb)%mat(-l:l,-l:l))
             densmat(is,iorb)%mat(-l:l,-l:l)=0.d0
           ENDDO
           ALLOCATE(D(-l:l,-l:l))
C All the computations can be performed in the new basis (using the field pr_orb%matn_rep)
           DO ik=1,nk
             DO iss=1,nsp
C
C Determination of the block indices :
C ------------------------------------
               IF(iss.LE.2) THEN
                is=iss
                is1=iss 
C If iss=1 (up), is=1 and is1=1 -> Description of the up/up block
C If iss=2 (down), is=2 and is1=2 -> Description of the dn/dn block
               ELSE 
                is=iss-2 
                is1=3-is
C If iss=3, is=1 (up) and is1=2 (down) -> Description of the up/dn block
C If iss=4, is=2 (down) and is1=1 (up) -> Description of the dn/up block
               ENDIF
C Only the k-points with inlcuded bands are considered for the projectors.
               IF(.NOT.kp(ik,is)%included) CYCLE
               IF(.NOT.kp(ik,is1)%included) CYCLE
               nbnd=kp(ik,is)%nb_top-kp(ik,is)%nb_bot+1
               nbbot=kp(ik,is)%nb_bot
               nbtop=kp(ik,is)%nb_top
C for the diagonal blocks, is=is1 ; thus nbtop and nbtop are the same for is and is1. 
C for the off-diagonal blocks (calculated only when SO is considered), 
C it was checked that nbtop(up)=nbtop(dn) and nbbot(up)=nbbot(dn) [in set_projections.f]
C thus nbtop and nbtop fit again for is and is1.
C
C Computation of the density matrix using the tetrahedron weights for the integration :
C -------------------------------------------------------------------------------------
               IF(tetr) THEN
C The representation of the projectors in the desired basis is used (field pr_orb%mat_rep)
                DO i=1,norm_radf(iorb)%n
                  ALLOCATE(mat(-l:l,nbbot:nbtop))
                  ALLOCATE(conj_mat(-l:l,nbbot:nbtop))
                  DO m=-l,l
                    mat(m,nbbot:nbtop)=pr_orb(iorb,ik,is)%matn_rep
     &                (m,nbbot:nbtop,i)*SQRT(kp(ik,is)%
     &                tetrweight(nbbot:nbtop))
                    conj_mat(m,nbbot:nbtop)=CONJG(
     =                pr_orb(iorb,ik,is1)%matn_rep(m,nbbot:nbtop,i))
     &                *SQRT(kp(ik,is1)%tetrweight(nbbot:nbtop))
C mat = Theta(iorb,ik,is)*sqrt(tetrahedron-weight(ik,is))
C conj_mat = conjugate[ Theta(iorb,ik,is1))*sqrt(tetrahedron-weight(ik,is1)) ]
                  ENDDO
                  CALL zgemm('N','T',(2*l+1),(2*l+1),nbnd,
     &              DCMPLX(1.D0,0.D0),mat,(2*l+1),conj_mat,(2*l+1),
     &              DCMPLX(0.D0,0.D0),D(-l:l,-l:l),(2*l+1))
                  DEALLOCATE(conj_mat,mat)
C D = mat*transpose(conj_mat) is a matrix of size (2*l+1)*(2*l+1)
                  densmat(iss,iorb)%mat(-l:l,-l:l)=
     =              densmat(iss,iorb)%mat(-l:l,-l:l)+D(-l:l,-l:l)
                ENDDO
C The summation over the |phi_j > basis is done.
C The summation over the k-points is done.
C densmat(iorb) is therefore the complete density matrix of the orbital iorb.
C
C Computation of the density matrix using a simple point integration :
C --------------------------------------------------------------------
               ELSE
                DO i=1,norm_radf(iorb)%n
                  ALLOCATE(mat(-l:l,nbbot:nbtop))
                  ALLOCATE(conj_mat(-l:l,nbbot:nbtop))
                  mat(-l:l,nbbot:nbtop)=
     =              pr_orb(iorb,ik,is)%matn_rep(-l:l,nbbot:nbtop,i)
                  conj_mat(-l:l,nbbot:nbtop)=
     =              pr_orb(iorb,ik,is1)%matn_rep(-l:l,nbbot:nbtop,i)
                  CALL zgemm('N','C',(2*l+1),(2*l+1),nbnd,
     &              DCMPLX(1.D0,0.D0),mat,(2*l+1),conj_mat,(2*l+1),
     &              DCMPLX(0.D0,0.D0),D(-l:l,-l:l),(2*l+1))
                  DEALLOCATE(mat,conj_mat)
C D = Theta(iorb,ik,is)*transpose(conjugate(Theta(iorb,ik,is))) is a matrix of size (2*l+1)*(2*l+1)
                  densmat(iss,iorb)%mat(-l:l,-l:l)=
     =              densmat(iss,iorb)%mat(-l:l,-l:l)
     +              +D(-l:l,-l:l)*kp(ik,is)%weight
C The weight used is a geometric factor asoociated to k and does not depend on the variable is. 
C That's why we merely multiply by the "weight" each term while summing over the k-points and the |phi_j> basis.
                ENDDO ! End of the |phi_j > basis loop
               ENDIF
             ENDDO    ! End of the iss loop  
           ENDDO      ! End of the ik loop
           DEALLOCATE(D)
C
          ENDIF       ! End of the basis if-then-else
C
        ENDDO         ! End of the iorb loop
C 
C
C ===============================
C Symmetrization to the full BZ :
C ===============================
        CALL symmetrize_mat(densmat,orb,norb)
C
C
C ============================================================================
C Application of the Rloc transformation to go back to the local coordinates :
C ============================================================================
       CALL rotdens_mat(densmat,orb,norb)
C
C
C =================================================================
C Printing the density matrices and the charge in the output file :
C =================================================================
        IF(ifprnt) THEN
         CALL printout(0)
         WRITE(buf,'(a)') '-------------------------------------'
         CALL printout(0)
         WRITE(buf,'(a)') 
     &     'Density Matrices for all the States of the System : '
         CALL printout(0)
C The density matrices and charge are printed for all the included orbitals
         DO iorb=1,norb
           CALL printout(0)
           isrt=orb(iorb)%sort
           l= orb(iorb)%l
C Description of the correlated orbital
           WRITE(buf,'(3(a,i3))')
     &       '  Sort = ',isrt,'  Atom = ',orb(iorb)%atom,
     &       '  and Orbital l = ',l
           CALL printout(0)
           q=0d0
C
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
           IF (l==0) THEN
C For a calculation spin-polarized with SO :
            IF (nsp==4) THEN
             WRITE(buf,'(2(2F12.6,2x))')
     &         densmat(1,iorb)%mat(1,1),
     &         densmat(3,iorb)%mat(1,1)
             CALL printout(0)
             WRITE(buf,'(2(2F12.6,2x))')
     &         densmat(4,iorb)%mat(1,1),
     &         densmat(2,iorb)%mat(1,1)
             CALL printout(0)
             q=q+densmat(1,iorb)%mat(1,1)+
     +         densmat(2,iorb)%mat(1,1)
C For a calculation spin-polarized without SO :
            ELSE IF (nsp==2) THEN
             WRITE(buf,'(2(2F12.6,2x))')
     &         densmat(1,iorb)%mat(1,1),DCMPLX(0.D0,0.D0)
             CALL printout(0)
             WRITE(buf,'(2(2F12.6,2x))')
     &         DCMPLX(0.D0,0.D0),densmat(2,iorb)%mat(1,1)
             CALL printout(0)
             q=q+densmat(1,iorb)%mat(1,1)+
     +         densmat(2,iorb)%mat(1,1)
C For a paramagnetic calculation without SO :
            ELSE 
             WRITE(buf,'(2F12.6,2x)') densmat(1,iorb)%mat(1,1)
             q=q+densmat(1,iorb)%mat(1,1)
             CALL printout(0)
            ENDIF
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
           ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the calculation is necessary spin-polarized with SO :
            WRITE(buf,'(a,a)') 'Writing the matrix as : ',
     &        '[ block 1 | block 2 ] with'
            CALL printout(0)
            WRITE(buf,'(a,a)') '                        ',
     &        '[ block 3 | block 4 ]'
            CALL printout(0)
C Printing the different blocks
            ALLOCATE(D(1,1:2*l+1))
            WRITE(buf,'(a,i2,a)') '   # For the block 1 :'
            CALL printout(0)
            DO m=1,2*l+1
              D(1,1:2*l+1)=densmat(1,iorb)%mat(m,1:(2*l+1))
              WRITE(buf,'(7(2F12.6),x)') D(:,:)
              CALL printout(0)
            ENDDO
            WRITE(buf,'(a,i2,a)') '   # For the block 2 :'
            CALL printout(0)
            DO m=1,2*l+1
              D(1,1:2*l+1)=
     &          densmat(1,iorb)%mat(m,2*l+2:2*(2*l+1))
              WRITE(buf,'(7(2F12.6),x)') D(:,:)
              CALL printout(0)
            ENDDO
            WRITE(buf,'(a,i2,a)') '   # For the block 3 :'
            CALL printout(0)
            DO m=2*l+2,2*(2*l+1)
              D(1,1:2*l+1)=
     &          densmat(1,iorb)%mat(m,1:(2*l+1))
              WRITE(buf,'(7(2F12.6),x)') D(:,:)
              CALL printout(0)
            ENDDO
            WRITE(buf,'(a,i2,a)') '   # For the block 4 :'
            CALL printout(0)
            DO m=2*l+2,2*(2*l+1)
              D(1,1:2*l+1)=
     &          densmat(1,iorb)%mat(m,2*l+2:2*(2*l+1))
              WRITE(buf,'(7(2F12.6),x)') D(:,:)
              CALL printout(0)
            ENDDO
            DEALLOCATE(D)
            DO m1=1,2*(2*l+1)
              q=q+densmat(1,iorb)%mat(m1,m1)
            ENDDO
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
           ELSE
            IF(nsp==4) THEN 
             WRITE(buf,'(a,a)') 'Writing the matrix as : ',
     &         '[ block up/up | block up/dn ] with'
             CALL printout(0)
             WRITE(buf,'(a,a)') '                        ',
     &         '[ block dn/up | block dn/dn ]'
             CALL printout(0)
            ELSEIF(nsp==2) THEN
             WRITE(buf,'(a,a)') 'Writing the matrix as : ',
     &         '[ block up/up |      0      ] with'
             CALL printout(0)
             WRITE(buf,'(a,a)') '                        ',
     &         '[      0      | block dn/dn ]'
             CALL printout(0)
            ENDIF
            DO iss=1,nsp
              IF(iss.LE.2) THEN
               is=iss
               is1=iss 
C If iss=1 (up), is=1 and is1=1 -> Description of the up/up block
C If iss=2 (down), is=2 and is1=2 -> Description of the down/down block
               IF (is==1) THEN
                IF ( ifSP.or.ifSO ) THEN
                 WRITE(buf,'(a)') '   # For the Up/Up block     :'  
                ENDIF
                CALL printout(0)
               ELSE
                WRITE(buf,'(a)') '   # For the Down/Down block :' 
                CALL printout(0)
               ENDIF
              ELSE 
               is=iss-2
               is1=3-is
C If iss=3, is=1 (up) and is1=2 (down) -> Description of the up/down block
C If iss=4, is=2 (down) and is1=1 (up) -> Description of the down/up block
               IF (is==1) THEN
                WRITE(buf,'(a)') '   # For the Up/Down block   :'  
                CALL printout(0)
               ELSE
                WRITE(buf,'(a)') '   # For the Down/Up block   :'  
                CALL printout(0)
               ENDIF
              ENDIF
              ALLOCATE(D(1,-l:l))
              DO m=-l,l
                D(1,-l:l)=densmat(iss,iorb)%mat(m,-l:l)
                WRITE(buf,'(7(2F12.6),x)') D(:,:)          
                CALL printout(0)
                IF(is==is1) q=q+densmat(iss,iorb)%mat(m,m)
              ENDDO
              DEALLOCATE(D)
            ENDDO
           ENDIF
C Displaying the charge q of the orbital
           CALL printout(0)
           WRITE(buf,'(a,f10.5)')'The charge of the orbital is : ',q
           CALL printout(0)
         ENDDO
        ENDIF 
C
C ==========================================================================
C Printing the total charge in the output file (only if only_corr =.FALSE.):
C ==========================================================================
      WRITE(buf,'(a,f11.5)')'TOTAL CHARGE = ',qtot
      CALL printout(1)
C 
      ENDIF    ! End of the .not.only_corr if-then-else
      RETURN
      END
       

          
          
             
      
