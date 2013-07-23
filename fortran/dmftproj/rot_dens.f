
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

      SUBROUTINE rotdens_mat(Dmat,orbit,norbit)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine applies to each density matrix in Dmat          %%
C %% the transformation to go from the global coordinates to the     %%
C %% local coordinates associated to the considered orbital.         %%
C %%                                                                 %%
C %% This version can be used for SO computations.                   %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definition of the variables :
C ----------------------------
      USE common_data
      USE projections
      USE symm
      USE reps
      IMPLICIT NONE
      INTEGER :: norbit
      TYPE(matrix), DIMENSION(nsp,norbit) :: Dmat
      COMPLEX(KIND=8),DIMENSION(:,:), ALLOCATABLE :: rot_dmat
      COMPLEX(KIND=8),DIMENSION(:,:), ALLOCATABLE :: tmp_mat
      COMPLEX(KIND=8):: ephase
      REAL(KIND=8):: factor
      TYPE(orbital), DIMENSION(norbit) :: orbit
      INTEGER :: iatom, isrt, iorb, is, is1, l, i, m
C
C
        DO iorb=1,norbit
          l=orbit(iorb)%l
          isrt=orbit(iorb)%sort
          iatom=orbit(iorb)%atom
C
          IF(ifSP.AND.ifSO) THEN
C In this case, the complete spinor rotation approach (matrices of size 2*(2*l+1) ) is used for rotloc.
           IF (l==0) THEN
C ------------------------------------------------------------------------------------------------------------
C For the s orbital, the spinor rotation matrix will be constructed directly from the Euler angles a,b and c :
C ------------------------------------------------------------------------------------------------------------
C Up/dn and Dn/up terms
            ALLOCATE(tmp_mat(1:2,1:2))
            ALLOCATE(rot_dmat(1:2,1:2))
            IF (rotloc(iatom)%timeinv) THEN
             factor=(rotloc(iatom)%a+rotloc(iatom)%g)/2.d0
             tmp_mat(2,1)=EXP(CMPLX(0.d0,factor))*
     &         DCOS(rotloc(iatom)%b/2.d0)
             tmp_mat(1,2)=-CONJG(tmp_mat(2,1))
C Up/dn and Dn/up terms
             factor=-(rotloc(iatom)%a-rotloc(iatom)%g)/2.d0
             tmp_mat(2,2)=-EXP(CMPLX(0.d0,factor))*
     &         DSIN(rotloc(iatom)%b/2.d0)
             tmp_mat(1,1)=CONJG(tmp_mat(2,2))
C definition of the total density matrix
             rot_dmat(1,1)=Dmat(1,iorb)%mat(1,1)
             rot_dmat(2,2)=Dmat(2,iorb)%mat(1,1)
             rot_dmat(1,2)=Dmat(3,iorb)%mat(1,1)
             rot_dmat(2,1)=Dmat(4,iorb)%mat(1,1)
C going to the local basis
             rot_dmat(1:2,1:2)=CONJG(MATMUl(
     &         rot_dmat(1:2,1:2),tmp_mat(1:2,1:2)))
             rot_dmat(1:2,1:2)=MATMUl(
     &         TRANSPOSE(tmp_mat(1:2,1:2)),
     &         rot_dmat(1:2,1:2))
            ELSE
             factor=(rotloc(iatom)%a+rotloc(iatom)%g)/2.d0
             tmp_mat(1,1)=EXP(CMPLX(0.d0,factor))*
     &         DCOS(rotloc(iatom)%b/2.d0)
             tmp_mat(2,2)=CONJG(tmp_mat(1,1))
C Up/dn and Dn/up terms
             factor=-(rotloc(iatom)%a-rotloc(iatom)%g)/2.d0
             tmp_mat(1,2)=EXP(CMPLX(0.d0,factor))*
     &         DSIN(rotloc(iatom)%b/2.d0)
             tmp_mat(2,1)=-CONJG(tmp_mat(1,2))
C definition of the total density matrix
             rot_dmat(1,1)=Dmat(1,iorb)%mat(1,1)
             rot_dmat(2,2)=Dmat(2,iorb)%mat(1,1)
             rot_dmat(1,2)=Dmat(3,iorb)%mat(1,1)
             rot_dmat(2,1)=Dmat(4,iorb)%mat(1,1)
C going to the local basis
             rot_dmat(1:2,1:2)=MATMUl(
     &         TRANSPOSE(CONJG(tmp_mat(1:2,1:2))),
     &         rot_dmat(1:2,1:2))
             rot_dmat(1:2,1:2)=MATMUl(
     &         rot_dmat(1:2,1:2),tmp_mat(1:2,1:2))
            ENDIF
            DEALLOCATE(tmp_mat)
C storing in Dmat
            Dmat(1,iorb)%mat(1,1)=rot_dmat(1,1)
            Dmat(2,iorb)%mat(1,1)=rot_dmat(2,2)
            Dmat(3,iorb)%mat(1,1)=rot_dmat(1,2)
            Dmat(4,iorb)%mat(1,1)=rot_dmat(2,1)
            DEALLOCATE(rot_dmat)
           ELSE
C -----------------------------------------------------------------------------------------------------
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) ) :
C -----------------------------------------------------------------------------------------------------
            IF (reptrans(l,isrt)%ifmixing) THEN
C We use the complete spin-space representation, so no trick on indices is necessary.
C
C Application of the operation inverse(Rloc).Dmat.(Rloc) :
C -------------------------------------------------------
             IF (rotloc(iatom)%timeinv) THEN
C In this case, the operators is antiunitary [ inverse(R)=transpose(R) ]
              Dmat(1,iorb)%mat(:,:)=CONJG(
     =          MATMUL(Dmat(1,iorb)%mat(:,:),
     &          rotloc(iatom)%rotrep(l)%mat(:,:) ))
              Dmat(1,iorb)%mat(:,:)=
     =          MATMUL(TRANSPOSE( rotloc(iatom)%
     &            rotrep(l)%mat(:,:) ),Dmat(1,iorb)%mat(:,:) )
C Dmat_{local} = inverse(Rloc) Dmat_{global}* Rloc*
C Dmat_{local} = transpose(Rloc) Dmat_{global}* Rloc*
             ELSE
C In this case, all the operators are unitary [ inverse(R)=transpose(conjugate(R)) ]
              Dmat(1,iorb)%mat(:,:)=
     =          MATMUL(Dmat(1,iorb)%mat(:,:),
     &          rotloc(iatom)%rotrep(l)%mat(:,:) )
              Dmat(1,iorb)%mat(:,:)=
     =          MATMUL(TRANSPOSE(CONJG( rotloc(iatom)%
     &            rotrep(l)%mat(:,:) )),Dmat(1,iorb)%mat(:,:) )
C Dmat_{local} = <x_local | x_global> Dmat_{global} <x_global | x_local> 
C Dmat_{local} = inverse(Rloc) Dmat_{global} Rloc
             ENDIF
C
            ELSE
C ----------------------------------------------------------------------------------------------
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only) :
C ----------------------------------------------------------------------------------------------
C definition of the total density matrix
             ALLOCATE(rot_dmat(1:2*(2*l+1),1:2*(2*l+1)))
             rot_dmat(1:(2*l+1),1:(2*l+1))=
     &         Dmat(1,iorb)%mat(-l:l,-l:l)
             rot_dmat(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))=
     &         Dmat(2,iorb)%mat(-l:l,-l:l)
             rot_dmat(1:(2*l+1),2*l+2:2*(2*l+1))=
     &         Dmat(3,iorb)%mat(-l:l,-l:l)
             rot_dmat(2*l+2:2*(2*l+1),1:(2*l+1))=
     &         Dmat(4,iorb)%mat(-l:l,-l:l)
             IF (rotloc(iatom)%timeinv) THEN
C In this case, the operator is antiunitary [ inverse(R)=transpose(R) ]
              rot_dmat(1:2*(2*l+1),1:2*(2*l+1))=CONJG(
     =          MATMUL(rot_dmat(1:2*(2*l+1),1:2*(2*l+1)),
     &          rotloc(iatom)%rotrep(l)
     &          %mat(1:2*(2*l+1),1:2*(2*l+1)) ))
              rot_dmat(1:2*(2*l+1),1:2*(2*l+1))=
     =          MATMUL(TRANSPOSE( rotloc(iatom)%
     &          rotrep(l)%mat(1:2*(2*l+1),1:2*(2*l+1)) ),
     &          rot_dmat(1:2*(2*l+1),1:2*(2*l+1)) )
C Dmat_{local} = inverse(Rloc) Dmat_{global}* Rloc*
C Dmat_{local} = transpose(Rloc) Dmat_{global}* Rloc*
             ELSE
C In this case, all the operators are unitary [ inverse(R)=transpose(conjugate(R)) ]
              rot_dmat(1:2*(2*l+1),1:2*(2*l+1))=
     =          MATMUL(rot_dmat(1:2*(2*l+1),1:2*(2*l+1)),
     &          rotloc(iatom)%rotrep(l)
     &          %mat(1:2*(2*l+1),1:2*(2*l+1)) )
              rot_dmat(1:2*(2*l+1),1:2*(2*l+1))=
     =          MATMUL(TRANSPOSE(CONJG( rotloc(iatom)%
     &          rotrep(l)%mat(1:2*(2*l+1),1:2*(2*l+1)) )),
     &          rot_dmat(1:2*(2*l+1),1:2*(2*l+1)) )
C Dmat_{local} = <x_local | x_global> Dmat_{global} <x_global | x_local> 
C Dmat_{local} = inverse(Rloc) Dmat_{global} Rloc
             ENDIF
C storing in dmat again
             Dmat(1,iorb)%mat(-l:l,-l:l)=
     &         rot_dmat(1:(2*l+1),1:(2*l+1))
             Dmat(2,iorb)%mat(-l:l,-l:l)=
     &         rot_dmat(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))
             Dmat(3,iorb)%mat(-l:l,-l:l)=
     &         rot_dmat(1:(2*l+1),2*l+2:2*(2*l+1))
             Dmat(4,iorb)%mat(-l:l,-l:l)=
     &         rot_dmat(2*l+2:2*(2*l+1),1:(2*l+1))
             DEALLOCATE(rot_dmat)
            ENDIF  ! End of the if mixing if-then-else
           ENDIF   ! End of the if "l=0" if-then-else
          ELSE
C ------------------------------------------------------------------------------
C The s-orbitals are a particular case of a "non-mixing" basis and is invariant.
C ------------------------------------------------------------------------------
           IF(l==0) CYCLE
C ----------------------------------------------------------------------------------------------
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only) :
C ----------------------------------------------------------------------------------------------
           ALLOCATE(rot_dmat(-l:l,-l:l))
           DO is=1,nsp
           rot_dmat=0.d0
C
C Application of the operation inverse(Rloc).Dmat.(Rloc) :
C -------------------------------------------------------
C In this case, (either a paramagnetic calculation or a spin-polarized one 
C but the symmetry operation does not change the magntization direction)
C all the operators are unitary [ inverse(R)=transpose(conjugate(R)) ]
              rot_dmat(-l:l,-l:l)=
     =          MATMUL(Dmat(is,iorb)%mat(-l:l,-l:l),
     &          rotloc(iatom)%rotrep(l)%mat(-l:l,-l:l) )
              rot_dmat(-l:l,-l:l)=
     =          MATMUL(TRANSPOSE(CONJG( rotloc(iatom)%
     &            rotrep(l)%mat(-l:l,-l:l) )),
     &            rot_dmat(-l:l,-l:l) )
C rotmat_{local} = <x_local | x_global> rotmat_{global} <x_global | x_local> 
C rotmat_{local} = inverse(Rloc) rotmat_{global} Rloc
C
C Storing the new value in Dmat : 
C -------------------------------
             Dmat(is,iorb)%mat(-l:l,-l:l)=rot_dmat(-l:l,-l:l)
           ENDDO
           DEALLOCATE(rot_dmat)
C
          ENDIF   ! End of the ifSO-ifSP if-then-else
        ENDDO     ! End of the iorb loop
C
        RETURN
        END





