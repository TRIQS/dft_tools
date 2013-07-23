
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

      SUBROUTINE timeinv_op(mat,lm,l,isrt)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine applies the time reversal operation to the      %%
C %% matrix mat which is associated to the l orbital of the atomic   %%
C %% isrt. (matrix size = lm) The matrix mat is assumed to already   %%
C %% be in the desired basis associated to isrt.                     %%
C %% The calculation done is :                                       %%
C %%       reptrans*T*conjg((inv(reptrans))*conjg(mat)               %%
C %%                                                                 %%
C %% If isrt=0, the matrix mat is assumed to be in the spherical     %%
C %% harmonics basis and no spin is considered. (lm = 2*l+1)         %%
C %% The calculation done is then : T*conjg(mat)                     %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE common_data
        USE reps
        IMPLICIT NONE
        INTEGER :: lm,l,isrt
        COMPLEX(KIND=8), DIMENSION(1:lm,1:lm) :: mat
        COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tinv
        COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmp_tinv
        COMPLEX(KIND=8), DIMENSION(-l:l,-l:l) :: tmat
        INTEGER :: m,n
C
C Definition of the complex conjugation operator in the spherical harmonic basis :
C --------------------------------------------------------------------------------
C
        tmat = CMPLX(0.d0,0.d0)
        DO m=-l,l
          tmat(m,-m)=(-1)**m
        END DO
C
C
C Calculation of the Time-reversal operator in the desired representation basis : 
C -------------------------------------------------------------------------------
C
        IF (isrt==0) THEN
C The case isrt=0 is a "default case" : 
C mat is in the spherical harmonic basis (without spinor representation)
         ALLOCATE(tinv(1:2*l+1,1:2*l+1))
         tinv(1:2*l+1,1:2*l+1)=tmat(-l:l,-l:l)
        ELSE
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) )
         IF (reptrans(l,isrt)%ifmixing) THEN
          ALLOCATE(tinv(1:2*(2*l+1),1:2*(2*l+1)))
          ALLOCATE(tmp_tinv(1:2*(2*l+1),1:2*(2*l+1)))
          tinv = CMPLX(0.d0,0.d0)
          tmp_tinv = CMPLX(0.d0,0.d0)
C Definition of the time-reversal operator as a spinor-operator (multiplication by -i.sigma_y)
          tinv(1:2*l+1,2*l+2:2*(2*l+1))=-tmat(-l:l,-l:l)
          tinv(2*l+2:2*(2*l+1),1:2*l+1)=tmat(-l:l,-l:l)
C The time reversal operator is put in the desired basis.
          tmp_tinv(1:2*(2*l+1),1:2*(2*l+1))=MATMUL(
     &      reptrans(l,isrt)%transmat(1:2*(2*l+1),1:2*(2*l+1)),
     &      tinv(1:2*(2*l+1),1:2*(2*l+1)))
          tinv(1:2*(2*l+1),1:2*(2*l+1))=MATMUL(
     &      tmp_tinv(1:2*(2*l+1),1:2*(2*l+1)),
     &      TRANSPOSE(reptrans(l,isrt)%transmat
     &           (1:2*(2*l+1),1:2*(2*l+1)) ) )
C the result tinv = (reptrans)*tinv*transpose(reptrans) 
C or tinv_{new_i} = <new_i|lm> tinv_{lm} (<lm|new_i>)*
C which is exactly the expression of the spinor operator in the new basis.
          DEALLOCATE(tmp_tinv)
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only)
         ELSE
          ALLOCATE(tinv(1:2*l+1,1:2*l+1))
          ALLOCATE(tmp_tinv(-l:l,-l:l))
          tinv = CMPLX(0.d0,0.d0)
          tmp_tinv = CMPLX(0.d0,0.d0)
C The time reversal operator is put in the desired basis.
          tmp_tinv(-l:l,-l:l)=MATMUL(
     &      reptrans(l,isrt)%transmat(-l:l,-l:l),
     &      tmat(-l:l,-l:l) )
          tinv(1:2*l+1,1:2*l+1)=MATMUL(
     &      tmp_tinv(-l:l,-l:l),TRANSPOSE(
     &      reptrans(l,isrt)%transmat(-l:l,-l:l)) )
          DEALLOCATE(tmp_tinv)
         END IF
C the result tinv = (reptrans)*tinv*transpose(reptrans) 
C or tinv_{new_i} = <new_i|lm> tinv_{lm} (<lm|new_i>)*
C which is exactly the expression of the operator in the new basis.
        END IF
C
C
C Multiplication of the matrix mat by the time reversal operator :
C ----------------------------------------------------------------
C
        mat(1:lm,1:lm) = MATMUL(
     &    tinv(1:lm,1:lm),CONJG(mat(1:lm,1:lm)) )
        DEALLOCATE(tinv)
C The multiplication is the product of tinv and (mat)*
C
        RETURN
        END



      SUBROUTINE add_timeinv(Dmat,orbit,norbit)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine calculates for each density matrix in Dmat      %%
C %% its image by the time-reversal operator and adds it to the      %%
C %% former one to get a time-symmetrized result.                    %%
C %%                                                                 %%
C %% This operation is done only if the computation is paramagnetic  %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
      USE common_data
      USE projections
      USE symm
      USE reps
      IMPLICIT NONE
      INTEGER :: norbit
      TYPE(matrix), DIMENSION(nsp,norbit) :: Dmat
      COMPLEX(KIND=8),DIMENSION(:,:,:), ALLOCATABLE :: rot_dmat
      COMPLEX(KIND=8),DIMENSION(:,:), ALLOCATABLE :: time_op
      COMPLEX(KIND=8),DIMENSION(:,:,:), ALLOCATABLE :: tmp_mat
      COMPLEX(KIND=8):: ephase
      TYPE(orbital), DIMENSION(norbit) :: orbit
      INTEGER :: isym, iorb, iatom, jorb, is, is1, l, i
      INTEGER :: isrt, jatom, imult, m
C
C
        DO iorb=1,norbit
          l=orbit(iorb)%l
          isrt=orbit(iorb)%sort
          iatom=orbit(iorb)%atom
C -----------------------------------------------------------------------------------
C The s-orbitals are a particular case of a "non-mixing" basis and are treated here :
C -----------------------------------------------------------------------------------
          IF(l==0) THEN
           IF (nsp==1) THEN
            Dmat(1,iorb)%mat(1,1) = 
     &        ( Dmat(1,iorb)%mat(1,1)+
     &        CONJG(Dmat(1,iorb)%mat(1,1)) )/2.d0
           ELSE
            ALLOCATE(tmp_mat(1,1,nsp))
            tmp_mat=0.d0
C Application of the time-reversal operation
C ------------------------------------------
            DO is=1,nsp
              is1=is+(-1)**(is+1)
C the time reversal operation transforms up/up -1- in dn/dn -2- and up/dn -3- in dn/up -4- (and vice versa)
              tmp_mat(1,1,is)=CONJG(Dmat(is1,iorb)%mat(1,1) )
              IF (is.gt.2) tmp_mat(1,1,is)=-tmp_mat(1,1,is)
C Off diagonal blocks are multiplied by (-1).
            ENDDO
C Symmetrization of Dmat :
C ------------------------
            DO is=1,nsp
              Dmat(is,iorb)%mat(1,1) = (Dmat(is,iorb)%mat(1,1)+
     &          tmp_mat(1,1,is) )/2.d0
            ENDDO
            DEALLOCATE(tmp_mat)
           ENDIF
C -----------------------------------------------------------------------------------------------------
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) ) :
C -----------------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C Calculation of the time-reversal operator :
C -------------------------------------------
           ALLOCATE(time_op(1:2*(2*l+1),1:2*(2*l+1)))
           time_op(:,:)=0.d0
           DO m=1,2*(2*l+1)
             time_op(m,m)=1.d0
           ENDDO
C time_op is Identity.
           CALL timeinv_op(time_op,2*(2*l+1),l,isrt)
C time_op is now the time-reversal operator in the desired basis ({new_i})
C
C Application of the time-reversal operation
C ------------------------------------------
           ALLOCATE(tmp_mat(1:2*(2*l+1),1:2*(2*l+1),1))
           tmp_mat(1:2*(2*l+1),1:2*(2*l+1),1)=
     =       MATMUL(Dmat(1,iorb)%mat(1:2*(2*l+1),1:2*(2*l+1)),
     &       TRANSPOSE(time_op(1:2*(2*l+1),1:2*(2*l+1)) ) )
           tmp_mat(1:2*(2*l+1),1:2*(2*l+1),1)=
     =       MATMUL(time_op(1:2*(2*l+1),1:2*(2*l+1)),
     &       CONJG(tmp_mat(1:2*(2*l+1),1:2*(2*l+1),1) ) )
C The operation performed is : time_op.conjugate(Dmat).transpose(conjugate(time_op))
C or in other words,           D(T)_{new_i} . Dmat* . D(inverse(T))*_{new_i}
C
C Symmetrization of Dmat :
C ------------------------       
           Dmat(1,iorb)%mat(1:2*(2*l+1),1:2*(2*l+1)) = 
     &       ( Dmat(1,iorb)%mat(1:2*(2*l+1),1:2*(2*l+1)) + 
     &       tmp_mat(1:2*(2*l+1),1:2*(2*l+1),1) )/2.d0
           DEALLOCATE(tmp_mat)
           DEALLOCATE(time_op)
C ----------------------------------------------------------------------------------------------
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only) :
C ----------------------------------------------------------------------------------------------
          ELSE
C Calculation of the time-reversal operator :
C -------------------------------------------
           ALLOCATE(time_op(-l:l,-l:l))
           time_op(:,:)=0.d0
           DO m=-l,l
             time_op(m,m)=1.d0
           ENDDO
C time_op is Identity.
           CALL timeinv_op(time_op,(2*l+1),l,isrt)
C time_op is now the time-reversal operator in the desired basis ({new_i})
C
           IF (nsp==1) THEN
C Application of the time-reversal operation and symmetrization :
C ---------------------------------------------------------------
            ALLOCATE(tmp_mat(-l:l,-l:l,1))
            tmp_mat(-l:l,-l:l,1)=
     =        MATMUL( Dmat(1,iorb)%mat(-l:l,-l:l),
     &        TRANSPOSE(time_op(-l:l,-l:l) ) )
            tmp_mat(-l:l,-l:l,1)=
     =        MATMUL(time_op(-l:l,-l:l),
     &        CONJG(tmp_mat(-l:l,-l:l,1)) )
C The operation performed is : time_op.conjugate(Dmat).transpose(conjugate(time_op))
C or in other words,           D(T)_{new_i} . Dmat* . D(inverse(T))*_{new_i}
            Dmat(1,iorb)%mat(-l:l,-l:l) = 
     &       ( Dmat(1,iorb)%mat(-l:l,-l:l) + 
     &       tmp_mat(-l:l,-l:l,1) )/2.d0
            DEALLOCATE(tmp_mat)
           ELSE
C Application of the time-reversal operation
C ------------------------------------------
            ALLOCATE(tmp_mat(-l:l,-l:l,nsp))
            DO is=1,nsp
              is1=is+(-1)**(is+1)
C the time reversal operation transforms up/up -1- in dn/dn -2- and up/dn -3- in dn/up -4 (and vice versa)
              tmp_mat(-l:l,-l:l,is)=
     =          MATMUL( Dmat(is1,iorb)%mat(-l:l,-l:l),
     &          TRANSPOSE( time_op(-l:l,-l:l) ) )
              tmp_mat(-l:l,-l:l,is)=
     =          MATMUL( time_op(-l:l,-l:l),
     &          CONJG( tmp_mat(-l:l,-l:l,is) ) )
C The operation performed is : time_op.conjugate(Dmat).transpose(conjugate(time_op))
C or in other words,           D(T)_{new_i} . Dmat* . D(inverse(T))*_{new_i}
              IF (is.gt.2) THEN 
               tmp_mat(-l:l,-l:l,is)=-tmp_mat(-l:l,-l:l,is)
              ENDIF
C Off diagonal terms are multiplied by (-1).
            ENDDO
C Symmetrization of Dmat :
C ------------------------
            DO is=1,nsp
              Dmat(is,iorb)%mat(-l:l,-l:l) = 
     &          (Dmat(is,iorb)%mat(-l:l,-l:l)+
     &          tmp_mat(-l:l,-l:l,is) )/2.d0
            ENDDO
            DEALLOCATE(tmp_mat)
           ENDIF
           DEALLOCATE(time_op)
C
          ENDIF     ! End of the type basis if-then-else
        ENDDO       ! End of the iorb loop
C         
        RETURN
        END




