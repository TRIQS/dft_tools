
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

        SUBROUTINE orthogonal_h(s1,ndim,inv)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine computes :                                      %%
C %%   - if inv = .FALSE. the square root of the Hermitian matrix s1 %%
C %%   - if inv = .TRUE. the inverse of the square root of the       %%
C %%        Hermitian matrix s1                                      %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE prnt
        IMPLICIT NONE
        INTEGER :: ndim, INFO, lm, lm1
        COMPLEX(KIND=8), DIMENSION(ndim) :: WORK
        COMPLEX(KIND=8), DIMENSION(ndim,ndim) :: s1
        INTEGER, DIMENSION(ndim,ndim) :: IPIV
        LOGICAL :: inv
C
C Calculation of S1^(1/2) or S1^(-1/2):
C -------------------------------------
        CALL sqrtm(s1,ndim,inv)
C The resulting matrix is stored in s1.
        RETURN
        END

        SUBROUTINE orthogonal_r(s2,ndim,inv)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine computes :                                      %%
C %%   - if inv = .FALSE. the square root of s1                      %%
C %%   - if inv = .TRUE. the inverse of the square root of s2        %%
C %% where s2 is a real symmetric matrix.                            %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE prnt
        IMPLICIT NONE
        INTEGER :: ndim, INFO, lm, lm1
        COMPLEX(KIND=8), DIMENSION(ndim) :: WORK
        COMPLEX(KIND=8), DIMENSION(ndim,ndim) :: s1
        REAL(KIND=8), DIMENSION(ndim,ndim) :: s2
        INTEGER, DIMENSION(ndim,ndim) :: IPIV
        LOGICAL :: inv
C
C Calculation of S2^(1/2) or S2^(-1/2):
C -------------------------------------
        s1=s2     
        CALL sqrtm(s1,ndim,inv)
        s2=REAL(s1)
C The resulting matrix is stored in s2.
        RETURN
        END


        SUBROUTINE sqrtm(cmat,m,inv)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine calculates the square root of a positively      %% 
C %% defined Hermitian matrix A=cmat using the decomposition         %%
C %%                       A=Z*D*Z^H                                 %% 
C %% where D is a diagonal matrix of eigenvalues of A,               %% 
C %%       Z is matrix of orthonormal eigenvectors of A,             %% 
C %%       Z^H is its Hermitian conjugate.                           %% 
C %% Then A^(1/2)=Z*D^(1/2)*Z^H.                                     %%
C %% Correction: the matrix A is allowed to be negatively defined.   %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        IMPLICIT NONE
        INTEGER :: m
        COMPLEX(KIND=8), DIMENSION(m,m):: cmat, D, D1 
        LOGICAL :: inv 
C Calculation of Z*D^(1/2):
C -------------------------
        CALL sqrt_eigenvec(cmat,D1,m,inv)
        WRITE(95,*) cmat
        WRITE(95,*) ' '
        WRITE(95,*) D1
        WRITE(95,*) ' '
C Calculation of A^(1/2)=Z*D^(1/2)*Z^H:
C -------------------------------------
        D=CONJG(cmat)
        call ZGEMM('N','T',m,m,m,DCMPLX(1.D0,0.D0),D1,
     &    m,D,m,DCMPLX(0.D0,0.D0),cmat,m)
C The resulting matrix is stored in cmat.
        RETURN
        END


        SUBROUTINE sqrt_eigenvec(cmat,D1,m,inv)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine computes :                                      %%
C %%   - if inv = .FALSE. Z*D^(1/2)                                  %%
C %%   - if inv = .TRUE.  Z*D^(-1/2)                                 %%
C %% where Z is a matrix of orthonormal eigenvectors of cmat and     %%
C %% D is the diagonal matrix of cmat's eigenvalues.                 %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE prnt
        IMPLICIT NONE
        LOGICAL :: inv, ifwrite
        INTEGER :: m, INFO, i, j
        INTEGER, PARAMETER :: nwork=40
C      
        COMPLEX(KIND=8), allocatable, DIMENSION(:) :: WORK     
        COMPLEX(KIND=8), DIMENSION(m,m) :: cmat, D1
        REAL(KIND=8), DIMENSION(m) :: W
        COMPLEX(KIND=8), DIMENSION(m) :: W_comp
        REAL(KIND=8), allocatable, DIMENSION(:) :: RWORK
C
C Finding the eigenvalues and the eigenvectors of cmat :
C ------------------------------------------------------
	ALLOCATE(rwork(3*m-2))
	ALLOCATE(work(2*m-1))
        CALL ZHEEV('V', 'U', m, cmat, m, W, WORK,2*m-1,RWORK,INFO) 
        IF (info.ne.0) THEN
         WRITE(buf,'(a)')
     &     'The subroutine zheev ends with info = ',info
         CALL printout(0)
         WRITE(buf,'(a)')'In sqrt_eigenvec, a pbm occurs in zheev.'
         CALL printout(0)
         WRITE(buf,'(a)')'END OF THE PRGM'
         CALL printout(0)
         STOP
        ENDIF       
C W contains the eigenvalues of cmat.
        W_comp=CMPLX(W,0d0)
C 
C Checking of the validity of the computation :
C ---------------------------------------------
        ifwrite=.FALSE.
        DO j=1,m
C The warning is written only once in the file case.outdmftpr
          IF (ifwrite) EXIT
C Checking if the eigenvalues are not negative.
          IF (W(j).lt.0.d0) THEN
           WRITE(buf,'(a,i2,a,a)')
     &       'WARNING : An eigenvalue (',j,') of the ',
     &       'overlap matrix is negative.'
           CALL printout(0)
           WRITE(buf,'(a,a)')'          The result ',
     &       'of the calculation may thus be wrong.'
           CALL printout(1)
           ifwrite=.TRUE.
          ENDIF       
          IF (ABS(W(j)).lt.1.d-12) THEN
           WRITE(buf,'(a,i2,a,a)')
     &       'WARNING : An eigenvalue (',j,') of the ',
     &       'overlap matrix is almost zero.'
           CALL printout(0)
           WRITE(buf,'(a,a)')'          The result ',
     &       'of the calculation may thus be wrong.'
           CALL printout(1)
           ifwrite=.TRUE.
          ENDIF       
        ENDDO
C
C Calculation of Z*D^(1/2) :
C --------------------------
C The result is stored in D1. 
        IF(.NOT.inv) THEN
         DO i=1,m
           DO j=1,m
             D1(i,j)=cmat(i,j)*SQRT(W_comp(j))
           ENDDO
         ENDDO   
        ELSE
C Calculation of Z*D^(-1/2) :
C --------------------------- 
C The result is stored in D1. 
         DO i=1,m
           DO j=1,m
             IF (ABS(W(j))==0.d0) THEN
              WRITE(buf,'(a,i2,a)')
     &          'An eigenvalue (',j,') of the ',
     &       'overlap matrix has the value 0.'
              CALL printout(0)
              WRITE(buf,'(a)')
     &          'The calculation can not be performed further.'
              CALL printout(0)
              CALL printout(0)
              WRITE(buf,'(a)')'END OF THE PRGM'
              CALL printout(0)
              STOP
             ENDIF
             D1(i,j)=cmat(i,j)/SQRT(W_comp(j))
           ENDDO
         ENDDO   
        ENDIF
C The resulting matrix is stored in D1 and cmat is now Z.
        RETURN
        END      
       
