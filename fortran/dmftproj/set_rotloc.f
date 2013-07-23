
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

       SUBROUTINE set_rotloc
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine sets up the Global->local coordinates           %%
C %% rotational matrices for each atom of the system.                %%
C %% These matrices will be used to create the projectors.           %%
C %% (They are the SR matrices defined in the tutorial file.)        %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
       USE common_data
       USE reps
       USE symm
       USE prnt
       IMPLICIT NONE
       COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmp_rot, spinrot
       REAL(KIND=8) :: alpha, beta, gama, factor
       INTEGER :: iatom, jatom, imu, isrt
       INTEGER :: is, is1, isym, l, lm
       INTEGER :: ind1, ind2, inof1, inof2
       COMPLEX(KIND=8) :: ephase
C
C ====================================================
C Multiplication by an S matrix for equivalent sites :
C ====================================================
C Up to now, rotloc is the rotloc matrix (from Global to local coordinates rotation : (rotloc)_ij = <x_global_i | x_local_j >)
C The matrix S to go from the representative atom of the sort to another one must be introduced. That's what is done here-after.
       DO isrt=1,nsort
         iatom=SUM(nmult(0:isrt-1))+1
         DO imu=1,nmult(isrt)
           jatom=iatom+imu-1
           DO isym=1,nsym
C If the symmetry operation isym transforms the representative atom iatom in the jatom, 
C the matrix rotloc is multiplied by the corresponding srot matrix, for each orbital number l.
C if R[isym](iatom) = jatom, rotloc is multiplied by R[isym] and Rloc is finally R[isym] X rotloc = <x_global|x_sym><x_sym|x_local>
             IF(srot(isym)%perm(iatom)==jatom) THEN
              WRITE(17,*) ' For jatom = ',jatom, ', isym =', isym 
              rotloc(jatom)%srotnum=isym
C Calculation of krotm and iprop.
              rotloc(jatom)%krotm(1:3,1:3)=
     =          MATMUL(srot(isym)%krotm(1:3,1:3), 
     &          rotloc(jatom)%krotm(1:3,1:3))
              rotloc(jatom)%iprop=rotloc(jatom)%iprop*
     *          srot(isym)%iprop
C Evaluation of the Euler angles of the final operation Rloc
              CALL euler(TRANSPOSE(rotloc(jatom)%krotm(1:3,1:3)), 
     &          alpha,beta,gama)
C According to Wien convention, euler takes in argument the transpose 
C of the matrix rotloc(jatom)%krotm to give a,b anc c of rotloc(jatom).
              rotloc(jatom)%a=alpha
              rotloc(jatom)%b=beta
              rotloc(jatom)%g=gama
C
C =============================================================================================================
C Calculation of the rotational matrices and evaluation of the fields timeinv and phase for the Rloc matrices :
C =============================================================================================================
              IF(ifSP.AND.ifSO) THEN
C No time reversal operation is applied to rotloc (alone). If a time reversal operation must be applied, 
C it comes from the symmetry operation R[isym]. That is why the field timeinv is the same as the one from srot.
               rotloc(jatom)%timeinv=srot(isym)%timeinv
               rotloc(jatom)%phase=0.d0
               DO l=1,lmax
                 ALLOCATE(tmp_rot(1:2*(2*l+1),1:2*(2*l+1)))
                 tmp_rot=0.d0
C Whatever the value of beta (0 or Pi), the spinor rotation matrix of isym is block-diagonal.
C because the time-reversal operation have been applied if necessary.
                 factor=srot(isym)%phase/2.d0
                 ephase=EXP(CMPLX(0.d0,factor))
C We remind that the field phase is (g-a) if beta=Pi. As a result, ephase = exp(+i(g-a)/2) = -exp(+i(alpha-gamma)/2)
C We remind that the field phase is (a+g) if beta=0. As a result, ephase = exp(+i(a+g)/2)=-exp(-i(alpha+gamma)/2)
C in good agreement with Wien conventions for the definition of this phase factor.
C Up/up block :
                 tmp_rot(1:2*l+1,1:2*l+1)=ephase*
     &              srot(isym)%rotl(-l:l,-l:l,l)
C Dn/dn block :
                 ephase=CONJG(ephase)
C now, ephase = exp(+i(a-g)/2) = -exp(-i(alpha-gamma)/2) if beta=Pi
C now, ephase = exp(-i(a+g)/2) = -exp(+i(alpha+gamma)/2) if beta=0
                 tmp_rot(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))=
     &             ephase*srot(isym)%rotl(-l:l,-l:l,l)
                 IF (rotloc(jatom)%timeinv) THEN
C In this case, the time reversal operator was applied to srot.
                  rotloc(jatom)%rotl(1:2*(2*l+1),1:2*(2*l+1),l)= 
     &              MATMUL(tmp_rot(1:2*(2*l+1),1:2*(2*l+1)),CONJG(
     &              rotloc(jatom)%rotl(1:2*(2*l+1),1:2*(2*l+1),l)))
C rotloc(jatom)%rotl now contains D(Rloc) = D(R[isym])*transpose[D(rotloc)].
                 ELSE
C In this case, no time reversal operator was applied to srot.
                  rotloc(jatom)%rotl(1:2*(2*l+1),1:2*(2*l+1),l)= 
     &              MATMUL(tmp_rot(1:2*(2*l+1),1:2*(2*l+1)),
     &              rotloc(jatom)%rotl(1:2*(2*l+1),1:2*(2*l+1),l))
C rotloc(jatom)%rotl now contains D(Rloc) = D(R[isym])*D(rotloc).
                 ENDIF
                 DEALLOCATE(tmp_rot)
               ENDDO
              ELSE
C Calculation of the rotational matrices associated to Rloc
               ALLOCATE(tmp_rot(1:2*lmax+1,1:2*lmax+1))
               DO l=1,lmax
C Use of the subroutine dmat to compute the rotational matrix 
C associated to the Rloc operation in a (2*l+1) space :
                 tmp_rot=0.d0
                 CALL dmat(l,rotloc(jatom)%a,rotloc(jatom)%b,
     &             rotloc(jatom)%g,
     &             REAL(rotloc(jatom)%iprop,KIND=8),tmp_rot,2*lmax+1)
                 rotloc(jatom)%rotl(-l:l,-l:l,l)=
     =             tmp_rot(1:2*l+1,1:2*l+1)
C rotloc(jatom)%rotl = table of the rotational matrices of the symmetry operation 
C for the different l orbital (from 1 to lmax), in the usual complex basis : dmat = D(R[isym])_l
C rotloc(jatom)%rotl = D(Rloc[jatom])_{lm}
               ENDDO
               DEALLOCATE(tmp_rot)
              ENDIF ! End of the "ifSO-ifSP" if-then-else
C
              EXIT
C Only one symmetry operation is necessary to be applied to R to get the complete rotloc matrix. 
C This EXIT enables to leave the loop as soon as a symmetry operation which transforms the representative atom in jatom is found.
             ENDIF ! End of the "perm" if-then-else
           ENDDO   ! End of the isym loop
C
C
C ===========================================================
C Computation of the rotational matrices in each sort basis :
C ===========================================================
           ALLOCATE(rotloc(jatom)%rotrep(lmax))
C
C Initialization of the rotloc(jatom)%rotrep field = D(Rloc)_{new_i}
C This field is a table of size lmax which contains the rotloc matrices 
C in the representation basis associated to each included orbital of the jatom.
           DO l=1,lmax
             ALLOCATE(rotloc(jatom)%rotrep(l)%mat(1,1))
             rotloc(jatom)%rotrep(l)%mat(1,1)=0.d0
           ENDDO
C
C Computation of the elements 'mat' in rotloc(jatom)%rotrep(l)
           DO l=1,lmax
C The considered orbital is not included, hence no computation
             IF (lsort(l,isrt)==0) cycle
C The considered orbital is included
             IF (ifSP.AND.ifSO) THEN
C In this case, the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) )
C --------------------------------------------------------------------------------------------------------------
              DEALLOCATE(rotloc(jatom)%rotrep(l)%mat)
              ALLOCATE(rotloc(jatom)%rotrep(l)%mat
     &          (1:2*(2*l+1),1:2*(2*l+1)))
              ALLOCATE(tmp_rot(1:2*(2*l+1),1:2*(2*l+1)))
C Computation of rotloc(jatom)%rotrep(l)%mat
              IF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the basis representation requires a complete spinor rotation approach too.
               IF(rotloc(jatom)%timeinv) THEN
                tmp_rot(1:2*(2*l+1),1:2*(2*l+1))=MATMUL(
     &           reptrans(l,isrt)%transmat(1:2*(2*l+1),1:2*(2*l+1)),
     &           rotloc(jatom)%rotl(1:2*(2*l+1),1:2*(2*l+1),l))
                rotloc(jatom)%rotrep(l)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =           MATMUL(tmp_rot(1:2*(2*l+1),1:2*(2*l+1)),
     &           TRANSPOSE(reptrans(l,isrt)%transmat
     &           (1:2*(2*l+1),1:2*(2*l+1))))
C Since the operation is antilinear, the field rotloc(jatom)%rotrep(l)%mat = (reptrans)*spinrot(l)*conjugate(inverse(reptrans))
C rotloc(jatom)%rotrep(l)%mat = D(Rloc)_{new_i} = <new_i|lm> D(Rloc)_{lm} [<lm|new_i>]^*
C which is exactly the expression of the spinor rotation matrix in the new basis.
               ELSE
                tmp_rot(1:2*(2*l+1),1:2*(2*l+1))=MATMUL(
     &           reptrans(l,isrt)%transmat(1:2*(2*l+1),1:2*(2*l+1)),
     &           rotloc(jatom)%rotl(1:2*(2*l+1),1:2*(2*l+1),l))
                rotloc(jatom)%rotrep(l)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =           MATMUL(tmp_rot(1:2*(2*l+1),1:2*(2*l+1)),
     &           TRANSPOSE(CONJG(reptrans(l,isrt)%transmat
     &           (1:2*(2*l+1),1:2*(2*l+1)))))
C Since the operation is linear, the field rotloc(jatom)%rotrep(l)%mat = (reptrans)*spinrot(l)*inverse(reptrans)
C rotloc(jatom)%rotrep(l)%mat = D(Rloc)_{new_i} = <new_i|lm> D(Rloc)_{lm} <lm|new_i>
C which is exactly the expression of the spinor rotation matrix in the new basis.
               ENDIF
              ELSE
C In this case, the basis representation is reduced to the up/up block and must be extended.
               ALLOCATE(spinrot(1:2*(2*l+1),1:2*(2*l+1)))
               spinrot(1:2*(2*l+1),1:2*(2*l+1))=0.d0
               spinrot(1:2*l+1,1:2*l+1)=
     &           reptrans(l,isrt)%transmat(-l:l,-l:l)
               spinrot(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))=
     &           reptrans(l,isrt)%transmat(-l:l,-l:l)
               IF(rotloc(jatom)%timeinv) THEN
                tmp_rot(1:2*(2*l+1),1:2*(2*l+1))=MATMUL(
     &           spinrot(1:2*(2*l+1),1:2*(2*l+1)),
     &           rotloc(jatom)%rotl(1:2*(2*l+1),1:2*(2*l+1),l))
                rotloc(jatom)%rotrep(l)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =           MATMUL(tmp_rot(1:2*(2*l+1),1:2*(2*l+1)),
     &           TRANSPOSE(spinrot(1:2*(2*l+1),1:2*(2*l+1))))
C Since the operation is antilinear, the field rotloc(jatom)%rotrep(l)%mat = (reptrans)*spinrot(l)*conjugate(inverse(reptrans))
C rotloc(jatom)%rotrep(l)%mat = D(Rloc)_{new_i} = <new_i|lm> D(Rloc)_{lm} [<lm|new_i>]^*
C which is exactly the expression of the spinor rotation matrix in the new basis.
               ELSE
                tmp_rot(1:2*(2*l+1),1:2*(2*l+1))=MATMUL(
     &           spinrot(1:2*(2*l+1),1:2*(2*l+1)),
     &           rotloc(jatom)%rotl(1:2*(2*l+1),1:2*(2*l+1),l))
                rotloc(jatom)%rotrep(l)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =           MATMUL(tmp_rot(1:2*(2*l+1),1:2*(2*l+1)),
     &           TRANSPOSE(CONJG(spinrot(1:2*(2*l+1),1:2*(2*l+1)))))
C Since the operation is linear, the field rotloc(jatom)%rotrep(l)%mat = (reptrans)*spinrot(l)*inverse(reptrans)
C rotloc(jatom)%rotrep(l)%mat = D(Rloc)_{new_i} = <new_i|lm> D(Rloc)_{lm} <lm|new_i>
C which is exactly the expression of the spinor rotation matrix in the new basis.
               ENDIF
               DEALLOCATE(spinrot)
              ENDIF  ! End of the if mixing if-then-else
              DEALLOCATE(tmp_rot)
C
             ELSE
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only)
C --------------------------------------------------------------------------------------------
              DEALLOCATE(rotloc(jatom)%rotrep(l)%mat)
              ALLOCATE(rotloc(jatom)%rotrep(l)%mat(-l:l,-l:l))
              ALLOCATE(tmp_rot(-l:l,-l:l))
C Computation of rotloc(jatom)%rotrep(l)%mat
              tmp_rot(-l:l,-l:l)=MATMUL(
     &          reptrans(l,isrt)%transmat(-l:l,-l:l),
     &          rotloc(jatom)%rotl(-l:l,-l:l,l))
              rotloc(jatom)%rotrep(l)%mat(-l:l,-l:l)=
     =          MATMUL(tmp_rot(-l:l,-l:l),
     &          TRANSPOSE(CONJG(reptrans(l,isrt)%transmat(-l:l,-l:l))))
C the field rotloc(jatom)%rotrep(l)%mat = (reptrans)*rotl*inverse(reptrans)
C rotloc(jatom)%rotrep(l)%mat = D(Rloc)_{new_i} = <new_i|lm> D(Rloc)_{lm} <lm|new_i>
C which is exactly the expression of the rotation matrix for the up/up block in the new basis.
              DEALLOCATE(tmp_rot)
             ENDIF
           ENDDO   ! End of the l loop
         ENDDO     ! End of the jatom loop
       ENDDO       ! End of the isrt loop
C
       RETURN
       END


       SUBROUTINE euler(Rot,a,b,c)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine calculates the Euler angles a, b and c of Rot.  %%
C %% The result are stored in a,b,c. (same as in SRC_lapwdm/euler.f) %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
        IMPLICIT NONE
        REAL(KIND=8) :: a,aa,b,bb,c,cc,zero,pi,y_norm,dot
        REAL(KIND=8), DIMENSION(3,3) :: Rot, Rot_temp
        REAL(KIND=8), DIMENSION(3) :: z,zz,y,yy,yyy,pom,x,xx
        INTEGER :: i,j
C Definition of the constants
        zero=0d0
        pi=ACOS(-1d0)
C Definition of Rot_temp=Id
        DO i=1,3
          DO j=1,3
            Rot_temp(i,j)=0
            IF (i.EQ.j) Rot_temp(i,i)=1
	  ENDDO
        ENDDO
C Initialization of y=e_y, z=e_z, yyy and zz
        DO j=1,3
          y(j)=Rot_temp(j,2)
          yyy(j)=Rot(j,2)
          z(j)=Rot_temp(j,3)
          zz(j)=Rot(j,3)
        ENDDO
C Calculation of yy
        CALL vecprod(z,zz,yy)
        y_norm=DSQRT(dot(yy,yy))
        IF (y_norm.lt.1d-10) THEN
C If yy=0, this implies that b is zero or pi
         IF (ABS(dot(y,yyy)).gt.1d0) THEN
          aa=dot(y,yyy)/ABS(dot(y,yyy))
          a=ACOS(aa)
         ELSE
          a=ACOS(dot(y,yyy))
         ENDIF
C
         IF (dot(z,zz).gt.zero) THEN
          c=zero
          b=zero
          IF (yyy(1).gt.zero) a=2*pi-a
         ELSE
          c=a
          a=zero
          b=pi
          IF (yyy(1).lt.zero) c=2*pi-c
         ENDIF
	ELSE
C If yy is not 0, then b belongs to ]0,pi[
         DO j=1,3
           yy(j)=yy(j)/y_norm
         ENDDO
C
         aa=dot(y,yy)
         bb=dot(z,zz)
         cc=dot(yy,yyy)
         IF (ABS(aa).gt.1d0) aa=aa/ABS(aa)
         IF (ABS(bb).gt.1d0) bb=bb/ABS(bb)
         IF (ABS(cc).gt.1d0) cc=cc/ABS(cc)
         b=ACOS(bb)
         a=ACOS(aa)
         c=ACOS(cc)
         IF (yy(1).gt.zero) a=2*pi-a
         CALL vecprod(yy,yyy,pom)
         IF (dot(pom,zz).lt.zero) c=2*pi-c
        ENDIF
C
	END


       SUBROUTINE vecprod(a,b,c)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine calculates the vector product of a and b.       %%
C %% The result is stored in c. (same as in SRC_lapwdm/euler.f)      %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(3) :: a,b,c
C
       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
C
       END

       REAL(KIND=8) FUNCTION dot(a,b)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This function calculates the scalar product of a and b.         %%
C %% The result is stored in dot. (same as in SRC_lapwdm/euler.f)    %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
       IMPLICIT NONE
       REAL(KIND=8) :: a,b
       INTEGER :: i
       dimension a(3),b(3)
       dot=0
       DO i=1,3
         dot=dot+a(i)*b(i)
       ENDDO
C
       END



