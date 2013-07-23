
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

        SUBROUTINE setsym
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine sets up the symmetry matrices of the structure  %%
C %% and the local rotation matrices for each atom of the system.    %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE common_data
        USE factorial
        USE file_names
        USE prnt
        USE reps
        USE symm
        IMPLICIT NONE
C
        COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmp_rot, spinrot
        COMPLEX(KIND=8), DIMENSION(:,:),ALLOCATABLE :: tmat
        COMPLEX(KIND=8), DIMENSION(:,:,:),ALLOCATABLE :: tmp_dmat
        REAL(KIND=8) :: factor
        INTEGER :: l, isym, mmax, nrefl, i, m, isrt, lms
        INTEGER :: lm, is, is1
        INTEGER :: iatom, imu, iatomref
        REAL(KIND=8) :: det
        REAL(KIND=8), DIMENSION(:),ALLOCATABLE :: bufreal
        COMPLEX(KIND=8), DIMENSION(:),ALLOCATABLE :: bufcomp
        COMPLEX(KIND=8), DIMENSION(:,:),ALLOCATABLE :: tmpcomp
        COMPLEX(KIND=8), DIMENSION(1:2,1:2) :: spmt

C
C
        WRITE(buf,'(a)')'======================================='
        CALL printout(0)
        WRITE(buf,'(a)')'Symmetry operations of the system'
        CALL printout(1)
C
C ===========================================
C Reading of the symmetry file case.dmftsym :
C ===========================================
        CALL setfact(170)
        READ(iusym,*)nsym
        WRITE(buf,'(a,i4)')'Number of Symmetries = ',nsym
        CALL printout(0)
        CALL printout(0)
C nsym = total number of symmetry operations for the structure
        lsym=lmax
        nlmsym=2*lsym+1
C lsym = maximal orbital number for the symmetry
C nlmsym = maximal size of the representation for the symmetry 
        ALLOCATE(srot(nsym))
        DO isym=1,nsym
          ALLOCATE(srot(isym)%perm(natom))
          READ(iusym,*)srot(isym)%perm 
        ENDDO
C srot = table of symop elements from to 1 to nsym.
C the field srot(isym)%perm = the table of permutation for the isym symmetry (table from 1 to natom)
C srot(isym)%perm(iatom) = R[isym](iatom) = image by R[isym] fo iatom
        WRITE(buf,'(a)')'Properties of the symmetry operations :'
        CALL printout(0)
        WRITE(buf,'(a)') '   alpha, beta, gamma are their Euler angles.'
        CALL printout(0)
        WRITE(buf,'(a)') '   iprop is the value of their determinant.'
        CALL printout(0)
        CALL printout(0)
        WRITE(buf,'(a)')' SYM.OP.  alpha      beta     gamma     iprop'
        CALL printout(0)
        DO isym=1,nsym
          READ(iusym,'()')
          READ(iusym,'()')
          READ(iusym,'(3(f6.1),i3)') srot(isym)%a, srot(isym)%b, 
     &        srot(isym)%g, srot(isym)%iprop
C Printing the matrices parameters in the file case.outdmftpr
          WRITE(buf,'(i5,3F10.1,5x,i3)')isym,
     &      srot(isym)%a,srot(isym)%b,srot(isym)%g,srot(isym)%iprop
          CALL printout(0)
          srot(isym)%a=srot(isym)%a/180d0*Pi
          srot(isym)%b=srot(isym)%b/180d0*Pi
          srot(isym)%g=srot(isym)%g/180d0*Pi
C the field srot(isym)%a is linked to the Euler precession angle (alpha)
C the field srot(isym)%b is linked to the Euler nutation angle (beta)
C the field srot(isym)%c is linked to the Euler intrinsic rotation angle (gamma)
C They are read in case.dmftsym in degree and are then transformed into radians 
C the field sort(isym)% iprop = value of the transformation determinant (1 or -1), 
C determines if there is an inversion in the transformation
          READ(iusym,*)(srot(isym)%krotm(1:3,i),i=1,3)
           srot(isym)%krotm(1:3,1:3)=
     &       TRANSPOSE(srot(isym)%krotm(1:3,1:3))
C the field srot(isym)%krotm = 3x3 matrices of rotation associated to the transformation (R[isym]).
C (without the global inversion). The matrix was multiplied by the value of iprop before being written in case.dmftsym.
C This reading line was chosen to be consistent with the writing line in rotmat_dmft (in SRC_lapw2)
        ENDDO
C
C =============================================================
C Determination of the properties for each symmetry operation :
C =============================================================
C
C Creation of the rotational matrices for each orbital :
C ------------------------------------------------------
        DO isym=1,nsym
          ALLOCATE(srot(isym)%rotl(-lsym:lsym,-lsym:lsym,lsym))
          srot(isym)%rotl=0.d0
          ALLOCATE(tmat(1:2*lsym+1,1:2*lsym+1))
          DO l=1,lsym
C Use of the subroutine dmat to compute the the rotational matrix 
C associated to the isym symmetry operation in a (2*l+1) space :
            CALL dmat(l,srot(isym)%a,srot(isym)%b,srot(isym)%g,
     &        REAL(srot(isym)%iprop,KIND=8),tmat,2*lsym+1)
            srot(isym)%rotl(-l:l,-l:l,l)=tmat(1:2*l+1,1:2*l+1)
C srot(isym)%rotl = table of the rotationnal matrices of the symmetry operation 
C for the different l orbital (from 1 to lsym), in the usual complex basis : dmat = D(R[isym])_l
C srot(isym)%rotl = D(R[isym])_{lm}
          ENDDO
          DEALLOCATE(tmat)
C
C
C Determination of the fields timeinv and phase (if SP+SO computations):
C ----------------------------------------------------------------------
C If the calculation is spin-polarized with spin-orbit, the magnetic spacegroup of the 
C system is of type III (black-and-white type). The operation must then be classified 
C according to their keeping the z-axis invariant or not.
C
C srot(isym)%timeinv = boolean indicating if a time reversal operation is required
          IF(ifSP.AND.ifSO) THEN
           det=srot(isym)%krotm(1,1)*srot(isym)%krotm(2,2)-
     -        srot(isym)%krotm(1,2)*srot(isym)%krotm(2,1)
C the value of det is cos(srot(isym)%b) even if the rotation is improper.
           IF(det < 0.0d0) THEN
            srot(isym)%timeinv=.TRUE.
C The direction of the magnetic moment is changed to its opposite ( srot(isym)%b=pi ),
C A time reversal operation is required.
            srot(isym)%phase=srot(isym)%g-srot(isym)%a
C In this case, we define a phase factor for the off-diagonal term (up/dn term) 
C which is srot(isym)%phase= g-a = 2pi+(alpha-gamma)
           ELSE
            srot(isym)%timeinv=.FALSE.
C The direction of the magnetic moment is unchanged ( srot(isym)%b=0 ), 
C no time reversal operation is required. 
            srot(isym)%phase=srot(isym)%a+srot(isym)%g
C In this case, we define a phase factor for the off-diagonal term (up/dn term) 
C which is srot(isym)%phase= a+g = 2pi-(alpha+gamma)
           ENDIF
          ELSE
C If the calculation is either spin-polarized without spin-orbit, or paramagnetic
C the magnetic spacegroup of the system is of type I (ordinary type). The operation 
C are thus merely applied.
           srot(isym)%timeinv=.FALSE.
           srot(isym)%phase=0.d0
          ENDIF    ! End of the ifSP if-then-else
C
C
C Computation of the rotational matrices in each sort basis :
C -----------------------------------------------------------
          ALLOCATE(srot(isym)%rotrep(lsym,nsort))
C
C Initialization of the srot(isym)%rotrep field
C This field is a table of size (lsym*nsort) which contains the rotation matrices 
C of isym in the representation basis associated to each included orbital of each atom. 
C srot(isym)%rotrep = D(R[isym])_{new_i}
          DO isrt=1,nsort
            DO l=1,lsym
              ALLOCATE(srot(isym)%rotrep(l,isrt)%mat(1,1))
              srot(isym)%rotrep(l,isrt)%mat(1,1)=0.d0
            ENDDO
          ENDDO
C
C Computation of the elements 'mat' in srot(isym)%rotrep(l,isrt)
          DO isrt=1,nsort
            IF (notinclude(isrt)) cycle
            DO l=1,lsym
C The considered orbital is not included, hence no computation
              IF (lsort(l,isrt)==0) cycle
C The considered orbital is included
              IF (reptrans(l,isrt)%ifmixing) THEN
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) )
C If this option is used, then ifSO=.TRUE. (because of the restriction in set_ang_trans.f)
C Moreover ifSP=.TRUE. (since ifSO => ifSP in this version) 
               DEALLOCATE(srot(isym)%rotrep(l,isrt)%mat)
               ALLOCATE(srot(isym)%rotrep(l,isrt)%mat
     &           (1:2*(2*l+1),1:2*(2*l+1)))
               ALLOCATE(tmp_rot(1:2*(2*l+1),1:2*(2*l+1)))
               ALLOCATE(spinrot(1:2*(2*l+1),1:2*(2*l+1)))
               spinrot=0.d0
C Computation of the full spinor rotation matrix associated to isym.
               CALL spinrotmat(spinrot,isym,l)
C Computation of srot(isym)%rotrep(l,isrt)%mat
               tmp_rot(1:2*(2*l+1),1:2*(2*l+1))=MATMUL(
     &           reptrans(l,isrt)%transmat(1:2*(2*l+1),1:2*(2*l+1)),
     &           spinrot(1:2*(2*l+1),1:2*(2*l+1)))
               srot(isym)%rotrep(l,isrt)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =           MATMUL(tmp_rot(1:2*(2*l+1),1:2*(2*l+1)),
     &           TRANSPOSE(CONJG(reptrans(l,isrt)%transmat
     &           (1:2*(2*l+1),1:2*(2*l+1)))))
C the field srot(isym)%rotrep(l,isrt)%mat = (reptrans)*spinrot(l)*inverse(reptrans) 
C or srot(isym)%rotrep = D(R[isym])_{new_i} = <new_i|lm> D(R[isym])_{lm} <lm|new_i>
C which is exactly the expression of the spinor rotation matrix in the new basis.
               DEALLOCATE(tmp_rot)
               DEALLOCATE(spinrot)
              ELSE
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only)
               DEALLOCATE(srot(isym)%rotrep(l,isrt)%mat)
               ALLOCATE(srot(isym)%rotrep(l,isrt)%mat(-l:l,-l:l))
               ALLOCATE(tmp_rot(-l:l,-l:l))
C Computation of srot(isym)%rotrep(l,isrt)%mat
               tmp_rot(-l:l,-l:l)=MATMUL(
     &           reptrans(l,isrt)%transmat(-l:l,-l:l),
     &           srot(isym)%rotl(-l:l,-l:l,l))
               srot(isym)%rotrep(l,isrt)%mat(-l:l,-l:l)=
     =           MATMUL(tmp_rot(-l:l,-l:l),
     &           TRANSPOSE(CONJG(reptrans(l,isrt)%transmat(-l:l,-l:l))))
C the field srot(isym)%rotrep(l,isrt)%mat = (reptrans)*rotl*inverse(reptrans)
C or srot(isym)%rotrep = D(R[isym])_{new_i} = <new_i|lm> D(R[isym])_{lm} <lm|new_i>
C which is exactly the expression of the rotation matrix for the up/up block in the new basis.
               DEALLOCATE(tmp_rot)
              ENDIF
            ENDDO   ! End of the l loop
          ENDDO     ! End of the isrt loop
        ENDDO       ! End of the isym loop
C
C
C =============================================================
C Printing the matrix parameters in the file fort.17 for test :
C =============================================================
        DO isym=1,nsym
          WRITE(17,'()')
          WRITE(17,'(a,i3)')' Sym. op.: ',isym
          DO i =1,3
            ALLOCATE(bufreal(3))
            bufreal(1:3)=srot(isym)%krotm(i,1:3)
            WRITE(17,'(3f10.4)') bufreal
            DEALLOCATE(bufreal)
          ENDDO
          WRITE(17,'(a,3f8.1,i4)')'a, b, g, iprop =',
     &      srot(isym)%a*180d0/Pi,srot(isym)%b*180d0/Pi,
     &      srot(isym)%g*180d0/Pi,srot(isym)%iprop
C Printing the data relative to SP option
          IF (ifSP) THEN
           WRITE(17,*)'If DIR. magn. mom. is inverted :'
     &       ,srot(isym)%timeinv
           WRITE(17,*)'phase = ',srot(isym)%phase
          ENDIF
C Printing the rotational matrices for each orbital number l.
          WRITE(17,'()')
          DO l=1,lsym
            WRITE(17,'(a,a,i2)')'Rotation matrix ',
     &        'D(R[isym])_{lm} for l = ',l
            DO m=-l,l
              ALLOCATE(bufcomp(-l:l))
              bufcomp(-l:l)=srot(isym)%rotl(m,-l:l,l)
              WRITE(17,'(7(2f7.3,x))') bufcomp
              DEALLOCATE(bufcomp)
            ENDDO
          ENDDO
C Printing the matrices rotrep(l,isrt)%mat
          WRITE(17,'()') 
          DO isrt=1,nsort
            IF (notinclude(isrt)) cycle
            DO l=1,lsym
              IF (lsort(l,isrt)==0) cycle
              WRITE(17,'(a,i2,a,i2)')'Representation for isrt = ',
     &          isrt,' and l= ',l 
              IF (reptrans(l,isrt)%ifmixing) THEN
               DO m=1,2*(2*l+1)
                 ALLOCATE(bufcomp(1:2*(2*l+1)))
                 bufcomp(1:2*(2*l+1))=
     &             srot(isym)%rotrep(l,isrt)%mat(m,1:2*(2*l+1))
                 WRITE(17,'(7(2f7.3,x))') bufcomp
                 DEALLOCATE(bufcomp)
               ENDDO
              ELSE
               DO m=-l,l
                 ALLOCATE(bufcomp(-l:l))
                 bufcomp(-l:l)=
     &             srot(isym)%rotrep(l,isrt)%mat(m,-l:l)
                 WRITE(17,'(7(2f7.3,x))') bufcomp
                  DEALLOCATE(bufcomp)
               ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C
C =================================================================================
C Applying time-reversal operator if the system is spin-polarized with Spin Orbit :
C =================================================================================
C
C If the calculation is spin-polarized with spin-orbit, the magnetic spacegroup of the compound 
C is of type III (black-and-white). The symmetry operations which reverse the z-axis must be 
C multiplied by the time-reversal operator.
C If spin-orbit is not taken into account, all the field timeinv are .FALSE. and no time-reversal
C is applied, since the magnetic spacegroup of the compound is of type I (ordinary).
        IF (ifSP) THEN
C The modification of srot(isym)%rotl is done for each isym
         DO isym=1,nsym
           DO l=1,lsym
             IF (srot(isym)%timeinv) THEN
C The field srot(isym)%rotl is multiplied by the time-reversal operator in the complex basis.
              ALLOCATE(tmpcomp(-l:l,-l:l))
              tmpcomp(-l:l,-l:l)=
     &          srot(isym)%rotl(-l:l,-l:l,l)
              CALL timeinv_op(tmpcomp,(2*l+1),l,0)
              srot(isym)%rotl(-l:l,-l:l,l)=tmpcomp(-l:l,-l:l)
              DEALLOCATE(tmpcomp)
C The field srot(isym)%phase must not be modified.
             END IF
           END DO
         END DO
C
C The other modification are done for each (isrt,l) included.
         DO isrt=1,nsort
           IF (notinclude(isrt)) cycle
           DO l=1,lsym
C The considered orbital is not included, hence no computation
             IF (lsort(l,isrt)==0) cycle
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) )
             IF (reptrans(l,isrt)%ifmixing) THEN
              DO isym=1,nsym
                IF (srot(isym)%timeinv) THEN
C The field srot(isym)%rotrep(l,isrt)%mat is multiplied by the time-reversal operator in the corresponding basis of isrt.
                 CALL timeinv_op(srot(isym)%rotrep(l,isrt)%mat,
     &             2*(2*l+1),l,isrt)
                END IF
              END DO    ! End of the isym loop
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only)
             ELSE
              DO isym=1,nsym
                IF (srot(isym)%timeinv) THEN
C The field srot(isym)%rotrep(l,isrt)%mat is multiplied by the time-reversal operator in the corresponding basis of isrt.
                 CALL timeinv_op(srot(isym)%rotrep(l,isrt)%mat,
     &             (2*l+1),l,isrt)
                END IF
              END DO    ! End of the isym loop
             END IF     ! End of the ifmixing if-then-else
           END DO       ! End of the l loop
         END DO         ! End of the isrt loop 
        END IF          ! End of the ifSP if-then-else
C
C
C ======================================================================
C Printing the time-reversal modification in the file fort.17 for test :
C ======================================================================
        IF (ifSP.AND.ifSO) THEN
         WRITE(17,'()')
         WRITE(17,'(a)') '---With time-reversal operation---' 
         WRITE(17,'()')
C Printing the srot(isym) operations if necessary :
         DO isym=1,nsym
           IF (srot(isym)%timeinv) THEN
            WRITE(17,'()')
            WRITE(17,'(a,i3)')' Sym. op.: ',isym
C Printing the new rotational matrices for each orbital number l.
            WRITE(17,'()') 
            DO l=1,lsym
              WRITE(17,'(a,a,i2)')'T*Rotation matrix ',
     &          'D(T.R[isym])_{lm} for l = ',l
              DO m=-l,l
                ALLOCATE(bufcomp(-l:l))
                bufcomp(-l:l)=srot(isym)%rotl(m,-l:l,l)
                WRITE(17,'(7(2f7.3,x))') bufcomp
                DEALLOCATE(bufcomp)
              ENDDO
            ENDDO
C Printing the new matrices rotrep(l,isrt)%mat
            WRITE(17,'()') 
            DO isrt=1,nsort
              IF (notinclude(isrt)) cycle
              DO l=1,lsym
                IF (lsort(l,isrt)==0) cycle
                 WRITE(17,'(a,i2,a,i2)')
     &             'Representation for isrt = ',isrt,' and l= ',l 
                IF (reptrans(l,isrt)%ifmixing) THEN
                 DO m=1,2*(2*l+1)
                   ALLOCATE(bufcomp(1:2*(2*l+1)))
                   bufcomp(1:2*(2*l+1))=
     &               srot(isym)%rotrep(l,isrt)%mat(m,1:2*(2*l+1))
                   WRITE(17,'(7(2f7.3,x))') bufcomp
                   DEALLOCATE(bufcomp)
                 ENDDO
                ELSE
                 DO m=-l,l
                   ALLOCATE(bufcomp(-l:l))
                   bufcomp(-l:l)=
     &               srot(isym)%rotrep(l,isrt)%mat(m,-l:l)
                   WRITE(17,'(7(2f7.3,x))') bufcomp
                     DEALLOCATE(bufcomp)
                 END DO
                END IF
              END DO
            END DO
           END IF
         ENDDO
        END IF
C
C
C ============================================================
C Creation of the global->local coordinate rotation matrices :
C ============================================================
        ALLOCATE(rotloc(natom))
        CALL printout(1)
        WRITE(buf,'(a)')'-------------------------------------'
        CALL printout(0)
        WRITE(buf,'(a)')'Global-to-local-coordinates rotations'
        CALL printout(1)
        WRITE(buf,'(a)')'Properties of the symmetry operations :'
        CALL printout(0)
        WRITE(buf,'(a)') '   alpha, beta, gamma are their Euler angles.'
        CALL printout(0)
        WRITE(buf,'(a)') '   iprop is the value of their determinant.'
        CALL printout(0)
        CALL printout(0)
        WRITE(buf,'(a)')'  SORT    alpha      beta     gamma     iprop'
        CALL printout(0)
        READ(iusym,'()')
        DO isrt=1,nsort
C Reading the data for the representative atom in case.dmftsym and printing them in case.outdmftpr :
C --------------------------------------------------------------------------------------------------
          iatomref=SUM(nmult(0:isrt-1))+1
          READ(iusym,'()')
          DO i=1,3
            ALLOCATE(bufreal(3))
            READ(iusym,*) bufreal
            rotloc(iatomref)%krotm(i,1:3)=bufreal(1:3)
            DEALLOCATE(bufreal)
          ENDDO
C the field rotloc(iatomref)%krotm = 3x3 matrices of rotation associated to the transformation Rloc
C Rloc = <x_global | x_local >. The matrix was not multiplied by the value of iprop before being 
C written in case.dmftsym (cf. SRC_lapw2/rotmat_dmft.f). 
C rotloc(iatomref)%krotm can thus be either a proper or an improper rotation (with inversion).
C This reading line was chosen to be consistent with the writing line in rotmat_dmft (in SRC_lapw2)
          READ(iusym,*)rotloc(iatomref)%a,rotloc(iatomref)%b,
     &      rotloc(iatomref)%g, rotloc(iatomref)%iprop
          WRITE(buf,'(i5,3F10.1,5x,i3)')isrt,
     &      rotloc(iatomref)%a, rotloc(iatomref)%b,
     &      rotloc(iatomref)%g, rotloc(iatomref)%iprop
          CALL printout(0)
          rotloc(iatomref)%a=rotloc(iatomref)%a/180d0*Pi
          rotloc(iatomref)%b=rotloc(iatomref)%b/180d0*Pi
          rotloc(iatomref)%g=rotloc(iatomref)%g/180d0*Pi
C the field rotloc%a is linked to the Euler precession angle (alpha)
C the field rotloc%b is linked to the Euler nutation angle (beta)
C the field rotloc%c is linked to the Euler intrinsic rotation angle (gamma)
C They are read in case.dmftsym and printed in case.outdmftpr in degree and are then transformed into radians 
C the field rotloc%iprop = value of the transformation determinant (should be 1 in almost all the cases), 
C determines if there is an inversion in the transformation from global to local basis.
          rotloc(iatomref)%krotm(1:3,1:3)=rotloc(iatomref)%iprop*
     &      rotloc(iatomref)%krotm(1:3,1:3)
C Now, the field rotloc(iatomref)%krotm described only the proper rotation associated to the transformation.
C
C Use of the subroutine dmat to compute the rotational matrix 
C associated to the rotloc(iatomref) operation in a (2*l+1) orbital space :
          ALLOCATE(tmat(1:2*lsym+1,1:2*lsym+1))
          ALLOCATE(tmp_dmat(1:2*lsym+1,1:2*lsym+1,1:lsym))
          DO l=1,lsym
            tmat=0.d0
            CALL dmat(l,rotloc(iatomref)%a,rotloc(iatomref)%b,
     &        rotloc(iatomref)%g,REAL(rotloc(iatomref)%iprop,KIND=8),
     &        tmat,2*lsym+1)
            tmp_dmat(1:2*l+1,1:2*l+1,l)=tmat(1:2*l+1,1:2*l+1)
C tmp_dmat = D(Rloc)_{lm}
          ENDDO
          DEALLOCATE(tmat)
C
C
C Storing the rotloc matrix and initializing the other fields for all equivalent atoms :
C --------------------------------------------------------------------------------------
C All the equivalent atoms will have the same rotloc description. These data 
C will be correctly redifined in the subroutine set_rotloc, where the action of the 
C symmetry operation which transforms the representative atom in the considered one 
C will be added.
          DO imu=1,nmult(isrt)
            iatom=SUM(nmult(0:isrt-1))+imu
            IF(ifSP.AND.ifSO) THEN
C In this case, we have to consider the spinor rotation matrix associated to rotloc 
C (the value of the Euler angle beta can be anything between 0 and Pi)
             ALLOCATE(rotloc(iatom)%rotl(1:2*(2*lsym+1),
     &         1:2*(2*lsym+1),lsym))
             rotloc(iatom)%rotl=0.d0
             DO l=1,lsym
C For each orbital (from l=0 to lsym)
C Calculation of the representation matrix of rotloc in the spin-space
C in agreement with Wien conventions used for the definition of spmt (in SRC_lapwdm/sym.f)
C Up/up and Dn/dn terms
               factor=(rotloc(iatomref)%a+rotloc(iatomref)%g)/2.d0
               spmt(1,1)=EXP(CMPLX(0.d0,factor))
     &           *DCOS(rotloc(iatomref)%b/2.d0)
               spmt(2,2)=CONJG(spmt(1,1))
C Up/dn and Dn/up terms
               factor=-(rotloc(iatomref)%a-rotloc(iatomref)%g)/2.d0
               spmt(1,2)=EXP(CMPLX(0.d0,factor))
     &           *DSIN(rotloc(iatomref)%b/2.d0)
               spmt(2,1)=-CONJG(spmt(1,2))
C Up/up block :
               rotloc(iatom)%rotl(1:2*l+1,1:2*l+1,l)=
     &           spmt(1,1)*tmp_dmat(1:2*l+1,1:2*l+1,l)
C Dn/dn block :
               rotloc(iatom)%rotl(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1),l)=
     &           spmt(2,2)*tmp_dmat(1:2*l+1,1:2*l+1,l)
C Up/dn block :
               rotloc(iatom)%rotl(1:2*l+1,2*l+2:2*(2*l+1),l)=
     &           spmt(1,2)*tmp_dmat(1:2*l+1,1:2*l+1,l)
C Dn/up block :
               rotloc(iatom)%rotl(2*l+2:2*(2*l+1),1:2*l+1,l)=
     &           spmt(2,1)*tmp_dmat(1:2*l+1,1:2*l+1,l)
C The fields rotloc(iatom)%rotl now contain D(rotloc)_{lm}xD(rotloc)_{1/2}
             ENDDO
            ELSE
C In this case, we can consider the spatial rotation matrix only  
C since each spin space is independent (paramagnetic or spin-polarized without SO computation)
             ALLOCATE(rotloc(iatom)%rotl(-lsym:lsym,-lsym:lsym,lsym))
             rotloc(iatom)%rotl=0.d0
             DO l=1,lsym
               rotloc(iatom)%rotl(-l:l,-l:l,l)=
     =           tmp_dmat(1:2*l+1,1:2*l+1,l)
C The fields rotloc(iatom)%rotl now contain D(rotloc)_{lm}
             ENDDO
            ENDIF
C The fields rotloc(iatom)%a,b and c will now contain the parameters linked to 
C the Euler angles of the local rotation rotloc.
            IF(imu.gt.1) THEN
             rotloc(iatom)%a=rotloc(iatomref)%a
             rotloc(iatom)%b=rotloc(iatomref)%b
             rotloc(iatom)%g=rotloc(iatomref)%g
             rotloc(iatom)%iprop=rotloc(iatomref)%iprop
             rotloc(iatom)%krotm(1:3,1:3)=
     =         rotloc(iatomref)%krotm(1:3,1:3)
            ENDIF
C The fields rotloc%phase, timeinv and srotnum are initialized to their 
C default value.
            rotloc(iatom)%phase=0.d0
            rotloc(iatom)%timeinv=.FALSE.
            rotloc(iatom)%srotnum=0
C the field rotloc(iatom)%srotnum and timeinv will be recalculated in set_rotloc.
          ENDDO
          DEALLOCATE(tmp_dmat)
        ENDDO     ! End of the isrt loop
C
C
C ====================================================================
C Printing the rotloc matrix parameters in the file fort.17 for test :
C ====================================================================
        DO isrt=1,nsort
          IF (notinclude(isrt)) cycle
          DO imu=1,nmult(isrt)
            iatom=SUM(nmult(0:isrt-1))+imu
            WRITE(17,'()')
            WRITE(17,'(2(a,i3))')' SORT ',isrt,' IMU= ',imu
            DO i=1,3
              ALLOCATE(bufreal(3))
              bufreal(1:3)=rotloc(iatom)%krotm(i,1:3)
              WRITE(17,'(3f10.4)') bufreal
              DEALLOCATE(bufreal)
            ENDDO
            WRITE(17,'(a,3f8.1,i4)')'a, b, g, iprop ==',
     &        rotloc(iatom)%a*180d0/Pi,rotloc(iatom)%b*180d0/Pi,
     &        rotloc(iatom)%g*180d0/Pi,rotloc(iatom)%iprop
C Printing the data relative to SP option
            IF (ifSP) THEN
              WRITE(17,*)'If DIR. magn. mom. is inverted :'
     &          ,rotloc(iatom)%timeinv
              WRITE(17,*)'phase = ',rotloc(iatom)%phase
            ENDIF
C Printing the rotloc matrices for each orbital number l.
            WRITE(17,'()')
            DO l=1,lsym
            WRITE(17,'(a,a,i2)')'Rotation matrix ',
     &        'D(R[isym])_{lm} for l = ',l 
            IF(ifSP.AND.ifSO) THEN
              DO m=1,2*(2*l+1)
                ALLOCATE(bufcomp(1:2*(2*l+1)))
                bufcomp(1:2*(2*l+1))=rotloc(iatom)%rotl(m,1:2*(2*l+1),l)
                WRITE(17,'(7(2f7.3,x))') bufcomp
                DEALLOCATE(bufcomp)
              ENDDO
             ELSE
              DO m=-l,l
                ALLOCATE(bufcomp(-l:l))
                bufcomp(-l:l)=rotloc(iatom)%rotl(m,-l:l,l)
                WRITE(17,'(7(2f7.3,x))') bufcomp
                DEALLOCATE(bufcomp)
              ENDDO
             ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C
C ==================================================================================
C Computation of the true local rotation matrices for each non representative atom : 
C ==================================================================================
        CALL set_rotloc 
C
C
C ====================================================================
C Printing the rotloc matrix parameters in the file fort.17 for test :
C ====================================================================
        DO isrt=1,nsort
          IF (notinclude(isrt)) cycle
          DO imu=1,nmult(isrt)
            iatom=SUM(nmult(0:isrt-1))+imu
            WRITE(17,'()')
            WRITE(17,'(2(a,i3))')' SORT ',isrt,' IMU= ',imu
            DO i=1,3
              ALLOCATE(bufreal(3))
              bufreal(1:3)=rotloc(iatom)%krotm(i,1:3)
              WRITE(17,'(3f10.4)') bufreal
              DEALLOCATE(bufreal)
            ENDDO
            WRITE(17,'(a,3f8.1,i4)')'a, b, g, iprop ==',
     &        rotloc(iatom)%a*180d0/Pi,rotloc(iatom)%b*180d0/Pi,
     &        rotloc(iatom)%g*180d0/Pi,rotloc(iatom)%iprop
C Printing the data relative to SP option
            IF (ifSP) THEN
              WRITE(17,*)'If DIR. magn. mom. is inverted :'
     &          ,rotloc(iatom)%timeinv
              WRITE(17,*)'phase = ',rotloc(iatom)%phase
            ENDIF
C Printing the rotloc matrices for each orbital number l.
            WRITE(17,'()')
            DO l=1,lsym
            WRITE(17,'(a,a,i2)')'Rotation matrix ',
     &        'D(R[isym])_{lm} for l = ',l 
            IF(ifSP.AND.ifSO) THEN
             DO m=1,2*(2*l+1)
               ALLOCATE(bufcomp(1:2*(2*l+1)))
               bufcomp(1:2*(2*l+1))=rotloc(iatom)%rotl(m,1:2*(2*l+1),l)
               WRITE(17,'(7(2f7.3,x))') bufcomp
               DEALLOCATE(bufcomp)
             ENDDO
            ELSE
             DO m=-l,l
               ALLOCATE(bufcomp(-l:l))
               bufcomp(-l:l)=rotloc(iatom)%rotl(m,-l:l,l)
               WRITE(17,'(7(2f7.3,x))') bufcomp
               DEALLOCATE(bufcomp)
             ENDDO
            ENDIF
            ENDDO
C Printing the matrices rotrep(l)%mat
            WRITE(17,'()')
            DO l=1,lsym
              IF (lsort(l,isrt)==0) cycle
              WRITE(17,'(a,i2)')'Representation for l= ',l 
              IF (ifSP.AND.ifSO) THEN
               DO m=1,2*(2*l+1)
                 ALLOCATE(bufcomp(1:2*(2*l+1)))
                 bufcomp(1:2*(2*l+1))=
     &             rotloc(iatom)%rotrep(l)%mat(m,1:2*(2*l+1))
                 WRITE(17,'(7(2f7.3,x))') bufcomp
                 DEALLOCATE(bufcomp)
               ENDDO
              ELSE
               DO m=-l,l
                 ALLOCATE(bufcomp(-l:l))
                 bufcomp(-l:l)=
     &             rotloc(iatom)%rotrep(l)%mat(m,-l:l)
                 WRITE(17,'(7(2f7.3,x))') bufcomp
                 DEALLOCATE(bufcomp)
               ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C                
        RETURN
        END


 

	Subroutine dmat(l,a,b,c,det,DD,length)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine computes the inverse of the matrix of the       %%
C %% representation of size (2*l+1) associated to the rotation       %%
C %% described by (a,b,c) angles in Euler description and with       %% 
C %% determinant det.                                                %%
C %% The obtained matrix is put in the variable DD.                  %% 
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER l,m,n,ifac,length
        COMPLEX*16 izero,imag, dd
        dimension DD(length,length)                       
	imag=(0d0,1d0)
	izero=(0d0,0d0)
	pi=acos(-1d0)

	do m=-l,l
	  do n=-l,l
	    call d_matrix(l,m,n,b,dm)
	    if (det.lt.-0.5) then
             dd(l+m+1,n+l+1)=(-1)**l*cdexp(imag*n*a)
     &        *cdexp(imag*m*c)*dm
	    else
              dd(l+m+1,n+l+1)=cdexp(imag*n*a)
     &          *cdexp(imag*m*c)*dm
	    end if
 3          format(2I3,2f10.6)
	  end do
	end do
	do j=1,2*l+1
	end do
 5      format(7(2f6.3,1X))

	end


	Subroutine d_matrix(l,m,n,b,dm)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine is called by the subroutine dmat to compute the %%
C %% the value of the coefficient dm.                                %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER l,m,n,t

 	sum=0d0
        
        f1=dfloat(ifac(l+m)*ifac(l-m))/  
     &     dfloat(ifac(l+n)*ifac(l-n))

	do t=0,2*l
          if ((l-m-t).ge.0.AND.(l-n-t).ge.0.AND.(t+n+m).ge.0) then
C general factor
           f2=dfloat(ifac(l+n)*ifac(l-n))/dfloat(ifac(l-m-t) 
     &       *ifac(m+n+t)*ifac(l-n-t)*ifac(t))
C factor with sin(b/2)
           if ((2*l-m-n-2*t).eq.0) then
            f3=1.
	   else
	    f3=(sin(b/2))**(2*l-m-n-2*t)
	   end if
C factor with cos(b/2)
	   if ((2*t+n+m).eq.0) then
	    f4=1.
	   else
	    f4=(cos(b/2))**(2*t+n+m)
	   end if
!	write(12,*)f1,f2,f3,f4
           sum=sum+(-1)**(l-m-t)*f2*f3*f4
          end if
	end do

	dm=sqrt(f1)*sum
	end


	Integer Function ifac(n)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine computes the factorial of the number n          %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
	if (n.eq.0) then
	 ifac=1
	else
	 ifac=1
	 do j=1,n
	   ifac=ifac*j
	 end do
	end if
	end


       SUBROUTINE spinrotmat(spinrot,isym,l)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine sets up the complete spinor rotation matrix     %%
C %% associated to the symmetry operation isym for the orbital l.    %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definition of the variables :
C -----------------------------
       USE common_data
       USE symm
       IMPLICIT NONE
       INTEGER :: l,isym
       COMPLEX(KIND=8) :: ephase, det
       REAL(KIND=8) :: factor
       COMPLEX(KIND=8), DIMENSION(1:2*(2*l+1),1:2*(2*l+1)) :: spinrot
       COMPLEX(KIND=8), DIMENSION(1:2,1:2) :: spmt
C
       spinrot=0.d0
C For a computation with spin polarized inputs :
       IF (ifSP) THEN
        IF (srot(isym)%timeinv) THEN
C In this case, the Euler angle Beta is Pi. The spinor rotation matrix is block-antidiagonal and
C the time reversal operation will be applied to keep the direction of the magnetization.
C Up/dn block :
         factor=srot(isym)%phase/2.d0
C We remind that the field phase is (g-a) in this case.
C as a result, ephase = exp(+i(g-a)/2) = -exp(+i(alpha-gamma)/2)
C in good agreement with Wien conventions for the definition of this phase factor.
         ephase=EXP(CMPLX(0.d0,factor))
         spinrot(1:2*l+1,2*l+2:2*(2*l+1))=
     =     ephase*srot(isym)%rotl(-l:l,-l:l,l)
C Dn/up block :
         ephase=-CONJG(ephase)
C now, ephase = -exp(+i(a-g)/2) = exp(-i(alpha-gamma)/2)
         spinrot(2*l+2:2*(2*l+1),1:2*l+1)=
     =     ephase*srot(isym)%rotl(-l:l,-l:l,l)
        ELSE
C In this case, the Euler angle Beta is 0. The spinor rotation matrix is block-diagonal and
C no time reversal operation will be applied.
C Up/up block :
         factor=srot(isym)%phase/2.d0
C We remind that the field phase is (a+g) in this case.
C as a result, ephase = exp(+i(a+g)/2)=-exp(-i(alpha+gamma)/2)
C in good agreement with Wien conventions for the definition of this phase factor.
         ephase=EXP(CMPLX(0.d0,factor))
         spinrot(1:2*l+1,1:2*l+1)=
     =     ephase*srot(isym)%rotl(-l:l,-l:l,l)
C Dn/dn block :
         ephase=CONJG(ephase)
C now, ephase = exp(-i(a+g)/2) = -exp(+i(alpha+gamma)/2)
         spinrot(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))=
     =     ephase*srot(isym)%rotl(-l:l,-l:l,l)
        ENDIF
       ELSE
C For a computation with paramagnetic treatment input files. (not used in this version)
C
C In this case, there is no restriction on the value of the Euler angle beta.
C The general definition of a spinor rotation matrix is used.
C
C Calculation of the representation matrix of isym in the spin-space
C in agreement with Wien conventions used for the definition of spmt (in SRC_lapwdm/sym.f)
C Up/up and Dn/dn terms
        factor=(srot(isym)%a+srot(isym)%g)/2.d0
        spmt(1,1)=EXP(CMPLX(0.d0,factor))*DCOS(srot(isym)%b/2.d0)
        spmt(2,2)=CONJG(spmt(1,1))
C Up/dn and Dn/up terms
        factor=-(srot(isym)%a-srot(isym)%g)/2.d0
        spmt(1,2)=EXP(CMPLX(0.d0,factor))*DSIN(srot(isym)%b/2.d0)
        spmt(2,1)=-CONJG(spmt(1,2))
C Up/up block :
        spinrot(1:2*l+1,1:2*l+1)=
     &    spmt(1,1)*srot(isym)%rotl(-l:l,-l:l,l)
C Dn/dn block :
        spinrot(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))=
     &    spmt(2,2)*srot(isym)%rotl(-l:l,-l:l,l)
C Up/dn block :
        spinrot(1:2*l+1,2*l+2:2*(2*l+1))=
     &    spmt(1,2)*srot(isym)%rotl(-l:l,-l:l,l)
C Dn/up block :
        spinrot(2*l+2:2*(2*l+1),1:2*l+1)=
     &    spmt(2,1)*srot(isym)%rotl(-l:l,-l:l,l)
       ENDIF
C
       RETURN
       END

