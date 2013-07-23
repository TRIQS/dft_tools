
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

      SUBROUTINE symmetrize_mat(Dmat,orbit,norbit)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine applies the symmetry operations to the          %%
C %% density matrices stored in Dmat and puts the resulting          %%
C %% density matrices into the local coordinate system.              %%
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
      COMPLEX(KIND=8),DIMENSION(:,:,:,:), ALLOCATABLE :: sym_dmat
      COMPLEX(KIND=8),DIMENSION(:,:), ALLOCATABLE :: tmp_mat
      COMPLEX(KIND=8):: ephase
      TYPE(orbital), DIMENSION(norbit) :: orbit
      INTEGER :: isym, iorb, iatom, jorb, is, is1, l, i, m
      INTEGER :: isrt, jatom, imult, asym
C
C =========================================
C Computation of the symmetrized matrices :
C =========================================
C
        iorb=1
C Initialization of the iorb index.
        DO WHILE (iorb.lt.(norbit+1)) 
C The use of the while-loop was motivated by the idea of studying 
C all the orbitals iorb associated to a same atomic sort isrt together.
C At the end, the index iorb is incremented by nmult(isrt) so that the 
C following studied orbitals are associated to another atomic sort. 
          l=orbit(iorb)%l
          isrt=orbit(iorb)%sort
C
C -----------------------------------------------------------------------------------
C The s-orbitals are a particular case of a "non-mixing" basis and are treated here :
C -----------------------------------------------------------------------------------
          IF (l==0) THEN
C The table sym_dmat will store the symmetrized value of the density matrices of Dmat
C associated to a same atomic sort isrt.
           ALLOCATE(sym_dmat(1,1,nsp,1:nmult(isrt)))
           sym_dmat=0.d0
C Its size is nmult(isrt) because symmetry operations can transform the representants 
C of a same atomic sort one into another.
C
C Loop on the representants of the atomic sort isrt
           DO imult=0,nmult(isrt)-1
             iatom=orbit(iorb+imult)%atom
C Loop on the symmetry operations of the system
             DO isym=1,nsym
               DO is=1,nsp
                 ALLOCATE(tmp_mat(1,1))
C If the calculation uses spin-polarized input files, the application of the symmetry operation
C depends on the field srot%timeinv.
                 IF(ifSP.AND.srot(isym)%timeinv) THEN
C In this case (spin-polarized computation), the symmetry operation is block-diagonal in spin-space but 
C the time reversal operator is included.
                  tmp_mat(1,1)=CONJG(Dmat(is,iorb+imult)%mat(1,1))
C because of the antiunitarity of the operator, the conjugate of Dmat must be use
                 ELSE 
                  tmp_mat(1,1)=Dmat(is,iorb+imult)%mat(1,1)
                 ENDIF
C
C Definition of the index where the transformed Dmat will be stored. [jorb = R[isym](iorb)] 
                 jorb=srot(isym)%perm(iatom)-iatom+(imult+1)
C
C Computation of the phase factors in the case of a SO computation :
C ------------------------------------------------------------------
C For up/up and dn/dn blocks, no phase factor is needed.
                 ephase=1.d0
C For the up/dn block, initialisation of the phase factor
                 IF(is==3) THEN
                  ephase=EXP(CMPLX(0d0,srot(isym)%phase))
C if srot%timeinv = .TRUE. , phase= g-a = 2pi+(alpha-gamma) and ephase = exp(+i(g-a)) = exp(+i(alpha-gamma)) 
C if srot%timeinv = .FALSE., phase= a+g = 2pi-(alpha+gamma) and ephase = exp(+i(a+g)) = exp(-i(alpha+gamma))
                 ENDIF
C For the dn/up block, initialisation of the phase factor
                 IF(is==4) THEN
                  ephase=EXP(CMPLX(0d0,-srot(isym)%phase))
C if srot%timeinv = .TRUE. , phase= g-a = 2pi+(alpha-gamma) and ephase = exp(-i(g-a)) = exp(-i(alpha-gamma)) 
C if srot%timeinv = .FALSE., phase= a+g = 2pi-(alpha+gamma) and ephase = exp(-i(a+g)) = exp(+i(alpha+gamma))
                 ENDIF
C
C Application of the symmetry operation which changes iorb in jorb=R[isym](iorb) :
C --------------------------------------------------------------------------------
C That's why the result is stored in the jorb section of sym_dmat.
                 sym_dmat(1,1,is,jorb)=
     =             sym_dmat(1,1,is,jorb)+tmp_mat(1,1)*ephase
                 DEALLOCATE(tmp_mat)
               ENDDO      ! End of the is loop
             ENDDO        ! End of the isym loop
           ENDDO          ! End of the imult loop
C
C Renormalization of the symmetrized density matrices :
C -----------------------------------------------------
           IF (nsym.gt.0) THEN
            DO imult=0,nmult(isrt)-1
              DO is=1,nsp
                Dmat(is,iorb+imult)%mat(1,1)=
     &            sym_dmat(1,1,is,imult+1)/nsym
              ENDDO
            ENDDO
           ENDIF
           DEALLOCATE(sym_dmat)
C Incrementation of the iorb index (for the while loop)
           iorb=iorb+nmult(isrt)
C
C -----------------------------------------------------------------------------------------------------
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) ) :
C -----------------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C The table sym_dmat will store the symmetrized value of the density matrices of Dmat
C associated to a same atomic sort isrt.
           ALLOCATE(sym_dmat(1:2*(2*l+1),1:2*(2*l+1),1,1:nmult(isrt)))
           sym_dmat=0.d0
C Its size is nmult(isrt) because symmetry operations can transform the representants 
C of a same atomic sort one into another.
C
C Loop on the representants of the atomic sort isrt
           DO imult=0,nmult(isrt)-1
             iatom=orbit(iorb+imult)%atom
C Loop on the symmetry operations of the system
             DO isym=1,nsym
               ALLOCATE(tmp_mat(1:2*(2*l+1),1:2*(2*l+1)))
C We use the complete spin-space representation, so no trick on indices is necessary.
               tmp_mat(1:2*(2*l+1),1:2*(2*l+1))=
     &           Dmat(1,iorb+imult)%mat(1:2*(2*l+1),1:2*(2*l+1))
C If the calculation is spin-polarized, the symmetry operator is antiunitary.
               IF(ifSP.AND.srot(isym)%timeinv) THEN
                tmp_mat(1:2*(2*l+1),1:2*(2*l+1))= 
     &            CONJG(tmp_mat(1:2*(2*l+1),1:2*(2*l+1)))
               ENDIF
C Definition of the index where the transformed Dmat will be stored. [jorb = R[isym](iorb)] 
               jorb=srot(isym)%perm(iatom)-iatom+(imult+1)
C Application of the symmetry operation :
C ---------------------------------------
C The transformation is :   srot%rotrep.tmpmat(iorb).inverse(sort%rotrep) = Dmat(jorb)
C or in other words, if R is a simple symmetry    D(R[isym]) tmpmat(iorb) D(inverse(R[isym]))  = Dmat(R[isym](iorb))
C                    if R is multiplied by Theta  D(R[isym]) tmpmat(iorb)* D(inverse(R[isym]))*  = Dmat(R[isym](iorb))
               tmp_mat(1:2*(2*l+1),1:2*(2*l+1))=
     =           MATMUL(tmp_mat(1:2*(2*l+1),1:2*(2*l+1)),
     &           TRANSPOSE(CONJG(srot(isym)%rotrep(l,isrt)
     %           %mat(1:2*(2*l+1),1:2*(2*l+1)) )) ) 
               sym_dmat(1:2*(2*l+1),1:2*(2*l+1),1,jorb)=
     =           sym_dmat(1:2*(2*l+1),1:2*(2*l+1),1,jorb)+ 
     +           MATMUL( srot(isym)%rotrep(l,isrt)
     %           %mat(1:2*(2*l+1),1:2*(2*l+1)) ,
     &           tmp_mat(1:2*(2*l+1),1:2*(2*l+1)))
               DEALLOCATE(tmp_mat)
             ENDDO       ! End of the isym loop        
           ENDDO         ! End of the imult loop
C Renormalization of the symmetrized density matrices :
C -----------------------------------------------------
           IF (nsym.gt.0) THEN
            DO imult=0,nmult(isrt)-1
              Dmat(1,iorb+imult)%mat(1:2*(2*l+1),1:2*(2*l+1))=
     =          sym_dmat(1:2*(2*l+1),1:2*(2*l+1),1,imult+1)/nsym
            ENDDO
           ENDIF
           DEALLOCATE(sym_dmat)
C Incrementation of the iorb index (for the while loop)
           iorb=iorb+nmult(isrt)
C
C ----------------------------------------------------------------------------------------------
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only) :
C ----------------------------------------------------------------------------------------------
          ELSE
C The table sym_dmat will store the symmetrized value of the density matrices of Dmat
C associated to a same atomic sort isrt.
           ALLOCATE(sym_dmat(-l:l,-l:l,nsp,1:nmult(isrt)))
           sym_dmat=0.d0
C Its size is nmult(isrt) because symmetry operations can transform the representants 
C of a same atomic sort one into another.
C 
C Loop on the representants of the atomic sort isrt
           DO imult=0,nmult(isrt)-1
             iatom=orbit(iorb+imult)%atom
C Loop on the symmetry operations of the system
             asym=0
             DO isym=1,nsym
               DO is=1,nsp
                 ALLOCATE(tmp_mat(-l:l,-l:l))
C If the calculation uses spin-polarized input files, the application of the symmetry operation
C depends on the field srot%timeinv.
                 IF(ifSP.AND.srot(isym)%timeinv) THEN
C In this case (spin-polarized computation), the symmetry operation is block-diagonal in spin-space but 
C the time reversal operatot is included.
                  tmp_mat(-l:l,-l:l)=CONJG(
     &               Dmat(is,iorb+imult)%mat(-l:l,-l:l))
C because of antiunitarity of the operator, the conjugate of Dmat must be use
                 ELSE 
                  tmp_mat(-l:l,-l:l)=
     &              Dmat(is,iorb+imult)%mat(-l:l,-l:l)
                 ENDIF
C
C Definition of the index where the transformed Dmat will be stored. [jorb = R[isym](iorb)] 
                 jorb=srot(isym)%perm(iatom)-iatom+(imult+1)
C
C Computation of the phase factors in the case of a SO computation :
C ------------------------------------------------------------------
C For up/up and dn/dn blocks, no phase factor is needed.
                 ephase=1.d0
C For the up/dn block, initialisation of the phase factor
                 IF(is==3) THEN
                  ephase=EXP(CMPLX(0d0,srot(isym)%phase))
C if srot%timeinv = .TRUE. , phase= g-a = 2pi+(alpha-gamma) and ephase = exp(+i(g-a)) = exp(+i(alpha-gamma)) 
C if srot%timeinv = .FALSE., phase= a+g = 2pi-(alpha+gamma) and ephase = exp(+i(a+g)) = exp(-i(alpha+gamma))
                 ENDIF
C For the dn/up block, initialisation of the phase factor
                 IF(is==4) THEN
                  ephase=EXP(CMPLX(0d0,-srot(isym)%phase))
C if srot%timeinv = .TRUE. , phase= g-a = 2pi+(alpha-gamma) and ephase = exp(-i(g-a)) = exp(-i(alpha-gamma)) 
C if srot%timeinv = .FALSE., phase= a+g = 2pi-(alpha+gamma) and ephase = exp(-i(a+g)) = exp(+i(alpha+gamma))
                 ENDIF
C
C Application of the symmetry operation which changes iorb in jorb :
C ------------------------------------------------------------------
C The transformation is :   srot%rotrep.tmpmat(iorb).inverse(sort%rotrep) = Dmat(jorb)
C or in other words, if R is a simple symmetry    D(R[isym]) tmpmat(iorb) D(inverse(R[isym]))  = Dmat(R[isym](iorb))
C                    if R is multiplied by T      D(R[isym]) tmpmat(iorb)* D(inverse(R[isym]))*  = Dmat(R[isym](iorb))
                 tmp_mat(-l:l,-l:l)=
     =             MATMUL(tmp_mat(-l:l,-l:l),
     &             TRANSPOSE(CONJG( srot(isym)
     &             %rotrep(l,isrt)%mat(-l:l,-l:l) )) ) 
                 sym_dmat(-l:l,-l:l,is,jorb)=
     =             sym_dmat(-l:l,-l:l,is,jorb)+ 
     +             MATMUL(srot(isym)%rotrep(l,isrt)%mat(-l:l,-l:l),
     &             tmp_mat(-l:l,-l:l) )*ephase
                 DEALLOCATE(tmp_mat)
               ENDDO       ! End of the is loop         
             ENDDO         ! End of the isym loop
           ENDDO           ! End of the imult loop
C
C Renormalization of the symmetrized density matrices :
C -----------------------------------------------------
           IF (nsym.gt.0) THEN
            DO imult=0,nmult(isrt)-1
              DO is=1,nsp
                Dmat(is,iorb+imult)%mat(-l:l,-l:l)=
     =            sym_dmat(-l:l,-l:l,is,imult+1)/(nsym-asym)
              ENDDO
            ENDDO
           ENDIF
           DEALLOCATE(sym_dmat)
C Incrementation of the iorb index (for the while loop)
           iorb=iorb+nmult(isrt)
          ENDIF     ! End of the type basis if-then-else
C
        ENDDO       ! End of the while(iorb) loop
C
C
C =============================================================
C Application of the time reversal operation if paramagnetism :
C =============================================================
C If the system is paramagnetic, the magnetic group of the system 
C is a type II Shubnikov group and time-reveral symmetry must be added 
C to achieve the complete symmetrization.
        IF (.not.ifSP) THEN
         CALL add_timeinv(Dmat,orbit,norbit)
        END IF
C
        RETURN
        END
