
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

        SUBROUTINE outqmc(elecn,qbbot)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine creates the output files :                      %% 
C %%   - case.ctqmcout, with all the informations necessary for a    %% 
C %%     CTQMC computation.                                          %%
C %%   - case.symqmc, describing all the symmetries of the system    %%
c %%     (necessary for CTQMC computation too).                      %%
C %%   - case.parproj which gives the partial charge projectors for  %%
C %%     orbitals and for all atoms.                                 %%
C %%   - case.sympar which contains the symmetry matrices for all    %%
C %%     the included orbitals (for partial charge analysis)         %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definition of the variables :
C -----------------------------
        USE almblm_data
        USE common_data
        USE file_names
        USE prnt
        USE reps
        USE symm
        USE projections
        IMPLICIT NONE
C
        INTEGER :: iorb, icrorb, irep, isrt
        INTEGER :: l, m, is, i1, i2, isym
        INTEGER :: ik, il, ib, ir, n
        INTEGER :: ind1, ind2, iatom
        INTEGER :: timeinvflag
        REAL(KIND=8) ::  qbbot, elecn, factor
        COMPLEX(KIND=8) ::  ephase
        COMPLEX(KIND=8),DIMENSION(:,:), ALLOCATABLE :: hk 
        COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: spinrot
        COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: densprint
        COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: time_op
C
C
C =======================================
C     Writing the file case.ctqmcout    :
C =======================================
C
        WRITE(buf,'(a)')'Writing the file case.ctqmcout...'
        CALL printout(0)
C
C Definition of 1 electron-Volt
        WRITE(ouctqmc,'(a)') '13.605698'     
C
C ---------------------------------------
C General informations about the system :
C ---------------------------------------
C
C Number of k-points in the I-BZ
        WRITE(ouctqmc,'(i6)') nk
C Definition of the spin-polarized flag ifSP
        IF (ifSP) THEN
         WRITE(ouctqmc,'(i6)') 1
        ELSE
         WRITE(ouctqmc,'(i6)') 0
        ENDIF
C Definition of the Spin-orbit flag ifSO
        IF (ifSO) THEN
         WRITE(ouctqmc,'(i6)') 1
        ELSE
         WRITE(ouctqmc,'(i6)') 0
        ENDIF
C The only possible combinations are :
C - (0,0) which stands for a paramagnetic computation without SO.
C - (1,0) which stands for computation without SO using spin-polarized input files.
C - (1,1) which stands for computation with SO using spin-polarized input files. 
C
C Writing the total charge below the lower limit of the energy window (variable "qbbot")
        WRITE(ouctqmc,*) qbbot
C Writing the total number of electrons (valence band+semicore) (variable "elecn"). 
C It is also the charge upto the Fermi level.
        WRITE(ouctqmc,*) elecn
C
C ------------------------------------------
C Description of all the included orbitals :
C ------------------------------------------ 
C
C Definition of the number of included orbitals "norb"
        WRITE(ouctqmc,'(i6)') norb
C Description of each orbital "iorb"
        DO iorb=1,norb
          IF (ifSO) THEN
          WRITE(ouctqmc,'(4(i6,x))') orb(iorb)%atom, orb(iorb)%sort, 
     &      orb(iorb)%l, 2*(2*orb(iorb)%l+1)
          ELSE
          WRITE(ouctqmc,'(4(i6,x))') orb(iorb)%atom, orb(iorb)%sort, 
     &      orb(iorb)%l, 2*orb(iorb)%l+1
          ENDIF
        ENDDO
C an orbital "iorb" is described by :
C  - the associated atom
C  - the corresponding atomic sort
C  - the considered orbital number l
C  - the size of the corresponding matrices : 2*l+1 without SO ; 2*(2*l+1) with SO 
C
C ----------------------------------------
C Description of the correlated orbitals :
C ---------------------------------------- 
C
C Definition of the number of correlated orbitals "ncrorb"
        WRITE(ouctqmc,'(i6)') ncrorb
C Description of each correlated orbital "icrorb"
        DO icrorb=1,ncrorb 
          l=crorb(icrorb)%l
          isrt=crorb(icrorb)%sort
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
          IF (l==0) THEN
C For the s-orbitals, the only irep possible is the matrix itself.
           IF (ifSP.AND.ifSO) THEN
C If SO is taken into account, spinor rotation matrix is considered. 
C The spin is not a good quantum number, so the whole representation is used. 
            WRITE(ouctqmc,'(6(i6,x))') crorb(icrorb)%atom, isrt, 0,
     &        2, crorb(icrorb)%ifSOat, 1
           ELSE
C Without SO, only the rotation matrix in orbital space is necessary.
            WRITE(ouctqmc,'(6(i6,x))') crorb(icrorb)%atom, isrt, 0,
     &        1, crorb(icrorb)%ifSOat, 1
           ENDIF
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
          ELSEIF(reptrans(l,isrt)%ifmixing) THEN
C In this case, the SO is necessary considered, spinor rotation matrices are used.
           IF (crorb(icrorb)%ifsplit) THEN
C If only an irep is correlated
            DO irep=1,reptrans(l,isrt)%nreps
              IF (crorb(icrorb)%correp(irep)) THEN
               WRITE(ouctqmc,'(6(i6,x))') crorb(icrorb)%atom, 
     &           isrt, l, reptrans(l,isrt)%dreps(irep),
     &           crorb(icrorb)%ifSOat, irep
              ENDIF
            ENDDO
           ELSE
C If no particular irep is correlated
            WRITE(ouctqmc,'(6(i6,x))') crorb(icrorb)%atom, isrt, l, 
     &       2*(2*l+1), crorb(icrorb)%ifSOat, 1
           ENDIF
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
          ELSE
           IF (ifSP.AND.ifSO) THEN 
C If SO is taken into account, spinor rotation matrices are considered.
C The spin is not a good quantum number, so the whole representation is used. 
            WRITE(ouctqmc,'(6(i6,x))') crorb(icrorb)%atom, isrt, l,
     &        2*(2*l+1), crorb(icrorb)%ifSOat, 1
           ELSE
C Without SO, only the rotation matrix in orbital space is necessary.
            IF (crorb(icrorb)%ifsplit) THEN
C If only an irep is correlated
             DO irep=1,reptrans(l,isrt)%nreps
               IF (crorb(icrorb)%correp(irep)) THEN
                WRITE(ouctqmc,'(6(i6,x))') crorb(icrorb)%atom,
     &            isrt, l, reptrans(l,isrt)%dreps(irep),
     &            crorb(icrorb)%ifSOat, irep
               ENDIF
             ENDDO
            ELSE
C If no particular irep is correlated
             WRITE(ouctqmc,'(6(i6,x))') crorb(icrorb)%atom, isrt, l, 
     &        (2*l+1), crorb(icrorb)%ifSOat, 1
            END IF    ! End of the ifsplit if-then-else
           END IF     ! End of the ifSO if-then-else
          END IF      ! End of the ifmixing if-then-else
        END DO        ! End of the icrorb loop
C an orbital "iorb" is described by :
C  - the associated atom
C  - the corresponding atomic sort
C  - the considered orbital number l
C  - the size of the "correlated" submatrix (can be the whole matrix) 
C  - the flag ifSOat which states that SO is considered for this orbital
C  - the number of the irep
C
C ------------------------------------------------------------------------------------------------
C Description of the global to local coordinates transformation Rloc for each correlated orbital :
C ------------------------------------------------------------------------------------------------
C Description of each transformation Rloc
        DO icrorb=1,ncrorb
          l=crorb(icrorb)%l
          isrt=crorb(icrorb)%sort
          iatom=crorb(icrorb)%atom
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
          IF(l==0) THEN
C For the s-orbitals, the only irep possible is the matrix itself.
           IF (ifSP.AND.ifSO) THEN 
C If SO is taken into account, spinor rotation matrix must be considered.
            ALLOCATE(spinrot(1:2,1:2))
            spinrot(:,:)=0.d0
C The spinor-rotation matrix is directly calculated from the Euler angles a,b and c.
            IF (rotloc(iatom)%timeinv) THEN
             factor=(rotloc(iatom)%a+rotloc(iatom)%g)/2.d0
             spinrot(2,1)=EXP(CMPLX(0.d0,factor))*
     &         DCOS(rotloc(iatom)%b/2.d0)
             spinrot(1,2)=-CONJG(spinrot(2,1))
C Up/dn and Dn/up terms
             factor=-(rotloc(iatom)%a-rotloc(iatom)%g)/2.d0
             spinrot(2,2)=-EXP(CMPLX(0.d0,factor))*
     &         DSIN(rotloc(iatom)%b/2.d0)
             spinrot(1,1)=CONJG(spinrot(2,2))
            ELSE
             factor=(rotloc(iatom)%a+rotloc(iatom)%g)/2.d0
             spinrot(1,1)=EXP(CMPLX(0.d0,factor))*
     &         DCOS(rotloc(iatom)%b/2.d0)
             spinrot(2,2)=CONJG(spinrot(1,1))
C Up/dn and Dn/up terms
             factor=-(rotloc(iatom)%a-rotloc(iatom)%g)/2.d0
             spinrot(1,2)=EXP(CMPLX(0.d0,factor))*
     &         DSIN(rotloc(iatom)%b/2.d0)
             spinrot(2,1)=-CONJG(spinrot(1,2))
            ENDIF
C Printing the transformation informations 
            DO m=1,2
             WRITE(ouctqmc,*) REAL(spinrot(m,1:2))
            ENDDO
            DO m=1,2
              WRITE(ouctqmc,*) AIMAG(spinrot(m,1:2))
            ENDDO
            DEALLOCATE(spinrot)
C Without SO, only the rotation matrix in orbital space is necessary.
           ELSE
C In this case, the Rloc matrix is merely identity.
            WRITE(ouctqmc,*) 1.d0
            WRITE(ouctqmc,*) 0.d0
           ENDIF
C Printing the time inversion flag if the calculation is spin-polarized (ifSP=1)
           IF (ifSP) THEN
            timeinvflag=0
            IF (rotloc(iatom)%timeinv) timeinvflag=1
            WRITE(ouctqmc,'(i6)') timeinvflag
           ENDIF
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the calculation is necessary spin-polarized with SO, spinor rotation matrices are used.
           IF(crorb(icrorb)%ifsplit) THEN
C If only an irep is correlated
            ind1=1
            DO irep=1,reptrans(l,isrt)%nreps
              IF(crorb(icrorb)%correp(irep)) THEN
               ind2=ind1+reptrans(l,isrt)%dreps(irep)-1
               DO m=ind1,ind2
                 WRITE(ouctqmc,*) REAL(rotloc(iatom)%
     &             rotrep(l)%mat(m,ind1:ind2))
               ENDDO
               DO m=ind1,ind2
                 WRITE(ouctqmc,*)AIMAG(rotloc(iatom)%
     &             rotrep(l)%mat(m,ind1:ind2))
               ENDDO
              ENDIF
              ind1=ind1+reptrans(l,isrt)%dreps(irep)
            ENDDO
           ELSE
C If no particular irep is correlated
            DO m=1,2*(2*l+1)
              WRITE(ouctqmc,*) REAL(rotloc(iatom)%
     &          rotrep(l)%mat(m,1:2*(2*l+1)))
            ENDDO
            DO m=1,2*(2*l+1)
              WRITE(ouctqmc,*) AIMAG(rotloc(iatom)%
     &          rotrep(l)%mat(m,1:2*(2*l+1)))
            ENDDO
           ENDIF
C Printing if the transformation included a time-reversal operation. 
           timeinvflag =0
           IF (rotloc(iatom)%timeinv) timeinvflag=1
           WRITE(ouctqmc,'(i6)') timeinvflag
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
          ELSE
           IF (ifSP.AND.ifSO) THEN 
C If SO is taken into account, spinor rotation matrices are considered.
C The spin is not a good quantum number, so the whole representation is used. 
            DO m=1,2*(2*l+1)
              WRITE(ouctqmc,*) REAL(rotloc(iatom)%
     &          rotrep(l)%mat(m,1:2*(2*l+1)))
            ENDDO
            DO m=1,2*(2*l+1)
              WRITE(ouctqmc,*) AIMAG(rotloc(iatom)%
     &          rotrep(l)%mat(m,1:2*(2*l+1)))
            ENDDO
           ELSE
C The calculation is either spin-polarized without SO or paramagnetic.
C The spin is a good quantum number and irep are possible.
            IF(crorb(icrorb)%ifsplit) THEN
C If only an irep is correlated
             ind1=-l
             DO irep=1,reptrans(l,isrt)%nreps
               IF(crorb(icrorb)%correp(irep)) THEN
                ind2=ind1+reptrans(l,isrt)%dreps(irep)-1
                DO m=ind1,ind2
                  WRITE(ouctqmc,*) REAL(rotloc(iatom)%
     &              rotrep(l)%mat(m,ind1:ind2))
                ENDDO
                DO m=ind1,ind2
                  WRITE(ouctqmc,*) AIMAG(rotloc(iatom)%
     &              rotrep(l)%mat(m,ind1:ind2))
                ENDDO
               ENDIF
               ind1=ind1+reptrans(l,isrt)%dreps(irep)
             ENDDO
            ELSE
C If no particular irep is correlated
             DO m=-l,l
               WRITE(ouctqmc,*) REAL(rotloc(iatom)%
     &           rotrep(l)%mat(m,-l:l))
             ENDDO
             DO m=-l,l
               WRITE(ouctqmc,*) AIMAG(rotloc(iatom)%
     &           rotrep(l)%mat(m,-l:l))
             ENDDO
            END IF   ! End of the ifsplit if-then-else
           END IF    ! End of the ifSO if-then-else
C Printing if the transformation included a time-reversal operation. 
           IF (ifSP) THEN
            timeinvflag =0
            IF (rotloc(iatom)%timeinv) timeinvflag=1
            WRITE(ouctqmc,'(i6)') timeinvflag
           END IF
          END IF     ! End of the ifmixing if-then-else
        END DO       ! End of the icrorb loop
C for each correlated orbital icrorb, the transformation Rloc is described by :
C  - the real part of the submatrix associated to "icrorb"
C  - the imaginary part of the submatrix associated to "icrorb"
C  - a flag which states if a time reversal operation is included in the transformation (if SP only ) 
C
C -------------------------------------------------------------------------------------------------------------
C Description of the transformation from complex harmonics to the basis associated to each correlated orbital :
C -------------------------------------------------------------------------------------------------------------
C
        DO icrorb=1,ncrorb
          l=crorb(icrorb)%l
          isrt=crorb(icrorb)%sort
C The transformation is printed only for the first (representative) atom of each sort.
          IF (crorb(icrorb)%first) THEN
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
           IF (l==0) THEN
C For the s-orbitals, the only basis considered is the complex basis.
             IF (ifSP.AND.ifSO) THEN
C If SO is taken into account, orbital+spin space is considered. 
              ALLOCATE(spinrot(1:2,1:2))
              spinrot(:,:)=0.d0
              spinrot(1,1)=1.d0
              spinrot(2,2)=1.d0
C Printing the number of irep in the considered basis and the size of each irep 
C in the considered basis (in this case, it's 1 and 2.)
              WRITE(ouctqmc,*) 1, 2
C Printing the transformation matrix
              DO m=1,2
                WRITE(ouctqmc,*) REAL(spinrot(m,1:2)) 
              ENDDO
              DO m=1,2
                WRITE(ouctqmc,*) AIMAG(spinrot(m,1:2)) 
              ENDDO
C 
              DEALLOCATE(spinrot)
             ELSE
C Without SO, only the matrix in orbital space is necessary.
C Printing the number of irep in the considered basis and the size of each irep 
C in the considered basis (in this case, it's 1 and 1.)
              WRITE(ouctqmc,'(2(i6,x))') 1, 1
C Printing the transformation matrix
              WRITE(ouctqmc,*) 1.d0 
              WRITE(ouctqmc,*) 0.d0 
             ENDIF
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
           ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the SO is necessary considered, spinor rotation matrices are used.
            IF (crorb(icrorb)%ifsplit) THEN
             WRITE(ouctqmc,'(15(i6,x))') reptrans(l,isrt)%nreps,
     &         reptrans(l,isrt)%dreps(1:reptrans(l,isrt)%nreps)
            ELSE
             WRITE(ouctqmc,'(2(i6,x))') 1,  2*(2*l+1)
            ENDIF
            DO m=1,2*(2*l+1)
              WRITE(ouctqmc,*) REAL(reptrans(l,isrt)%
     &          transmat(m,1:2*(2*l+1)) )
            ENDDO
            DO m=1,2*(2*l+1)
              WRITE(ouctqmc,*) AIMAG(reptrans(l,isrt)%
     &          transmat(m,1:2*(2*l+1)) )
            ENDDO
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
           ELSE
            IF (ifSP.AND.ifSO) THEN
C If SO is taken into account, orbital+spin space is considered.
            ALLOCATE(spinrot(1:2*(2*l+1),1:2*(2*l+1)))
            spinrot(:,:)=0.d0
            spinrot(1:2*l+1,1:2*l+1)= 
     &        reptrans(l,isrt)%transmat(-l:l,-l:l)
            spinrot(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))= 
     &        reptrans(l,isrt)%transmat(-l:l,-l:l)
C Printing the number of irep in the considered basis and the size of each irep 
C in the considered basis (in this case, it's 1 and 2*(2*l+1).)
              WRITE(ouctqmc,*) 1, 2*(2*l+1)
C Printing the transformation matrix
              DO m=1,2*(2*l+1)
                WRITE(ouctqmc,*) REAL(spinrot(m,1:2*(2*l+1)) ) 
              ENDDO
              DO m=1,2*(2*l+1)
                WRITE(ouctqmc,*) AIMAG(spinrot(m,1:2*(2*l+1)) ) 
              ENDDO
C 
              DEALLOCATE(spinrot)
C Without SO, only the rotation matrix in orbital space is necessary.
            ELSE
             IF (crorb(icrorb)%ifsplit) THEN
              WRITE(ouctqmc,'(8(i6,x))') reptrans(l,isrt)%nreps,
     &          reptrans(l,isrt)%dreps(1:reptrans(l,isrt)%nreps)
             ELSE
              WRITE(ouctqmc,'(2(i6,x))') 1, 2*l+1
             ENDIF
             DO m=-l,l
               WRITE(ouctqmc,*) REAL(reptrans(l,isrt)%
     &           transmat(m,-l:l) )
             ENDDO
             DO m=-l,l
               WRITE(ouctqmc,*) AIMAG(reptrans(l,isrt)%
     &           transmat(m,-l:l) )
             END DO
            END IF    ! End of the ifSO if-then-else
           END IF     ! End of the ifmixing if-then-else
          END IF      ! End of the iffirst if-then-else
        END DO        ! End of the icrorb loop
C for each correlated orbital icrorb, the basis transformation is described by :
C  - the number of irep in the new basis
C  - the dimension of each irep in the new basis
C  - the real part of the basis transformation
C  - the imaginary part of the basis transformation
C
C -------------------------------------------------------------------------
C Description of the number of bands in the energy window at each k_point :
C -------------------------------------------------------------------------
C
         DO is=1,ns
C If SO is considered, the number of up and dn bands are the same. 
           IF ((ifSP.AND.ifSO).and.(is.eq.2)) cycle
C Printing the number of included bands in the window for each k_point
           DO ik=1,nk 
             WRITE(ouctqmc,'(i6)')
     &         ABS(kp(ik,is)%nb_top-kp(ik,is)%nb_bot+1)
           ENDDO
         ENDDO
C
C -----------------------------------------------------------
C Description of the projectors for the correlated orbitals :
C -----------------------------------------------------------
        DO ik=1,nk
          DO icrorb=1,ncrorb
            l=crorb(icrorb)%l
            isrt=crorb(icrorb)%sort
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
            IF (l==0) THEN
C For the s-orbitals, the only irep possible is the matrix itself.
C if the calculation is spin-polarized, up and dn projectors are written the one above the other (orbital+spin-space of size 2)
             DO is=1,ns
               WRITE(ouctqmc,*) 
     &           REAL(pr_crorb(icrorb,ik,is)%mat_rep(1,
     &           kp(ik,is)%nb_bot:kp(ik,is)%nb_top))
             ENDDO
             DO is=1,ns 
               WRITE(ouctqmc,*) 
     &           AIMAG(pr_crorb(icrorb,ik,is)%mat_rep(1,
     &           kp(ik,is)%nb_bot:kp(ik,is)%nb_top))
             ENDDO
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
            ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the calculation is necessary spin-polarized with SO, spinor rotation matrices are used.
             IF(crorb(icrorb)%ifsplit) THEN
C If only 1 irep is correlated
              ind1=1
              DO irep=1,reptrans(l,isrt)%nreps
                IF(crorb(icrorb)%correp(irep)) THEN
                 ind2=ind1+reptrans(l,isrt)%dreps(irep)-1
                 DO m=ind1,ind2
                   WRITE(ouctqmc,*) 
     &               REAL(pr_crorb(icrorb,ik,1)%mat_rep(m,
     &               kp(ik,1)%nb_bot:kp(ik,1)%nb_top)) 
                 ENDDO
                 DO m=ind1,ind2
                   WRITE(ouctqmc,*) 
     &               AIMAG(pr_crorb(icrorb,ik,1)%mat_rep(m,
     &               kp(ik,1)%nb_bot:kp(ik,1)%nb_top))
                 ENDDO
                ENDIF
                ind1=ind1+reptrans(l,isrt)%dreps(irep)
              ENDDO
             ELSE
C If no particular irep is correlated
              DO m=1,2*(2*l+1)
                WRITE(ouctqmc,*) 
     &            REAL(pr_crorb(icrorb,ik,1)%mat_rep(m,
     &            kp(ik,1)%nb_bot:kp(ik,1)%nb_top)) 
              ENDDO
              DO m=1,2*(2*l+1)
                WRITE(ouctqmc,*) 
     &            AIMAG(pr_crorb(icrorb,ik,1)%mat_rep(m,
     &            kp(ik,1)%nb_bot:kp(ik,1)%nb_top))
              ENDDO
             ENDIF
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
            ELSE
             IF ((.not.(ifSP.AND.ifSO)).AND.crorb(icrorb)%ifsplit) THEN
C If only 1 irep is correlated (case without SO)
              ind1=-l
              DO irep=1,reptrans(l,isrt)%nreps
                IF(crorb(icrorb)%correp(irep)) THEN
                 ind2=ind1+reptrans(l,isrt)%dreps(irep)-1
                 DO is=1,ns
                   DO m=ind1,ind2
                     WRITE(ouctqmc,*) 
     &                 REAL(pr_crorb(icrorb,ik,is)%mat_rep(m,
     &                 kp(ik,is)%nb_bot:kp(ik,is)%nb_top)) 
                   ENDDO
                 ENDDO
                 DO is=1,ns
                   DO m=ind1,ind2
                     WRITE(ouctqmc,*) 
     &                 AIMAG(pr_crorb(icrorb,ik,is)%mat_rep(m,
     &                 kp(ik,is)%nb_bot:kp(ik,is)%nb_top))
                   ENDDO
                 ENDDO
                ENDIF
                ind1=ind1+reptrans(l,isrt)%dreps(irep)
              ENDDO
             ELSE
C If no particular irep is correlated (case with and without SO)
              DO is=1,ns
                DO m=-l,l
                  WRITE(ouctqmc,*) 
     &              REAL(pr_crorb(icrorb,ik,is)%mat_rep(m,
     &              kp(ik,is)%nb_bot:kp(ik,is)%nb_top)) 
                ENDDO
              ENDDO
              DO is=1,ns
                DO m=-l,l
                  WRITE(ouctqmc,*) 
     &              AIMAG(pr_crorb(icrorb,ik,is)%mat_rep(m,
     &              kp(ik,is)%nb_bot:kp(ik,is)%nb_top))
                ENDDO
              ENDDO 
             END IF   ! End of the ifsplit if-then-else
            END IF    ! End of the ifmixing if-then-else
          END DO      ! End of the icrorb loop
        END DO        ! End of the ik loop
C for each k-point and each correlated orbital, the corresponding projector is described by :
C  - the real part of the "correlated" submatrix
C  - the imaginary part of the "correlated" submatrix 
C
C ----------------------------------------------------------------------------
C Description of the weight of each k_point for the simple point integration :
C ----------------------------------------------------------------------------
        DO ik=1,nk
          WRITE(ouctqmc,*) kp(ik,1)%weight
C This is a geometrical factor independent of the spin value.
        ENDDO
C
C -----------------------------------------------------
C Description of the Hamiltonian H(k) at each k_point :
C -----------------------------------------------------
        DO is=1,ns
          DO ik=1,nk
C If the calculation is spin-polarized with SO, the numbers for up and dn bands are the same. 
            IF ((ifSP.AND.ifSO).AND.(is.eq.2)) cycle
            DO ib=kp(ik,is)%nb_bot,kp(ik,is)%nb_top
              WRITE(ouctqmc,*) kp(ik,is)%eband(ib)
            ENDDO
          ENDDO
        ENDDO
C for each spin value is and each k-point, 
C  - the energies of the band with spin is at point k
C
C
C ======================================
C     Writing the file case.symqmc     :
C ======================================
C
        WRITE(buf,'(a)')'Writing the file case.symqmc...'
        CALL printout(0)
C
C ----------------------------------------------------------
C Description of the general informations about the system :
C ----------------------------------------------------------
        WRITE(ousymqmc,'(2(i6,x))') nsym, natom
C nysm is the total number of symmetries in the system 
C natom is the total number of atom in the unit cell
C
C -------------------------------------------------------------------
C Description of the permutation matrix associated to each symmetry :
C -------------------------------------------------------------------
        DO isym=1,nsym
          WRITE(ousymqmc,'(100(i6,x))') srot(isym)%perm(1:natom)
        ENDDO
C
C ------------------------------------------------------------
C Description of the time-reversal property for each symmetry :
C ------------------------------------------------------------
        IF (ifSP) THEN
         ALLOCATE(timeflag(nsym))
         timeflag=0 
         DO isym=1,nsym
           IF (srot(isym)%timeinv) timeflag(isym)= 1
         ENDDO
         WRITE(ousymqmc,'(100(i6,x))') timeflag(1:nsym)
         DEALLOCATE(timeflag)
C When the calculation is spin-polarized (with SO), a flag which states 
C if a time reversal operation is included in the transformation isym is written.
        ENDIF
C
C -----------------------------------------------------------------------------------------
C Description of the representation matrices of each symmetry for each correlated orbital :
C -----------------------------------------------------------------------------------------
        DO isym=1,nsym
          DO icrorb=1,ncrorb
            isrt=crorb(icrorb)%sort
            l=crorb(icrorb)%l
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
            IF(l==0) THEN
C For the s-orbitals, the only irep possible is the matrix itself.
             IF (ifSP.AND.ifSO) THEN 
C If SO is taken into account, spinor rotation matrix is considered. 
              ALLOCATE(spinrot(1:2,1:2))
              spinrot(:,:)=0.d0
              IF (srot(isym)%timeinv) THEN
C In this case, the Euler angle Beta is Pi. The spinor rotation matrix is block-diagonal, 
C since the time-reversal operator is included in the definition of the transformation.
C
               factor=srot(isym)%phase/2.d0
C We remind that the field phase is (g-a) if beta=Pi. 
C Up/up and Dn/dn terms
               spinrot(1,1)=EXP(CMPLX(0.d0,factor))
               spinrot(2,2)=CONJG(spinrot(1,1))
C spinrot(1,1) = -exp(+i(alpha-gamma)/2) ; spinrot(2,2) = -exp(-i(alpha-gamma)/2)
C in good agreement with Wien conventions for the definition of this phase factor.
              ELSE
C In this case, the Euler angle Beta is 0. The spinor rotation matrix is block-diagonal, 
C No time-reversal treatment was applied to the transformation.
C
               factor=srot(isym)%phase/2.d0
C We remind that the field phase is (a+g) if beta=0. 
C Up/up and Dn/dn terms
               spinrot(1,1)=EXP(CMPLX(0.d0,factor))
               spinrot(2,2)=CONJG(spinrot(1,1))
C the field phase is 2pi-(alpha+gamma) in this case.
C spinrot(1,1) = -exp(-i(alpha+gamma)/2) ; spinrot(2,2) = -exp(i(alpha-gamma)/2)
C in good agreement with Wien conventions for the definition of this phase factor.
              END IF
C Printing the transformation informations 
              DO m=1,2
               WRITE(ousymqmc,*) REAL(spinrot(m,1:2))
              ENDDO
              DO m=1,2
                WRITE(ousymqmc,*) AIMAG(spinrot(m,1:2))
              ENDDO
C
              DEALLOCATE(spinrot)
C Without SO, only the rotation matrix in orbital space is necessary.
             ELSE
C In this case, the Rloc matrix is merely identity.
              WRITE(ousymqmc,*) 1.d0
              WRITE(ousymqmc,*) 0.d0
             ENDIF
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
            ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the SO is necessary considered, spinor rotation matrices are used.
             IF(crorb(icrorb)%ifsplit) THEN
C If only an irep is correlated
              ind1=1
              DO irep=1,reptrans(l,isrt)%nreps
                IF(crorb(icrorb)%correp(irep)) THEN
                 ind2=ind1+reptrans(l,isrt)%dreps(irep)-1
                 DO m=ind1,ind2
                   WRITE(ousymqmc,*) REAL(srot(isym)%
     &               rotrep(l,isrt)%mat(m,ind1:ind2))
                 ENDDO
                 DO m=ind1,ind2
                   WRITE(ousymqmc,*)AIMAG(srot(isym)%
     &               rotrep(l,isrt)%mat(m,ind1:ind2))
                 ENDDO
                ENDIF
                ind1=ind1+reptrans(l,isrt)%dreps(irep)
              ENDDO
             ELSE
C If no particular irep is correlated
              DO m=1,2*(2*l+1)
                WRITE(ousymqmc,*) REAL(srot(isym)%
     &            rotrep(l,isrt)%mat(m,1:2*(2*l+1)))
              ENDDO
              DO m=1,2*(2*l+1)
                WRITE(ousymqmc,*) AIMAG(srot(isym)%
     &            rotrep(l,isrt)%mat(m,1:2*(2*l+1)))
              ENDDO
             ENDIF
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
            ELSE
             IF (ifSP.AND.ifSO) THEN 
C If the calculation is spin-polarized with SO, spinor rotation matrices are considered.
C The spin is not a good quantum number, so the whole representation is used. 
              ALLOCATE(spinrot(1:2*(2*l+1),1:2*(2*l+1)))
              spinrot(:,:)=0.d0
              IF (srot(isym)%timeinv) THEN
C In this case, the Euler angle Beta is Pi. The spinor rotation matrix is block-diagonal, 
C since the time-reversal operator is included in the definition of the transformation.
C
               factor=srot(isym)%phase/2.d0
C We remind that the field phase is (g-a) in this case.
C Up/up block :
               ephase=EXP(CMPLX(0.d0,factor))
C As a result, ephase = -exp(i(alpha-gamma)/2)
               spinrot(1:2*l+1,1:2*l+1)=
     =           ephase*srot(isym)%rotrep(l,isrt)%mat(-l:l,-l:l)
C Dn/dn block :
               ephase=CONJG(ephase)
C Now, ephase = -exp(-i(alpha-gamma)/2)
               spinrot(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))=
     =           ephase*srot(isym)%rotrep(l,isrt)%mat(-l:l,-l:l)
              ELSE
C In this case, the Euler angle Beta is 0. The spinor rotation matrix is block-diagonal, 
C No time-reversal treatment was applied to the transformation.
C
               factor=srot(isym)%phase/2.d0
C We remind that the field phase is (a+g) in this case.
C Up/up block :
               ephase=EXP(CMPLX(0.d0,factor))
C As a result, ephase = -exp(-i(alpha+gamma)/2)
               spinrot(1:2*l+1,1:2*l+1)=
     =           ephase*srot(isym)%rotrep(l,isrt)%mat(-l:l,-l:l)
C Dn/dn block :
               ephase=CONJG(ephase)
C Now, ephase = -exp(i(alpha+gamma)/2)
               spinrot(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))=
     =           ephase*srot(isym)%rotrep(l,isrt)%mat(-l:l,-l:l)
              END IF
C Printing the transformation informations 
              DO m=1,2*(2*l+1)
                WRITE(ousymqmc,*) REAL(spinrot(m,1:2*(2*l+1)) )
              ENDDO
              DO m=1,2*(2*l+1)
                WRITE(ousymqmc,*) AIMAG(spinrot(m,1:2*(2*l+1)) )
              ENDDO
C
              DEALLOCATE(spinrot)
C In the other cases (spin-polarized without SO or paramagnetic), only the rotation matrix in orbital space is necessary.
             ELSE
              IF(crorb(icrorb)%ifsplit) THEN
C If only an irep is correlated
               ind1=-l
               DO irep=1,reptrans(l,isrt)%nreps
                 IF(crorb(icrorb)%correp(irep)) THEN
                  ind2=ind1+reptrans(l,isrt)%dreps(irep)-1
                  DO m=ind1,ind2
                    WRITE(ousymqmc,*) REAL(srot(isym)%
     &                rotrep(l,isrt)%mat(m,ind1:ind2))
                  ENDDO
                  DO m=ind1,ind2
                    WRITE(ousymqmc,*) AIMAG(srot(isym)%
     &                rotrep(l,isrt)%mat(m,ind1:ind2))
                  ENDDO
                 ENDIF
                 ind1=ind1+reptrans(l,isrt)%dreps(irep)
               ENDDO
              ELSE
C If no particular irep is correlated
               DO m=-l,l
                 WRITE(ousymqmc,*) REAL(srot(isym)%
     &             rotrep(l,isrt)%mat(m,-l:l))
               ENDDO
               DO m=-l,l
                 WRITE(ousymqmc,*) AIMAG(srot(isym)%
     &             rotrep(l,isrt)%mat(m,-l:l))
               ENDDO
              END IF   ! End of the ifsplit if-then-else
             END IF    ! End of the ifSO if-then-else
            END IF     ! End of the ifmixing if-then-else
          END DO       ! End of the icrorb loop
        END DO         ! End of the isym loop
C for each symmetry and each correlated orbital icrorb, is described :
C  - the real part of the submatrix associated to "isym" for the "icrorb" 
C  - the imaginary part of the submatrix associated to "isym" for the "icrorb"
C
C -----------------------------------------------------------------------
C Description of the time reversal operator for each correlated orbital :
C -----------------------------------------------------------------------
C This description occurs only if the computation is paramagnetic. (ifSP=0)
        IF (.not.ifSP) THEN
         DO icrorb=1,ncrorb
           l=crorb(icrorb)%l
           isrt=crorb(icrorb)%sort
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
           IF (l==0) THEN
C For the s-orbitals, the only basis considered is the complex basis.
            WRITE(ousymqmc,*) 1.d0 
            WRITE(ousymqmc,*) 0.d0 
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
           ELSE
C Calculation of the time-reversal operator
            ALLOCATE(time_op(-l:l,-l:l))
            time_op(:,:)=0.d0
            DO m=-l,l
              time_op(m,m)=1.d0
            ENDDO
C time_op is Identity.
            CALL timeinv_op(time_op,(2*l+1),l,isrt)
C time_op is now the time-reversal operator in the desired basis ({new_i})
C
            IF(crorb(icrorb)%ifsplit) THEN
C If only an irep is correlated
             ind1=-l
             DO irep=1,reptrans(l,isrt)%nreps
               IF(crorb(icrorb)%correp(irep)) THEN
                ind2=ind1+reptrans(l,isrt)%dreps(irep)-1
                DO m=ind1,ind2
                  WRITE(ousymqmc,*) REAL(time_op(m,ind1:ind2))
                ENDDO
                DO m=ind1,ind2
                  WRITE(ousymqmc,*) AIMAG(time_op(m,ind1:ind2))
                ENDDO
               ENDIF
               ind1=ind1+reptrans(l,isrt)%dreps(irep)
             ENDDO
            ELSE
C If no particular irep is correlated
             DO m=-l,l
               WRITE(ousymqmc,*) REAL(time_op(m,-l:l))
             ENDDO
             DO m=-l,l
               WRITE(ousymqmc,*) AIMAG(time_op(m,-l:l))
             ENDDO
            END IF   ! End of the ifsplit if-then-else
            DEALLOCATE(time_op)
           END IF    ! End of the l if-then-else
         END DO       ! End of the icrorb loop
        END IF        ! End of the ifsp if-then-else
C for each correlated orbital icrorb, the time-reversal operator is described by :
C  - its real part 
C  - its imaginary part 
C
C
C ===============================
C Writing the file case.parproj :
C ===============================
C
        WRITE(buf,'(a)')'Writing the file case.parproj...'
        CALL printout(0)
C
C ----------------------------------------------------------------
C Description of the size of the basis for each included orbital :
C ----------------------------------------------------------------
        DO iorb=1,norb
          WRITE(oupartial,'(3(i6))') norm_radf(iorb)%n
        ENDDO
C There is not more than 1 LO for each orbital (hence n < 4 )
C 
C The following descriptions are made for each included orbital.
        DO iorb=1,norb
          l=orb(iorb)%l
          isrt=orb(iorb)%sort
C ------------------------------------
C Description of the Theta projector :
C ------------------------------------
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
          IF (l==0) THEN
           DO ik=1,nk
             DO ir=1,norm_radf(iorb)%n
               DO is=1,ns
                 WRITE(oupartial,*) 
     &             REAL(pr_orb(iorb,ik,is)%matn_rep(1,
     &             kp(ik,is)%nb_bot:kp(ik,is)%nb_top,ir))
               ENDDO
               DO is=1,ns 
                 WRITE(oupartial,*) 
     &             AIMAG(pr_orb(iorb,ik,is)%matn_rep(1,
     &             kp(ik,is)%nb_bot:kp(ik,is)%nb_top,ir))
               ENDDO
             ENDDO    ! End of the ir loop
           ENDDO      ! End of the ik loop
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the calculation is necessary spin-polarized with SO, spinor rotation matrices are used.
           DO ik=1,nk
             DO ir=1,norm_radf(iorb)%n
               DO m=1,2*(2*l+1)
                 WRITE(oupartial,*) 
     &             REAL(pr_orb(iorb,ik,1)%matn_rep(m,
     &             kp(ik,1)%nb_bot:kp(ik,1)%nb_top,ir)) 
               ENDDO
               DO m=1,2*(2*l+1)
                 WRITE(oupartial,*) 
     &             AIMAG(pr_orb(iorb,ik,1)%matn_rep(m,
     &             kp(ik,1)%nb_bot:kp(ik,1)%nb_top,ir))
               ENDDO
             ENDDO    ! End of the ir loop
           ENDDO      ! End of the ik loop
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
          ELSE
           DO ik=1,nk
             DO ir=1,norm_radf(iorb)%n
               DO is=1,ns
                 DO m=-l,l
                   WRITE(oupartial,*) 
     &               REAL(pr_orb(iorb,ik,is)%matn_rep(m,
     &               kp(ik,is)%nb_bot:kp(ik,is)%nb_top,ir)) 
                 ENDDO
               ENDDO  ! End of the is loop
               DO is=1,ns
                 DO m=-l,l
                   WRITE(oupartial,*) 
     &               AIMAG(pr_orb(iorb,ik,is)%matn_rep(m,
     &               kp(ik,is)%nb_bot:kp(ik,is)%nb_top,ir))
                 ENDDO
               ENDDO  ! End of the is loop 
             ENDDO    ! End of the ir loop
           ENDDO      ! End of the ik loop
          ENDIF       ! End of the ifmixing if-then-else
C for each included orbital, for each k-point and each |phi_j> elmt, 
C the corresponding Thetaprojector is described by :
C  - the real part of the matrix
C  - the imaginary part of the matrix 
C
C -------------------------------------------------------------------------------
C Description of the density matrices below the lower limit e_bot of the window :
C -------------------------------------------------------------------------------
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
          IF (l==0) THEN
C With SO, the density matrix is printed for the complete spin+orbital subspace. 
C (in this case, this means a matrix of size 2)
           IF (ifSP.AND.ifSO) THEN
            ALLOCATE(densprint(2,2))
            densprint(1,1)=densmat(1,iorb)%mat(1,1)
            densprint(2,2)=densmat(2,iorb)%mat(1,1)
            densprint(1,2)=densmat(3,iorb)%mat(1,1)
            densprint(2,1)=densmat(4,iorb)%mat(1,1)
            DO m=1,2
              WRITE(oupartial,*) REAL(densprint(m,1:2))
            ENDDO
            DO m=1,2
              WRITE(oupartial,*) AIMAG(densprint(m,1:2))
            ENDDO
            DEALLOCATE(densprint)
C In the other cases the density matrix is printed for each spin subspace independently.
C (in this case, this means two matrices of size 1)
           ELSE
            DO is=1,ns
              WRITE(oupartial,*) REAL(densmat(is,iorb)%mat(1,1))
            ENDDO
            DO is=1,ns
              WRITE(oupartial,*) AIMAG(densmat(is,iorb)%mat(1,1))
            ENDDO
           ENDIF  ! End of the ifSO if-then-else
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the calculation is necessary spin-polarized with SO, spinor rotation matrices are used.
C As a result, the density matrix is printed for the complete spin+orbital subspace. 
C (this means a matrix of size 2*(2*l+1))
           DO m=1,2*(2*l+1)
             WRITE(oupartial,*) 
     &         REAL(densmat(1,iorb)%mat(m,1:2*(2*l+1)))
           ENDDO
           DO m=1,2*(2*l+1)
             WRITE(oupartial,*) 
     &         AIMAG(densmat(1,iorb)%mat(m,1:2*(2*l+1)))
           ENDDO
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
          ELSE
C With SO, the density matrix is printed for the complete spin+orbital subspace. 
C (this means a matrix of size 2*(2*l+1))
           IF (ifSP.AND.ifSO) THEN
            ALLOCATE(densprint(1:2*(2*l+1),1:2*(2*l+1)))
            densprint(1:2*l+1,1:2*l+1)=
     &        densmat(1,iorb)%mat(-l:l,-l:l)
            densprint(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))=
     &        densmat(2,iorb)%mat(-l:l,-l:l)
            densprint(1:2*l+1,2*l+2:2*(2*l+1))=
     &        densmat(3,iorb)%mat(-l:l,-l:l)
            densprint(2*l+2:2*(2*l+1),1:2*l+1)=
     &        densmat(4,iorb)%mat(-l:l,-l:l)
            DO m=1,2*(2*l+1)
              WRITE(oupartial,*) REAL(densprint(m,1:2*(2*l+1)))
            ENDDO
            DO m=1,2*(2*l+1)
              WRITE(oupartial,*) AIMAG(densprint(m,1:2*(2*l+1)))
            ENDDO
            DEALLOCATE(densprint)
C In the other cases, the density matrix is printed for each spin subspace independently.
C (this means two matrices of size 2*l+1)
           ELSE
            DO is=1,ns
              DO m=-l,l
                WRITE(oupartial,*) 
     &            REAL(densmat(is,iorb)%mat(m,-l:l))
              ENDDO
            ENDDO
            DO is=1,ns
              DO m=-l,l
                WRITE(oupartial,*) 
     &            AIMAG(densmat(is,iorb)%mat(m,-l:l))
             ENDDO
            ENDDO
           ENDIF  ! End of the ifSO if-then-else
          ENDIF       ! End of the ifmixing if-then-else
C for each included orbital, the corresponding density matrix is described by :
C  - the real part of the matrix
C  - the imaginary part of the matrix 
C
C --------------------------------------------------------------------
C Description of the global to local coordinates transformation Rloc :
C --------------------------------------------------------------------
C Description of each transformation Rloc
          iatom=orb(iorb)%atom
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
          IF(l==0) THEN
C For the s-orbitals, the only irep possible is the matrix itself.
           IF (ifSP.AND.ifSO) THEN 
C If SO is taken into account, spinor rotation matrix must be considered.
            ALLOCATE(spinrot(1:2,1:2))
            spinrot(:,:)=0.d0
C The spinor-rotation matrix is directly calculated from the Euler angles a,b and c.
            IF (rotloc(iatom)%timeinv) THEN
             factor=(rotloc(iatom)%a+rotloc(iatom)%g)/2.d0
             spinrot(2,1)=EXP(CMPLX(0.d0,factor))*
     &         DCOS(rotloc(iatom)%b/2.d0)
             spinrot(1,2)=-CONJG(spinrot(2,1))
C Up/dn and Dn/up terms
             factor=-(rotloc(iatom)%a-rotloc(iatom)%g)/2.d0
             spinrot(2,2)=-EXP(CMPLX(0.d0,factor))*
     &         DSIN(rotloc(iatom)%b/2.d0)
             spinrot(1,1)=CONJG(spinrot(2,2))
            ELSE
             factor=(rotloc(iatom)%a+rotloc(iatom)%g)/2.d0
             spinrot(1,1)=EXP(CMPLX(0.d0,factor))*
     &         DCOS(rotloc(iatom)%b/2.d0)
             spinrot(2,2)=CONJG(spinrot(1,1))
C Up/dn and Dn/up terms
             factor=-(rotloc(iatom)%a-rotloc(iatom)%g)/2.d0
             spinrot(1,2)=EXP(CMPLX(0.d0,factor))*
     &         DSIN(rotloc(iatom)%b/2.d0)
             spinrot(2,1)=-CONJG(spinrot(1,2))
            ENDIF
C Printing the transformation informations 
            DO m=1,2
             WRITE(oupartial,*) REAL(spinrot(m,1:2))
            ENDDO
            DO m=1,2
              WRITE(oupartial,*) AIMAG(spinrot(m,1:2))
            ENDDO
            DEALLOCATE(spinrot)
C Without SO, only the rotation matrix in orbital space is necessary.
           ELSE
C In this case, the Rloc matrix is merely identity.
            WRITE(oupartial,*) 1.d0
            WRITE(oupartial,*) 0.d0
           ENDIF
C Printing the time inversion flag if the calculation is spin-polarized (ifSP=1)
           IF (ifSP) THEN
            timeinvflag=0
            IF (rotloc(iatom)%timeinv) timeinvflag=1
            WRITE(oupartial,'(i6)') timeinvflag
           ENDIF
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the calculation is necessary spin-polarized with SO, spinor rotation matrices are used.
           DO m=1,2*(2*l+1)
             WRITE(oupartial,*) REAL(rotloc(iatom)%
     &         rotrep(l)%mat(m,1:2*(2*l+1)))
           ENDDO
           DO m=1,2*(2*l+1)
             WRITE(oupartial,*) AIMAG(rotloc(iatom)%
     &         rotrep(l)%mat(m,1:2*(2*l+1)))
           ENDDO
C Printing if the transforamtion included a time-reversal operation. 
           timeinvflag =0
           IF (rotloc(iatom)%timeinv) timeinvflag=1
           WRITE(oupartial,'(i6)') timeinvflag
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
          ELSE
           IF (ifSP.AND.ifSO) THEN 
C If SO is taken into account, spinor rotation matrices are considered.
C The spin is not a good quantum number, so the whole representation is used. 
            DO m=1,2*(2*l+1)
              WRITE(oupartial,*) REAL(rotloc(iatom)%
     &          rotrep(l)%mat(m,1:2*(2*l+1)))
            ENDDO
            DO m=1,2*(2*l+1)
              WRITE(oupartial,*) AIMAG(rotloc(iatom)%
     &          rotrep(l)%mat(m,1:2*(2*l+1)))
            ENDDO
           ELSE
C The calculation is either spin-polarized without SO or paramagnetic.
C The spin is a good quantum number and irep are possible.
            DO m=-l,l
              WRITE(oupartial,*) REAL(rotloc(iatom)%
     &          rotrep(l)%mat(m,-l:l))
            ENDDO
            DO m=-l,l
              WRITE(oupartial,*) AIMAG(rotloc(iatom)%
     &          rotrep(l)%mat(m,-l:l))
            ENDDO
           END IF    ! End of the ifSO if-then-else
C Printing if the transformation included a time-reversal operation. 
           IF (ifSP) THEN
            timeinvflag =0
            IF (rotloc(iatom)%timeinv) timeinvflag=1
            WRITE(oupartial,'(i6)') timeinvflag
           END IF
          END IF     ! End of the ifmixing if-then-else
C for each included orbital iorb, the transformation Rloc is described by :
C  - the real part of the submatrix associated to "iorb"
C  - the imaginary part of the submatrix associated to "iorb"
C  - a flag which states if a time reversal operation is included in the transformation (if SP only ) 
C
        END DO       ! End of the iorb loop
C 
C
C ==============================
C Writing the file case.sympar :
C ==============================
C
        WRITE(buf,'(a)')'Writing the file case.sympar...'
        CALL printout(0)
C
C ----------------------------------------------------------
C Description of the general informations about the system :
C ----------------------------------------------------------
        WRITE(ousympar,'(2(i6,x))') nsym, natom
C nysm is the total number of symmetries of the system 
C natom is the total number of atom in the unit cell
C
C -------------------------------------------------------------------
C Description of the permutation matrix associated to each symmetry :
C -------------------------------------------------------------------
        DO is=1,nsym
          WRITE(ousympar,'(100(i6,x))') srot(is)%perm
        ENDDO
C
C ------------------------------------------------------------
C Description of the time-reversal property for each symmetry :
C ------------------------------------------------------------
        IF (ifSP) THEN
         ALLOCATE(timeflag(nsym))
         timeflag=0 
         DO isym=1,nsym
           IF (srot(isym)%timeinv) timeflag(isym)= 1
         ENDDO
         WRITE(ousympar,'(100(i6,x))') timeflag(1:nsym)
         DEALLOCATE(timeflag)
C When the calculation is spin-polarized (with SO), a flag which states 
C if a time reversal operation is included in the transformation isym is written.
        ENDIF
C
C -------------------------------------------------------------------------------------------
C Description of the representation matrices of each symmetry for all the included orbitals :
C -------------------------------------------------------------------------------------------
        DO isym=1,nsym
          DO iorb=1,norb
            isrt=orb(iorb)%sort
            l=orb(iorb)%l
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
            IF(l==0) THEN
C For the s-orbitals, the only irep possible is the matrix itself.
             IF (ifSP.AND.ifSO) THEN 
C If SO is taken into account, spinor rotation matrix is considered. 
              ALLOCATE(spinrot(1:2,1:2))
              spinrot(:,:)=0.d0
              IF (srot(isym)%timeinv) THEN
C In this case, the Euler angle Beta is Pi. The spinor rotation matrix is block-diagonal, 
C since the time-reversal operator is included in the definition of the transformation.
C
               factor=srot(isym)%phase/2.d0
C We remind that the field phase is (g-a) if beta=Pi. 
C Up/up and Dn/dn terms
               spinrot(1,1)=EXP(CMPLX(0.d0,factor))
               spinrot(2,2)=CONJG(spinrot(1,1))
C spinrot(1,1) = -exp(i(alpha-gamma)/2) ; spinrot(2,2) = -exp(-i(alpha-gamma)/2)
C in good agreement with Wien conventions for the definition of this phase factor.
              ELSE
C In this case, the Euler angle Beta is 0. The spinor rotation matrix is block-diagonal, 
C No time-reversal treatment was applied to the transformation.
C
               factor=srot(isym)%phase/2.d0
C We remind that the field phase is (a+g) if beta=0. 
C Up/up and Dn/dn terms
               spinrot(1,1)=EXP(CMPLX(0.d0,factor))
               spinrot(2,2)=CONJG(spinrot(1,1))
C spinrot(1,1) = -exp(-i(alpha+gamma)/2) ; spinrot(2,2) = -exp(i(alpha-gamma)/2)
C in good agreement with Wien conventions for the definition of this phase factor.
              END IF
C Printing the transformation informations 
              DO m=1,2
               WRITE(ousympar,*) REAL(spinrot(m,1:2))
              ENDDO
              DO m=1,2
                WRITE(ousympar,*) AIMAG(spinrot(m,1:2))
              ENDDO
C
              DEALLOCATE(spinrot)
C Without SO, only the rotation matrix in orbital space is necessary.
             ELSE
C In this case, the Rloc matrix is merely identity.
              WRITE(ousympar,*) 1.d0
              WRITE(ousympar,*) 0.d0
             ENDIF
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
            ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the SO is necessary considered, spinor rotation matrices are used.
             DO m=1,2*(2*l+1)
               WRITE(ousympar,*) REAL(srot(isym)%
     &           rotrep(l,isrt)%mat(m,1:2*(2*l+1)))
             ENDDO
             DO m=1,2*(2*l+1)
               WRITE(ousympar,*) AIMAG(srot(isym)%
     &           rotrep(l,isrt)%mat(m,1:2*(2*l+1)))
             ENDDO
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
            ELSE
             IF (ifSP.AND.ifSO) THEN 
C If the calculation is spin-polarized with SO, spinor rotation matrices are considered.
C The spin is not a good quantum number, so the whole representation is used. 
              ALLOCATE(spinrot(1:2*(2*l+1),1:2*(2*l+1)))
              spinrot(:,:)=0.d0
              IF (srot(isym)%timeinv) THEN
C In this case, the Euler angle Beta is Pi. The spinor rotation matrix is block-diagonal, 
C since the time-reversal operator is included in the definition of the transformation.
C
               factor=srot(isym)%phase/2.d0
C We remind that the field phase is (g-a) in this case.
C Up/up block :
               ephase=EXP(CMPLX(0.d0,factor))
C AS a result, ephase = -exp(i(alpha-gamma)/2)
               spinrot(1:2*l+1,1:2*l+1)=
     =           ephase*srot(isym)%rotrep(l,isrt)%mat(-l:l,-l:l)
C Dn/dn block :
               ephase=CONJG(ephase)
C Now, ephase = -exp(-i(alpha-gamma)/2)
               spinrot(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))=
     =           ephase*srot(isym)%rotrep(l,isrt)%mat(-l:l,-l:l)
              ELSE
C In this case, the Euler angle Beta is 0. The spinor rotation matrix is block-diagonal, 
C No time-reversal treatment was applied to the transformation.
C
               factor=srot(isym)%phase/2.d0
C We remind that the field phase is 2pi-(alpha+gamma) in this case.
C Up/up block :
               ephase=EXP(CMPLX(0.d0,factor))
C As a result, ephase = -exp(-i(alpha+gamma)/2)
               spinrot(1:2*l+1,1:2*l+1)=
     =           ephase*srot(isym)%rotrep(l,isrt)%mat(-l:l,-l:l)
C Dn/dn block :
               ephase=CONJG(ephase)
C Now, ephase = -exp(i(alpha+gamma)/2)
               spinrot(2*l+2:2*(2*l+1),2*l+2:2*(2*l+1))=
     =           ephase*srot(isym)%rotrep(l,isrt)%mat(-l:l,-l:l)
              END IF
C Printing the transformation informations 
              DO m=1,2*(2*l+1)
                WRITE(ousympar,*) REAL(spinrot(m,1:2*(2*l+1)) )
              ENDDO
              DO m=1,2*(2*l+1)
                WRITE(ousympar,*) AIMAG(spinrot(m,1:2*(2*l+1)) )
              ENDDO
C
              DEALLOCATE(spinrot)
C In the other cases (spin-polarized without SO or paramagnetic), only the rotation matrix in orbital space is necessary.
             ELSE
              DO m=-l,l
                WRITE(ousympar,*) REAL(srot(isym)%
     &            rotrep(l,isrt)%mat(m,-l:l))
              ENDDO
              DO m=-l,l
                WRITE(ousympar,*) AIMAG(srot(isym)%
     &            rotrep(l,isrt)%mat(m,-l:l))
              ENDDO
             END IF    ! End of the ifSO if-then-else
            END IF     ! End of the ifmixing if-then-else
          END DO       ! End of the iorb loop
        END DO         ! End of the isym loop
C for each symmetry and each included orbital iorb, is described :
C  - the real part of the matrix associated to "isym" for the "iorb" 
C  - the imaginary part of the matrix associated to "isym" for the "iorb"
C
C ---------------------------------------------------------------------
C Description of the time reversal operator for each included orbital :
C ---------------------------------------------------------------------
C This description occurs only if the computation is paramagnetic (ifSP=0)
        IF (.not.ifSP) THEN
         DO iorb=1,norb
           l=orb(iorb)%l
           isrt=orb(iorb)%sort
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
           IF (l==0) THEN
C For the s-orbitals, the only basis considered is the complex basis.
            WRITE(ousympar,*) 1.d0 
            WRITE(ousympar,*) 0.d0 
C
C If the basis representation can be reduce to the up/up block (basis without "mixing").
C --------------------------------------------------------------------------------------
           ELSE
C Calculation of the time-reversal operator
            ALLOCATE(time_op(-l:l,-l:l))
            time_op(:,:)=0.d0
            DO m=-l,l
              time_op(m,m)=1.d0
            ENDDO
C time_op is Identity.
            CALL timeinv_op(time_op,(2*l+1),l,isrt)
C time_op is now the time-reversal operator in the desired basis ({new_i})
C
            DO m=-l,l
              WRITE(ousympar,*) REAL(time_op(m,-l:l))
            ENDDO
            DO m=-l,l
              WRITE(ousympar,*) AIMAG(time_op(m,-l:l))
            ENDDO
            DEALLOCATE(time_op)
           END IF    ! End of the l if-then-else
         END DO       ! End of the icrorb loop
        END IF        ! End of the ifsp if-then-else
C for each included orbital iorb, the time-reversal operator is described by :
C  - its real part 
C  - its imaginary part 
C
        RETURN
        END
        
        
           
