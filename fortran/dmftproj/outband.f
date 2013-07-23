
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

        SUBROUTINE outband
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine creates the output file case.outband, with all  %%
C %% the informations necessary for the computation of the spectral  %%
C %% function of the system.                                         %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definition of the variables :
C -----------------------------
        USE almblm_data
        USE bands
        USE common_data
        USE file_names
        USE prnt
        USE projections
        USE reps
        IMPLICIT NONE
C
        INTEGER :: iorb, icrorb, irep, isrt
        INTEGER :: l, m, is, i1, i2, i
        INTEGER :: ik, il, ib, ir, n
        INTEGER :: ind1, ind2, iatom
C
        WRITE(buf,'(a)')'Writing the file case.outband...'
        CALL printout(0)
C
C ======================================
C Informations about the chosen k-path :
C ======================================
C
C Number of k-points along the chosen k-path
        WRITE(ouband,'(i6)') nkband
C Description of the number of bands in the energy window at each k_point
C
        DO is=1,ns
C If SO is considered, the number of up and dn bands are the same. 
          IF ((ifSP.AND.ifSO).and.(is.eq.2)) cycle
          DO ik=1,nk 
            WRITE(ouband,'(i6)')
     &        ABS(kp(ik,is)%nb_top-kp(ik,is)%nb_bot+1)
          ENDDO   ! End of the ik loop
        ENDDO     ! End of the is loop
C for each k-point, the number of band included in the energy window is written.
C ===========================================================
C Description of the projectors for the correlated orbitals :
C ===========================================================
        DO ik=1,nk
          DO icrorb=1,ncrorb
            l=crorb(icrorb)%l
            isrt=crorb(icrorb)%sort
C 
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
            IF (l==0) THEN
C For the s-orbitals, the only irep possible is the matrix itself.
             DO is=1,ns
               WRITE(ouband,*) 
     &           REAL(pr_crorb(icrorb,ik,is)%mat_rep(1,
     &           kp(ik,is)%nb_bot:kp(ik,is)%nb_top))
             ENDDO
             DO is=1,ns 
               WRITE(ouband,*) 
     &           AIMAG(pr_crorb(icrorb,ik,is)%mat_rep(1,
     &           kp(ik,is)%nb_bot:kp(ik,is)%nb_top))
             ENDDO
C
C If the basis representation needs a complete spinor rotation approach (basis with "mixing" ).
C ---------------------------------------------------------------------------------------------
            ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C In this case, the SO is necessary considered, spinor rotation matrices are used.
             IF(crorb(icrorb)%ifsplit) THEN
C If only 1 irep is correlated
              ind1=1
              DO irep=1,reptrans(l,isrt)%nreps
                IF(crorb(icrorb)%correp(irep)) THEN
                 ind2=ind1+reptrans(l,isrt)%dreps(irep)-1
                 DO m=ind1,ind2
                   WRITE(ouband,*) 
     &               REAL(pr_crorb(icrorb,ik,1)%mat_rep(m,
     &               kp(ik,1)%nb_bot:kp(ik,1)%nb_top)) 
                 ENDDO
                 DO m=ind1,ind2
                   WRITE(ouband,*) 
     &               AIMAG(pr_crorb(icrorb,ik,1)%mat_rep(m,
     &               kp(ik,1)%nb_bot:kp(ik,1)%nb_top))
                 ENDDO
                ENDIF
                ind1=ind1+reptrans(l,isrt)%dreps(irep)
              ENDDO
             ELSE
C If no particular irep is correlated
              DO m=1,2*(2*l+1)
                WRITE(ouband,*) 
     &            REAL(pr_crorb(icrorb,ik,1)%mat_rep(m,
     &            kp(ik,1)%nb_bot:kp(ik,1)%nb_top)) 
              ENDDO
              DO m=1,2*(2*l+1)
                WRITE(ouband,*) 
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
                     WRITE(ouband,*) 
     &                 REAL(pr_crorb(icrorb,ik,is)%mat_rep(m,
     &                 kp(ik,is)%nb_bot:kp(ik,is)%nb_top)) 
                   ENDDO
                 ENDDO
                 DO is=1,ns
                   DO m=ind1,ind2
                     WRITE(ouband,*) 
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
                  WRITE(ouband,*) 
     &              REAL(pr_crorb(icrorb,ik,is)%mat_rep(m,
     &              kp(ik,is)%nb_bot:kp(ik,is)%nb_top)) 
                ENDDO
              ENDDO
              DO is=1,ns
                DO m=-l,l
                  WRITE(ouband,*) 
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
C ======================================================
C Description of the Hamiltonian H(k) at each k_point :
C ======================================================
        DO is=1,ns
          DO ik=1,nk
C If SO is considered, the numbers of up and dn bands are the same. 
            IF (ifSO.and.is.eq.2) cycle
            DO ib=kp(ik,is)%nb_bot,kp(ik,is)%nb_top
              WRITE(ouband,*) kp(ik,is)%eband(ib)
            ENDDO
          ENDDO   ! End of the ik loop
        ENDDO     ! End of the is loop
C for each spin value is and each k-point, 
C  - the energies of the band with spin is at point k
C
C ================================================================
C Description of the size of the basis for each included orbital :
C ================================================================
        DO iorb=1,norb
          WRITE(ouband,'(3(i6))') norm_radf(iorb)%n
        ENDDO
C There is not more than 1 LO for each orbital (hence n < 4 )
C 
C ====================================
C Description of the Theta projector :
C ====================================
        DO iorb=1,norb
          l=orb(iorb)%l
          isrt=orb(iorb)%sort
C
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
          IF (l==0) THEN
           DO ik=1,nk
             DO ir=1,norm_radf(iorb)%n
               DO is=1,ns
                 WRITE(ouband,*) 
     &             REAL(pr_orb(iorb,ik,is)%matn_rep(1,
     &             kp(ik,is)%nb_bot:kp(ik,is)%nb_top,ir))
               ENDDO
               DO is=1,ns 
                 WRITE(ouband,*) 
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
                 WRITE(ouband,*) 
     &             REAL(pr_orb(iorb,ik,1)%matn_rep(m,
     &             kp(ik,1)%nb_bot:kp(ik,1)%nb_top,ir)) 
               ENDDO
               DO m=1,2*(2*l+1)
                 WRITE(ouband,*) 
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
                   WRITE(ouband,*) 
     &               REAL(pr_orb(iorb,ik,is)%matn_rep(m,
     &               kp(ik,is)%nb_bot:kp(ik,is)%nb_top,ir)) 
                 ENDDO
               ENDDO  ! End of the is loop
               DO is=1,ns
                 DO m=-l,l
                   WRITE(ouband,*) 
     &               AIMAG(pr_orb(iorb,ik,is)%matn_rep(m,
     &               kp(ik,is)%nb_bot:kp(ik,is)%nb_top,ir))
                 ENDDO
               ENDDO  ! End of the is loop 
             ENDDO    ! End of the ir loop
           ENDDO      ! End of the ik loop
          ENDIF       ! End of the ifmixing if-then-else
        ENDDO         ! End of the iorb loop
C for each included orbital, for each k-point and each |phi_j> elmt, 
C the corresponding Thetaprojector is described by :
C  - the real part of the matrix
C  - the imaginary part of the matrix 
C
C =============================
C Description of the k-labels :
C =============================
        DO i=1,nlab
          WRITE(ouband,'(2i6,a)') i,labels(i)%pos,labels(i)%kname
        ENDDO
C for each label, are written :
C  - the number of the corresponding k-point in the k-path
C  - the name associated to this label
C
       RETURN
        END
        
        
           
