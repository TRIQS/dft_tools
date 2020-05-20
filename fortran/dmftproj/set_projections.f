
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

        SUBROUTINE set_projections(e1,e2)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine sets up the projection matrices in the energy   %%
C %% window [e1,e2]. If proj_mode is 1 or 2 then e1 and e2 are not   %%
C %% energies, but directly used as band indices.                    %%
C %% Two types of projection can be defined :                        %%
C %%   - The projectors <u_orb|ik,ib,is> for the correlated orbital  %%
C %%     only. (orb = iatom,is,m)                                    %%
C %%     (They are stored in the table pr_crorb)                     %%
C %%   - The Theta projectors <theta_orb|k,ib> for all the orbitals  %%
C %%     (They are stored in the table pr_orb)                       %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
C
        USE almblm_data
        USE common_data
        USE prnt
        USE projections
        USE reps
        USE symm       
        IMPLICIT NONE
C
        REAL(KIND=8) :: e1, e2
        INTEGER :: iorb, icrorb, ik, is, ib, m, l, lm, nbbot, nbtop
        INTEGER :: isrt, n, ilo, iatom, i, imu, jatom, jorb,isym, jcrorb
        LOGICAL :: included,param
        COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: coeff
        COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmp_mat
        COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmp_matbis
        COMPLEX(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: tmp_matn
C
C
C
        WRITE(buf,'(a)')'Creation of the projectors...'
        CALL printout(0)
C
C
C ======================================================================
C Selection of the bands which lie in the chosen energy window [e1;e2]
C or if proj_mode = [1,2] e1 and e2 are directly used as band indices:       
C ======================================================================
C
        
        IF(proj_mode.gt.0) THEN
C The same number of bands are included at each k-point
         kp(:,:)%included=.TRUE.
C Use directly e1 and e2 as band indices. Set to nbmin or
C nbmax if too small or too large
         DO is=1,ns
           DO ik=1,nk
            IF(e1 > kp(ik,is)%nbmin) THEN
             kp(ik,is)%nb_bot=INT(e1)
            ELSE
             kp(ik,is)%nb_bot=kp(ik,is)%nbmin
            ENDIF
            IF(e2 < kp(ik,is)%nbmax) THEN
             kp(ik,is)%nb_top=INT(e2)
            ELSE
             kp(ik,is)%nb_top=kp(ik,is)%nbmax
            ENDIF
           ENDDO
         ENDDO
        ELSE
         kp(:,:)%included=.FALSE.
C the field kp%included = boolean which is .TRUE. when there is at least one band 
C at this k-point whose energy eignevalue is in the energy window.
         DO is=1,ns
           DO ik=1,nk
             included=.FALSE.
             DO ib=kp(ik,is)%nbmin,kp(ik,is)%nbmax
               IF(.NOT.included.AND.kp(ik,is)%eband(ib) > e1.AND.
     &          kp(ik,is)%eband(ib).LE.e2) THEN
C If the energy eigenvalue E of the band ib at the k-point ik is such that e1 < E =< e2, 
C then all the band with ib1>ib must be "included" in the computation and kp%nb_bot is initialized at the value ib.
                included=.TRUE.
                kp(ik,is)%nb_bot=ib
               ELSEIF(included.AND.kp(ik,is)%eband(ib) > e2) THEN
C If the energy eigenvalue E of the current band ib at the k-point ik is such that E > e2 and all the previous 
C band are "included", then the field kp%included = .TRUE. and kp%nb_top = ib-1 (the index of the previous band)
                kp(ik,is)%nb_top=ib-1
                kp(ik,is)%included=.TRUE.
                EXIT
C The loop on the band ib is stopped, since all the bands after ib have an energy > that of ib. 
               ELSEIF(ib==kp(ik,is)%nbmax.AND.kp(ik,is)%eband(ib)
     &          > e1.AND.kp(ik,is)%eband(ib).LE.e2) THEN
C If the energy eigenvalue E of the last band ib=kp%nbmax at the k-point ik is such that e1 < E =< e2 and all the 
C previous bands are "included", then the band ib must be "included" and kp%nb_bot is initialized at the value kp%nbmax.
                kp(ik,is)%nb_top=ib
                kp(ik,is)%included=.TRUE.
               ENDIF
C If the eigenvalues of the bands at the k-point ik are < e1 and included=.FALSE. 
C of if the eigenvalues of the bands at the k-point ik are in the energy window [e1,e2] and included=.TRUE.,
C nothing is done...
             ENDDO    ! End of the ib loop
C If all the eigenvalues of the bands at the k-point ik are not in the window, 
C then kp%included remains at the value .FALSE. and the field kp%nb_top and kp%nb_bot are set to 0
             IF (.not.kp(ik,is)%included) THEN
              kp(ik,is)%nb_bot=0
              kp(ik,is)%nb_top=0 
             ENDIF
           ENDDO      ! End of the ik loop
         ENDDO        ! End of the is loop
        ENDIF
C ---------------------------------------------------------------------------------------
C Checking of the input files if spin-polarized inputs and SO is taken into account:
C There should not be any difference between up and dn limits for each k-point.
C Printing a Warning if this is not the case. 
C -------------------
C
        IF (ifSP.AND.ifSO) THEN
         param=.TRUE.
         DO ik=1,nk
           param=param.AND.(kp(ik,1)%included.eqv.kp(ik,2)%included)
           param=param.AND.(kp(ik,1)%nb_bot==kp(ik,2)%nb_bot)
           param=param.AND.(kp(ik,1)%nb_top==kp(ik,2)%nb_top)
           IF (.not.param) EXIT
C For a valid compoutation, the same k-points must be included for up and dn states, 
C and the upper and lower limits must be the same in both case.
         ENDDO
         IF (.not.param) THEN
          WRITE(buf,'(a,a)')'A Spin-orbit computation for this',
     &      ' compound is not possible with these input files.'
          CALL printout(0)
          WRITE(buf,'(a)')'END OF THE PRGM'
          CALL printout(0)
          STOP
         ENDIF
        ENDIF
C ---------------------------------------------------------------------------------------
C
C
C ==================================================================
C Orthonormalization of the radial wave functions for each orbital :
C ==================================================================
C
C This step is essential for setting the Theta projectors.
        IF(.NOT.ALLOCATED(norm_radf)) THEN
         ALLOCATE(norm_radf(norb))
C norm_radf is a table of "ortfunc" elements, its size ranges from 1 to norb.
         DO iorb=1,norb
           l=orb(iorb)%l
           isrt=orb(iorb)%sort
           norm_radf(iorb)%n=nLO(l,isrt)+2
           n=norm_radf(iorb)%n
           ALLOCATE(norm_radf(iorb)%s12(n,n,ns))
C norm_radf%n = size of the matrix s12
C norm_radf%s12 = matrix of size n*n (one for spin up, one for spin down, if necessary)
           DO is=1,ns
             norm_radf(iorb)%s12(1:n,1:n,is)=0d0
             norm_radf(iorb)%s12(1,1,is)=1d0
             norm_radf(iorb)%s12(2,2,is)=u_dot_norm(l,isrt,is)
C Initialization of the matrix norm_radf%s12 for each orbital (l,isrt).
C We remind tha it is assumed that nLO has the value 0 or 1 only !!
             DO ilo=1,nLO(l,isrt)
               norm_radf(iorb)%s12(2+ilo,2+ilo,is)=1d0
               norm_radf(iorb)%s12(2+ilo,1,is)=
     =           ovl_LO_u(ilo,l,isrt,is)
               norm_radf(iorb)%s12(1,2+ilo,is)=
     =           ovl_LO_u(ilo,l,isrt,is)
               norm_radf(iorb)%s12(2+ilo,2,is)=
     =           ovl_LO_udot(ilo,l,isrt,is)
               norm_radf(iorb)%s12(2,2+ilo,is)=
     =           ovl_LO_udot(ilo,l,isrt,is)
             ENDDO
C Computation of the square root of norm_radf:
             CALL orthogonal_r(norm_radf(iorb)%
     &         s12(1:n,1:n,is),n,.FALSE.)
C the field norm_radf%s12 is finally the C matrix described in the tutorial (or in equation (3.63) in my thesis) 
           ENDDO
         ENDDO
        ENDIF
C
C =====================================
C Creation of the projection matrices :
C =====================================
C
        IF(.NOT.ALLOCATED(pr_orb)) THEN
         ALLOCATE(pr_crorb(ncrorb,nk,ns))
         ALLOCATE(pr_orb(norb,nk,ns))
        ENDIF
C pr_crorb = table of "proj_mat" elements for the correlated orbitals (size from 1 to ncrorb, from 1 to nk, from 1 to ns)
C pr_orb = table of "proj_mat_n" elements for all the orbitals (size from 1 to norb, from 1 to nk, from 1 to ns)
        DO is=1,ns
          DO ik=1,nk
C Only the k-points with inlcuded bands are considered for the projectors.
            IF(.NOT.kp(ik,is)%included) CYCLE
C ------------------------------------------------
C Wannier Projectors for the correlated orbitals :
C ------------------------------------------------
            DO icrorb=1,ncrorb
              l=crorb(icrorb)%l
              iatom=crorb(icrorb)%atom
              isrt=crorb(icrorb)%sort
C Case of l=0 :
C -------------
              IF (l==0) THEN
               IF(ALLOCATED(pr_crorb(icrorb,ik,is)%mat)) THEN
                DEALLOCATE(pr_crorb(icrorb,ik,is)%mat)
               ENDIF
               ALLOCATE(pr_crorb(icrorb,ik,is)%
     %           mat(1,kp(ik,is)%nb_bot:kp(ik,is)%nb_top))
C pr_crorb%mat = the projection matrix with 1 line and (nb_top-nb_bot) columns
               DO ib=kp(ik,is)%nb_bot,kp(ik,is)%nb_top
                 pr_crorb(icrorb,ik,is)%mat(1,ib)=
     =             kp(ik,is)%Alm(1,iatom,ib)
                 DO ilo=1,nLO(l,isrt)
                   pr_crorb(icrorb,ik,is)%mat(1,ib)=
     =               pr_crorb(icrorb,ik,is)%mat(1,ib)+
     +               kp(ik,is)%Clm(ilo,1,iatom,ib)*
     *               ovl_LO_u(ilo,l,isrt,is)
                 ENDDO ! End of the ilo loop
               ENDDO   ! End of the ib loop
C prcrorb(icrorb,ik,is)%mat(1,ib)= <ul1(icrorb,1,is)|psi(is,ik,ib)> = Alm+Clm*ovl_LO_u
C 
C Case of any other l :
C ---------------------
              ELSE
               lm=l*l
C Since the correlated orbital is the l orbital, the elements range from l*l+1 to (l+1)^2
C the sum from 0 to (l-1) of m (from -l to l) is l^2. 
               IF(ALLOCATED(pr_crorb(icrorb,ik,is)%mat)) THEN
                DEALLOCATE(pr_crorb(icrorb,ik,is)%mat)
               ENDIF
               ALLOCATE(pr_crorb(icrorb,ik,is)%
     %           mat(-l:l,kp(ik,is)%nb_bot:kp(ik,is)%nb_top))
C pr_crorb%mat = the projection matrix with (2*l+1) lines and (nb_top-nb_bot) columns
               DO m=-l,l
                 lm=lm+1
                 DO ib=kp(ik,is)%nb_bot,kp(ik,is)%nb_top
                   pr_crorb(icrorb,ik,is)%mat(m,ib)=
     =               kp(ik,is)%Alm(lm,iatom,ib)
                   DO ilo=1,nLO(l,isrt)
                     pr_crorb(icrorb,ik,is)%mat(m,ib)=
     =                 pr_crorb(icrorb,ik,is)%mat(m,ib)+
     +                 kp(ik,is)%Clm(ilo,lm,iatom,ib)*
     *                 ovl_LO_u(ilo,l,isrt,is)
                   ENDDO ! End of the ilo loop
                 ENDDO   ! End of the ib loop
               ENDDO     ! End of the m loop
C prcrorb(icrorb,ik,is)%mat(m,ib)= <ul1(icrorb,m,is)|psi(is,ik,ib)> = Alm+Clm*ovl_LO_u
              ENDIF      ! End of the if l=0 if-then-else
            ENDDO        ! End of the icrorb loop
C
C ---------------------------------------
C Theta Projectors for all the orbitals :
C ---------------------------------------
            DO iorb=1,norb
              l=orb(iorb)%l
              n=norm_radf(iorb)%n
              iatom=orb(iorb)%atom
C Case of l=0 :
C -------------
              IF (l==0) THEN
               IF(ALLOCATED(pr_orb(iorb,ik,is)%matn)) THEN 
                DEALLOCATE(pr_orb(iorb,ik,is)%matn)
               ENDIF
               ALLOCATE(pr_orb(iorb,ik,is)%
     %           matn(1,kp(ik,is)%nb_bot:kp(ik,is)%nb_top,n))
               ALLOCATE(coeff(1:n))
C pr_orb%matn = the projection matrix with 1 line and (nb_top-nb_bot) columns for the n (size of s12) coefficients
C coeff = table of size n which will contain the decomposition of the Bloch state |psi_ik,ib,is> 
C as in equation 22 of the tutorial (Alm, Blm, and Clm ) 
               DO ib=kp(ik,is)%nb_bot,kp(ik,is)%nb_top
                 coeff(1)=kp(ik,is)%Alm(1,iatom,ib)
                 coeff(2)=kp(ik,is)%Blm(1,iatom,ib)
                 coeff(3:n)=kp(ik,is)%Clm(1:n-2,1,iatom,ib)
                 coeff=MATMUL(coeff,norm_radf(iorb)%s12(1:n,1:n,is))
C coeff = coefficients c_(j,lm) of the decomposition of the state |psi> in the orthogonalized basis |phi_j> 
C as defined in the tutorial (equation 25)
                 pr_orb(iorb,ik,is)%matn(1,ib,1:n)=coeff(1:n)
               ENDDO
               DEALLOCATE(coeff)
C pr_orb(iorb,ik,is)%matn(m,ib,1:n) is then the Theta projector as defined in equation 26 of the tutorial.
C 
C Case of any other l :
C ---------------------
              ELSE
               lm=l*l
C As the orbital is the l orbital, the elements range from l*l+1 to (l+1)^2
C the sum from 0 to (l-1) of m (from -l to l) is l^2. 
               IF(ALLOCATED(pr_orb(iorb,ik,is)%matn)) THEN 
                DEALLOCATE(pr_orb(iorb,ik,is)%matn)
               ENDIF
               ALLOCATE(pr_orb(iorb,ik,is)%
     %           matn(-l:l,kp(ik,is)%nb_bot:kp(ik,is)%nb_top,n))
               ALLOCATE(coeff(1:n))
C pr_orb%matn = the projection matrix with (2*l+1) lines and (nb_top-nb_bot) columns for the n (size of s12) coefficients
C coeff = table of size n which will contain the decomposition of the Bloch state |psi_ik,ib,is> 
C as in equation 22 of the tutorial (Alm, Blm, and Clm ) 
               DO m=-l,l
                 lm=lm+1
                 DO ib=kp(ik,is)%nb_bot,kp(ik,is)%nb_top
                   coeff(1)=kp(ik,is)%Alm(lm,iatom,ib)
                   coeff(2)=kp(ik,is)%Blm(lm,iatom,ib)
                   coeff(3:n)=kp(ik,is)%Clm(1:n-2,lm,iatom,ib)
                   coeff=MATMUL(coeff,
     &               norm_radf(iorb)%s12(1:n,1:n,is))
C coeff = coefficients c_(j,lm) of the decomposition of the state |psi> in the orthogonalized basis |phi_j> 
C as defined in the tutorial (equation 25)
                   pr_orb(iorb,ik,is)%matn(m,ib,1:n)=coeff(1:n)
                ENDDO
               ENDDO  ! End of the m loop
               DEALLOCATE(coeff)
C pr_orb(iorb,ik,is)%matn(m,ib,1:n) is then the Theta projector as defined in equation 26 of the tutorial.
              ENDIF      ! End of the if l=0 if-then-else
            ENDDO    ! End of the iorb loop
C 
          ENDDO      ! End of the loop on ik
        ENDDO        ! End of the loop on is
C
C
C ==========================================================================
C Multiplication of the projection matrices by the local rotation matrices :
C ==========================================================================
C
C ------------------------------------------------
C Wannier Projectors for the correlated orbitals :
C ------------------------------------------------
C
        DO jcrorb=1,ncrorb
          jatom=crorb(jcrorb)%atom
          isrt=crorb(jcrorb)%sort
          l=crorb(jcrorb)%l
C
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
          IF (l==0) THEN
C For the s orbital, no multiplication is needed, since the matrix representation of any rotation 
C (and thus Rloc) is always 1.
           DO ik=1,nk
             DO is=1,ns
C Only the k-points with inlcuded bands are considered for the projectors.
               IF(.NOT.kp(ik,is)%included) CYCLE
               nbtop=kp(ik,is)%nb_top
               nbbot=kp(ik,is)%nb_bot
               IF(ALLOCATED(pr_crorb(jcrorb,ik,is)%mat_rep)) THEN
                DEALLOCATE(pr_crorb(jcrorb,ik,is)%mat_rep)
               ENDIF
               ALLOCATE(pr_crorb(jcrorb,ik,is)
     &           %mat_rep(1,nbbot:nbtop))
               pr_crorb(jcrorb,ik,is)%mat_rep(1,nbbot:nbtop)=
     =           pr_crorb(jcrorb,ik,is)%mat(1,nbbot:nbtop)
C As a result, prcrorb%matrep = prcrorb%mat
             ENDDO
           ENDDO
C
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) )
C ---------------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C If this option is used, then ifSO=.TRUE. (because of the restriction in set_ang_trans.f)
C Moreover ifSP=.TRUE. (since ifSO => ifSP in this version) 
C As a result, we know that nb_bot(up)=nb_bot(dn) and nb_top(up)=nb_top(dn)
           DO ik=1,nk
C Only the k-points with inlcuded bands are considered for the projectors.
             IF(.NOT.kp(ik,1)%included) CYCLE
             nbbot=kp(ik,1)%nb_bot
             nbtop=kp(ik,1)%nb_top
C In this case, the projection matrix will be stored in prcrorb%matrep with is=1.
             IF(ALLOCATED(pr_crorb(jcrorb,ik,1)%mat_rep)) THEN
              DEALLOCATE(pr_crorb(jcrorb,ik,1)%mat_rep)
             ENDIF
             ALLOCATE(pr_crorb(jcrorb,ik,1)%
     %         mat_rep(1:2*(2*l+1),nbbot:nbtop))
C The element prcrorb%matrep for is=2 is set to 0, since all the matrix will be stored in the matrix matrep for is=1 
             IF(.not.ALLOCATED(pr_crorb(jcrorb,ik,2)%mat_rep)) THEN
              ALLOCATE(pr_crorb(jcrorb,ik,2)%mat_rep(1,1))
              pr_crorb(jcrorb,ik,2)%mat_rep(1,1)=0.d0
             ENDIF 
C Creation of a matrix tmp_mat which "concatenates" up and dn parts of pr_crorb.
             ALLOCATE(tmp_mat(1:2*(2*l+1),nbbot:nbtop))
             tmp_mat(1:(2*l+1),nbbot:nbtop)=
     =         pr_crorb(jcrorb,ik,1)%mat(-l:l,nbbot:nbtop)         
C The first (2l+1) lines are the spin-up part of the projection matrix prcrorb%mat.
C
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if there is no dn part of pr_orb.
C -------------------------
C
             IF(.not.ifSP) THEN
              WRITE(buf,'(a,a,i2,a)')'The projectors on ',
     &          'the dn states are required for isrt = ',isrt,
     &          ' but there is no spin-polarized input files.'
              CALL printout(0)
              WRITE(buf,'(a)')'END OF THE PRGM'
              CALL printout(0)
              STOP
             ENDIF
C ---------------------------------------------------------------------------------------
C
C The last (2l+1) lines are the spin-dn part of the projection matrix prcrorb%mat.
             tmp_mat((2*l+2):2*(2*l+1),nbbot:nbtop)=
     =         pr_crorb(jcrorb,ik,2)%mat(-l:l,nbbot:nbtop)
C
C Multiplication by the local rotation matrix ; Up and dn parts are treated independently
C since in lapw2 (-alm) the coefficients Alm, Blm and Clm were calculated in the local frame 
C but without taking into account the spinor-rotation matrix.
             ALLOCATE(tmp_matbis(1:(2*l+1),nbbot:nbtop))
             tmp_matbis(1:(2*l+1),nbbot:nbtop)=
     =         tmp_mat(1:(2*l+1),nbbot:nbtop)
             CALL rot_projectmat(tmp_matbis,
     &         l,nbbot,nbtop,jatom,isrt)
             tmp_mat(1:(2*l+1),nbbot:nbtop)=
     =         tmp_matbis(1:(2*l+1),nbbot:nbtop)
             tmp_matbis(1:(2*l+1),nbbot:nbtop)=
     =         tmp_mat(2*l+2:2*(2*l+1),nbbot:nbtop)
             CALL rot_projectmat(tmp_matbis,
     &         l,nbbot,nbtop,jatom,isrt)
             tmp_mat(2*l+2:2*(2*l+1),nbbot:nbtop)=
     =         tmp_matbis(1:(2*l+1),nbbot:nbtop)
             DEALLOCATE(tmp_matbis)
C
C Putting pr_crorb in the desired basis associated to (l,isrt)
C 
             pr_crorb(jcrorb,ik,1)%mat_rep(1:2*(2*l+1),nbbot:nbtop)=
     =         MATMUL(reptrans(l,isrt)%transmat
     &         (1:2*(2*l+1),1:2*(2*l+1)),
     &         tmp_mat(1:2*(2*l+1),nbbot:nbtop))
C pr_crorb%mat_rep = proj_{new_i} = reptrans*proj_{lm} = <new_i|lm>*proj_{lm}
             DEALLOCATE(tmp_mat)
           ENDDO    ! End of the ik loop
C
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only)
C --------------------------------------------------------------------------------------------
          ELSE
           DO ik=1,nk
             DO is=1,ns
C Only the k-points with inlcuded bands are considered for the projectors.
               IF(.NOT.kp(ik,is)%included) CYCLE
C In this case, nb_top(up) and nb_bot(up) can differ from nb_top(dn) and nb_bot(dn)
               nbbot=kp(ik,is)%nb_bot
               nbtop=kp(ik,is)%nb_top
               IF(ALLOCATED(pr_crorb(jcrorb,ik,is)%mat_rep)) THEN
                DEALLOCATE(pr_crorb(jcrorb,ik,is)%mat_rep)
               ENDIF
               ALLOCATE(pr_crorb(jcrorb,ik,is)
     &           %mat_rep(-l:l,nbbot:nbtop))
               pr_crorb(jcrorb,ik,is)%mat_rep(-l:l,nbbot:nbtop)=
     =           pr_crorb(jcrorb,ik,is)%mat(-l:l,nbbot:nbtop)
C
C Multiplication by the local rotation matrix 
C since in lapw2 (-alm) the coefficients Alm, Blm and Clm were calculated in the local frame 
               CALL rot_projectmat(pr_crorb(jcrorb,ik,is)
     &           %mat_rep(-l:l,nbbot:nbtop),l,nbbot,nbtop,jatom,isrt)
C
C Putting pr_crorb in the desired basis associated to (l,isrt)
               pr_crorb(jcrorb,ik,is)%mat_rep(-l:l,nbbot:nbtop)=
     =           MATMUL(reptrans(l,isrt)%transmat(-l:l,-l:l),
     &           pr_crorb(jcrorb,ik,is)%mat_rep(-l:l,nbbot:nbtop))
C pr_crorb%mat_rep = proj_{new_i} = reptrans*proj_{lm} = <new_i|lm>*proj_{lm}
             ENDDO  ! End of the is loop
           ENDDO    ! End of the ik loop
          ENDIF     ! End of the if mixing if-then-else
        ENDDO       ! End of the jcrorb loop
C
C ---------------------------------------
C Theta Projectors for all the orbitals :
C ---------------------------------------
C
        DO jorb=1,norb
          jatom=orb(jorb)%atom
          isrt=orb(jorb)%sort
          n=norm_radf(jorb)%n
          l=orb(jorb)%l
C
C The case l=0 is a particular case of "non-mixing" basis.
C --------------------------------------------------------
          IF (l==0) THEN
C For the s orbital, no multiplication is needed, since the matrix representation of any rotation 
C (and therefore Rloc) is always 1.
           DO ik=1,nk
             DO is=1,ns
C Only the k-points with inlcuded bands are considered for the projectors.
               IF(.NOT.kp(ik,is)%included) CYCLE
               nbtop=kp(ik,is)%nb_top
               nbbot=kp(ik,is)%nb_bot
               IF(ALLOCATED(pr_orb(jorb,ik,is)%matn_rep)) THEN
                DEALLOCATE(pr_orb(jorb,ik,is)%matn_rep)
               ENDIF
               ALLOCATE(pr_orb(jorb,ik,is)%matn_rep
     &           (1,nbbot:nbtop,1:n))
               pr_orb(jorb,ik,is)%matn_rep(1,nbbot:nbtop,1:n)=
     =           pr_orb(jorb,ik,is)%matn(1,nbbot:nbtop,1:n)
C As a result, prorb%matnrep = prorb%matn
             ENDDO
           ENDDO
C
C If the basis representation needs a complete spinor rotation approach (matrices of size 2*(2*l+1) )
C ---------------------------------------------------------------------------------------------------
          ELSEIF (reptrans(l,isrt)%ifmixing) THEN
C If this option is used, then ifSO=.TRUE. (restriction in set_ang_trans.f)
C Moreover ifSP=.TRUE. (since ifSO => ifSP) 
C As a result, we know that nb_bot(up)=nb_bot(dn) and nb_top(up)=nb_top(dn)
           DO ik=1,nk
C Only the k-points with inlcuded bands are considered for the projectors.
             IF(.NOT.kp(ik,1)%included) CYCLE
             nbbot=kp(ik,1)%nb_bot
             nbtop=kp(ik,1)%nb_top
C In this case, the projection matrix will be stored in prorb%matnrep with is=1.
             IF(ALLOCATED(pr_orb(jorb,ik,1)%matn_rep)) THEN
              DEALLOCATE(pr_orb(jorb,ik,1)%matn_rep)
             ENDIF
             ALLOCATE(pr_orb(jorb,ik,1)%
     %         matn_rep(1:2*(2*l+1),nbbot:nbtop,1:n))
C The element prorb%matnrep for is=2 is set to 0, since all the matrix will be stored in the matrix matnrep for is=1 
             IF(.not.ALLOCATED(pr_orb(jorb,ik,2)%matn_rep)) THEN
              ALLOCATE(pr_orb(jorb,ik,2)%matn_rep(1,1,1))
              pr_orb(jorb,ik,2)%matn_rep(1,1,1)=0.d0
             ENDIF 
C Creation of a matrix tmp_matn which "concatenates" up and dn parts of pr_orb
             ALLOCATE(tmp_matn(1:2*(2*l+1),nbbot:nbtop,1:n))
             tmp_matn(1:(2*l+1),nbbot:nbtop,1:n)=
     =         pr_orb(jorb,ik,1)%matn(-l:l,nbbot:nbtop,1:n)
C The first (2l+1) lines are the spin-up part of the projection matrix prorb%matn.
C
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if there is no dn part of pr_orb.
C -------------------------
C
             IF(.not.ifSP) THEN
              WRITE(buf,'(a,a,i2,a)')'The projectors on ',
     &          'the down states are required for isrt = ',isrt,
     &          ' but there is no spin-polarized input files.'
              CALL printout(0)
              WRITE(buf,'(a)')'END OF THE PRGM'
              CALL printout(0)
              STOP
             ENDIF
C ---------------------------------------------------------------------------------------
C
C The last (2l+1) lines are the spin-dn part of the projection matrix prorb%matn.
             tmp_matn(2*l+2:2*(2*l+1),nbbot:nbtop,1:n)=
     =         pr_orb(jorb,ik,2)%matn(-l:l,nbbot:nbtop,1:n)
C
             DO i=1,n
C Multiplication by the local rotation matrix ; Up and dn parts are treated independently
C since in lapw2 (-alm) the coefficients Alm, Blm and Clm were calculated in the local frame 
C but without taking into account the spinor-rotation matrix.
             ALLOCATE(tmp_matbis(1:(2*l+1),nbbot:nbtop))
             tmp_matbis(1:(2*l+1),nbbot:nbtop)=
     =         tmp_matn(1:(2*l+1),nbbot:nbtop,i)
             CALL rot_projectmat(tmp_matbis,
     &         l,nbbot,nbtop,jatom,isrt)
             tmp_matn(1:(2*l+1),nbbot:nbtop,i)=
     =         tmp_matbis(1:(2*l+1),nbbot:nbtop)
             tmp_matbis(1:(2*l+1),nbbot:nbtop)=
     =         tmp_matn(2*l+2:2*(2*l+1),nbbot:nbtop,i)
             CALL rot_projectmat(tmp_matbis,
     &         l,nbbot,nbtop,jatom,isrt)
             tmp_matn(2*l+2:2*(2*l+1),nbbot:nbtop,i)=
     =         tmp_matbis(1:(2*l+1),nbbot:nbtop)
             DEALLOCATE(tmp_matbis)
C Putting pr_orb in the desired basis associated to (l,isrt)
               pr_orb(jorb,ik,1)%matn_rep
     &           (1:2*(2*l+1),nbbot:nbtop,i)=
     =           MATMUL(reptrans(l,isrt)%
     &           transmat(1:2*(2*l+1),1:2*(2*l+1)),
     &           tmp_matn(1:2*(2*l+1),nbbot:nbtop,i))
C pr_orb%matn_rep = proj_{new_i} = reptrans*proj_{lm} = <new_i|lm>*proj_{lm}
             ENDDO   ! End of the i-loop
             DEALLOCATE(tmp_matn)
           ENDDO       ! End of the ik loop
C
C If the basis representation can be reduce to the up/up block (matrices of size (2*l+1) only)
C --------------------------------------------------------------------------------------------
          ELSE
           DO ik=1,nk
             DO is=1,ns
C Only the k-points with inlcuded bands are considered for the projectors.
               IF(.NOT.kp(ik,is)%included) CYCLE
C In this case, nb_top(up) and nb_bot(up) can differ from nb_top(dn) and nb_bot(dn)
               nbbot=kp(ik,is)%nb_bot
               nbtop=kp(ik,is)%nb_top
               IF(ALLOCATED(pr_orb(jorb,ik,is)%matn_rep)) THEN
                DEALLOCATE(pr_orb(jorb,ik,is)%matn_rep)
               ENDIF
               ALLOCATE(pr_orb(jorb,ik,is)%
     &           matn_rep(-l:l,nbbot:nbtop,1:n))
               pr_orb(jorb,ik,is)%matn_rep(-l:l,nbbot:nbtop,1:n)=
     =           pr_orb(jorb,ik,is)%matn(-l:l,nbbot:nbtop,1:n)
C
               DO i=1,n
C Multiplication by the local rotation matrix 
C since in lapw2 (-alm) the coefficients Alm, Blm and Clm were calculated in the local frame 
                 CALL rot_projectmat(pr_orb(jorb,ik,is)
     &             %matn_rep(-l:l,nbbot:nbtop,i),
     &             l,nbbot,nbtop,jatom,isrt)
C Putting pr_orb in the desired basis associated to (l,isrt)
                 pr_orb(jorb,ik,is)%matn_rep(-l:l,nbbot:nbtop,i)=
     =             MATMUL(reptrans(l,isrt)%transmat(-l:l,-l:l),
     &             pr_orb(jorb,ik,is)%matn_rep(-l:l,nbbot:nbtop,i))
C pr_orb%matn_rep = proj_{new_i} = reptrans*proj_{lm} = <new_i|lm>*proj_{lm}
               ENDDO   ! End of the i loop
             ENDDO     ! End of the is loop
           ENDDO       ! End of the ik loop
          ENDIF        ! End of the if mixing if-then-else
        ENDDO          ! End of the jorb loop
C
C
C =============================================================================
C Printing the projectors with k-points 1 and nk in the file fort.18 for test :
C =============================================================================
        DO icrorb=1,ncrorb
          iatom=crorb(icrorb)%atom
          isrt=crorb(icrorb)%sort
          l=crorb(icrorb)%l
          WRITE(18,'()')
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
        DO iorb=1,norb
          iatom=orb(iorb)%atom
          isrt=orb(iorb)%sort
          l=orb(iorb)%l
          n=norm_radf(iorb)%n
          WRITE(18,'()')
          WRITE(18,'(a,i4)') 'iorb = ', iorb
          WRITE(18,'(a,i4,a,i4)') 'isrt = ', isrt, ' l = ', l
          IF (l==0) THEN
          WRITE(18,'(a,i4)') 'ik = ', 1
           DO i=1,n
             WRITE(18,'(i4)') i
             DO ib = kp(1,1)%nb_bot,kp(1,1)%nb_top
               WRITE(18,*) pr_orb(iorb,1,1)%matn_rep(:,ib,i)
                IF (ifSP)
     &         WRITE(18,*) pr_orb(iorb,1,2)%matn_rep(:,ib,i)
               WRITE(18,'()')
             ENDDO
           ENDDO
           WRITE(18,'(a,i4)') 'ik = ', nk
           DO i=1,n
             WRITE(18,'(i4)') i
             DO ib = kp(nk,1)%nb_bot,kp(nk,1)%nb_top
               WRITE(18,*) pr_orb(iorb,nk,1)%matn_rep(:,ib,i)
               IF (ifSP)
     &         WRITE(18,*) pr_orb(iorb,nk,2)%matn_rep(:,ib,i)
               WRITE(18,'()')
             ENDDO
           ENDDO
          ELSEIF(reptrans(l,isrt)%ifmixing) THEN
           DO i=1,n
             WRITE(18,'(i4)') i
             DO ib = kp(1,1)%nb_bot,kp(1,1)%nb_top
               WRITE(18,*) pr_orb(iorb,1,1)%matn_rep(:,ib,i)
               WRITE(18,'()')
             ENDDO
           ENDDO
           WRITE(18,'(a,i4)') 'ik = ', nk
           DO i=1,n
             WRITE(18,'(i4)') i
             DO ib = kp(nk,1)%nb_bot,kp(nk,1)%nb_top
               WRITE(18,*) pr_orb(iorb,nk,1)%matn_rep(:,ib,i)
               WRITE(18,'()')
             ENDDO
           ENDDO
          ELSE
           DO i=1,n
             WRITE(18,'(i4)') i
             DO ib = kp(1,1)%nb_bot,kp(1,1)%nb_top
               WRITE(18,*) pr_orb(iorb,1,1)%matn_rep(:,ib,i)
                IF (ifSP)
     &         WRITE(18,*) pr_orb(iorb,1,2)%matn_rep(:,ib,i)
               WRITE(18,'()')
             ENDDO
           ENDDO
           WRITE(18,'(a,i4)') 'ik = ', nk
           DO i=1,n
             WRITE(18,'(i4)') i
             DO ib = kp(nk,1)%nb_bot,kp(nk,1)%nb_top
               WRITE(18,*) pr_orb(iorb,nk,1)%matn_rep(:,ib,i)
               IF (ifSP)
     &         WRITE(18,*) pr_orb(iorb,nk,2)%matn_rep(:,ib,i)
               WRITE(18,'()')
             ENDDO
           ENDDO
          ENDIF
        ENDDO
C
        RETURN
        END
       
