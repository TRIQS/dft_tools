
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

        PROGRAM dmftproj
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This prgm computes projections to a local (correlated) set of   %%
C %% orbitals from the set of eigenfunctions obtained with Wien2k.   %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE almblm_data
        USE common_data
        USE file_names
        USE prnt
        USE symm
        USE reps
        IMPLICIT NONE
C
        REAL(KIND=8) :: e_win, e_sum, elecn, qtot, qdum
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Alm_sum, Qlm_sum
        COMPLEX(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: occ_mat
        COMPLEX(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: occ_mat_sym
C
        COMPLEX(KIND=8) :: coff
        COMPLEX(KIND=8),DIMENSION(-3:3,-3:3) :: tmpmat
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: lnreps
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: correps
        INTEGER :: isrt, ie, l, m, isym, jatom
        INTEGER :: lm, ik, ilo, ib, iatom, imu, io
        INTEGER :: idum, i1, i2
        INTEGER :: m1, m2, lm1, lm2
        INTEGER :: is, irep, nbrep
        INTEGER :: iorb, icrorb, nmaxrep
        INTEGER :: paramflag, lcorr
        LOGICAL :: ifcorr
        REAL(KIND=8) :: fdum, rtetr
        REAL(KIND=8),PARAMETER :: Elarge=1d6
        character(len=120) line
C ================================
C Processing of the command line :
C ================================
        CALL readcomline
C ====================================================
C Initialization of the variable ns (number of spin) :
C ====================================================
C If the computation uses spin-polarized input files, ns=2
        ns=1
        IF(ifSP) ns=2
C ===================================
C Opening of the input/output files :
C ===================================
        CALL openfiles
C =========================================
C Reading of the input file case.indmftpr :
C =========================================
        READ(iuinp,*)nsort
C nsort = number of sorts of atom
        ALLOCATE(nmult(0:nsort))
        nmult(0)=0
        READ(iuinp,*)nmult(1:nsort)
C nmult = multiplicity for each sort of atom, table from 1 to nsort
        natom=SUM(nmult(1:nsort))
C natom = total number of atoms in the unit cell
        ALLOCATE(isort(natom))
        iatom=0
        DO isrt=1,nsort
          DO imu=1,nmult(isrt)
            iatom=iatom+1
            isort(iatom)=isrt
          ENDDO
        ENDDO
C isort = table of correspondance iatom -> isort (from 1 to natom)
        READ(iuinp,*)lmax
C lmax = maximal orbital number l for all the atoms
        IF(ifSO) THEN
         nlm=(lmax+1)*(lmax+1)*2
        ELSE
         nlm=(lmax+1)*(lmax+1)
        ENDIF
C nlm = maximal number of matrix elements for an l-orbital
C only doubled when SO because of the up and down independent parts... 
        ALLOCATE(lsort(0:lmax,nsort))
        ALLOCATE(defbasis(nsort))
        ALLOCATE(lnreps(0:lmax,nsort))
        IF(.not.ifSO) THEN
C Spin is a good quantum number and ireps are considered in orbital space only.
         ALLOCATE(correps(2*lmax+1,0:lmax,nsort))
        ELSE
C Spin is not a good quantum number anymore (possibility of basis which mixes up and dn states)
C the ireps are considered in spin+orbital space.
         ALLOCATE(correps(2*(2*lmax+1),0:lmax,nsort))
        ENDIF
        ALLOCATE(ifSOflag(nsort))
        DO isrt=1,nsort
          READ(iuinp,*) defbasis(isrt)%typebasis
          IF (defbasis(isrt)%typebasis(1:8)=='fromfile') THEN
           READ(iuinp,*) defbasis(isrt)%sourcefile
          ELSE
           defbasis(isrt)%sourcefile = 'null'
          ENDIF
C defbasis = table of correspondance isort -> "basistrans" element, table from 1 to nsort
C defbasis(isrt)%typebasis = "cubic", "complex" or "fromfile"
C defbasis(isrt)%sourcefile = the name of the file to read if typebasis="fromfile"
          READ(iuinp,*)lsort(0:lmax,isrt)
          READ(iuinp,*)lnreps(0:lmax,isrt)
C ifcorr is a flag who states if the atomic sort isrt has correlated orbitals.
          ifcorr=.FALSE.
          DO l=0,lmax
            IF (lsort(l,isrt)==2) THEN
             ifcorr=.TRUE.
C If lnreps(l,isrt)=1, the treatment is the same as a 0 value.
C because if the number of irep is 1, this irep will be the correlated one.
C
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if the number of irep is not correct.
C -------------------------
C
             IF (ifSO) THEN 
C With SO, the number of ireps must not exceed 2*(2*l+1).
              IF(lnreps(l,isrt).gt.(2*(2*l+1))) THEN
               WRITE(buf,'(a,a,i2,a,i2,a)')' The number of ireps ',
     &           'considered for l=',l,' and isrt=',isrt,
     &           ' is not possible.'
               CALL printout(0)
               WRITE(buf,'(a)')'END OF THE PRGM'
               CALL printout(0)
               STOP
              ENDIF
             ELSE
C Without SO, the number of ireps must not exceed (2*l+1).
              IF(lnreps(l,isrt).gt.(2*l+1)) THEN
               WRITE(buf,'(a,a,i2,a,i2,a)')' The number of ireps ',
     &           'considered for l=',l,' and isrt=',isrt,
     &           ' is not possible.'
               CALL printout(0)
               WRITE(buf,'(a)')'END OF THE PRGM'
               CALL printout(0)
               STOP
              ENDIF
             ENDIF
C ---------------------------------------------------------------------------------------
C
C The description of the different ireps is considered only if there are more than 1 irep.
C that is to say if lnreps(l,isrt)=2, 3,...
             IF(lnreps(l,isrt)>0) THEN
              READ(iuinp,'(14i1)') correps(1:lnreps(l,isrt),l,isrt)
             ENDIF
            ENDIF
          ENDDO
C The ifSO_flag is read only if there is a correlated orbital for the sort isrt.
          IF (ifcorr) THEN
           READ(iuinp,'(i1)') ifSOflag(isrt)
          ENDIF
        ENDDO
C lsort = index for each orbital (0 : not include / 1 : include / 2 : correlated), table from 0 to lmax, from 1 to nsort
C lnreps = number of irreducible representations for each orbital, table from 0 to lmax, from 1 to nsort (temporary variables)
C correps = index for each irreducible representations of the correlated orbital, table from 1 to lnreps(l,isrt), from 0 to lmax, from 1 to nsort (temporary variable)
C ifSOflag = table of correspondance isort -> optionSO (1 or 0). Only used for isort with correlated orbitals

        READ(iuinp,'(A)',iostat=io) line
C Try reading the energies/bandindices and the proj_mode
        READ(line,*,iostat=io) e_bot, e_top, proj_mode
C If it fails we know that we are dealing with an older version of the indmftpr file
C with only 2 values on the window line. proj_mode = 0.
        IF(io.ne.0) THEN
         proj_mode = 0
         READ(line,*,iostat=io) e_bot, e_top
         IF(io.ne.0) THEN
          WRITE(buf,'(a,a)')' The energy window line',
     &     ' is ill-defined.'
          CALL printout(0)
          STOP
          WRITE(buf,'(a)')'END OF THE PRGM'
         ENDIF
        ENDIF

C ---------------------------------------------------------------------------------------
C proj_mode:
C 0: use energy window for projection
C 1: use all band indices present in the given energy window
C    (same number of bands at all kpoints)
C 2: use given band indices (same number of bands at all kpoints)
C ---------------------------------------------------------------------------------------

C ---------------------------------------------------------------------------------------
C e_bot, e_top : lower/upper energy limits of window (used in mode 0)
C b_bot, b_top : lower/upper band index of window (used in mode 2)
C In mode 1 e_bot/e_top are provided in the input file and then
C translated into b_bot/b_top
C ---------------------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if the energy/band window or proj_mode is not well-defined.
C ---------------------------------------------------------------------------------------
C
        IF((proj_mode.lt.0) .or. (proj_mode.gt.2)) THEN
         WRITE(buf,'(a,a)')' The energy window mode (3rd value)',
     &     ' must be 0,1 or 2.'
         CALL printout(0)
         WRITE(buf,'(a)')'END OF THE PRGM'
         CALL printout(0)
         STOP
        ENDIF
        
        IF(proj_mode==0) THEN
                b_bot = 0
                b_top = 0
        ELSEIF(proj_mode==1) THEN
                b_bot = 1e3
                b_top = 1
        ELSEIF(proj_mode==2) THEN
                b_bot = INT(e_bot)
                b_top = INT(e_top)
                e_bot = 0.0
                e_top = 0.0
        ENDIF

        IF((proj_mode.lt.2) .and. (e_bot.gt.e_top)) THEN
         WRITE(buf,'(a,a)')' The energy window ',
     &     ' is ill-defined.'
         CALL printout(0)
         WRITE(buf,'(a)')'END OF THE PRGM'
         CALL printout(0)
         STOP
        ENDIF

        IF((proj_mode==2) .and. (b_bot.gt.b_top)) THEN
          WRITE(buf,'(a,a)')' The k-point index window ',
     &     ' is ill-defined.'
          CALL printout(0)
          WRITE(buf,'(a)')'END OF THE PRGM'
          CALL printout(0)
          STOP
        ENDIF
        
C Writing in the output file case.outdmftpr the previous informations :
C =====================================================================
        WRITE(buf,'(a,a)')'Welcome in DMFTPROJ: ',
     &    'PROJECTION TO LOCALIZED BASIS'
        CALL printout(1)
        WRITE(buf,'(a,a)')'This prgm will build',
     &    ' the Wannier projectors to the'
        CALL printout(0)
        WRITE(buf,'(a,a)')'localized orbitals of an atom',
     &    ' onto which DMFT will be applied.'
        CALL printout(1)
        WRITE(buf,'(a)')'You are performing a computation'
        CALL printout(0)
C Spin orbit option
        IF(ifSO) THEN
         WRITE(buf,'(a)')'in which Spin-Orbit is included.'
        ELSE
         WRITE(buf,'(a)')'without Spin-Orbit.'
        ENDIF
        CALL printout(0)
C Spin polarized option
        IF(ifSP) THEN
         WRITE(buf,'(a)')'using Spin-Polarized Wien2k input files.'
        ELSE
         WRITE(buf,'(a)')'using Paramagnetic Wien2k input files.'
        ENDIF
        CALL printout(0)
        IF (ifSO.AND.(.not.ifSP)) THEN
         WRITE(buf,'(a,a)')'You must use Spin-Polarized input files',
     &     ' to perform Spin-Orbit computation, with this version.'
         CALL printout(0)
         WRITE(buf,'(a)')'END OF THE PRGM'
         CALL printout(0)
         STOP
        ENDIF
C Printing nsort, nmult
        WRITE(buf,'(a)')'======================================='
        CALL printout(0)
        WRITE(buf,'(a,i3)')'Sorts of atoms = ',nsort
        CALL printout(0)
        WRITE(buf,'(a,50i2)')'Equivalent sites per each sort:',
     &          nmult(1:nsort)
        CALL printout(1)
C
        norb=0
        ncrorb=0
        ALLOCATE(notinclude(1:nsort))
        DO isrt=1,nsort
          WRITE(buf,'(a)')'-------------------------------------'
          CALL printout(0)
          WRITE(buf,'(a,i2,a)')'For the sort ',isrt,' :'
          CALL printout(0)
          notinclude(isrt)=.TRUE.
C Printing the name of the included orbitals for each sort
          DO l=0,lmax
            IF(lsort(l,isrt).NE.0) THEN
             WRITE(buf,'(a,i2,a)')'The orbital l=',l,' is included.'
             CALL printout(0)
             norb=norb+nmult(isrt)
             notinclude(isrt)=.FALSE.
            ENDIF
          ENDDO
C The variable notinclude(isrt) is a boolean which precises whether the sort isrt
C is considered in the pbm. (whether there is at least one lsort(l,isrt) not 0.)
          IF (notinclude(isrt)) THEN 
           WRITE(buf,'(a)')'No orbital is included.'
           CALL printout(0)
           CALL printout(0)
           cycle
C If no orbital of isrt is included, they can't be correlated orbitals.
          END IF
          CALL printout(0)
C Determination of the total number of correlated orbitals for each sort
          DO l=0,lmax
            IF(lsort(l,isrt)==2) THEN
             ncrorb=ncrorb+nmult(isrt)
            ENDIF ! End of the lsort=2 if-then-else
          ENDDO   ! End of the l loop
        ENDDO     ! End of the isrt loop
C norb = total number of included orbitals in the system
C ncrorb = total number of correlated orbitals in the system
C
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if no orbital is included.
C -------------------------
C
        IF (norb==0) THEN
         WRITE(buf,'(a,a)')'You must include at least one orbital.'
         CALL printout(0)
         WRITE(buf,'(a)')'END OF THE PRGM'
         CALL printout(0)
         STOP
        ENDIF
C ---------------------------------------------------------------------------------------
C
C ===========================================================================================
C Initialization of the "orbital-type" tables orb and crorb, tables of size norb and ncrorb :
C ===========================================================================================
        ALLOCATE(orb(norb),crorb(ncrorb))
        iorb=0
        icrorb=0    
        DO isrt=1,nsort
          IF (notinclude(isrt)) cycle
          DO l=0,lmax
            IF(lsort(l,isrt).NE.0) THEN
C -------------------------------
C For all the included orbitals :
C -------------------------------
             DO imu=1,nmult(isrt)
               iatom=SUM(nmult(0:isrt-1))+imu
               iorb=iorb+1
               orb(iorb)%atom=iatom
C the field orb%atom = number of the atom when classified in the order (isort,imult)
               orb(iorb)%sort=isrt
C the field orb%sort = sort of the associated atom
               orb(iorb)%l=l
C the field orb%l = the orbital number l
               IF(imu==1) THEN
                orb(iorb)%first=.TRUE.
               ELSE
                orb(iorb)%first=.FALSE.
               ENDIF
C the field orb%first = boolean (if first_atom of the sort isort or not) 
               IF(lnreps(l,isrt).NE.0) THEN
                orb(iorb)%ifsplit=.TRUE.
               ELSE
                orb(iorb)%ifsplit=.FALSE.
               ENDIF
C the field orb%ifsplit = boolean (if ireps are used or not)
             ENDDO
C
             IF(lsort(l,isrt)==2) THEN
C ---------------------------------
C For all the correlated orbitals :
C ---------------------------------
              DO imu=1,nmult(isrt)
                iatom=SUM(nmult(0:isrt-1))+imu
                icrorb=icrorb+1
                crorb(icrorb)%atom=iatom
C the field crorb%atom = number of the atom when classified in the order (isort,imult)
                crorb(icrorb)%sort=isrt
C the field crorb%sort = sort of the associated atom
                crorb(icrorb)%l=l
C the field crorb%l = the orbital number l
                IF(imu==1) THEN
                 crorb(icrorb)%first=.TRUE.
                ELSE
                 crorb(icrorb)%first=.FALSE.
                ENDIF
C the field orb%first = boolean (if first_atom of the sort isort or not) 
                IF(lnreps(l,isrt).NE.0) THEN
                 crorb(icrorb)%ifsplit=.TRUE.
                 ALLOCATE(crorb(icrorb)%correp(lnreps(l,isrt)))
                 crorb(icrorb)%correp=.FALSE.
                 DO irep=1,lnreps(l,isrt)
                   IF(correps(irep,l,isrt)==1) 
     &             crorb(icrorb)%correp(irep)=.TRUE.
                 ENDDO
C the field crorb%correp is defined only when crorb%ifsplit= true
C the field orb%correp = boolean table of size lnreps(l,isrt) : True if the ireps is correlated, False otherwise  
                ELSE
                 crorb(icrorb)%ifsplit=.FALSE.
                ENDIF
C the field orb%ifsplit = boolean (if ireps are used or not)
                IF (ifSOflag(isrt)==1) THEN
                 crorb(icrorb)%ifSOat=1
                ELSE
                 crorb(icrorb)%ifSOat=0
                ENDIF
C the field crorb%ifSOflag = boolean (if SO are used or not)
              ENDDO
             ENDIF   ! End of the lsort=2 if-then-else
            ENDIF    ! End of the lsort>0 if-then-else
          ENDDO      ! End of the l loop
        ENDDO        ! End of the isrt loop
C
C =======================================================================================
C Reading of the transformation matrices from the complex to the required angular basis :
C =======================================================================================
        CALL set_ang_trans
C 
C ======================================================================================
C Comparing data about correlated ireps and the description of transformation matrices :
C ======================================================================================
C
        CALL printout(0)
        CALL printout(0)
        WRITE(buf,'(a)')'======================================='
        CALL printout(0)
        WRITE(buf,'(a)')'Precisions about correlated orbitals.' 
        CALL printout(0)
        CALL printout(0)
        DO isrt=1,nsort
          IF (notinclude(isrt)) cycle
          WRITE(buf,'(a)')'-------------------------------------'
          CALL printout(0)
          WRITE(buf,'(a,i2,a)')'For the sort ',isrt,' :'
          CALL printout(0)
          lcorr=0
          DO l=0,lmax
C Only correlated orbital l of isrt are considered here.
            IF (lsort(l,isrt)==2) THEN
             lcorr=lcorr+1
C If the whole orbital is correlated (lnreps=0 in this case)
             IF (lnreps(l,isrt)==0) THEN
              WRITE(buf,'(a,i2,a)')'The whole orbital l=',l,
     &          ' is included as correlated.'
              CALL printout(0)
C If only one particular irep of the orbital is correlated
             ELSE
C
C For a computation without spin-orbit or a computation with SO and with a basis which mixes up and dn states.
C ------------------------------------------------------------------------------------------------------------
              IF ((.not.ifSO).OR.
     &          (ifSO.AND.(l.NE.0).AND.reptrans(l,isrt)%ifmixing))
     &          THEN
C without SO, the case l=0 can not occur since lnreps(0,isrt)=0.
C
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if the data about ireps are conflicting.
C -------------------------
C
               IF (lnreps(l,isrt).NE.reptrans(l,isrt)%nreps) THEN
                WRITE(buf,'(a,a,i2,a)')
     &            'The number of ireps considered ',
     &            'for the orbital l= ', l ,' is wrong.' 
                CALL printout(0)
                WRITE(buf,'(a)')'END OF THE PRGM'
                CALL printout(0)
                STOP
C ---------------------------------------------------------------------------------------
C
C Writing in the output file case.outdmftpr the irep considered as correlated.
               ELSE
                nbrep=0
                DO irep=1,lnreps(l,isrt)
                  IF (correps(irep,l,isrt)==1) THEN
                   WRITE(buf,'(a,i2,a,i2,a)') 
     &              'The irep ',irep,' of orbital l= ', l,
     &               ' is considered as correlated.'
                   CALL printout(0)
                   nbrep=nbrep+1
                  ENDIF
                ENDDO
C ---------------------------------------------------------------------------------------
C Printing a Warning if more than one irep for one value of l is considered. 
C -------------------
C
                IF (nbrep.gt.1) THEN
                 CALL printout(0)
                 WRITE(buf,'(a,a)') 'WARNING : ',
     &             'more than 1 irep is included as correlated.'
                 CALL printout(0)
                 WRITE(buf,'(a,a,a)') '          ',
     &             'The calculation may not be correct ',
     &             'in this case.'
                 CALL printout(1)
                ENDIF
               ENDIF    ! End of the data-conflict if-then-else
C
C For a computation with spin-orbit with basis which doesn't mix up and dn states.
C --------------------------------------------------------------------------------
              ELSE
               WRITE(buf,'(a,i2,a)')'The whole orbital l=',l,
     &           ' is included as correlated.'
               CALL printout(0)
               WRITE(buf,'(a,a)')'because this computation ',
     &           'includes Spin-Orbit coupling.'
               CALL printout(0)
              ENDIF      ! End of the ifSo if-then-else
             ENDIF       ! End of the lnreps=0 if-then-else
            ENDIF        ! End of the lsort=2 if-then-else
C In the case of no correlated orbitals are considered for the atomic sort isrt :
          ENDDO      ! End of the l loop
          IF (lcorr==0) THEN 
            WRITE(buf,'(a,a)')'No orbital is included as correlated.'
            CALL printout(0)
          ENDIF    ! End of the lcorr=0 if-then-else
        ENDDO        ! End of the isrt loop
        CALL printout(0)
        DEALLOCATE(lnreps,correps)
C lnreps and correps can not be used anymore...
C
C ==================================
C Setting of the symmetry matrices :
C ==================================
        CALL setsym
C
C =========================================================================================
C Reading of the Wien2k informations in the case.almblm file (generated by x lapw2 -almd) : 
C =========================================================================================
C
        CALL printout(0)
        CALL printout(0)
        WRITE(buf,'(a)')'======================================='
        CALL printout(0)
        CALL printout(0)
        WRITE(buf,'(a,a)')'Reading of the file ',almblm_file
        CALL printout(0)
C Reading of the klist_band file if the computation if band oriented (option -band)
        IF(ifBAND) CALL read_k_list
        DO is=1,ns
C If the computation is spin-polarized, there are two differents file (up and down)
          IF(is==2) THEN
           CLOSE(iualmblm)
           OPEN(iualmblm,file=almblm_file_sp2,status='old')
           WRITE(buf,'(a,a)')'Reading of the file ',almblm_file_sp2
           CALL printout(0)
          ENDIF
C -------------------------------------------------------------
C Reading of the general informations in the case.almblm file :
C -------------------------------------------------------------
          READ(iualmblm,*)elecn
          READ(iualmblm,*)nk
          READ(iualmblm,*)nloat
C elecn = total number of semicore+valence electrons in the system
C nk = total number of k_points
C nloat = maximal number of LO (local orbitals in LAPW expansion)
          IF(ifBAND) THEN
           IF (is==1) READ(iuinp,*)eferm
           READ(iualmblm,*)
          ELSE
           READ(iualmblm,*)eferm
          ENDIF
C eferm = fermi level (if the computation is band-oriented, it is read in case.indmftpr)
          IF(is==1) THEN
           ALLOCATE(kp(nk,ns),u_dot_norm(0:lmax,nsort,ns))
           ALLOCATE(ovl_LO_u(nloat,0:lmax,nsort,ns))
           ALLOCATE(ovl_LO_udot(nloat,0:lmax,nsort,ns))
           ALLOCATE(nLO(0:lmax,nsort))
          ENDIF
          nLO=0
          DO isrt=1,nsort
C Beginning of the loop on the sort of atoms (isort)
            
            DO l=0,lmax
              READ(iualmblm,*)u_dot_norm(l,isrt,is)
              READ(iualmblm,*)nLO(l,isrt)
C
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if nLO is more than 1.
C -------------------------
C
              IF (nLO(l,isrt) > 1) THEN
               WRITE(buf,'(a,a)')'The current version of DMFTproj ',
     &           ' cannot be used with more than 1 LO orbital by atom. '
               CALL printout(0)
               WRITE(buf,'(a,i2,a,i2)')
     &           ' This is not the case for the orbital l= ',l,
     &           ' of the atomic sort ',isrt
               CALL printout(0)
               WRITE(buf,'(a)')'END OF THE PRGM'
               CALL printout(0)
               STOP
              ENDIF
C ---------------------------------------------------------------------------------------
C
C It is assumed in the following that nLO is 0 or 1.
              DO ilo=1,nLO(l,isrt)
                READ(iualmblm,*)ovl_LO_u(ilo,l,isrt,is),
     &            ovl_LO_udot(ilo,l,isrt,is)
              ENDDO
            ENDDO
C kp = table of "kp_data" elements. It ranges from 1 to nk and from 1 to ns.
C u_dot_norm(isort,l) = norm <u_dotl1|u_dotl1> for the orbital
C nLO(isort,l) = number of LO (local orbitals) for each orbital of each sort (its value is assumed to be 0 or 1)
C ovl_LO_u(isort, l) = overlap element <ul2|ul1> for the LO orbitals 
C ovl_LO_udot(isort, l) = overlap element <ul2|u-dotl1> for the LO orbitals
C These informations are relative to the basis set for the atomic eigenstates (LAPW-APW expansion)
C
C --------------------------------------------------------------
C For each kpoints and isrt, the "kp_data" elements are filled :
C --------------------------------------------------------------
            DO ik=1,nk
              READ(iualmblm,'()')
              READ(iualmblm,'()')
              READ(iualmblm,*)idum,kp(ik,is)%nbmin,kp(ik,is)%nbmax
C idum = useless variable in case.almblm
C kp(ik,is)%nbmin = index of the lowest band
C kp(ik,is)%nbmzx = index of the uppest band 
              IF(.NOT.ALLOCATED(kp(ik,is)%Alm)) THEN
               ALLOCATE(kp(ik,is)%eband(kp(ik,is)
     &           %nbmin:kp(ik,is)%nbmax))
               ALLOCATE(kp(ik,is)%Alm(nlm,natom,
     &           kp(ik,is)%nbmin:kp(ik,is)%nbmax))
               ALLOCATE(kp(ik,is)%Blm(nlm,natom,
     &           kp(ik,is)%nbmin:kp(ik,is)%nbmax))
               ALLOCATE(kp(ik,is)%Clm(nloat,nlm,natom,
     &           kp(ik,is)%nbmin:kp(ik,is)%nbmax))
               ALLOCATE(kp(ik,is)%tetrweight(kp(ik,is)%nbmin:
     &           kp(ik,is)%nbmax))
              ENDIF
              DO ib=kp(ik,is)%nbmin,kp(ik,is)%nbmax
                READ(iualmblm,*)rtetr,kp(ik,is)%eband(ib)
                kp(ik,is)%tetrweight(ib)=CMPLX(rtetr,0d0)
              ENDDO
C rtetr = tetrahedron weights of the band ib at this kpoint 
C the field kp(ik,is)%eband(ib) = eigenvalues of the ib band at this kpoint
C the field kp(ik,is)%tetrweight(ib) = the tetrahedron weights are set as complex number to avoid problems with SQRT(tetrweight)
              kp(ik,is)%weight=REAL(kp(ik,is)%tetrweight
     &          (kp(ik,is)%nbmin))
C the field kp(ik,is)%weight = value of the tetrahedron weight of the lowest band (fully occupied) at this kpoint -> "a geometric factor" 
              kp(ik,is)%eband=kp(ik,is)%eband-eferm
C the eigenvalues kp(ik,is)%eband are shifted with respect to the fermi level.
C 
C Reading of the Alm, Blm and Clm coefficient
              DO imu=1,nmult(isrt)
                iatom=SUM(nmult(0:isrt-1))+imu
                READ(iualmblm,'()')
                READ(iualmblm,*)idum
                DO ib=kp(ik,is)%nbmin,kp(ik,is)%nbmax
                  lm=0
                  DO l=0,lmax
                    DO m=-l,l
                      lm=lm+1
                      READ(iualmblm,*)kp(ik,is)%Alm(lm,iatom,ib),
     &                              kp(ik,is)%Blm(lm,iatom,ib)
                      DO ilo=1,nLO(l,isrt)
                        READ(iualmblm,*)kp(ik,is)%Clm(ilo,lm,iatom,ib)
                      ENDDO
                    ENDDO    ! End of the m loop
                  ENDDO      ! End of the l loop
                ENDDO        ! End of the ib loop
              ENDDO          ! End of the imu loop
C the field kp(ik,is)%Alm = coefficient A_(lm,ib,iatom)(ik,is) as defined in equation (2.34) of my thesis (equation (??) of the tutorial)
C the field kp(ik,is)%Blm = coefficient B_(lm,ib,iatom)(ik,is) as defined in equation (2.34) of my thesis (equation (??) of the tutorial)
C the field kp(ik,is)%Clm = coefficient C_(ilo,lm,ib,iatom)(ik,is) as defined in equation (2.34) of my thesis (equation (??) of the tutorial)
C Their explicit expression depends of the representation (LAPW or APW). They enable to compute the projectors.
C These values are given for all the orbitals (even those which are not included in the study)
            ENDDO      ! End of the loop on kp
          ENDDO        ! End of the loop on isort
        ENDDO          ! End of the loop on ns (spin)
C End of reading the case.almblm.file
C Printing in the file case.outdmftpr the fermi level (in Rydberg)
        CALL printout(0)
        WRITE(buf,'(a,f10.5,a)')'The value of the Fermi Energy is ',
     &    eferm,' Ry.'
        CALL printout(0)
        WRITE(buf,'(a,a)')'All the considered energies are now given ',
     &    'with respect to this value. (E_Fermi is now 0 Ry)'
        CALL printout(1)
C
C ---------------------------------------------------------------------------------------
C If proj_mode=1 find now the lowest and highes band index
C ---------------------------------------------------------------------------------------
C
        IF(proj_mode==1) THEN
         DO is=1,ns
           DO ik=1,nk
             DO ib=kp(ik,is)%nbmin,kp(ik,is)%nbmax
               IF(kp(ik,is)%eband(ib) > e_bot.AND.
     &          kp(ik,is)%eband(ib).LE.e_top) THEN
                IF(ib.gt.b_top) THEN
                 b_top = ib
                ENDIF
                IF(ib.lt.b_bot) THEN
                 b_bot = ib
                ENDIF             
               ENDIF
             ENDDO    ! End of the ib loop
           ENDDO      ! End of the ik loop
         ENDDO        ! End of the is loop
         e_top = 0.0
         e_bot = 0.0
        ENDIF

C ---------------------------------------------------------------------------------------
C Printing the size of the Energy window
C ---------------------------------------------------------------------------------------

        CALL printout(0)
        IF(proj_mode==0) THEN
         WRITE(buf,'(2(a,f10.5),a)')
     &     'The Eigenstates are projected in an energy window from ',
     &        e_bot,' Ry  to  ',e_top,' Ry around the Fermi level.'
        ELSEIF(proj_mode==1) THEN
         WRITE(buf,'(a,2(a,i3),a)')
     &     'The Eigenstates are projected for the band indices from ',
     &     'band Nr. ', b_bot,' to  ',b_top,'.'
        ELSEIF(proj_mode==2) THEN         
         WRITE(buf,'(a,2(a,i3),a)')
     &     'The Eigenstates are projected for the band indices from ',
     &     'band Nr. ',b_bot,' to  ',b_top,'.'
        ENDIF
        CALL printout(1)
C
C ==============================================================
C Computation of the density matrices up to the Fermi level Ef :
C ==============================================================
C
        WRITE(buf,'(a)')'======================================='
        CALL printout(0)
        WRITE(buf,'(a,a)')'Computation of the Occupancies ',
     &    'and Density Matrices up to E_Fermi'
        CALL printout(1)
C ----------------------------------------
C Setting up the projections for all bands
C ----------------------------------------
        IF(proj_mode==0) THEN
         CALL set_projections(-Elarge,Elarge)
        ELSE
         CALL set_projections(1d0,1d6)
        ENDIF


C Elarge is an energy variable equal to 1.d6 Rydberg (very large !!!)
C
C ---------------------------------------------------------
C Computation of the density matrices and the total charges
C ---------------------------------------------------------
C 
        IF(.NOT.ifBAND) CALL density(.TRUE.,.FALSE.,qdum,.TRUE.)
C For the integration, tetrahedron weights are used. 
C The computation is performed for all the included orbitals
C and the density matrices are printed in the file case.outdmftpr
C qdum is the total charge density. (unused variable) 
C
C The calculation of Wannier projectors is performed only if correlated orbitals are included.
        IF(ncrorb.NE.0) THEN
C 
C ==========================================================================
C Computation of the charge below the lower limit e_bot/b_bot of the window :
C ==========================================================================
C
         WRITE(buf,'(a)')'======================================='
         CALL printout(0)
         IF(proj_mode==0) THEN
          WRITE(buf,'(a,a,f10.5,a)')'Computation of the total ',
     &     'Charge below the lower limit of the energy window :',
     &     e_bot,' Ry'  
         ELSE
          WRITE(buf,'(a,a,i3)')'Computation of the total ',
     &     'Charge below the lower band index Nr. ', b_bot  
         ENDIF
         CALL printout(1)
C
C ----------------------------------------
C Setting up the projections for all bands
C ----------------------------------------
         IF(proj_mode==0) THEN
          CALL set_projections(-Elarge,e_bot)
         ELSE
C set_projections expects REAL(8)
          CALL set_projections(1d0,REAL(b_bot-1,8))
         ENDIF
C
C ---------------------------------------------------------
C Computation of the density matrices and the total charges
C ---------------------------------------------------------
C
         IF(.NOT.ifBAND) CALL density(.FAlSE.,.FALSE.,qtot,.FALSE.)
C A simple point integration is used. 
C The computation is performed for all the included orbitals.
C qtot is the total charge density below e_bot/b_bot. 
C Nothing will be printed in the file case.outdmftpr apart from the total charge qtot.
C
C
C ============================================================
C Computation of the Wannier projectors in the energy window :
C ============================================================
C
         WRITE(buf,'(a)')'======================================='
         CALL printout(0)
         IF(proj_mode==0) THEN
          WRITE(buf,'(a,a,a,f10.5,a,f10.5,a)')'Computation of the ',
     &     'Occupancies and Density Matrices in the desired ',
     &     'energy window [ ',e_bot,'; ',e_top,']' 
         ELSE
          WRITE(buf,'(a,a,a,i3,a,i3,a)')'Computation of the ',
     &     'Occupancies and Density Matrices in the desired ',
     &     'band range [ ',b_bot,'; ',b_top,']' 
         ENDIF
         CALL printout(1)
C
C ----------------------------------------
C Setting up the projections for all bands
C ----------------------------------------
         IF(proj_mode==0) THEN
          CALL set_projections(e_bot,e_top)
         ELSE
          CALL set_projections(REAL(b_bot,8),REAL(b_top,8))
         ENDIF
C
C ------------------------------------------------------------------------------
C Orthonormalization of the projectors for correlated orbitals P(icrorb,ik,is) :
C ------------------------------------------------------------------------------
         IF(ifSO) THEN
C In this case, up and dn states must be orthogonalized together 
C because the spin is not a good quantum number anymore.
          CALL orthogonal_wannier_SO
         ELSE
C In this case, up and dn states can be orthogonalized separately 
          CALL orthogonal_wannier
         ENDIF
C 
C ---------------------------------------------------------
C Computation of the density matrices and the total charges
C ---------------------------------------------------------
C Tetrahedron weights are used, the computation are done for correlated orbitals only and are printed in the outputfile. 
         IF(.NOT.ifBAND) CALL density(.TRUE.,.TRUE.,qdum,.TRUE.)
C For the integration, tetrahedron weights are used. 
C The computation is performed for the correlated orbitals only
C and the density matrices are printed in the file case.outdmftpr
C qdum is the total charge density in the energy window. (unused variable) 
C
C
C Writing the output files for DMFT computations :
C ------------------------------------------------ 
         IF(.NOT.ifBAND) THEN
          CALL outqmc(elecn,qtot)
          CALL outbwin
         ELSE
          CALL outband
         ENDIF
        ENDIF
C End of the prgm
        CALL printout(0)
        WRITE(buf,'(a)')'END OF THE PRGM'
        CALL printout(0)
C
        END

        



           
