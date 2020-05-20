
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

C--------------------
C MODULE almblm_data
C--------------------
       MODULE almblm_data
        INTEGER :: nk, nloat
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: nLO
        REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: u_dot_norm
        REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: ovl_LO_u
        REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: ovl_LO_udot
        TYPE kp_data
          LOGICAL :: included
          INTEGER :: nb_bot, nb_top
          INTEGER :: nbmin,nbmax
          REAL(KIND=8) :: weight
          COMPLEX(KIND=8), DIMENSION(:),  ALLOCATABLE :: tetrweight
          REAL(KIND=8),DIMENSION(:),  ALLOCATABLE :: eband
          COMPLEX(KIND=8),DIMENSION(:,:,:),  ALLOCATABLE :: Alm, Blm
          COMPLEX(KIND=8),DIMENSION(:,:,:,:),  ALLOCATABLE :: Clm
        ENDTYPE
        TYPE(kp_data), DIMENSION(:,:), ALLOCATABLE :: kp
       ENDMODULE almblm_data
C
C--------------
C MODULE bands
C--------------
       MODULE bands
        INTEGER :: nlab, nkband
        TYPE label
          CHARACTER(len=20) :: kname
          INTEGER :: pos
        ENDTYPE
        TYPE(label), DIMENSION(:), ALLOCATABLE :: labels
       ENDMODULE
C
C--------------------
C MODULE common_data
C--------------------
       MODULE common_data
C 11/03/10 : Modification of the fullpath for myDMFTproj-2
C        CHARACTER(len=*), PARAMETER :: wien_path=
C     &    '/workpmc/martins/DMFTprojectors/newDMFTproj'
        CHARACTER(len=250) :: wien_path
        INTEGER :: natom, nsort, lmax, nlm, ns, nsp
        INTEGER, DIMENSION(:), ALLOCATABLE :: isort
        INTEGER, DIMENSION(:), ALLOCATABLE :: nmult
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: lsort
        INTEGER, DIMENSION(:), ALLOCATABLE :: ifSOflag
        INTEGER, DIMENSION(:), ALLOCATABLE :: timeflag
        INTEGER :: b_bot, b_top, proj_mode
        LOGICAL :: ifSO, ifSP, ifBAND
        LOGICAL, DIMENSION(:), ALLOCATABLE :: notinclude
        REAL(KIND=8) :: eferm
        REAL(KIND=8) :: e_bot, e_top

        REAL(KIND=8), PARAMETER :: PI=3.1415926535898d0
C New type structure basistrans
        TYPE deftrans
          CHARACTER(len=8) :: typebasis
C The size of typebasis is limited to 8 characters !
          CHARACTER(len=25) :: sourcefile
C The size of sourcefile is limited to 25 characters !
        ENDTYPE
        TYPE(deftrans), DIMENSION(:), ALLOCATABLE :: defbasis
C Type structure orbital
        TYPE orbital 
          INTEGER :: atom
          INTEGER :: sort
          INTEGER :: l
          LOGICAL :: first
          LOGICAL :: ifsplit
          INTEGER :: ifSOat  
          LOGICAL,DIMENSION(:), ALLOCATABLE :: correp
        ENDTYPE
        TYPE(orbital), DIMENSION(:), ALLOCATABLE :: orb, crorb
        INTEGER :: norb, ncrorb
       ENDMODULE common_data
C
C------------------
C MODULE factorial
C------------------
       MODULE factorial
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fac
        INTEGER :: nfctrl
       CONTAINS
        SUBROUTINE setfact(n)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine sets the factorial array                        %%
C %%        FAC(I+1) = I! for I=0,...,N-1                            %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        IMPLICIT NONE
        INTEGER :: n, i
C
        nfctrl=n
        ALLOCATE(fac(nfctrl))
C     I! = FAC(I+1)
        fac(1)=1.0d00
        DO i=1,nfctrl-1
          fac(i+1)=i*fac(i)
        ENDDO
        RETURN
        END SUBROUTINE setfact
       END MODULE factorial
C
C-------------------
C MODULE file_names
C-------------------
       MODULE file_names
        INTEGER :: iudef, iuinp, iusym, iualmblm, iumatfile, iuradwf
        INTEGER :: iuklist
        INTEGER :: ouproj, ouprn, ouctqmc, oupartial,ousymqmc, ousympar
        INTEGER :: ouband, oubwinup, oubwindn, oubwin
        INTEGER :: outw2kpath
        CHARACTER(len=25) :: jobname
        CHARACTER(len=35) :: inp_file, sym_file, almblm_file 
        CHARACTER(len=35) :: almblm_file_sp2
        CHARACTER(len=35) :: radwf_file, radwf_file_sp2
        CHARACTER(len=35) :: prn_file, ctqmc_file, partial_file
        CHARACTER(len=35) :: klist_file
        CHARACTER(len=35) :: symqmc_file, sympar_file, outband_file 
        CHARACTER(len=35) :: oubwin_file, oubwinup_file, oubwindn_file
        CHARACTER(len=8), PARAMETER :: inp_ext='indmftpr'
        CHARACTER(len=7), PARAMETER :: sym_ext='dmftsym'
        CHARACTER(len=6), PARAMETER :: almblm_ext='almblm'
        CHARACTER(len=8), PARAMETER :: almblmup_ext='almblmup'
        CHARACTER(len=8), PARAMETER :: almblmdn_ext='almblmdn'
        CHARACTER(len=9), PARAMETER :: prn_ext='outdmftpr'
        CHARACTER(len=8), PARAMETER :: ctqmc_ext='ctqmcout'
        CHARACTER(len=7), PARAMETER :: partial_ext='parproj'
        CHARACTER(len=6), PARAMETER :: symqmc_ext='symqmc'
        CHARACTER(len=6), PARAMETER :: sympar_ext='sympar'
        CHARACTER(len=7), PARAMETER :: radwfup_ext='radwfup'
        CHARACTER(len=7), PARAMETER :: radwfdn_ext='radwfdn'
        CHARACTER(len=10), PARAMETER :: klist_ext='klist_band'
        CHARACTER(len=7), PARAMETER :: outband_ext='outband'
        CHARACTER(len=6), PARAMETER :: oubwin_ext='oubwin'
        CHARACTER(len=8), PARAMETER :: oubwinup_ext='oubwinup'
        CHARACTER(len=8), PARAMETER :: oubwindn_ext='oubwindn'
       CONTAINS
        SUBROUTINE set_file_name(filename,exten)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine sets the file name                              %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        IMPLICIT NONE
        CHARACTER(len=*) :: filename, exten
        INTEGER :: i1, i2, i
        i1=LEN_TRIM(jobname)
        i2=LEN(exten)
        i=i1+i2+1
        IF(LEN(filename) < i) THEN
         WRITE(*,'(i3,3a)')
     &     i,' characters required for the $case.',exten,
     &     ' filename, too long'
         STOP
        ENDIF
        filename=' '
        filename(1:i)=jobname(1:i1)//'.'//exten(1:i2)
        END SUBROUTINE set_file_name
C
        SUBROUTINE openfiles
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine opens the input and output units for dmftproj   %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE common_data, ONLY: ifSP, ifSO, ifBAND, wien_path
        IMPLICIT NONE
        CHARACTER(len=120) :: buf
        INTEGER :: i1, i2
C initialize input/output channels
        CALL setchannels
C Get working directory name:
        CALL system('pwd > dir_name.tmp')
        OPEN(outw2kpath,file='dir_name.tmp',status='old')
        READ(outw2kpath,'(a)')buf
        CLOSE(outw2kpath,status='delete')
        i1=INDEX(buf,'/',.TRUE.)
        i2=LEN_TRIM(buf)
        jobname(1:i2-i1)=buf(i1+1:i2)
        jobname(i2-i1+1:)=' '         
C Construct file names
        CALL set_file_name(inp_file,inp_ext)
        CALL set_file_name(sym_file,sym_ext)
        IF(.NOT.ifSP) THEN
         CALL set_file_name(almblm_file,almblm_ext)
        ELSE
         CALL set_file_name(almblm_file,almblmup_ext)
         CALL set_file_name(almblm_file_sp2,almblmdn_ext)
        ENDIF
        CALL set_file_name(prn_file,prn_ext)
        CALL set_file_name(ctqmc_file,ctqmc_ext)
        CALL set_file_name(partial_file,partial_ext)
        CALL set_file_name(symqmc_file,symqmc_ext)
        CALL set_file_name(sympar_file,sympar_ext)
        IF(ifSP.AND.ifSO) THEN
         CALL set_file_name(radwf_file,radwfup_ext)
         CALL set_file_name(radwf_file_sp2,radwfdn_ext)
        ENDIF
        IF(ifBAND) THEN
         CALL set_file_name(klist_file,klist_ext)
         CALL set_file_name(outband_file,outband_ext)
        ENDIF
        IF(ifSP) THEN
         CALL set_file_name(oubwinup_file,oubwinup_ext)
         CALL set_file_name(oubwindn_file,oubwindn_ext)
        ELSE
         CALL set_file_name(oubwin_file,oubwin_ext)
        ENDIF
C Open units      
        OPEN(iuinp,file=inp_file,status='old')
        OPEN(iusym,file=sym_file,status='old')
        OPEN(iualmblm,file=almblm_file,status='old')
        OPEN(ouprn,file=prn_file)
        OPEN(ouctqmc,file=ctqmc_file)
        OPEN(oupartial,file=partial_file)
        OPEN(ousymqmc,file=symqmc_file)
        OPEN(ousympar,file=sympar_file)
        IF(ifBAND) THEN
         OPEN(iuklist,file=klist_file,status='old')
         OPEN(ouband,file=outband_file)
        ENDIF
        IF(ifSP) THEN
          OPEN(oubwinup,file=oubwinup_file)
          OPEN(oubwindn,file=oubwindn_file)
        ELSE
          OPEN(oubwin,file=oubwin_file)
        ENDIF
C
C Set path to Wien2k
        CALL system('echo $WIENROOT > path_wienroot.tmp')
        OPEN(outw2kpath,file='path_wienroot.tmp',status='old')
        READ(outw2kpath,'(a)')wien_path
        CLOSE(outw2kpath,status='delete')
C
        RETURN
        END SUBROUTINE
C
        SUBROUTINE setchannels
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine opens the input and output channels             %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE common_data, ONLY: ifSP
        IMPLICIT NONE
C Channels
C input
        iudef=5        ! def-file
        iuinp=7        ! input data
        iusym=8        ! symmetries
        iualmblm=9        ! almblm matrices from Wien
        iumatfile=15        !transformation matrices between different angular basises 
        iuradwf=16        !radial mesh and wave functions
        iuklist=20        !bands
C output
        ouprn=10        ! print-out file
        ouproj=11      ! projection matrices and other data for DMFT run
        ouctqmc=12      ! output for ctqmc
        oupartial=13    ! output for partial charges projectors for ctqmc
        ousymqmc=14      ! output for permutations and rotation matrices
        ousympar=19      ! output for permutations and rotation matrices
                           ! for partial charges analisis
        ouband=21        ! bands
        IF(ifSP) THEN
         oubwinup=22        ! included bands information for lapw2(up)
         oubwindn=23        ! included bands information for lapw2(dn)
        ELSE
         oubwin=22        ! included bands information for lapw2
        ENDIF
C
        RETURN
        END SUBROUTINE
C
C
       ENDMODULE file_names
C
       MODULE prnt
        CHARACTER(len=250) :: buf
        CONTAINS
        SUBROUTINE printout(newline)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine prints the string in buf to the screen          %%
C %% and to the output file and renitializes buf                     %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE file_names
        IMPLICIT NONE 
        INTEGER :: newline, i
        i=LEN_TRIM(buf)
        WRITE(ouprn,'(a)')buf(1:i)
        WRITE(*,'(a)')buf(1:i)
        buf=' '
        IF(newline==1) THEN
         WRITE(ouprn,'(/)')
         WRITE(*,'(/)')
        ENDIF
        RETURN
        END subroutine 
       ENDMODULE prnt
C
C--------------------
C MODULE projections
C--------------------
       MODULE projections
        TYPE proj_mat 
          COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: mat
          COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: mat_rep
        ENDTYPE
        TYPE(proj_mat), DIMENSION(:,:,:), ALLOCATABLE :: pr_crorb
C
        TYPE proj_mat_n 
          COMPLEX(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: matn
          COMPLEX(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: matn_rep
        ENDTYPE
        TYPE(proj_mat_n), DIMENSION(:,:,:), ALLOCATABLE :: pr_orb
C
        TYPE ortfunc 
          INTEGER :: n
          REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: s12
        ENDTYPE
        TYPE(ortfunc), DIMENSION(:),  ALLOCATABLE :: norm_radf
       ENDMODULE projections
C
C-------------
C MODULE reps
C-------------
       MODULE reps
        TYPE ang_bas
          INTEGER :: nreps
          INTEGER, DIMENSION(:), ALLOCATABLE :: dreps
          LOGICAL :: ifmixing
          COMPLEX(KIND=8),DIMENSION(:,:), ALLOCATABLE :: transmat
        ENDTYPE
        TYPE(ang_bas), DIMENSION(:,:), ALLOCATABLE :: reptrans
       ENDMODULE
C
C-------------
C MODULE symm
C-------------
       MODULE symm
        TYPE matrix
          COMPLEX(KIND=8),DIMENSION(:,:),ALLOCATABLE :: mat
        ENDTYPE
        TYPE symop
          LOGICAL :: timeinv
          INTEGER, DIMENSION(:), ALLOCATABLE :: perm
          INTEGER :: iprop
          REAL(KIND=8) :: a, b, g
          REAL(KIND=8) :: phase
          REAL(KIND=8) :: krotm(3,3)
          COMPLEX(KIND=8),DIMENSION(:,:,:),ALLOCATABLE ::rotl
          TYPE(matrix),DIMENSION(:,:),ALLOCATABLE ::rotrep
        ENDTYPE
        TYPE symoploc
          LOGICAL :: timeinv
          INTEGER :: iprop
          INTEGER :: srotnum
          REAL(KIND=8) :: a, b, g
          REAL(KIND=8) :: phase
          REAL(KIND=8) :: krotm(3,3)
          COMPLEX(KIND=8),DIMENSION(:,:,:),ALLOCATABLE ::rotl
          TYPE(matrix),DIMENSION(:),ALLOCATABLE ::rotrep
        ENDTYPE
        INTEGER :: nsym
        INTEGER :: lsym, nlmsym
        TYPE(symop), DIMENSION(:), ALLOCATABLE :: srot
        TYPE(symoploc), DIMENSION(:), ALLOCATABLE ::  rotloc
        TYPE(matrix), DIMENSION(:,:), ALLOCATABLE :: densmat
        TYPE(matrix), DIMENSION(:,:), ALLOCATABLE :: crdensmat
       END MODULE symm

