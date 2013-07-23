
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

        SUBROUTINE set_ang_trans
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine sets up the matrices for transformation between %%
C %% the default complex spherical harmonics used in Wien2k and an   %%
C %% angular basis chosen, for each orbital of each atom.            %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE common_data
        USE file_names
        USE reps
        USE prnt
        IMPLICIT NONE
        CHARACTER(len=150) :: fullpath
        CHARACTER(len=250) :: buf1
        CHARACTER(len=25) :: basis_file
        CHARACTER(len=1) :: repsign
        INTEGER, DIMENSION(2*(2*lmax+1)) :: degrep
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rtrans,itrans
        INTEGER :: m, l, m1, irep, isrt, ind, ind1, ind2
        COMPLEX(KIND=8),DIMENSION(:,:), ALLOCATABLE :: tempmat
        LOGICAL :: flag
C
C
        WRITE(buf,'(a)')'======================================='
        CALL printout(0)
        WRITE(buf,'(a)')'Basis representation for each sort.'
        CALL printout(0)
        CALL printout(0)
C =================================
C Creation of the reptrans matrix :
C =================================
C 
C For the s-electrons : no transformation is necessary (it's always the scalar 1)
        ALLOCATE(reptrans(1:lmax,1:nsort))
C Definition of the size of reptrans (size lmax*nsort)
C Each element of this table is an "ang_bas" element, which will be defined below.
        DO isrt=1,nsort
C -----------------------------------------------
C Case of a representation in the complex basis :
C -----------------------------------------------
          IF (defbasis(isrt)%typebasis(1:7)=='complex') THEN
           DO l=1,lmax
             IF (lsort(l,isrt)==0) THEN
C The considered orbital is not included, all the fields are set up to default value.
              reptrans(l,isrt)%nreps=1
              ALLOCATE(reptrans(l,isrt)%dreps(1))
              ALLOCATE(reptrans(l,isrt)%transmat(1,1))
              reptrans(l,isrt)%transmat=0d0
              reptrans(l,isrt)%dreps(1)=0
              reptrans(l,isrt)%ifmixing=.FALSE.
             ELSE
C The considered orbital is included.
              reptrans(l,isrt)%nreps=1
              ALLOCATE(reptrans(l,isrt)%dreps(1))
              ALLOCATE(reptrans(l,isrt)%transmat(-l:l,-l:l))
              reptrans(l,isrt)%transmat=0d0
              reptrans(l,isrt)%dreps(1)=2*l+1
              reptrans(l,isrt)%ifmixing=.FALSE.
              DO m=-l,l
                reptrans(l,isrt)%transmat(m,m)=1d0
              ENDDO
C In this case, the transformation matrix is just the Identity (hence 1 irep).
C Spin up and Spin down states are not mixed in the basis representation. 
             ENDIF
           ENDDO
C ---------------------------------------------
C Case of a representation in the cubic basis :
C ---------------------------------------------
          ELSEIF (defbasis(isrt)%typebasis(1:5)=='cubic') THEN
           DO l=1,lmax
             IF (lsort(l,isrt)==0) THEN
C The considered orbital is not included, all the fields are set up to default value.
              reptrans(l,isrt)%nreps=1
              ALLOCATE(reptrans(l,isrt)%dreps(1))
              ALLOCATE(reptrans(l,isrt)%transmat(1,1))
              reptrans(l,isrt)%transmat=0d0
              reptrans(l,isrt)%dreps(1)=0
              reptrans(l,isrt)%ifmixing=.FALSE.
             ELSE
C The considered orbital is included.
C The cubic basis is described in the format transpose(P) where P is the usual matrix 
C of the eigenvectors of a matrix D ( D.P=Delta.P with Delta diagonal or P=<lm|new_i>). 
C In other words, each line of the file describes the coefficient of the "new basis vector"
C in the basis { |l,-l,up>,...|l,l,up>,|l,-l,dn>,...|l,l,dn> }.
C The transformation matrices are stored in the directory SRC_templates, the variable "fullpath" 
C must be updated if this prgm is copied.
              ALLOCATE(reptrans(l,isrt)%transmat(-l:l,-l:l))
              ALLOCATE(rtrans(-l:l))
              ALLOCATE(itrans(-l:l))
C              write(*,*)fullpath
              IF (l==1) CALL 
     &          set_harm_file(fullpath,'case.cf_p_cubic')
C standard cubic representation of p electrons : px,py,pz 
              IF (l==2) CALL 
     &          set_harm_file(fullpath,'case.cf_d_eg_t2g')
C standard cubic representation of d-electrons : dz2, dx2-y2, dxy, dxz,dyz (Wien-convention for the phase)
              IF (l==3) CALL 
     &          set_harm_file(fullpath,'case.cf_f_mm2')
C mm2 representation of the f electrons (standard definition with complex coefficients)
C
C Reading of the file
              OPEN(iumatfile,file=fullpath,status='old') 
              ind=-l  
              irep=0
              DO m=-l,l
                READ(iumatfile,'(a)')buf1
                READ(buf1(1:1),'(a)')repsign
                IF(repsign=='*') THEN
C Finding the different ireps in the new basis (a "*" means the end of an irep)
                 irep=irep+1
                 degrep(irep)=m-ind+1
                 ind=m+1
                ENDIF
                READ(buf1(2:250),*)(rtrans(m1),itrans(m1),m1=-l,l)
C The line of the file is stored in the column of reptrans, which is temporarly "P".
                reptrans(l,isrt)%transmat(-l:l,m)= 
     &            CMPLX(rtrans(-l:l),itrans(-l:l))
              ENDDO
              reptrans(l,isrt)%transmat(-l:l,-l:l)=
     =          TRANSPOSE(CONJG(reptrans(l,isrt)%transmat(-l:l,-l:l)))
C reptrans%transmat = inverse(P) = <new_i|lm>, the transformation matrix from complex basis to the cubic one.
C ( inverse(P) is the decomposition of the complex basis in the new basis...) 
              reptrans(l,isrt)%nreps=irep
              ALLOCATE(reptrans(l,isrt)%dreps(irep))
              reptrans(l,isrt)%dreps(1:irep)=degrep(1:irep)
              reptrans(l,isrt)%ifmixing=.FALSE.
C reptrans%nreps = the total number of ireps in the cubic basis
C reptrans%dreps = table of the size of the different ireps
C reptrans%ifmixing = .FALSE. because Spin up and Spin down states are not mixed in the basis representation.
              CLOSE(iumatfile)
              DEALLOCATE(rtrans)
              DEALLOCATE(itrans)
             ENDIF 
           ENDDO
C ---------------------------------------------------------
C Case of a representation defined in an added input file :
C ---------------------------------------------------------
          ELSEIF (defbasis(isrt)%typebasis(1:8)=='fromfile') THEN
           basis_file=defbasis(isrt)%sourcefile
           OPEN(iumatfile,file=basis_file,status='old')
           DO l=1,lmax
             IF (lsort(l,isrt)==0) THEN
C The considered orbital is not included, all the fields are set up to default value.
              reptrans(l,isrt)%nreps=1
              ALLOCATE(reptrans(l,isrt)%dreps(1))
              ALLOCATE(reptrans(l,isrt)%transmat(1,1))
              reptrans(l,isrt)%transmat=0d0
              reptrans(l,isrt)%dreps(1)=0
             ELSE
C The considered orbital is included.
C The new basis is described in the format transpose(P) where P is the usual matrix 
C of the eigenvectors of a matrix D ( D.P=Delta.P with Delta diagonal or P=<lm|new_i>). 
C In other words, each line of the file describes the coefficient of the "new basis vector"
C in the basis { |l,-l,up>,...|l,l,up>,|l,-l,dn>,...|l,l,dn> }.
C The transformation matrices are stored in the directory SRC_templates, the variable "fullpath" 
C must be updated if this prgm is copied.
              ind=1  
              irep=0
              ALLOCATE(tempmat(1:2*(2*l+1),1:2*(2*l+1)))
              ALLOCATE(rtrans(1:2*(2*l+1)))
              ALLOCATE(itrans(1:2*(2*l+1)))
C 
C Reading of the file
              DO m=1,2*(2*l+1)
                READ(iumatfile,'(a)')buf1
                READ(buf1(1:1),'(a)')repsign
                IF(repsign=='*') THEN
C Finding the different ireps in the new basis (a "*" means the end of an irep)
                 irep=irep+1
                 degrep(irep)=m-ind+1
                 ind=m+1
                ENDIF
                READ(buf1(2:250),*)(rtrans(m1),itrans(m1),
     &            m1=1,2*(2*l+1))
                tempmat(1:2*(2*l+1),m)=
     =            CMPLX(rtrans(1:2*(2*l+1)),itrans(1:2*(2*l+1)))
C The lines of the read matrix are stored in the column of tempmat, which is then P.
              ENDDO
C
C Determination if the basis mixes Spin up and Spin down states 
              flag=.TRUE.
              ind1=1
              ind2=1
C The "do while" loop stops when flag=FALSE or i=2*(l+1)
              DO WHILE (flag.AND.(ind1.lt.2*(l+1)))
                flag=flag.AND.
     &          (tempmat((2*l+1)+ind1,(2*l+1)+ind2)==tempmat(ind1,ind2))
                flag=flag.AND.(tempmat((2*l+1)+ind1,ind2)==0.d0)
                flag=flag.AND.(tempmat(ind1,(2*l+1)+ind2)==0.d0)           
                IF (ind2==(2*l+1)) THEN 
                 ind1=ind1+1
                 ind2=1
                ELSE
                 ind2=ind2+1
                END IF 
              ENDDO
              IF (flag) THEN
C If flag=TRUE (then i=2*l+2), the tempmat matrix is block diagonal in spin with 
C the condition block up/up = block down/down.
C The Spin up and Spin down states are not mixed in the basis representation.
               reptrans(l,isrt)%ifmixing=.FALSE.
C reptrans%ifmixing = .FALSE. because Spin up and Spin down states are not mixed in the basis representation.
C
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if the basis description is not correct.
C -------------------------
C
               IF (SUM(degrep(1:irep/2)).ne.(2*l+1)) THEN
                WRITE(buf,'(a,a,i2,a,i2,a)')'The basis description ',
     &            'for isrt = ',isrt,' and l = ',l,' is not recognized.'
                CALL printout(0)
                WRITE(buf,'(a,a)')'Check the structure of the file ',
     &            defbasis(isrt)%sourcefile
                CALL printout(0)
                WRITE(buf,'(a)')'END OF THE PRGM'
                CALL printout(0)
                STOP
               END IF
C ---------------------------------------------------------------------------------------
C
               ALLOCATE(reptrans(l,isrt)%transmat(-l:l,-l:l))
               reptrans(l,isrt)%transmat(-l:l,-l:l)=
     =           tempmat(1:(2*l+1),1:(2*l+1))
               reptrans(l,isrt)%transmat(-l:l,-l:l)=
     =           TRANSPOSE(CONJG(reptrans(l,isrt)%transmat(-l:l,-l:l)))
C The up/up block is enough to describe the transformation (as for cubic or complex bases)
C reptrans%transmat = inverse(P) = <new_i|lm> 
C inverse(P) is indeed the decomposition of the complex basis in the new basis. 
               reptrans(l,isrt)%nreps=irep/2
               ALLOCATE(reptrans(l,isrt)%dreps(reptrans(l,isrt)%nreps))
               reptrans(l,isrt)%dreps(1:reptrans(l,isrt)%nreps)=
     =           degrep(1:reptrans(l,isrt)%nreps)
C reptrans%nreps = the number of ireps in the desired basis for up spin
C reptrans%dreps = table of the size of the different ireps for up spin
              ELSE
C If flag=FALSE, either the tempmat matrix either mixes Spin up and Spin down states 
C or the representation basis for Spin up and Spin down states differ.
C In this case, it is not possible to reduce the description only to the up/up block.
C The whole tempmat matrix is necessary.
C
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if the basis description is not correct.
C -------------------------
C
               IF (SUM(degrep(1:irep)).ne.(2*(2*l+1))) THEN
                WRITE(buf,'(a,a,i2,a,i2,a)')'The basis description ',
     &            'for isrt = ',isrt,' and l = ',l,' is not recognized.'
                CALL printout(0)
                WRITE(buf,'(a,a)')'Check the structure of the file ',
     &            defbasis(isrt)%sourcefile
                CALL printout(0)
                WRITE(buf,'(a)')'END OF THE PRGM'
                CALL printout(0)
                STOP
               END IF
C ---------------------------------------------------------------------------------------
C
               reptrans(l,isrt)%ifmixing=.TRUE.
C reptrans%ifmixing = .TRUE. because Spin up and Spin down states are mixed in the basis representation.
               ALLOCATE(reptrans(l,isrt)%transmat
     &           (1:2*(2*l+1),1:2*(2*l+1)))
               reptrans(l,isrt)%transmat(1:2*(2*l+1),1:2*(2*l+1))=
     =           tempmat(1:2*(2*l+1),1:2*(2*l+1))
               reptrans(l,isrt)%transmat(1:2*(2*l+1),1:2*(2*l+1))=
     =           TRANSPOSE(CONJG(reptrans(l,isrt)%transmat
     &           (1:2*(2*l+1),1:2*(2*l+1))))
C In this case, reptrans%transmat is a square matrix which ranges from 1 to 2*(2*l+1).
C reptrans%transmat = inverse(P) =  <new_i|lm> 
C inverse(P) is indeed the decomposition of the complex basis in the new basis.
               reptrans(l,isrt)%nreps=irep
               ALLOCATE(reptrans(l,isrt)%dreps(irep))
               reptrans(l,isrt)%dreps(1:irep)=degrep(1:irep)
C reptrans%nreps = the total number of ireps in the desired basis
C reptrans%dreps = table of the size of the different ireps
C 
C Restriction for simplicity in the following (and for physical reasons) : 
C a basis with ifmixing=.TRUE. is allowed only if the computation includes SO.
               IF (.not.ifSO) THEN
                WRITE(buf,'(a,a,i2,a,i2,a)')'The basis description ',
     &            'for isrt = ',isrt,' and l = ',l,
     &            ' mixes up and down states.'
                CALL printout(0)
                WRITE(buf,'(a,a)')'This option can not ',
     &            'be used in a computation without Spin-Orbit.'
                CALL printout(0)
                WRITE(buf,'(a,a)')'Modify the structure of the file ',
     &            defbasis(isrt)%sourcefile
                CALL printout(0)
                WRITE(buf,'(a)')'END OF THE PRGM'
                CALL printout(0)
                STOP
               END IF
              END IF
              DEALLOCATE(tempmat)
              DEALLOCATE(rtrans)
              DEALLOCATE(itrans)
             ENDIF
           ENDDO
           CLOSE(iumatfile)
C ----------------------------------------------
C Case of a wrong definition in the input file :
C ----------------------------------------------
          ELSE
C
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if the file has not the expected structure.
C -------------------------
C
           WRITE(buf,'(a,i2,a)')'The basis description for isrt = ',
     &       isrt,' is not recognized.'
           CALL printout(0)
           WRITE(buf,'(a)')'END OF THE PRGM'
           CALL printout(0)
           STOP
          ENDIF
C ---------------------------------------------------------------------------------------
C
        ENDDO
C
C
C ===============================================
C Printing the basis representation information :
C ===============================================
C
        DO isrt=1,nsort
          IF (notinclude(isrt)) cycle
          CALL printout(0)
          WRITE(buf,'(a)')'-------------------------------------'
          CALL printout(0)
          WRITE(buf,'(a,i2,a)')'For the sort ',isrt,' :'
          CALL printout(0)
          IF (defbasis(isrt)%typebasis(1:7)=='complex') THEN
C -----------------------------------------------
C Case of a representation in the complex basis :
C -----------------------------------------------
           WRITE(buf,'(a,i2,a)')'The atomic sort', isrt,
     &       ' is studied in the complex basis representation.'
           CALL printout(0)
           CALL printout(0)
          ELSEIF (defbasis(isrt)%typebasis(1:5)=='cubic') THEN
C ---------------------------------------------
C Case of a representation in the cubic basis :
C ---------------------------------------------
           WRITE(buf,'(a,i2,a)')'The atomic sort', isrt,
     &       ' is studied in the cubic basis representation.'
           CALL printout(0)
           CALL printout(0)
           DO l=0,lmax
C The considered orbital is not included.
             IF (lsort(l,isrt)==0) cycle
C Case of the s-electrons
             IF (l==0) THEN
              WRITE(buf,'(a,a,(F12.6))')'The basis for s-orbital ',
     &          'is still',1.d0
              CALL printout(0)
             ELSE
C Case of the other orbitals
              WRITE(buf,'(a,i2,a,a,a)')'The basis for orbital l=',l,
     &         ' has the following properties :'
              CALL printout(0)
              WRITE(buf,'(a,i2)')' - number of ireps : ',
     &          reptrans(l,isrt)%nreps
              CALL printout(0)
              WRITE(buf,'(a,14(i2,x))')' - degree of each ireps : ',
     &          reptrans(l,isrt)%dreps(1:reptrans(l,isrt)%nreps)
              CALL printout(0)
              WRITE(buf,'(a,a,a)')'The transformation matrix is block',
     &          ' diagonal in the spin-space. The up/up and down/down',
     &          ' blocks are the same and defined as :'
              CALL printout(0)
C The transformation matrix "P = <lm|new_i>" is displayed.
              DO m=-l,l
                WRITE(buf,'(7(2F12.6),x)')
     &            CONJG(reptrans(l,isrt)%transmat(-l:l,m))
                CALL printout(0)
              ENDDO
              CALL printout(0)
             ENDIF
           ENDDO 
           CALL printout(0)
          ELSE
C ---------------------------------------------------------
C Case of a representation defined in an added input file :
C ---------------------------------------------------------
           WRITE(buf,'(a,i2,a,a,a)')'The atomic sort', isrt,
     &       ' is studied in the basis representation', 
     &       ' defined in the file ',
     &       defbasis(isrt)%sourcefile
           CALL printout(0)
           CALL printout(0)
           DO l=0,lmax
C The considered orbital is not included.
             IF (lsort(l,isrt)==0) cycle
C Case of the s-electrons
             IF (l==0) THEN
              WRITE(buf,'(a,a,(F12.6))')'The basis for s-orbital ',
     &          'is still',1.d0
              CALL printout(0)
              CALL printout(0)
             ELSE
C Case of the other orbitals
              WRITE(buf,'(a,i2,a)')'The basis for orbital l=',l,
     &         ' has the following properties :'
              CALL printout(0)
              WRITE(buf,'(a,i2)')' - number of ireps : ',
     &          reptrans(l,isrt)%nreps
              CALL printout(0)
              WRITE(buf,'(a,14(i2,x))')' - degree of each ireps : ',
     &          reptrans(l,isrt)%dreps(1:reptrans(l,isrt)%nreps)
              CALL printout(0)
              IF (reptrans(l,isrt)%ifmixing) THEN
C If the whole matrix description is necessary.
               WRITE(buf,'(a,a)')'The transformation matrix mixes',
     &           ' up and down states in the spin-space'
               CALL printout(0)
               WRITE(buf,'(a,a)') ' and is defined as : ',
     &           '[ block 1 | block 2 ] with'
               CALL printout(0)
               WRITE(buf,'(a,a)') '                     ',
     &           '[ block 3 | block 4 ]'
               CALL printout(0)
C The transformation matrix "P = <lm|new_i>" is displayed.
               WRITE(buf,'(a,i2,a)') 'For the block 1 :'
               CALL printout(0)
               DO m=1,2*l+1
                 WRITE(buf,'(7(2F12.6),x)') 
     &             CONJG(reptrans(l,isrt)%transmat(1:(2*l+1),m))
                 CALL printout(0)
               ENDDO
               WRITE(buf,'(a,i2,a)') 'For the block 2 :'
               CALL printout(0)
               DO m=1,2*l+1
                 WRITE(buf,'(7(2F12.6),x)') 
     &             CONJG(reptrans(l,isrt)%transmat(2*l+2:2*(2*l+1),m))
                 CALL printout(0)
               ENDDO
               WRITE(buf,'(a,i2,a)') 'For the block 3 :'
               CALL printout(0)
               DO m=2*l+2,2*(2*l+1)
                 WRITE(buf,'(7(2F12.6),x)') 
     &             CONJG(reptrans(l,isrt)%transmat(1:(2*l+1),m))
                 CALL printout(0)
               ENDDO
               WRITE(buf,'(a,i2,a)') 'For the block 4 :'
               CALL printout(0)
               DO m=2*l+2,2*(2*l+1)
                 WRITE(buf,'(7(2F12.6),x)') 
     &             CONJG(reptrans(l,isrt)%
     &             transmat(2*l+2:2*(2*l+1),m))
                 CALL printout(0)
               ENDDO
              ELSE
C If the matrix description can be reduced to its up/up block.
               WRITE(buf,'(a,a,a)')'The transformation matrix is block',
     &           ' diagonal in the spin-space. The up/up and down/down',
     &           ' blocks are the same and defined as :'
               CALL printout(0)
C The transformation matrix "P = <lm|new_i>" is displayed.
               DO m=-l,l
                 WRITE(buf,'(7(2F12.6),x)')
     &             CONJG(reptrans(l,isrt)%transmat(-l:l,m))
                 CALL printout(0)
               ENDDO
              ENDIF   ! End of the ifmixing if-then-else
              CALL printout(0)
             ENDIF    ! End of the l if-then-else
           ENDDO      ! End of the l loop
           CALL printout(0)
          ENDIF       ! End of the basis description if-then-else
        ENDDO         ! End of the isrt loop
C
        RETURN
        END        
        

        SUBROUTINE set_harm_file(fullpath,filename)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine sets the fullpath variable                      %%
C %% Be careful, wien_path is defined in modules.f !!!               %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
         USE common_data, ONLY : wien_path
         USE prnt
         IMPLICIT NONE
         CHARACTER(len=*) :: filename, fullpath
         CHARACTER(len=*), PARAMETER :: dir='SRC_templates'
         INTEGER :: i1, i2, i, i3
C
         i1=LEN_TRIM(wien_path)
         i2=LEN(dir)
         i3=LEN(filename)
         i=i1+i2+i3+2
         IF(LEN(fullpath) < i) THEN
           WRITE(buf,'(a)')
     &     'Characters required for the basis transformation ',
     &     ' filename is too long.'
         CALL printout(0)
         WRITE(buf,'(a)')'END OF THE PRGM'
         CALL printout(0)
         STOP
           STOP
         ENDIF
         fullpath=' '
         fullpath(1:i)=wien_path(1:i1)//'/'//dir//'/'//filename(1:i3)
        END SUBROUTINE set_harm_file
        


