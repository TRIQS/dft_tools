
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

        SUBROUTINE readcomline
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine reads and process the command line options      %%
C %% (Only -so, -sp and -band are the possible ones).                %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
        USE common_data
        USE prnt
        IMPLICIT NONE
        CHARACTER(len=100) :: buf1
        CHARACTER(len=100), DIMENSION(:), ALLOCATABLE :: flags
        INTEGER :: i1, iargc, iarg
        LOGICAL :: ifError
C Process the command line :
C ---------------------------
        iarg=iargc()
        ALLOCATE(flags(iarg))
        ifSP=.FALSE.
        ifSO=.FALSE.
        ifBAND=.FALSE.
        ifError=.FALSE.
        DO i1=1,iarg
          CALL getarg(i1,buf1)
          READ(buf1,*)flags(i1)
          flags(i1)=ADJUSTL(flags(i1))
          CALL makelowcase(flags(i1))
          SELECT CASE(flags(i1)(1:5))
            CASE('-sp  ')
              ifSP=.TRUE.
            CASE('-so  ')
              ifSO=.TRUE.
            CASE('-band')
              ifBAND=.TRUE.
            CASE DEFAULT
              ifError=.TRUE.
              EXIT
          END SELECT
        ENDDO
C ---------------------------------------------------------------------------------------
C Interruption of the prgm if the options are not recognized.
C -------------------------
C
        IF (ifError) THEN
         WRITE(6,'(a,a,a)')'Command line option: ',flags(i1)(1:5),
     &     ' is not recognized.'
         WRITE(6,'(a)')'END OF THE PRGM'
         STOP
        ENDIF
C ---------------------------------------------------------------------------------------
        DEALLOCATE(flags)
C
        RETURN
        END

      SUBROUTINE makelowcase(string)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine modifies the input string into low case         %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
      IMPLICIT NONE
      CHARACTER*26 :: upperabc
      CHARACTER*26 :: lowabc
      CHARACTER* (*) string
      INTEGER :: i,k
      PARAMETER(upperabc='ABCDEFHGIJKLMNOPQRSTUVWXYZ')
      PARAMETER(lowabc='abcdefhgijklmnopqrstuvwxyz')
      DO i=1,len(string)
        DO k=1,26
         IF(string(i:i)==upperabc(k:k)) string(i:i)=lowabc(k:k)
        ENDDO
      ENDDO
      RETURN
      END


