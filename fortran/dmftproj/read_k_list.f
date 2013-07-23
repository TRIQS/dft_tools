
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

        SUBROUTINE read_k_list
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine reads the labels of high-symmetry points        %%
C %% along the k-path chosen for plotting the k-resolved spectral    %%
C %% function.                                                       %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ---------------------------
        USE bands
        USE file_names
        IMPLICIT NONE
        CHARACTER(len=100) :: buf
        INTEGER :: ilab, pos, i
C
C ========================================================================
C Determination of the total number of labels and the number of k-points :
C ========================================================================
        buf=' '
        nlab=0
	pos=0
C nlab will count the number of labels met.
C pos will count the number of lines 
C (which is also the number of k-points along the k-path)
        DO WHILE (buf(1:3).NE.'END') 
          READ(iuklist,'(a)') buf
          pos=pos+1
          IF(buf(1:1).NE.' ') THEN
           nlab=nlab+1
          ENDIF
        ENDDO
C pos is now the number of line of the file.
C nlab is the number of labels, includind the label 'END'
C
C nkband = number of k-points along the k-path.
        nkband=pos-1
C The label 'END' must not be taken into account
        nlab=nlab-1
C The last line of the file case.klist_band contains "END".
C So the while loop can have an end too.
C
C =============================
C Determination of the labels :
C =============================
        ALLOCATE(labels(nlab))
C The file case.klist_band is read again.
        REWIND(iuklist)
        ilab=0
        DO pos=1,nkband
          READ(iuklist,'(a)') buf
          IF(buf(1:1).NE.' ') THEN
           ilab=ilab+1
           labels(ilab)%pos=pos
C labels(ilab)%pos is the number of the corresponding k-point 
           i=INDEX(buf,' ')
C determination of the size of buf
C (index is a function which finds the index of ' ' in buf)
           labels(ilab)%kname=' '
           labels(ilab)%kname(1:i)=buf(1:i)
C labels(ilab)%kname is the corresponding label
          ENDIF
        ENDDO
C ======================================
C Printing the labels read for testing :
C ======================================
        WRITE(*,*) nkband
        WRITE(*,*)'nlab = ', nlab
        DO i=1,nlab
          WRITE(*,*) i, labels(i)%pos, labels(i)%kname
        ENDDO
C
        RETURN
        END
        
