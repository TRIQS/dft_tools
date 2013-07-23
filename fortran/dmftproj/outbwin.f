
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

        SUBROUTINE outbwin
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine creates the output file case.oubwin             %% 
C %% which contains all the informations for the charge density      %%
C %% self-consistency.                                               %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definition of the variables :
C ----------------------------
        USE almblm_data
        USE common_data
        USE file_names
        USE prnt
        IMPLICIT NONE
        INTEGER :: is, ik, ou
C
        WRITE(buf,'(a)')'Writing the file case.outbwin...'
        CALL printout(0)
C
        DO is=1,ns
C ====================================
C Definition of the file case.oubwin :
C ====================================
C If the computations is spin-polarized, the output file is divided 
C in two files : case.oubwinup and case.oubwindn
          IF(ifSP.AND.is==1) THEN 
            ou=oubwinup
          ELSEIF(ifSP.AND.is==2) THEN
            ou=oubwindn
          ELSE
            ou=oubwin
          ENDIF
C =======================================
C General informations about the system :
C =======================================
C
C Number of k-points in the I-BZ
          WRITE(ou,'(i6)') nk
C Definition of the Spin-orbit flag ifSO
          IF(ifSO) THEN
           WRITE(ou,'(i6)') 1
          ELSE
           WRITE(ou,'(i6)') 0 
          ENDIF
C ====================================================
C Description of the main properties of each k-point :
C ====================================================
          DO ik=1,nk
C Description of the if-included flag
            IF(kp(ik,is)%included) THEN
             WRITE(ou,'(i6)') 1
            ELSE
             WRITE(ou,'(i6)') 0
            ENDIF
            IF(kp(ik,is)%included) THEN
C Range of bands included at each k-point
             WRITE(ou,'(2(i6))') kp(ik,is)%nb_bot,kp(ik,is)%nb_top
C Weight associated to each k-point (for the simple point integration)
             WRITE(ou,*) kp(ik,is)%weight
            ENDIF
          ENDDO   ! End of the ik loop
        ENDDO     ! End of the is loop
C
        RETURN
        END 
        
        

