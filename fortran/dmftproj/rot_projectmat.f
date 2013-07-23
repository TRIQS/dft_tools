
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

       SUBROUTINE rot_projectmat(mat,l,bottom,top,jatom,isrt)
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine makes the transformation from local to global   %%
C %% frame coordinates for the matrices mat in agreement with        %%
C %% the atom j considered.                                          %%
C %%                                                                 %%
C %% mat SHOULD BE IN THE COMPLEX SPHERICAL HARMONICS BASIS.         %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definiton of the variables :
C ----------------------------
       USE almblm_data, ONLY : nk
       USE common_data
       USE symm
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: l, bottom, top, jatom, isrt
       COMPLEX(KIND=8), DIMENSION(-l:l,bottom:top) :: mat
       COMPLEX(KIND=8), DIMENSION(-l:l,bottom:top) :: mattmp
       COMPLEX(KIND=8), DIMENSION(1:2*l+1,1:2*l+1) :: rot_dmat
       INTEGER :: is, ik, isym, lm, lms, ind1, ind2, m
C
       DO m=-l,l
         mattmp(m,bottom:top)= mat(m,bottom:top)
       END DO
C mat is the projector in the local frame (spherical harmonic basis).
C
C The subroutine lapw2 has actually made the computation in the local frame 
C BUT with considering the up and the dn elements in the global frame (no rotation in spin-space), 
C That's why we have to make the computation only in the spin-space to put entirely the matrix mat in the global frame.
C Moreover, no time-reversal symmetry should be taken into account, since the true "rotloc" matrix is considered in lapw2 (-alm). 
C
C The transformation is thus simply achieved by performing the multiplication by rotloc = <x_global | x_local >
C (use of the subroutine dmat)
       rot_dmat=0.d0
       CALL dmat(l,rotloc(jatom)%a,rotloc(jatom)%b,
     &   rotloc(jatom)%g,
     &   REAL(rotloc(jatom)%iprop,KIND=8),rot_dmat,2*l+1)
C Performing the rotation
        mattmp(-l:l,bottom:top)=
     =    MATMUL(rot_dmat(1:2*l+1,1:2*l+1),
     &    mattmp(-l:l,bottom:top))
C The variable mattmp is then the projector in the global frame (spherical harmonic basis).
C The resulting matrix is stored in mat.
        mat(-l:l,bottom:top)=mattmp(-l:l,bottom:top) 
C
        RETURN
        END

