
##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2019 by A. D. N. James, M. Zingl and M. Aichhorn
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################

from types import *
import numpy
from triqs_dft_tools.converters.converter_tools import *
import os.path
from locale import atof


class readElkfiles:
    """
    Read in the Elk output files for the elk_converter.py routines.
    """

    def __init__(self):
        self.dft_file = 'PROJ.OUT'
        self.band_file = 'BAND.OUT'
        self.eval_file = 'EIGVAL.OUT'
        self.efermi_file = 'EFERMI.OUT'
        self.kp_file = 'KPOINTS.OUT'
        self.geom_file='GEOMETRY.OUT'

    def read_elk_file(self, filename, to_replace):
        #this function is almost identical to read_fortran_file, but it removes words after ':'
        """
        Returns a generator that yields all numbers in the Fortran file as float, with possible replacements.

        Parameters
        ----------
        filename : string
                   Name of Fortran-produced file.
        to_replace : dict of str:str
                     Dictionary defining old_char:new_char.

        Returns
        -------
        string
            The next number in file.

        """
        #import os.path
        #import string
        import locale
        from locale import atof
        if not(os.path.exists(filename)):
            raise IOError("File %s does not exist." % filename)
        for line in open(filename, 'r'):
            for old, new in to_replace.items():
                line = line.replace(old, new)
            # removes the substring after ':' character - if it exists
            if ":" in line:
              line = self.split_string(line)
            if "(" in line:
              line = self.split_string3(line)
            for x in line.split():
                yield atof(x)

    def read_elk_file2(self, filename, to_replace):
        #this function is almost identical to read_fortran_file, but it removes words after '('
        """
        Returns a generator that yields all numbers in the Fortran file as float, with possible replacements.

        Parameters
        ----------
        filename : string
                   Name of Fortran-produced file.
        to_replace : dict of str:str
                     Dictionary defining old_char:new_char.

        Returns
        -------
        string
            The next number in file.

        """
        import os.path
        import string
        if not(os.path.exists(filename)):
            raise IOError("File %s does not exist." % filename)
        for line in open(filename, 'r'):
            for old, new in to_replace.items():
                line = line.replace(old, new)
            # removes the substring after '(' character - if it exists
            if "(" in line:
              line = self.split_string2(line,'(')
            # removes the substring after '+' character - if it exists
            if "+" in line:
              line = self.split_string2(line,'+')
            # removes the substring after '|' character - if it exists
            if "|" in line:
              line = self.split_string2(line,'|')
            #only include lines which have multiple characters
            if(len(line)>1):
              yield line.split()

    def split_string(self,line):
        """
        This removes the excess information after ':' of the string read in from the file
        """
        temp = line.split(':')
        line = temp[0]
        return line

    def split_string2(self,line,str):
        """
        This removes the excess information after 'str' of the string read in from the file
        """
        temp = line.split(str)
        line = temp[0]
        return line

    def split_string3(self,line):
        """
        This removes the excess information after '(' of the string read in from the file
        """
        temp = line.split('(')
        line = temp[0]
        return line


    def read_proj(self,dft_file):#,things_to_set):
      """
      This function reads the contents of PROJ.OUT and returns them for general use.
      """
      #read in PROJ.OUT file
      R = self.read_elk_file( dft_file, self.fortran_to_replace)
      try:
        #arrays of entry names
        gen_info_entries = ['nproj', 'n_k', 'spinpol', 'SO','natm']
        proj_entries = ['sort', 'natom', 'l', 'dim']
        at_entries = ['spatom','atom']
        #Read in the No. of projectors, no. of k-points, spin and spinorb info
        gen_info = {name: int(val) for name, val in zip(gen_info_entries, R)}
        #mpi.report(gen_info)
        #read in the information of all the projectors
        #initialized variables
        icorr=0
        lmax=3 #maximum l value that will be used
        n_shells=0
        n_inequiv_shells=0
        #local arrays
        neqatom=[]#numpy.zeros([n_shells], int)
        proj=[]
        shells=[]#numpy.zeros([n_shells], int)
        corr_shells=[]#numpy.zeros([n_shells], int)
        prjtype=[]
        wan=[]
        proj_info=[]#numpy.zeros([n_shells], int)
        T=[]
        basis=[]
        inequiv_to_corr=[]
        corr_to_inequiv=[]
        ind=[]
        proj_idx=[]
        at=[]
        idxlm=[]
        #assign global variables
        n_k=gen_info['n_k']
        SO=gen_info['SO']
        #loop over projector flags
        for ip in range(gen_info['nproj']):
          #projector index
          proj_idx.append(int(next(R)))
          proj.append({name: int(val) for name, val in zip(proj_entries, R)})
          na=proj[ip]['natom'] # integer which reduces when reading in atom info
          while(na>0):
            neqatom.append(int(next(R)))
            na-=neqatom[n_inequiv_shells]
            #appends the mapping array to index of inequivalent atom in shells
            inequiv_to_corr.append(icorr)
            for ia in range(neqatom[n_inequiv_shells]):
              corr_to_inequiv.append(n_inequiv_shells)
              at.append({name: int(val) for name, val in zip(at_entries, R)})
              shells.append(proj[ip].copy())
              shells[n_shells].update({'atom':at[n_shells]['atom'], 'spatom':at[n_shells]['spatom']})
              corr_shells.append(shells[n_shells].copy())
              n_orb=2*shells[n_shells]['l']+1
              #lm submatrix indices
              idxlm.append(numpy.zeros(2*lmax+1, dtype=int))
              nrep=proj[ip]['dim']
              for i in range(nrep):
                idxlm[n_shells][i]=next(R)-1
              ind.append(idxlm[n_shells][0:nrep])
              #determine basis type and transformation matrix
              basis.append(int(next(R)))
              #determine whether which basis the projectors where generated in
              #spherical harmonics
              T.append(numpy.zeros([n_orb, n_orb], dtype=complex))
              #Elk generated unitary basis
              if (basis[n_shells]==2):
                #reads the transformation matrix
                for i in range(n_orb):
                  for j in range(n_orb):
                    T[n_shells][i, j] = next(R)
                for i in range(n_orb):
                  for j in range(n_orb):
                    T[n_shells][i, j] += 1j * next(R)
           #transformation matrix from spherical harmonics to cubic basis
              elif (basis[n_shells]==1):
                T[n_shells]=self.determine_T(shells[n_shells]['l'])
              else:
                for i in range(n_orb):
                  T[n_shells][i,i] = 1.0
           #index for the next inequivalent atom (+1 to index is not needed as this is incorporated in
           #neqatom[ish]
              n_shells+=1
           #increase the inequiv_to_corr value
            icorr+=neqatom[n_inequiv_shells]
           #increase the numer of inequivalent atoms
            n_inequiv_shells+=1
           #end reading the file if read the last projector
          if (ip+1==gen_info['nproj']):
            break
        #determine the irreducible representation
        dim_reps=[]
        n_reps=[]
        irep=[]
        for ish in range(n_inequiv_shells):
          isheq=inequiv_to_corr[ish]
          [n_reps,dim_reps,irep] = self.determine_rep(isheq,ish,corr_shells,basis,ind,n_reps,dim_reps,irep)
        for ish in range(n_shells):
          ishin=corr_to_inequiv[ish]
          corr_shells[ish].update({'SO':SO, 'irep':irep[ishin]})

      except StopIteration:  # a more explicit error if the file is corrupted.
          #raise "Elk_converter : reading PROJ.OUT file failed!"
          raise IOError("Elk_converter : reading PROJ.OUT file failed!")

      R.close()
      #output desired information
      return (gen_info,n_shells,n_inequiv_shells,corr_to_inequiv,inequiv_to_corr,corr_shells,n_reps,dim_reps,ind,basis,T)

    def determine_T(self,l):
        """
        Current version calculates the transformation matrix T to convert the inputs from spherical
        harmonics to the cubic basis (as used in TRIQS and Wien2k).
        This routine is currently very similar to spherical_to_cubic in the TRIQS library.
        This routine can be extended to include other unitary rotation matrices
        """
        from math import sqrt
        size=2*l+1
        T = numpy.zeros((size,size),dtype=complex)
        if(l == 0):
          T[0,0]= 1.0
        elif(l == 1):
          T[0,0] = 1.0/sqrt(2);   T[0,2] = -1.0/sqrt(2)
          T[1,0] = 1j/sqrt(2);    T[1,2] = 1j/sqrt(2)
          T[2,1] = 1.0
        elif l == 2:
          T[0,2] = 1.0
          T[1,0] = 1.0/sqrt(2);   T[1,4] = 1.0/sqrt(2)
          T[2,0] =-1.0/sqrt(2);   T[2,4] = 1.0/sqrt(2)
          T[3,1] = 1.0/sqrt(2);   T[3,3] =-1.0/sqrt(2)
          T[4,1] = 1.0/sqrt(2);   T[4,3] = 1.0/sqrt(2)
        elif l == 3:
        #This needs to be checked
          T[0,0] = 1.0/sqrt(2);    T[0,6] = -1.0/sqrt(2)
          T[1,1] = 1.0/sqrt(2);    T[1,5] = 1.0/sqrt(2)
          T[2,2] = 1.0/sqrt(2);    T[2,4] = -1.0/sqrt(2)
          T[3,3] = 1.0
          T[4,2] = 1j/sqrt(2);   T[4,4] = 1j/sqrt(2)
          T[5,1] = 1j/sqrt(2);   T[5,5] = -1j/sqrt(2)
          T[6,0] = 1j/sqrt(2);   T[6,6] = 1j/sqrt(2)
        else: raise ValueError("determine: implemented only for l=0,1,2,3")
        return numpy.matrix(T)


    def determine_rep(self,ish,ishin,corr_shells,basis,ind,n_reps,dim_reps,irep):
      """
      Determines the irreducible representation used for projection calculation.
      Only for Cubic harmonics at the moment
      """
      #if all of the orbitals were used to construct the projectors
      rep=corr_shells[ish]['dim']
      if(rep==2*corr_shells[ish]['l']+1):
         n_reps.append(1)
         dim_reps.append([corr_shells[ish]['dim']])
         irep.append(1)
      #if a subset of obitals were used for generating the projectors (l>1)
      #if Elk generated T matrix is used
      elif((basis[ish]==2)&(corr_shells[ish]['l']>1)):
         if (corr_shells[ish]['l']==2):
            n_reps.append(2)
            dim_reps.append([3, 2])
            irep.append(dim_reps[ishin].index(rep)+1)
         elif (corr_shells[ish]['l']==3):
            n_reps.append(3)
            dim_reps.append([1, 3, 3])
      #cubic T matrix
      elif((basis[ish]==1)&(corr_shells[ish]['l']>1)):
         if (corr_shells[ish]['l']==2):
            n_reps.append(2)
            dim_reps.append([2, 3])
            irep.append(dim_reps[ishin].index(rep)+1)
         elif (corr_shells[ish]['l']==3):
            n_reps.append(3)
            dim_reps.append([2, 2, 3])
      #determine the dim_reps from the lm indices in ind
            if(rep==3):
              irep.append(3)
            else:
              for i in range(2):
                #look for ind[i]==0
                if(ind[i]==0):
                  irep.append(1)
                  break
                #default to irep=2
                elif(i==1):
                  irep.append(2)

      else:
         raise ValueError("Elk_converter (determine_rep) : irreducible representations were not found!")

      return (n_reps,dim_reps,irep)

    def read_projector(self,shell,n_spin_blocks,ish,proj_mat,ind,T,basis,filext):
        """
        This function reads the contents of WANPROJ_L**_S**_A****.OUT
        """
        import string
        l = str(shell[ish]['l']).zfill(2)
        s = str(shell[ish]['sort']).zfill(2)
        a = str(shell[ish]['spatom']).zfill(4)
        #construct file name from atoms info
        wan_shells_file = 'WANPROJ_L'+l+'_S'+s+'_A'+a+filext
        R = self.read_elk_file(wan_shells_file, self.fortran_to_replace)
        try:
        #no. of kpts and number of orbital
          gen_entries = ['n_k', 'lmmax','dim_rep']
          gen = {name: int(val) for name, val in zip(gen_entries, R)}
          #projector lm size
          dim=gen['lmmax']
          #no. of lm used to construct projector
          dim_rep=gen['dim_rep']
          lat=[]
          n_k=gen['n_k']
          n_orbitals = numpy.zeros([n_k, n_spin_blocks], int)
          band_window = [None for isp in range(n_spin_blocks)]
          for isp in range(n_spin_blocks):
            band_window[isp] = numpy.zeros([n_k, 2], dtype=int)
          for ik in range(0,n_k):
          #k-point index and correlated band window indices
             lat_entries =  ['ik','vklx','vkly','vklz']
             lat.append({name: int(val) for name, val in zip(lat_entries, R)})
             #loop over each spin
             for isp in range(0,n_spin_blocks):
               #band window indices
               proj_entries = ['isp', 'ist_min','ist_max']
               proj_dim={name: int(val) for name, val in zip(proj_entries, R)}
               if(isp+1!=proj_dim['isp']):
                 raise IOError("Elk_converter (",wan_shells_file,") : reading spin projecotrs failed!")
                 return
               #setting band window arrays
               n_orbitals[ik,isp]=proj_dim['ist_max']-proj_dim['ist_min']+1
               band_window[isp][ik, 0] = proj_dim['ist_min']
               band_window[isp][ik, 1] = proj_dim['ist_max']
               #define temporary matrix for reading in the projectors
               mat = numpy.zeros([dim, n_orbitals[ik,isp]], complex)
               # Real part
               for j in range(dim):
                  for i in range(n_orbitals[ik,isp]):
                     mat[j, i] = next(R)
               # Imag part:
               for j in range(dim):
                  for i in range(n_orbitals[ik,isp]):
                     mat[j, i] += 1j * next(R)
               #rotate projectors into basis T
               if(basis[ish]!=0):
                  mat[:,:]=numpy.matmul(T[ish],mat[:,:])

               #put desired projector subset into master array
               proj_mat[ik,isp,ish,0:dim_rep,0:n_orbitals[ik,isp]]=mat[ind[ish],0:n_orbitals[ik,isp]]
               #delete temporary index
               del mat
          #resize array length of band indices to maximum number used
          proj_mat=proj_mat[:,:,:,:,:numpy.max(n_orbitals)]

        except StopIteration:  # a more explicit error if the file is corrupted.
            raise IOError("Elk_converter (read projectors): reading file failed!")
        R.close()
        return(n_orbitals,band_window,dim_rep,proj_mat)

    def read_eig(self,filext=".OUT"):
        """
        This function reads the contents of EIGVAL.OUT and EFERMI.OUT
        """
        import string
        #construct file name from atoms info
        R = self.read_elk_file( self.efermi_file, self.fortran_to_replace)
        try:
          efermi=next(R)
        except StopIteration:  # a more explicit error if the file is corrupted.
            raise IOError("Elk_converter (EFERMI.OUT): reading file failed!")
        eval_file = 'EIGVAL'+filext
        R = self.read_elk_file( eval_file, self.fortran_to_replace)
        try:
        #no. of kpts
          n_k=int(next(R))
        #no. second variation states (no. bands)
          nstsv=int(next(R))
          en=[]
          occ=[]
          kp2=[]
          for ik in range(0,n_k):
            #k index and corresponing lattice coordinates
            k_entries = ['ik', 'vklx','vkly','vklz']
            kp2.append([{name: val for name, val in zip(k_entries, R)}])
            en.append([])
            occ.append([])
            for ist in range(0,nstsv):
            #stores the band index, energy eigenvalues and occupancies
              en_entries = ['ist', 'hopping','occ']
              read_en={name: val for name, val in zip(en_entries, R)}
            #Energy eigenvalues and occupations
              en[ik].append(read_en['hopping']-efermi)
              occ[ik].append(read_en['occ'])

        except StopIteration:  # a more explicit error if the file is corrupted.
            raise IOError("Elk_converter (EIGVAL.OUT): reading file failed!")
        R.close()
        return(en,occ,nstsv)

    def read_kpoints(self,filext=".OUT"):
        """
        This function reads the contents of KPOINTS.OUT
        """
        import string
        kp_file = 'KPOINTS'+filext
        R = self.read_elk_file( kp_file, self.fortran_to_replace)
        try:
        #no. of kpts
          n_k=int(next(R))
          kp=[]
        #reads in the k index, lattice vectors, weights and nmat for each kpt
          #array for bz weights
          bz_weights = numpy.ones([n_k], float) / float(n_k)
          #array for lattice vectors
          vkl = numpy.ones([n_k,3], float)
          for ik in range(n_k):
          #k-grid info
            k_entries = ['ik', 'vklx','vkly','vklz', 'bz_weights', 'nmat']
            kp.append({name: val for name, val in zip(k_entries, R)})
          #convert to integers
            kp[ik]['ik']=int(kp[ik]['ik'])
          #read in BZ weights
            bz_weights[ik]=kp[ik]['bz_weights']
            j=0
            for i in range(1,4):
              vkl[ik,j]=kp[ik][k_entries[i]]
              j=j+1
        except StopIteration:  # a more explicit error if the file is corrupted.
            raise IOError("Elk_converter (KPOINTS.OUT): reading file failed!")
        R.close()
        return(bz_weights,vkl)

    def readsym(self):
        """
        Read in the (crystal) symmetries in lattice coordinates
        """
        dft_file='SYMCRYS.OUT'
        R = self.read_elk_file2( dft_file, self.fortran_to_replace)
        try:
          symmat=[]
          spinmat=[]
          tr=[]
          #read the number of crystal symmetries
          x = next(R)
          nsym = int(atof(x[0]))
          #set up symmetry matrices
          for isym in range(nsym):
            symmat.append(numpy.zeros([3, 3], float))
            spinmat.append(numpy.zeros([3, 3], float))
            tr.append(numpy.zeros([3], float))
          #read each symmetry
          for isym in range(nsym):
            #read the symmetry index and check it
            x = next(R)
            isymm = int(x[3])
            if(isym+1!=isymm):
              raise IOError("Elk_converter : reading symmetries failed!")
            #next read has no useful infomation
            x = next(R)
            #read in the translation vector used for symmetries
            x = next(R)
            for i in range(3):
              tr[isym][i]=atof(x[i])
            #read in the spatial symmetries
            #string with no useful information
            x = next(R)
            #read in the spatial symmetry
            for i in range(3):
              x = next(R)
              for j in range(3):
                symmat[isym][i,j]=int(atof(x[j]))
            #string with no useful information
            x = next(R)
            #read in the spin symmetries
            for i in range(3):
              x = next(R)
              for j in range(3):
                spinmat[isym][i,j]=int(atof(x[j]))

        except StopIteration:  # a more explicit error if the file is corrupted.
          raise IOError("Elk_converter : reading SYMCRYS.OUT file failed!")
        R.close()
        return nsym,spinmat,symmat,tr

    def readlat(self):
        """
        Read in information about the lattice.
        """
        dft_file='LATTICE.OUT'
        R = self.read_elk_file2( dft_file, self.fortran_to_replace)
        try:
          amat = numpy.zeros([3, 3], float)
          amatinv = numpy.zeros([3, 3], float)
          bmat = numpy.zeros([3, 3], float)
          bmatinv = numpy.zeros([3, 3], float)
          #real space lattice matrices
          #cycling through information which is not needed
          for i in range(4):
            x = next(R)
          #reading in the lattice vectors as matrix
          for i in range(3):
            x = next(R)
            for j in range(3):
              amat[i,j] = atof(x[j])
          #cycling through information which is not needed
          x = next(R)
          #reading in the inverse lattice matrix
          for i in range(3):
            x = next(R)
            for j in range(3):
              amatinv[i,j] = atof(x[j])
          #read in cell volume (for transport)
          x = next(R)
          cell_vol = atof(x[-1])    
          #cycling through information which is not needed
          for i in range(4):
            x = next(R)
          #reading in the reciprocal lattice vectors as matrix
          for i in range(3):
            x = next(R)
            for j in range(3):
              bmat[i,j] = atof(x[j])
          #cycling through information which is not needed
          x = next(R)
          #reading in the inverse reciprocal lattice matrix
          for i in range(3):
            x = next(R)
            for j in range(3):
              bmatinv[i,j] = atof(x[j])

        except StopIteration:  # a more explicit error if the file is corrupted.
          raise IOError("Elk_converter : reading PROJ.OUT file failed!")
        R.close()
        return amat,amatinv,bmat,bmatinv,cell_vol

    def read_geometry(self):
      """
      This function reads the contents of GEOMETRY.OUT
      """
      import string
      #read in the Geometery file
      dft_file='GEOMETRY.OUT'
      #atom positions in lattice coordinates
      R = self.read_elk_file2( dft_file, self.fortran_to_replace)
      nspecies=0
      na=[]
      atpos_entries = ['kx', 'ky','kz','Bkx', 'Bky', 'Bkz']
      atpos=[]
      try:
         i=0
         #cycle through file until "atoms" is found
         while(i<1000):
           x = next(R)
           if "atoms" in x:
             break
           elif(i==1000):
             raise IOError("Elk_converter : Could not find 'atoms' in GEOMETRY.OUT!")
           i+=1
         #read in the number of species
         x = next(R)
         ns = int(atof(x[0]))
         #loop over species
         for js in range(ns):
           #Skip species name
           x = next(R)
           #read in the number of atoms per species
           x = next(R)
           na.append(int(atof(x[0])))
           #loop over atomss pre species
           atpos.append([])
           for ia in range(na[js]):
             atpos[js].append(numpy.zeros(6, float))
             x = next(R)
             for j in range(6):
               atpos[js][ia][j]=atof(x[j])
      except StopIteration:  # a more explicit error if the file is corrupted.
        raise IOError("Elk_converter : Could not open GEOMETRY.OUT!")
      R.close()
      return ns, na, atpos

#commented out for now - unsure this will produce DFT+DMFT PDOS
##band character dependent calculations
#    def read_bc(self,fileext):
#        """
#        Read in the ELK generated band characters from BC.OUT
#        """

#        #import string
#        file = 'BC'+fileext
#        R = self.read_elk_file(file, self.fortran_to_replace)
#        try:
#        #no. of kpts and number of orbital
#          gen_entries = ['maxlm', 'nspinor','natmtot','nstsv','nkpt','irep']
#          gen = {name: int(val) for name, val in zip(gen_entries, R)}
#          #projector lm size
#          #check the read in information complies with previous read in data
#          nspinor=self.SP+1
#          if(gen['nspinor'] != nspinor):
#            mpi.report('HDF file nspinor = %s'%nspinor)
#            mpi.report('BC.OUT nspinor = %s'%gen['nspinor'])
#            raise IOError("Elk_converter (",file,") : reading nspinor failed!")
#            return
#          if(gen['natmtot'] != self.n_atoms):
#            raise IOError("Elk_converter (",file,") : reading no. of atoms failed!")
#            return
#          if(gen['nstsv'] != self.nstsv):
#            raise IOError("Elk_converter (",file,") : reading all states failed!")
#            return
#          if(gen['nkpt'] != self.n_k):
#            raise IOError("Elk_converter (",file,") : reading kpoints failed failed!")
#            return
#          if(gen['irep'] == 0):
#            raise IOError("Elk_converter (",file,") : Band characters are in spherical hamonics, may have issues with the PDOS!")
#            return

#          dim=gen['maxlm']
#          lmax=numpy.sqrt(dim)-1
#          bc = numpy.zeros([dim,nspinor,self.n_atoms,self.nstsv,self.n_k], float)

#          for ik in range(0,self.n_k):
#            for iatom in range(0,self.n_atoms):
#              for ispn in range(0,nspinor):
#                entry =  ['ispn','ias','is','ia','ik']
#                ent = {name: int(val) for name, val in zip(entry, R)}
#          #k-point index and correlated band window indices
#          #check read in values
#                if(ent['ispn'] != ispn+1):
#                  raise IOError("Elk_converter (",file,") : reading ispn failed!")
#                  return
#                if(ent['ias'] != iatom+1):
#                  raise IOError("Elk_converter (",file,") : reading iatom failed!")
#                  return
#                if(ent['ik'] != ik+1):
#                  raise IOError("Elk_converter (",file,") : reading ik failed!")
#                  return

#                for ist in range(self.nstsv):
#                  for lm in range(dim):
#                     bc[lm,ispn,iatom,ist,ik] = next(R)

#        except StopIteration:  # a more explicit error if the file is corrupted.
#            raise IOError("Elk_converter (read BC.OUT): reading file failed!")
#        R.close()
#        return(bc,dim)

