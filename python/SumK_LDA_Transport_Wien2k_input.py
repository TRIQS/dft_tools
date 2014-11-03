
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
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
################################################################################

#=======================================================================================================================
# #################################################################
# Code for Transport/Optic calculations based on SumK_LDA... class
# by Xiaoyu Deng <xiaoyu.deng@gmail.com>
# The code read in files needed for transport/Optic calculations from Wien outputs.
# including:  symmetry, velocity, lattice constants.
# The HDF convention is not adopted here since momentum file from Wien output is usually quite large
# and it is not necessary to keep in in the HDF file.
# #################################################################
#=======================================================================================================================


import numpy, sys, os.path

def Read_Fortran_File (filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exists" % filename
    for line in open(filename, 'r') :
	for x in line.replace('D', 'E').split() : 
	    yield string.atof(x)

class Velocity_k:
    """momentum matrix for a single k points.
    """
    def __init__(self):
        self.kp = [0.0, 0.0, 0.0]
        self.bandwin = [999, 0]
        # here a matrix for 3 vels. use list since the size of vel is to be determined.
        self.vel = []


class Velocities:
    """Class containig the velocities.
       Provides container for velocities 
       as a function of k and method to read them from case.pmat Wien2k file
       use as ClassInstance.vks[ik].vel[iband][jband][ix]
    """

    def __init__(self, wiencase, spinbl=""):
        if not(os.path.exists(wiencase + ".pmat" + spinbl)) : raise IOError, "File %s does not exists" % wiencase + ".pmat"

        # expected format:
        # k  nu nu2      (k denotes k-point, nu1 starting band index, nu2 ending band index
        # kpos           ( read in if needed.
        # (real, im)     (of velocity for all nu <=nu_prime within nu1 and nu2)
        #   ...             ...
        # k  nu nu2
        #   ...
        f = open(wiencase + ".pmat" + spinbl)
        self.vks = []

        while 1:
            try:
                s = f.readline()
                if (s == ""):
                    break               
            except:
                break
            try:
                vec = Velocity_k()
                [k, nu1, nu2] = [int (x) for x in s.strip().split()]
                vec.bandwin[0] = nu1
                vec.bandwin[1] = nu2

                vec.kp = f.readline().strip().split()                          
                dim = vec.bandwin[1] - vec.bandwin[0] + 1
                shape = (dim, dim, 3)    
                vxyz = numpy.zeros(shape, dtype=complex)
                for nu in xrange(dim):
                    for nu_prime in xrange(nu, dim):
                        for i in xrange(3):
                            s = f.readline().strip("\n ()").split(',')
                            vxyz[nu][nu_prime][i] = float(s[0]) + float(s[1]) * 1j
                            if(nu_prime != nu):
                                vxyz[nu_prime][nu][i] = vxyz[nu][nu_prime][i].conjugate()
                                                
                vec.vel = vxyz
                self.vks.append(vec)
             
            except IOError: 
                print("Reading case.pmat error. Wrong format?\n ")
                raise
        f.close()    


    def getvel(self, k):
        # should return array at a given k. use as [iband][jband][i]
        return self.vks[k].vel
    
    def plot(self):
        f = open("velband.dat", "w")
        bandid = numpy.array([vk.bandwin for vk in self.vks]).flatten()
        minb = bandid.min()
        maxb = bandid.max()
        for ib in range(minb, maxb + 1):
            ik = 0
            for vk in self.vks:
                if(vk.bandwin[0] <= ib and vk.bandwin[1] >= ib):
                    f.write(str(ik) + "  ")
                    for i in range(3):
                        f.write(str(vk.vel[ib - vk.bandwin[0]][ib - vk.bandwin[0]][i].real) + "  ")
                    f.write("\n") 
                ik = ik + 1
            f.write("&\n")            
        f.close()   
                
    
class SGsymmetry():
    """ read in symmetry of space group from wiencase.outputs other than wiencase.struct since 
        in this file there are also symmetry operations in xyz coordinates..
    """
    def __init__(self, wiencase):
        structfile = wiencase + ".outputs"
        self.nsymm = 1
        self.symm = []
        self.tau = []
        self.bravaismatrix = numpy.zeros((3, 3), dtype=numpy.float_)
        self.symmcartesian = []
        self.taucartesian = []
        with open(structfile, "r") as f:
            f.readline()
            f.readline()# bravais matrix
            for i in range(3):            
                line = f.readline().strip().split()
                self.bravaismatrix[i, :] = numpy.array([numpy.float(item) for item in line])[:]
            print "bravais matrix", self.bravaismatrix            
            
            while 1:
                s = f.readline().strip(" ").split()
                try:
                    if(s[0] == "PGBSYM:"):
                        self.nsymm = int(s[-1])
                        break
                except:
                    continue
            
            f.readline()
            f.readline()
            for i in range(self.nsymm):
                f.readline()
                ## read symmcartesian
                symmt = numpy.zeros((3, 3), dtype=float)
                taut = numpy.zeros((3), dtype=float)
                for ir in range(3):
                    s = f.readline().strip().split()
                    for ic in range(3):
                        symmt[ir, ic] = float(s[ic]) 
                s = f.readline().strip().split()
                for ir in range(3):
                    taut[ir] = float(s[ir])
                    
                self.symmcartesian.append(symmt)
                self.taucartesian.append(taut)
                
                ##read symm
                symmt = numpy.zeros((3, 3), dtype=float)
                taut = numpy.zeros((3), dtype=float)
                for ir in range(3):
                    s = f.readline().strip().split()
                    for ic in range(3):
                        symmt[ir, ic] = float(s[ic]) 
                    taut[ir] = float(s[3])
                    
                self.symm.append(symmt)
                self.tau.append(taut)
                
                f.readline()
            # end
            f.close()
            print "Read wiencase.outputs done!"   

    def checksymmxyz(self):
        '''  This is to check symm in cartesian coordinates and in primitive cell lattice
            For details of symm, one should check wien/SRC_symmetry/, latsym.f, pglsym.f,pgbsym.f
            In general, if a lattice is orthorhombic, then symmcartesian is the same as symmprimitive
            (at most different by a transpose (or inversion)). If a lattice is not orthorhombic, these two symms could 
            be related by bravais matrix. 
            One special case is CXZ lattice. In Wien CXZ lattice contains two cases: orthorhombic, and 
            monoclinic (with only gamma not equal to 90). For CXZ lattice monoclinic, symmcartesian is also the same as
            symmprimitive (or with transpose (or inversion)). 
        '''
        for i in range(self.nsymm):
            mat = self.symmcartesian[i].transpose()  # according to wien2k, why?
            bm = self.bravaismatrix
            bminv = numpy.linalg.inv(bm)#.transpose()
            res = numpy.dot(bminv, numpy.dot(mat, bm)) - self.symm[i]
            print i, res
        
    
    def size(self):
        return self.nsymm

def cellvolume(latticetype, latticeconstants, latticeangle):
    """ calculate cell volume: volumecc conventional cell, volumepc, primitive cell.
    """
    for i in range(3):
        latticeangle[i] *= 1.0 / 180 * numpy.pi
    a = latticeconstants[0]
    b = latticeconstants[1]
    c = latticeconstants[2]
    c_al = numpy.cos(latticeangle[0])
    c_be = numpy.cos(latticeangle[1])
    c_ga = numpy.cos(latticeangle[2])
    volumecc = a * b * c * numpy.sqrt(1 + 2 * c_al * c_be * c_ga - c_al ** 2 - c_be * 82 - c_ga ** 2)
    
    det = {"P":1,
         "F":4,
         "B":2,
         "R":3,
         "H":1,
         "CXY":2,
         "CYZ":2,
         "CXZ":2
         }
    volumepc = volumecc / det[latticetype]
    
    return volumecc, volumepc

    
class WienStruct():
    """ parsing Wien Struct file
    """          
    def __init__(self, wiencase):
        structfile = wiencase + ".struct"
        
        with open(structfile, "r") as infile:
            print "read in Wien case file %s" % structfile
            infile.readline()#title
            
            tem = infile.readline() #lattice
            self.latticetype = tem[0:10].split()[0]
            self.ineqvsite = int(tem[27:30])
            try:
                self.sgrnumber = int(tem[30:33])
            except:
                self.sgrnumber = None
            try:
                self.sgrlabel = tem[34:38]
            except:
                self.sgrlabel = None
            
            print self.latticetype, self.ineqvsite, self.sgrnumber, self.sgrlabel
            infile.readline()
            
            tem = infile.readline() # lattice constants
            self.latticeconstants = [float(tem[0:10]), float(tem[10:20]), float(tem[20:30])]
            self.latticeangle = [float(tem[30:40]), float(tem[40:50]), float(tem[50:60])]
            print "Cell"
            print self.latticeconstants[:]
            print self.latticeangle[:]
            
            self.positions = []
            self.atomsymbols = []
            self.multi = []
            self.atomnumbers = []
            self.locrotmatrix = []
            for isite in range(self.ineqvsite):
                tem = infile.readline()
                positions = []
                positions.append([float(tem[12:22]), float(tem[25:35]), float(tem[38:48])])
                
                tem = infile.readline()
                multi = int(tem[15:17])
                self.multi.append(multi)
                
                for im in range(multi - 1):
                    tem = infile.readline()
                    positions.append([float(tem[12:22]), float(tem[25:35]), float(tem[38:48])])
                #print positions
                self.positions.append(positions)
                #print self.positions
                tem = infile.readline().strip(" ").split()
                self.atomsymbols.append(tem[0])
                self.atomnumbers.append(float(tem[-1]))
                
                mat = numpy.zeros((3, 3), dtype=numpy.float_)
                tem = infile.readline().strip(" ").split()
                #print tem[-3:],mat[0,:]
                mat[0, :] = tem[-3:]
                tem = infile.readline().strip(" ").split()
                mat[1, :] = tem[-3:]
                tem = infile.readline().strip(" ").split()
                mat[2, :] = tem[-3:]
                self.locrotmatrix.append(mat)
                
            for ia in range(len(self.atomsymbols)):
                print "atom symbol :  %s   atom number: %d   atom multi:  %d" % (self.atomsymbols[ia], self.atomnumbers[ia], self.multi[ia])
                print "positions:"
                for im in self.positions[ia]:
                    print im[:]
                                     
                
            # symmetry with lattice vector
            tem = infile.readline().strip(" ").split()
            self.symm = []
            self.tau = []
            self.Nsymm = int(tem[0])
            for isymm in range(self.Nsymm):
                symmt = numpy.zeros((3, 3), dtype=float)
                taut = numpy.zeros((3), dtype=float) 
                for ir in range(3):
                    s = infile.readline()
                    for ic in range(3):
                        symmt[ir][ic] = float(s[ic * 2:ic * 2 + 2]) 
                    taut[ir] = float(s[7:17])
                self.symm.append(symmt)
                self.tau.append(taut)
                infile.readline()
                
            #############
            print "Read in %s.struct done!" % wiencase
            
            
            
        ## convential Cell Volume and primitive Cell. In bohr^3 unit
        self.VolumeCC, self.VolumePC = cellvolume(self.latticetype, self.latticeconstants, self.latticeangle)



    def readSGsymm(self, wiencase):
        """ read in symmetry of space group from wiencase.outputs other than wiencase.struct since 
        in this file there are also symmetry operations in xyz coordinates..
        """
        structfile = wiencase + ".outputs"
        self.nsymm = 1
        self.symm = []
        self.tau = []
        # note bravaismatrix is not accurate enough in wiencase.outputs file. Just use it for test.
        self.bravaismatrix = numpy.zeros((3, 3), dtype=numpy.float_)
        self.symmcartesian = []
        self.taucartesian = []
        with open(structfile, "r") as f:
            f.readline()
            f.readline()# bravais matrix
            for i in range(3):            
                line = f.readline().strip().split()
                self.bravaismatrix[i, :] = numpy.array([numpy.float(item) for item in line])[:]
            print "bravais matrix", self.bravaismatrix            
            
            while 1:
                try:
                    s = f.readline().strip(" ").split()
                    if(s[0] == "PGBSYM:"):
                        self.nsymm = int(s[-1])
                        break
                except:
                    assert "Error in read case.outputs"
            
            #f.readline()
            #f.readline()
            for i in range(self.nsymm):
                while 1:
                    s = f.readline().strip().split()
                    if s[0] == "Symmetry":
                        break
                  
                ## read symmcartesian
                symmt = numpy.zeros((3, 3), dtype=float)
                taut = numpy.zeros((3), dtype=float)
                for ir in range(3):
                    s = f.readline().strip().split()
                    for ic in range(3):
                        symmt[ir, ic] = float(s[ic]) 
                s = f.readline().strip().split()
                for ir in range(3):
                    taut[ir] = float(s[ir])
                    
                self.symmcartesian.append(symmt)
                self.taucartesian.append(taut)
                
                ##read symm
                symmt = numpy.zeros((3, 3), dtype=float)
                taut = numpy.zeros((3), dtype=float)
                for ir in range(3):
                    s = f.readline().strip().split()
                    for ic in range(3):
                        symmt[ir, ic] = numpy.float(s[ic]) 
                    taut[ir] = numpy.float(s[3])
                    
                self.symm.append(symmt)
                self.tau.append(taut)
                
            # end
            f.close()

    def checksymmxyz(self):
        '''  This is to check symm in cartesian coordinates and in primitive cell lattice
            For details of symm, should check wien/SRC_symmetry/, latsym.f, pglsym.f,pgbsym.f
            In general, if a lattice is orthorhombic, then symmcartesian is the same as symmprimitive
            (at most different by a transpose (or inversion)). If a lattice is not orthorhombic, these two symms could 
            be related by bravais matrix. 
            One special case is CXZ lattice. In Wien CXZ lattice contains two cases: orthorhombic, and 
            monoclinic (with only gamma not equal to 90). For CXZ lattice monoclinic, symmcartesian is also the same as
            symmprimitive (or with transpose (or inversion)). 
        '''
        ortho = numpy.abs(numpy.array(self.latticeangle) - 90.0).sum() <= 1e-6
        for i in range(self.nsymm):
            mat = self.symmcartesian[i].transpose()  # according to wien2k, why?
            res = mat
            if (not ortho) and (self.latticetype != "CXZ") :
                bm = self.bravaismatrix
                bminv = numpy.linalg.inv(bm)#.transpose()
                res = numpy.dot(bminv, numpy.dot(mat, bm))
            res -= self.symm[i]
            print i, numpy.abs(res).sum()

        ## primitive cell vectors.+

