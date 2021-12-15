
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
import os.path
from locale import atof


class ElkConverterTools:
    """
    Conversion Tools required to help covert Elk outputs into the TRIQS format.
    """

    def __init__(self):
        pass

    def rotaxang(self,rot):
      """
      This routine determines the axis of rotation vector (v) and the angle of rotation (th).
      If R corresponds to an improper rotation then only the proper part is used and the determinant
      is set to -1. The rotation convention follows the "right-hand rule". See Elk's rotaxang
      routine.
      """
      eps=1E-8
      v=numpy.zeros([3], numpy.float_)
    # find the determinant
      det=numpy.linalg.det(rot)
      if (abs(det-1.0)<eps):
        det=1.0
      elif (abs(det+1.0)<eps):
        det=-1.0
      else:
        raise "sym_converter : Invalid rotation matrix!"
    # proper rotation matrix
      rotp=det*rot
      v[0]=(rotp[1,2]-rotp[2,1])/2.0
      v[1]=(rotp[2,0]-rotp[0,2])/2.0
      v[2]=(rotp[0,1]-rotp[1,0])/2.0
      t1=numpy.sqrt(numpy.dot(v,v))
      t2=(rotp[0,0]+rotp[1,1]+rotp[2,2]-1.0)/2.0
      if (abs(abs(t2)-1.0)>eps):
      # theta not equal to 0 or pi
        th=-numpy.arctan2(t1,t2)
        v[:]=v[:]/t1
      else:
      # special case of sin(th)=0
        if(t2>eps):
      # zero angle: axis arbitrary
          th=0.0
          v[:]=1.0/numpy.sqrt(3.0)
        else:
      # rotation by pi
          th=numpy.pi
          if((rotp[0,0]>=rotp[1,1])&(rotp[0,0]>=rotp[2,2])):
            if(rotp[0,0]<(-1.0+eps)):
              mpi.report(rotp[0,0],-1.0+eps)
              raise "sym_converter : Invalid rotation matrix!"
            v[0]=numpy.sqrt(abs(rotp[0,0]+1.0)/2.0)
            v[1]=(rotp[1,0]+rotp[0,1])/(4.0*v[0])
            v[2]=(rotp[2,0]+rotp[0,2])/(4.0*v[0])
          elif((rotp[1,1]>=rotp[0,0])&(rotp[1,1]>=rotp[2,2])):
            if(rotp[1,1]<(-1.0+eps)):
              mpi.report(rotp[1,1],-1.0+eps)
              raise "sym_converter : Invalid rotation matrix!"
            v[1]=numpy.sqrt(abs(rotp[1,1]+1.0)/2.0)
            v[2]=(rotp[2,1]+rotp[1,2])/(4.0*v[1])
            v[0]=(rotp[0,1]+rotp[1,0])/(4.0*v[1])
          else:
            if(rotp[2,2]<(-1.0+eps)):
              mpi.report(rotp[2,2],-1.0+eps)
              raise "sym_converter : Invalid rotation matrix!"
            v[2]=numpy.sqrt(abs(rotp[2,2]+1.0)/2.0)
            v[0]=(rotp[0,2]+rotp[2,0])/(4.0*v[2])
            v[1]=(rotp[1,2]+rotp[2,1])/(4.0*v[2])
      # return -theta and v. -theta is returned as TRIQS does not rotate
      # the observable (such as the density matrix) which is done in Elk
      return v,-th

    def axangsu2(self,v,th):
      """
      Calculate the rotation SU(2) matrix - see Elk's axangsu2 routine.
      """
      su2=numpy.zeros([2,2], numpy.complex_)
      t1=numpy.sqrt(numpy.dot(v,v))
      if(t1<1E-8):
          raise "sym_converter : zero length axis vector!"
      # normalise the vector
      t1=1.0/t1
      x=v[0]*t1; y=v[1]*t1; z=v[2]*t1
      #calculate the SU(2) matrix
      cs=numpy.cos(0.5*th)
      sn=numpy.sin(0.5*th)
      su2[0,0]=cs-z*sn*1j
      su2[0,1]=-y*sn-x*sn*1j
      su2[1,0]=y*sn-x*sn*1j
      su2[1,1]=cs+z*sn*1j
      #return the SU(2) matrix
      return su2

    def v3frac(self,v,eps):
       """
       This finds the fractional part of 3-vector v components. This uses the
       same method as in Elk (version 6.2.8) r3fac subroutine.
       """
       v[0]=v[0]-numpy.floor(v[0])
       if(v[0] < 0): v[0]+=1
       if((1-v[0]) < eps): v[0]=0
       if(v[0] < eps): v[0]=0
       v[1]=v[1]-numpy.floor(v[1])
       if(v[1] < 0): v[1]+=1
       if((1-v[1]) < eps): v[1]=0
       if(v[1] < eps): v[1]=0
       v[2]=v[2]-numpy.floor(v[2])
       if(v[2] < 0): v[2]+=1
       if((1-v[2]) < eps): v[2]=0
       if(v[2] < eps): v[2]=0
       return v

    def gen_perm(self,nsym,ns,na,natmtot,symmat,tr,atpos,epslat=1E-6):
        """
        Generate the atom permutations per symmetry.
        """
        perm=[]
        iea=[]
        for isym in range(nsym):
          iea.append(numpy.zeros([natmtot,ns], numpy.int_))
       #loop over species
          for js in range(ns):
       #loop over species atoms
            v=numpy.zeros([3,na[js]], numpy.float_)
            v2=numpy.zeros(3, numpy.float_)
            for ia in range(na[js]):
              v[:,ia]=self.v3frac(atpos[js][ia][0:3],epslat)
            for ia in range(na[js]):
              v2[:]=numpy.matmul(symmat[isym][:,:],(atpos[js][ia][0:3]+tr[isym]))
              v2[:]=self.v3frac(v2,epslat)
              for ja in range(na[js]):
                t1=sum(abs(v[:,ja]-v2[:])) #check
                if(t1 < epslat):
                  iea[isym][ja,js]=ia
                  break
        #put iea into perm format
        for isym in range(nsym):
          perm.append([])
          ja=0
          prv_atms=0
          for js in range(ns):
            for ia in range(na[js]):
              perm[isym].append(iea[isym][ia,js]+prv_atms+1)
              ja+=1
            prv_atms+=na[js]
        #output perm
        return perm

    def symlat_to_complex_harmonics(self,nsym,n_shells,symlat,shells):
        """
        This calculates the Elk (crystal) symmetries in complex spherical harmonics
        This follows the methodology used in Elk's rotzflm, ylmrot and ylmroty routines.
        """
        #need SciPy routines to get Euler angles - need version 1.4+
        #from scipy.spatial.transform import Rotation as R
        symmat=[]
        rot=numpy.identity(3, numpy.float_)
        angi=numpy.zeros(3, numpy.float_)
        #loop over symmetries
        for isym in range(nsym):
          symmat.append([])
          for ish in range(n_shells):
            l=shells[ish]['l']
            symmat[isym].append(numpy.zeros([2*l+1, 2*l+1], numpy.complex_))
            #get determinant
            det=numpy.linalg.det(symlat[isym])
            p=1
            #p is -1 for improper symmetries
            if(det<0.0): p=-1
            rot[:,:]=p*symlat[isym][:,:]
            #r=R.from_matrix(rot)
            #get the y-convention Euler angles as used by Elk.
            #ang=r.as_euler('zyz')
            ang=self.zyz_euler(rot)
            #Elk uses inverse rotations, i.e. the function is being rotated, not the spherical harmonics
            #TRIQS rotates the spherical harmonics instead
            angi[0]=ang[0]
            angi[1]=ang[1]
            angi[2]=ang[2]
            #calculate the symmetry in the complex spherical harmonic basis.
            d = self.ylmrot(p,angi,l)
            symmat[isym][ish][:,:] = d[:,:]
        #return the complex spherical harmonic
        return symmat

    def zyz_euler(self,rot):
        """
        This calculates the Euler angles of matrix rot in the y-convention.
        See Elk's roteuler routine.
        This will be made redundent when TRIQS uses scipy version 1.4+
        """
        eps=1E-8
        pi=numpy.pi
        ang=numpy.zeros(3, numpy.float_)
        #get the Euler angles
        if((abs(rot[2,0])>eps) or (abs(rot[2,1])>eps)):
          ang[0]=numpy.arctan2(rot[2,1],rot[2,0])
          if(abs(rot[2,0])>abs(rot[2,1])):
            ang[1]=numpy.arctan2(rot[2,0]/numpy.cos(ang[0]),rot[2,2])
          else:
            ang[1]=numpy.arctan2(rot[2,1]/numpy.sin(ang[0]),rot[2,2])
          ang[2]=numpy.arctan2(rot[1,2],-rot[0,2])
        else:
          ang[0]=numpy.arctan2(rot[0,1],rot[0,0])
          if(rot[2,2]>0.0):
            ang[1]=0.0
            ang[2]=0.0
          else:
            ang[1]=pi
            ang[2]=pi
        #return Euler angles
        return ang


    def ylmrot(self,p,angi,l):
        """
        calculates the rotation matrix in complex spherical harmonics for l.
        THIS HAS ONLY BEEN TESTED FOR l=2.
        """
        d=numpy.identity(2*l+1, numpy.complex_)
        # generate the rotation matrix about the y-axis
        dy=self.ylmroty(angi[1],l)
        # apply inversion to odd l values if required
        if(p==-1):
          if(l % 2.0 != 0):
            dy*=-1
        # rotation by alpha and gamma
        for m1 in range(-l,l+1,1):
          lm1=l+m1
          for m2 in range(-l,l+1,1):
            lm2=l+m2
            t1=-m1*angi[0]-m2*angi[2]
            d[lm1,lm2]=dy[lm1,lm2]*(numpy.cos(t1)+1j*numpy.sin(t1))
        #return the rotation matrix
        return d

    def ylmroty(self,beta,l):
        """
        returns the rotation matrix around the y-axis with angle beta.
        This uses the same real matrix formual as in Elk - see Elk's manual for ylmroty description
        """
        #import the factorial function - needed for later versions of scipy (needs testing)
        from scipy import special as spec
        #calculates the rotation matrix in complex spherical harmonics for l
        dy=numpy.identity(2*l+1, numpy.float_)
        #sine and cosine of beta
        cb=numpy.cos(beta/2.0)
        sb=numpy.sin(beta/2.0)
        # generate the rotaion operator for m-components of input l
        for m1 in range(-l,l+1,1):
          for m2 in range(-l,l+1,1):
            sm=0.0
            minlm=numpy.amin([l+m1, l-m2]) + 1
            for k in range(minlm):
              if(((l+m1-k)>=0) and ((l-m2-k)>=0) and ((m2-m1+k)>=0)):
                j=2*(l-k)+m1-m2
                if(j==0):
                  t1=1.0
                else:
                  t1=cb**j
                j=2*k+m2-m1
                if(j!=0):
                  t1=t1*sb**j
                t2=t1/(spec.factorial(k)*spec.factorial(l+m1-k)*spec.factorial(l-m2-k)*spec.factorial(m2-m1+k))
                if(k % 2.0 != 0):
                  t2=-t2
                sm+=t2
            t1=numpy.sqrt(spec.factorial(l+m1)*spec.factorial(l-m1)*spec.factorial(l+m2)*spec.factorial(l-m2))
            dy[m1+l,m2+l]=t1*sm
        #return y-rotation matrix
        return dy

