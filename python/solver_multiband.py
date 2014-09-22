from itertools import product
from hamiltonians import *

class LocalProblem():

    def __init__(self,spin_names,orb_names,orb_hyb,h_loc_type,**h_loc_params):

        self.spin_names = spin_names
        self.orb_names = orb_names
        self.orb_hyb = orb_hyb
        self.h_loc_type = h_loc_type
        self.h_loc_params = h_loc_params

        self.gf_struct = self.get_gf_struct()
        self.h_loc = self.get_h_loc()

    # Set block structure of GF
    def get_gf_struct(self):
        gf_struct = {} 
        if self.orb_hyb: # outer blocks are spin blocks
            for sn in self.spin_names:
                gf_struct[sn] = [int(i) for i in self.orb_names]
        else: # outer blocks are spin-orbital blocks
            for sn, an in product(self.spin_names,self.orb_names):
                gf_struct[sn+'_'+an] = [0]

    	return gf_struct

    # Pick desired Hamiltonian
    def get_h_loc(self):

        if self.h_loc_type == "slater":
            return h_loc_slater(self.spin_names,self.orb_names,self.orb_hyb,**self.h_loc_params)     # h_loc_params must include U_matrix, and optionally H_dump
        elif self.h_loc_type == "kanamori":
            return h_loc_kanamori(self.spin_names,self.orb_names,self.orb_hyb,**self.h_loc_params)   # h_loc_params must include U, Uprime, J_hund, and optionally H_dump
        elif self.h_loc_type == "density":
            return h_loc_density(self.spin_names,self.orb_names,self.orb_hyb,**self.h_loc_params)    # h_loc_params must include U, Uprime, and optionally H_dump
        elif self.h_loc_type == "other":
            return self.h_loc_params["h_loc"]                                                # user provides h_loc with argument h_loc
        else:
            raise RuntimeError("Hamiltonian type not implemented.")
