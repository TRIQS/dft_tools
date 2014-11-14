    def lattice_gf(self, ik, mu, iw_or_w="iw", beta=40, broadening, mesh=None, with_Sigma=True, with_dc=True): 
        """Calculates the lattice Green function from the LDA hopping and the self energy at k-point number ik
           and chemical potential mu."""

        ntoi = self.spin_names_to_ind[self.SO]
        spn = self.spin_block_names[self.SO]
        if (iw_or_w != "iw") and (iw_or_w != "w"): 
            raise ValueError, "lattice_gf: implemented only for Re/Im frequency functions."

        # Are we including Sigma?
        if with_Sigma:
            if with_dc: sigma_inc_dc = self.add_dc()
            Sigma_imp = getattr(self,"Sigma_imp_"+iw_or_w)
            if iw_or_w == "iw": beta = Sigma_imp[0].mesh.beta   # override beta if Sigma_iw is present
        else:
            if (iw_or_w == "w") and (mesh is None):
                raise ValueError, "Give the mesh=(om_min,om_max,n_points) for the lattice GfReFreq." 

        # Check if G_upfold is present
        set_up_G_upfold = False                       # Assume not
        if not hasattr(self,"G_upfold_"+iw_or_w)):
            set_up_G_upfold = True                    # Need to create G_upfold_iw
            setattr(self,"G_upfold_"+iw_or_w,0)       # Set temporarily to zero
        else:                                         # Check that existing GF is consistent
            GFsize = [ gf.N1 for bname,gf in getattr(self,"G_upfold_"+iw_or_w)]
            unchangedsize = all( [ self.n_orbitals[ik,ntoi[spn[isp]]]==GFsize[isp]
                                   for isp in range(self.n_spin_blocks[self.SO]) ] )
            if not unchangedsize: set_up_G_upfold = True
            if (iw_or_w != "iw") and (self.G_upfold_iw.mesh.beta != beta): # additional check for ImFreq
            set_up_G_upfold = True
        G_upfold = getattr(self,"G_upfold_"+iw_or_w) 

        # Set up G_upfold
        if set_up_G_upfold:
            block_structure = [ range(self.n_orbitals[ik,ntoi[sp]]) for sp in spn ]
            gf_struct = [ (spn[isp], block_structure[isp]) for isp in range(self.n_spin_blocks[self.SO]) ]
            block_ind_list = [block for block,inner in gf_struct]
            if with_Sigma:
                glist = lambda : [ GfImFreq(indices=inner,mesh=Sigma_imp[0].mesh) for block,inner in gf_struct]
            else:
                if iw_or_w == "iw":
                    glist = lambda : [ GfImFreq(indices=inner,beta=beta) for block,inner in gf_struct]
                elif iw_or_w == "w":
                    glist = lambda : [ GfImFreq(indices=inner,window=(mesh[0],mesh[1]),n_points=mesh[2]) for block,inner in gf_struct]
            G_upfold = BlockGf(name_list = block_ind_list, block_list = glist(), make_copies = False)
            G_upfold.zero()

        if iw_or_w == "iw":
            G_upfold << iOmega_n
        elif iw_or_w == "w":
            G_upfold << Omega + 1j*broadening

        idmat = [numpy.identity(self.n_orbitals[ik,ntoi[sp]],numpy.complex_) for sp in spn]
        M = copy.deepcopy(idmat)
        for ibl in range(self.n_spin_blocks[self.SO]):
            ind = ntoi[spn[ibl]]
            n_orb = self.n_orbitals[ik,ind]
            M[ibl] = self.hopping[ik,ind,0:n_orb,0:n_orb] - (idmat[ibl]*mu) - (idmat[ibl] * self.h_field * (1-2*ibl))
        G_upfold -= M

        if with_Sigma:
            for icrsh in range(self.n_corr_shells):
                for bname,gf in G_upfold: gf -= self.upfold(ik,icrsh,bname,sigma_inc_dc[icrsh][bname],gf)

        G_upfold.invert()

        return G_upfold
