#from triqs_dft_tools.sumk_dft import *
from sumk_dft import *
#from triqs_dft_tools.converters.wien2k_converter import *
from converters.vasp_converter import *
#from pytriqs.applications.impurity_solvers.hubbard_I.hubbard_solver import Solver
from hf_solver import Solver
import shutil

class TestSumkDFT(SumkDFT):
# calculate and save occupancy matrix in the Bloch basis for VASP charge denity recalculation
    def calc_density_correction(self, filename='GAMMA', dm_type='wien2k'):
        r"""
        Calculates the charge density correction and stores it into a file.
        
        The charge density correction is needed for charge-self-consistent DFT+DMFT calculations.
        It represents a density matrix of the interacting system defined in Bloch basis
        and it is calculated from the sum over Matsubara frequecies of the full GF,

        ..math:: N_{\nu\nu'}(k) = \sum_{i\omega_{n}} G_{\nu\nu'}(k, i\omega_{n})
        
        The density matrix for every `k`-point is stored into a file.

        Parameters
        ----------
        filename : string
                   Name of the file to store the charge density correction.

        Returns
        -------
        (deltaN, dens) : tuple
                         Returns a tuple containing the density matrix `deltaN` and
                         the corresponing total charge `dens`.

        """
        assert type(filename) == StringType, "calc_density_correction: filename has to be a string!"
        assert dm_type in ('vasp', 'wien2k'), "'type' must be either 'vasp' or 'wienk'"

        ntoi = self.spin_names_to_ind[self.SO]
        spn = self.spin_block_names[self.SO]
        dens = {sp: 0.0 for sp in spn}

# Fetch Fermi weights and energy window band indices
        if dm_type == 'vasp':
            fermi_weights = 0
            band_window = 0
            if mpi.is_master_node():
                ar = HDFArchive(self.hdf_file,'r')
                fermi_weights = ar['dft_misc_input']['dft_fermi_weights']
                band_window = ar['dft_misc_input']['band_window']
                del ar
            fermi_weights = mpi.bcast(fermi_weights)
            band_window = mpi.bcast(band_window)

# Convert Fermi weights to a density matrix
            dens_mat_dft = {}
            for sp in spn:
                dens_mat_dft[sp] = [fermi_weights[ik, ntoi[sp], :].astype(numpy.complex_) for ik in xrange(self.n_k)]

        # Set up deltaN:
        deltaN = {}
        for sp in spn:
            deltaN[sp] = [numpy.zeros([self.n_orbitals[ik,ntoi[sp]],self.n_orbitals[ik,ntoi[sp]]], numpy.complex_) for ik in range(self.n_k)]

        ikarray = numpy.array(range(self.n_k))
        for ik in mpi.slice_array(ikarray):
            G_latt_iw = self.lattice_gf(ik = ik, mu = self.chemical_potential, iw_or_w = "iw")
            for bname,gf in G_latt_iw:
                deltaN[bname][ik][:, :] = G_latt_iw[bname].density()
                dens[bname] += self.bz_weights[ik] * G_latt_iw[bname].total_density()
                if dm_type == 'vasp':
# In 'vasp'-mode subtract the DFT density matrix
                    diag_inds = numpy.diag_indices(self.n_orbitals[ik, ntoi[bname]])
                    deltaN[bname][ik][diag_inds] -= dens_mat_dft[bname][ik]
                    dens[bname] -= self.bz_weights[ik] * dens_mat_dft[bname][ik].sum().real

        # mpi reduce:
        for bname in deltaN:
            for ik in range(self.n_k):
                deltaN[bname][ik] = mpi.all_reduce(mpi.world, deltaN[bname][ik], lambda x,y : x+y)
            dens[bname] = mpi.all_reduce(mpi.world, dens[bname], lambda x,y : x+y)
        mpi.barrier()

        # now save to file:
        if dm_type == 'wien2k':
            if mpi.is_master_node():
                if self.SP == 0:
                    f = open(filename,'w')
                else:
                    f = open(filename+'up','w')
                    f1 = open(filename+'dn','w')
                # write chemical potential (in Rydberg):
                f.write("%.14f\n"%(self.chemical_potential/self.energy_unit))
                if self.SP != 0: f1.write("%.14f\n"%(self.chemical_potential/self.energy_unit))
                # write beta in rydberg-1
                f.write("%.14f\n"%(G_latt_iw.mesh.beta*self.energy_unit))
                if self.SP != 0: f1.write("%.14f\n"%(G_latt_iw.mesh.beta*self.energy_unit))

                if self.SP == 0: # no spin-polarization

                    for ik in range(self.n_k):
                        f.write("%s\n"%self.n_orbitals[ik,0])
                        for inu in range(self.n_orbitals[ik,0]):
                            for imu in range(self.n_orbitals[ik,0]):
                                valre = (deltaN['up'][ik][inu,imu].real + deltaN['down'][ik][inu,imu].real) / 2.0
                                valim = (deltaN['up'][ik][inu,imu].imag + deltaN['down'][ik][inu,imu].imag) / 2.0
                                f.write("%.14f  %.14f "%(valre,valim))
                            f.write("\n")
                        f.write("\n")
                    f.close()

                elif self.SP == 1: # with spin-polarization

                    # dict of filename: (spin index, block_name)
                    if self.SO == 0: to_write = {f: (0, 'up'), f1: (1, 'down')}
                    if self.SO == 1: to_write = {f: (0, 'ud'), f1: (0, 'ud')}
                    for fout in to_write.iterkeys():
                        isp, sp = to_write[fout]
                        for ik in range(self.n_k):
                            fout.write("%s\n"%self.n_orbitals[ik,isp])
                            for inu in range(self.n_orbitals[ik,isp]):
                                for imu in range(self.n_orbitals[ik,isp]):
                                    fout.write("%.14f  %.14f "%(deltaN[sp][ik][inu,imu].real,deltaN[sp][ik][inu,imu].imag))
                                fout.write("\n")
                            fout.write("\n")
                        fout.close()
        elif dm_type == 'vasp':
#            assert self.SP == 0, "Spin-polarized density matrix is not implemented"

            if mpi.is_master_node():
                with open(filename, 'w') as f:
                    f.write(" %i  -1  ! Number of k-points, default number of bands\n"%(self.n_k))
                    for sp in spn:
                        for ik in xrange(self.n_k):
                            ib1 = band_window[0][ik, 0]
                            ib2 = band_window[0][ik, 1]
                            f.write(" %i  %i  %i\n"%(ik + 1, ib1, ib2))
                            for inu in xrange(self.n_orbitals[ik, 0]):
                                for imu in xrange(self.n_orbitals[ik, 0]):
                                    if self.SP == 0:
                                        valre = (deltaN['up'][ik][inu, imu].real + deltaN['down'][ik][inu, imu].real) / 2.0
                                        valim = (deltaN['up'][ik][inu, imu].imag + deltaN['down'][ik][inu, imu].imag) / 2.0
                                    else:
                                        valre = deltaN[sp][ik][inu, imu].real
                                        valim = deltaN[sp][ik][inu, imu].imag

                                    f.write(" %.14f  %.14f"%(valre, valim))
                                f.write("\n")
        else:
            raise NotImplementedError("Unknown density matrix type: '%s'"%(dm_type))

        return deltaN, dens

    def calc_hamiltonian_correction(self, filename='GAMMA'):
        r"""
        Calculates the charge density correction and stores it into a file.
        
        The charge density correction is needed for charge-self-consistent DFT+DMFT calculations.
        It represents a density matrix of the interacting system defined in Bloch basis
        and it is calculated from the sum over Matsubara frequecies of the full GF,

        ..math:: N_{\nu\nu'}(k) = \sum_{i\omega_{n}} G_{\nu\nu'}(k, i\omega_{n})
        
        The density matrix for every `k`-point is stored into a file.

        Parameters
        ----------
        filename : string
                   Name of the file to store the charge density correction.

        Returns
        -------
        (deltaN, dens) : tuple
                         Returns a tuple containing the density matrix `deltaN` and
                         the corresponing total charge `dens`.

        """
        assert type(filename) == StringType, "calc_density_correction: filename has to be a string!"

        ntoi = self.spin_names_to_ind[self.SO]
        spn = self.spin_block_names[self.SO]
        dens = {sp: 0.0 for sp in spn}

# Fetch Fermi weights and energy window band indices
        fermi_weights = 0
        band_window = 0
        if mpi.is_master_node():
            ar = HDFArchive(self.hdf_file,'r')
            fermi_weights = ar['dft_misc_input']['dft_fermi_weights']
            band_window = ar['dft_misc_input']['band_window']
            del ar
        fermi_weights = mpi.bcast(fermi_weights)
        band_window = mpi.bcast(band_window)

        # Set up deltaH:
        deltaH = {}
        for sp in spn:
            deltaH[sp] = [numpy.zeros([self.n_orbitals[ik,ntoi[sp]],self.n_orbitals[ik,ntoi[sp]]], numpy.complex_) for ik in range(self.n_k)]

        ikarray = numpy.array(range(self.n_k))
        for ik in mpi.slice_array(ikarray):
            sigma_minus_dc = [s.copy() for s in self.Sigma_imp_iw]
            sigma_minus_dc = self.add_dc('iw')

            beta = self.Sigma_imp_iw[0].mesh.beta   # override beta if Sigma_iw is present
            n_iw = len(self.Sigma_imp_iw[0].mesh)

            block_structure = [ range(self.n_orbitals[ik,ntoi[sp]]) for sp in spn ]
            gf_struct = [ (spn[isp], block_structure[isp]) for isp in range(self.n_spin_blocks[self.SO]) ]
            block_ind_list = [block for block,inner in gf_struct]
            glist = lambda : [ GfImFreq(indices=inner,beta=beta,n_points=n_iw) for block,inner in gf_struct]
            G_latt = BlockGf(name_list = block_ind_list, block_list = glist(), make_copies = False)
            G_latt.zero()

            for icrsh in range(self.n_corr_shells):
                for bname, gf in G_latt:
                    gf += self.upfold(ik,icrsh,bname,sigma_minus_dc[icrsh][bname],gf)

            for sp in spn:
                deltaH[sp][ik][:, :] = G_latt[sp](0) # Any Matsubara frequency will do
#            G_latt_iw = self.lattice_gf(ik = ik, mu = self.chemical_potential, iw_or_w = "iw")
#            for bname,gf in G_latt_iw:
#                deltaN[bname][ik][:, :] = G_latt_iw[bname].density()
#                dens[bname] += self.bz_weights[ik] * G_latt_iw[bname].total_density()
#                if dm_type == 'vasp':
## In 'vasp'-mode subtract the DFT density matrix
#                    diag_inds = numpy.diag_indices(self.n_orbitals[ik, ntoi[bname]])
#                    deltaN[bname][ik][diag_inds] -= dens_mat_dft[bname][ik]
#                    dens[bname] -= self.bz_weights[ik] * dens_mat_dft[bname][ik].sum().real

        # mpi reduce:
        for bname in deltaH:
            for ik in range(self.n_k):
                deltaH[bname][ik] = mpi.all_reduce(mpi.world, deltaH[bname][ik], lambda x,y : x+y)
        mpi.barrier()

        # now save to file:
        if mpi.is_master_node():
            with open(filename, 'w') as f:
                f.write("H  %i  -1  ! Number of k-points, default number of bands\n"%(self.n_k))
                for sp in spn:
                    for ik in xrange(self.n_k):
                        ib1 = band_window[0][ik, 0]
                        ib2 = band_window[0][ik, 1]
                        f.write(" %i  %i  %i\n"%(ik + 1, ib1, ib2))
                        for inu in xrange(self.n_orbitals[ik, 0]):
                            for imu in xrange(self.n_orbitals[ik, 0]):
                                if self.SP == 0:
                                    valre = (deltaH['up'][ik][inu, imu].real + deltaH['down'][ik][inu, imu].real) / 2.0
                                    valim = (deltaH['up'][ik][inu, imu].imag + deltaH['down'][ik][inu, imu].imag) / 2.0
                                else:
                                    valre = deltaH[sp][ik][inu, imu].real
                                    valim = deltaH[sp][ik][inu, imu].imag

                                f.write(" %.14f  %.14f"%(valre, valim))
                            f.write("\n")
        return deltaH

def dmft_cycle():
    lda_filename = 'vasp'
    beta = 400
    U_int = 4.00
    J_hund = 0.70
    Loops =  1                       # Number of DMFT sc-loops
    Mix = 1.0                        # Mixing factor in QMC
    DC_type = 0                      # 0...FLL, 1...Held, 2... AMF, 3...Lichtenstein
    DC_Mix = 1.0                     # 1.0 ... all from imp; 0.0 ... all from Gloc   
    useBlocs = False                 # use bloc structure from LDA input
    useMatrix = True                 # use the U matrix calculated from Slater coefficients instead of (U+2J, U, U-J)
    chemical_potential_init=-0.0      # initial chemical potential

    use_dudarev = True

    HDFfilename = lda_filename+'.h5'

# Convert DMFT input:
# Can be commented after the first run
    Converter = VaspConverter(filename=lda_filename)
    Converter.convert_dft_input()

#check if there are previous runs:
    previous_runs = 0
    previous_present = False

    if mpi.is_master_node():
        ar = HDFArchive(HDFfilename,'a')
        if 'iterations' in ar:
            previous_present = True
            previous_runs = ar['iterations']
        else:
            previous_runs = 0
            previous_present = False
        del ar

    mpi.barrier()
    previous_runs    = mpi.bcast(previous_runs)
    previous_present = mpi.bcast(previous_present)

# Init the SumK class
#SK=SumkDFT(hdf_file=lda_filename+'.h5',use_dft_blocks=False)
    SK=TestSumkDFT(hdf_file=lda_filename+'.h5',use_dft_blocks=False)

    Norb = SK.corr_shells[0]['dim']
    l    = SK.corr_shells[0]['l']

# Init the Hubbard-I solver:
    S = Solver(beta = beta, l = l, dudarev=use_dudarev)

# DEBUG
#SK.put_Sigma(Sigma_imp=[S.Sigma])
#dH = SK.calc_hamiltonian_correction(filename='GAMMA')
# END DEBUG

    chemical_potential=chemical_potential_init
# load previous data: old self-energy, chemical potential, DC correction
    if (previous_present):
        mpi.report("Using stored data for initialisation")
        if (mpi.is_master_node()):
            ar = HDFArchive(HDFfilename,'a')
            S.Sigma <<= ar['SigmaF']
            del ar
            things_to_load=['chemical_potential','dc_imp']
            old_data=SK.load(things_to_load)
            chemical_potential=old_data[0]
            SK.dc_imp=old_data[1]
        S.Sigma = mpi.bcast(S.Sigma)
        chemical_potential=mpi.bcast(chemical_potential)
        SK.dc_imp=mpi.bcast(SK.dc_imp)


# DMFT loop:
    for Iteration_Number in range(1,Loops+1):
        
            itn = Iteration_Number + previous_runs
           
            # put Sigma into the SumK class:
            S.Sigma.zero() # !!!!
            SK.put_Sigma(Sigma_imp = [ S.Sigma ])

            # Compute the SumK, possibly fixing mu by dichotomy
            if SK.density_required and (Iteration_Number > 1):
                chemical_potential = SK.calc_mu( precision = 0.001 )
            else:
                mpi.report("No adjustment of chemical potential\nTotal density  = %.3f"%SK.total_density(mu=chemical_potential))

            # Density:
            S.G <<= SK.extract_G_loc()[0]
            mpi.report("Total charge of Gloc : %.6f"%S.G.total_density())

            # calculated DC at the first run to have reasonable initial non-interacting atomic level positions
            if ((Iteration_Number==1)and(previous_present==False)):
                if use_dudarev:
                    dc_value_init = 0.0
                else:
                    dc_value_init=U_int/2.0
                dm=S.G.density()
                SK.calc_dc( dm, U_interact = U_int, J_hund = J_hund, orb = 0, use_dc_formula = DC_type, use_dc_value=dc_value_init)

            # calculate non-interacting atomic level positions:
#        eal = SK.eff_atomic_levels()[0]
#        S.set_atomic_levels( eal = eal )

            # solve it:
            corr_energy, dft_dc = S.solve(U_int = U_int, J_hund = J_hund, verbosity = 1)

            # Now mix Sigma and G:
            if ((itn>1)or(previous_present)):
                if (mpi.is_master_node()and (Mix<1.0)):
                    ar = HDFArchive(HDFfilename,'r')
                    mpi.report("Mixing Sigma and G with factor %s"%Mix)
                    if ('SigmaF' in ar):
                        S.Sigma <<= Mix * S.Sigma + (1.0-Mix) * ar['SigmaF']
                    if ('GF' in ar):
                        S.G <<= Mix * S.G + (1.0-Mix) * ar['GF']
                    del ar
                S.G = mpi.bcast(S.G)
                S.Sigma = mpi.bcast(S.Sigma)
            
            # after the Solver has finished, set new double counting: 
#            dm = S.G.density()
#            if not use_dudarev:
#                SK.calc_dc( dm, U_interact = U_int, J_hund = J_hund, orb = 0, use_dc_formula = DC_type )

            # correlation energy calculations:
            correnerg = 0.5 * (S.G * S.Sigma).total_density()
            mpi.report("Corr. energy = %s"%correnerg)

            # store the impurity self-energy, GF as well as correlation energy in h5
            if (mpi.is_master_node()):
                ar = HDFArchive(HDFfilename,'a')
                ar['iterations'] = itn
                ar['chemical_cotential%s'%itn] = chemical_potential
                ar['SigmaF'] = S.Sigma
                ar['GF'] = S.G
                ar['correnerg%s'%itn] = correnerg
                ar['DCenerg%s'%itn] = SK.dc_energ
                del ar

            #Save essential SumkDFT data:
            things_to_save=['chemical_potential','dc_energ','dc_imp']
            SK.save(things_to_save)

            # print out occupancy matrix of Ce 4f
#            mpi.report("Orbital densities of impurity Green function:")
#            for s in dm:
#                mpi.report("Block %s: "%s)
#                for ii in range(len(dm[s])):
#                    str = ''
#                    for jj in range(len(dm[s])):
#                        if (dm[s][ii,jj].real>0):
#                            str += "   %.4f"%(dm[s][ii,jj].real)
#                        else:
#                            str += "  %.4f"%(dm[s][ii,jj].real)
#                    mpi.report(str)
            mpi.report("Total charge of impurity problem : %.6f"%S.G.total_density())


# find exact chemical potential
#if (SK.density_required):
#    SK.chemical_potential = SK.calc_mu( precision = 0.000001 )
#    SK.chemical_potential = SK.calc_mu( precision = 0.01 )

#dN, d = SK.calc_density_correction(filename='GAMMA', dm_type='vasp')
#mpi.report("Trace of Density Matrix: %s"%d)
    mpi.report("Storing Hamiltonian correction GAMMA...")
    SK.put_Sigma(Sigma_imp=[S.Sigma])
    dH = SK.calc_hamiltonian_correction(filename='GAMMA')
#    shutil.copyfile('GAMMA', 'it%i.GAMMA'%(itn))

# store correlation energy contribution to be read by Wien2ki and then included to DFT+DMFT total energy
    if (mpi.is_master_node()):
        ar = HDFArchive(HDFfilename)
        itn = ar['iterations']
        correnerg = ar['correnerg%s'%itn]
        DCenerg = ar['DCenerg%s'%itn]
        del ar
        correnerg -= DCenerg[0]
        f=open(lda_filename+'.qdmft','a')
        f.write("%.16f\n"%correnerg)
        f.close()

    return corr_energy, dft_dc

if __name__ == '__main__':
    dmft_cycle()
