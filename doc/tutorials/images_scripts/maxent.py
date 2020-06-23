from triqs.gf import *
from h5 import *
from triqs_maxent import *

filename = 'nio'

ar = HDFArchive(filename+'.h5','a')
if 'iteration_count' in ar['DMFT_results']:
    iteration_offset = ar['DMFT_results']['iteration_count']+1
    G_latt = ar['DMFT_results']['Iterations']['G_latt_orb_it'+str(iteration_offset-1)]


tm = TauMaxEnt(cost_function='bryan', probability='normal')

print((G_latt['up'][0,0]))
t2g_orbs = [0,1,3]
eg_orbs = [2,4]
op_orbs = [5,6,7]

orbs = [t2g_orbs, eg_orbs, op_orbs]
#orbs = [t2g_orbs]

for orb in orbs:

    print('\n'+str(orb[0])+'\n')

    gf = 0*G_latt['up'][0,0]
    for iO in orb:
        gf = gf + G_latt['up'][iO,iO]
    tm.set_G_iw(gf)
    tm.omega =LinearOmegaMesh(omega_min=-20, omega_max=20, n_points=201)
    tm.alpha_mesh = LogAlphaMesh(alpha_min=0.01, alpha_max=20000, n_points=60)

    tm.set_error(1.e-3)
    result=tm.run()
    result.get_A_out('LineFitAnalyzer')
    
    if 'iteration_count' in ar['DMFT_results']:
        iteration_offset = ar['DMFT_results']['iteration_count']+1
        for oo in orb:
            ar['DMFT_results']['Iterations']['G_latt_orb_w_o'+str(oo)+'_it'+str(iteration_offset-1)] = result.analyzer_results['LineFitAnalyzer']['A_out']
        ar['DMFT_results']['Iterations']['w_it'+str(iteration_offset-1)] = result.omega


    # you may be interested in the details of the line analyzer:
    # from triqs.plot.mpl_interface import oplot
    #plt.figure(2)
    #result.analyzer_results['LineFitAnalyzer'].plot_linefit()
    #plt.savefig('ana'+str(orb[0])+'.pdf',fmt='pdf')

del ar
