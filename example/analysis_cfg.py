import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.stats import rankdata
from scipy.linalg.lapack import dsyevd
from scipy.sparse.linalg import eigs
from optparse import OptionParser


import soft_spot_project.sim_lib as sim
from soft_spot_project.init_model import InitializeModel
from soft_spot_project.visualization_tools import render_snapshot,render_field
from soft_spot_project.modes_tools import get_sparse_hessian,participation_ratio
import soft_spot_project.lib_bond_order as lbo


from softspot.softspot import SoftSpot



def read_andrea_cfg(file_name):
	# read snapshot
	with open(file_name) as f:
	    first_line = f.readline().split()
	    second_line = f.readline().split()

	L = np.float64(second_line)
	data = np.loadtxt(file_name,skiprows=2)

	diameter = data[:,0]
	pos = data[:,1:]
	pos += L/2.

	return L,diameter,pos



def smooth_field_over_neighbors(field):
	# takes the mean over interacting neighbors
	smooth_field = np.zeros_like(field)
	for i in range(sim.system.npart):
		j = sim.pairwise.pair_list_ij[:sim.pairwise.pair_list_ni[i],i]
		smooth_field[i] = np.mean(field[j])
	return smooth_field





desc=""" stuff """
parser = OptionParser(description=desc, usage='Usage: %prog (<infile>) <options>',
    version='%prog version 1.0', epilog='Author: David Richard')



parser.add_option("-f","--namefile", dest="namefile", action="store")

(options, args) = parser.parse_args()





# read configuration file
file_name = options.namefile
L,diameter,pos = read_andrea_cfg(file_name)
box_size = L*np.ones(pos.shape[1])



# inititalize system (takes positions, diameters, box geometry and local strain (zero for you))
sim = InitializeModel(model_name='swap_ipl',r=pos,diameter=diameter,box_size=box_size,strain=0.0).sim
sim.system.print_system_info()



# compute "everything" (information needed for the forces and hessian), has to be called before any mode analysis
sim.pairwise.compute_poly_everything()

print('U/N=',sim.system.thermo_pot/sim.system.npart)
print('Pressure (excess)=',sim.system.thermo_pre)
print('Stress (xy)=',sim.system.thermo_sigma)
print('Typical force grad =',sim.system.typical_grad)




# compute hexatic bond order parameter
bond_order_hexatic = lbo.bond_order.compute_2d_bop(sim.system.r,sim.pairwise.pair_list_ni,sim.pairwise.pair_list_ij+1,sim.system.boxl,sim.system.box_delrx)
# compute Tanaka bond metric
bond_order_theta = lbo.bond_order.compute_2d_theta(sim.system.r,sim.system.d,0.2,sim.pairwise.pair_list_ni,sim.pairwise.pair_list_ij+1,sim.system.boxl,sim.system.box_delrx)



# compute Hessian
hessian = get_sparse_hessian()



# partial diagonalization
n_modes = 50
n_modes += sim.system.ndim # add d to take into account the translational modes
vals, vecs = eigs(hessian, k=n_modes, which='LM', sigma=0, tol=1e-6)
# Throw away infinitesimal imaginary part and sort.
vals  = np.real(vals)
vecs = np.real(vecs)
idx = vals.argsort()
vals   = vals[idx]
vecs  = vecs[:,idx]



# full diagonalization (get all harmonic eigenmodes)
# dense_h = hessian.todense()
# st = time.time()
# print('start dsyev(d)')
# vals, vecs, out = dsyevd(dense_h)
# print('elapsed time (s):',time.time()-st)




# init the SoftSpot object with the current hessian
softspot = SoftSpot(ndim=sim.system.ndim, npart=sim.system.npart, hessian=hessian)



# plt.rcParams['text.usetex'] = True
# plt.rcParams['text.latex.preamble'] = r"\usepackage{cmbright} \usepackage{amsmath}"
plt.rcParams['axes.linewidth'] = 0.5 # set the value globally

fig = plt.gcf()
fig.set_facecolor('white') 
fig.set_size_inches(10,7)


ax1 = plt.subplot(2,3,1,aspect='equal')
plt.title('snapshot',fontsize=9)
ax2 = plt.subplot(2,3,2,aspect='equal')
plt.title('smooth Tanaka bond order',fontsize=9)
ax3 = plt.subplot(2,3,3,aspect='equal')
plt.title('rank (smooth Tanaka bond order)',fontsize=9)
ax4 = plt.subplot(2,3,4)
plt.title('participation vs. frequency',fontsize=9)
plt.ylabel(r'$e_p$',fontsize=9)
plt.xlabel(r'$\omega=\sqrt{\kappa}$',fontsize=9)
plt.yscale('log')
plt.xscale('log')
ax5 = plt.subplot(2,3,5,aspect='equal')
plt.title('example mapping',fontsize=9)
ax6 = plt.subplot(2,3,6,aspect='equal')
plt.title('modes detection',fontsize=9)

render_snapshot(ax1,pos,diameter,box_size,color='dimgray',edgecolor='black',alpha=0.75)
render_snapshot(ax2,pos,diameter,box_size,color=smooth_field_over_neighbors(bond_order_theta),edgecolor='black',alpha=1,cmap='cividis',cbar=True)
render_snapshot(ax3,pos,diameter,box_size,color=rankdata(smooth_field_over_neighbors(bond_order_theta))/sim.system.npart,edgecolor='black',alpha=1,cmap='cividis',cbar=True)

render_snapshot(ax5,pos,diameter,box_size,color='dimgray',edgecolor='black',alpha=0.75)
render_snapshot(ax6,pos,diameter,box_size,color='dimgray',edgecolor='black',alpha=0.75)


# example PHM modes
psi = vecs[:,2] # defined an harmonic mode as a starting solution
render_field(ax5,pos,psi,color='skyblue',scale=0.02)
# map psi onto pi
result_cg = softspot.find(psi, mode='cg', options={'tol': 1e-6})
pi = result_cg['pi']
render_field(ax5,pos,pi,color='crimson',scale=0.02)



pi_modes = []
pi_kappa = []
ep_psi = []
ep_pi  = []

count_pi = 0

for k in range(2,len(vals)):
	result_cg = softspot.find(vecs[:,k], mode='cg', options={'tol': 1e-6})
	pi = result_cg['pi']
	kappa_pi = result_cg['kappa']
	ep_psi.append(participation_ratio(vecs[:,k]))

	already_found = False
	for m,mode in enumerate(pi_modes):
		dot = np.abs(np.dot(pi,mode))
		if(dot>0.999):
			already_found=True
			break

	if (already_found==False):
		pi_modes.append(pi)
		pi_kappa.append(kappa_pi)
		ep_pi.append(participation_ratio(pi))
		count_pi+=1
		print('#',k,'pi count ',count_pi,'kappa(pi)',kappa_pi)


pi_modes  = np.array(pi_modes)
pi_kappa  = np.array(pi_kappa)
psi_kappa = vals[2:]
ep_psi = np.array(ep_psi)
ep_pi  = np.array(ep_pi)

# sort according to stiffness
idx = pi_kappa.argsort()
pi_modes = pi_modes[idx]
pi_kappa = pi_kappa[idx]
ep_pi = ep_pi[idx]


for k,mode in enumerate(pi_modes[::-1]): # display the softest at the end
	render_field(ax6,pos,mode,color=cm.magma(1.*k/count_pi),scale=0.02, remove_fraction=0.95)


# participation versus mode frequency
ax4.plot(np.sqrt(psi_kappa),ep_psi,'k.',label='harmonic')
ax4.plot(np.sqrt(pi_kappa),ep_pi,'ro',markerfacecolor='none',label='pseudo harmonic')
ax4.legend()



plt.tight_layout(pad=0.2)
plt.savefig('mode_analysis.jpg', dpi=500)





