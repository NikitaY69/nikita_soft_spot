import numpy as np
import sys


from soft_spot_project import sim_lib as sim




class InitializeModel():

    ntypes = {
        'swap_ipl': 'poly',
    }


    def __init__(self, model_name, r, diameter, box_size, strain=0.0, options=None):
        """
        r: [0, L]
        """
        self.model_name = model_name
        self.options = options

        npart = r.shape[0] 
        ndim  = r.shape[1]

        sim.system.ndim  = ndim
        sim.system.npart = npart


        if self.ntypes[model_name] == 'poly':
            sim.system.ntypes =  npart
        else:
            sim.system.ntypes = self.ntypes[model_name]

        sim.system.init_system_arrays()
        sim.system.boxl[:]    = box_size
        sim.system.boxl2[:]   = box_size/2
        sim.system.box_volume = np.prod(box_size)
        sim.system.box_strain = strain
        sim.system.box_delrx  = strain*box_size[1]
        sim.system.box_rho    = 1.0*npart/np.prod(box_size)

        sim.system.r = np.array(r.T)
        sim.system.d = np.array(diameter)


        self.sim = sim

        eval('self.set_%s()' % model_name)




    def set_swap_ipl(self):
        rcut_buffer = 1.05
        rcut_one    = 1.25
        epsilon     = 0.2

        sim.cell.frc_cell_rc_min = rcut_buffer * rcut_one * sim.system.d.max()
        sim.cell.init_cell()

        sim.pairwise.init_pairwise_arrays()
        sim.potentials.init_poly_force_field('swap_ipl',[1.,rcut_one,epsilon])




