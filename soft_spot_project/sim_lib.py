from __future__ import print_function, absolute_import, division
from . import _sim_lib
import f90wrap.runtime
import logging

class Kinds(f90wrap.runtime.FortranModule):
    """
    Module kinds
    
    
    Defined at src.f95 lines 1-25
    
    """
    @f90wrap.runtime.register_class("sim_lib.const_vec")
    class const_vec(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=const_vec)
        
        
        Defined at src.f95 lines 11-12
        
        """
        def __init__(self, handle=None):
            """
            self = Const_Vec()
            
            
            Defined at src.f95 lines 11-12
            
            
            Returns
            -------
            this : Const_Vec
            	Object to be constructed
            
            
            Automatically generated constructor for const_vec
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _sim_lib.f90wrap_const_vec_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Const_Vec
            
            
            Defined at src.f95 lines 11-12
            
            Parameters
            ----------
            this : Const_Vec
            	Object to be destructed
            
            
            Automatically generated destructor for const_vec
            """
            if self._alloc:
                _sim_lib.f90wrap_const_vec_finalise(this=self._handle)
        
        @property
        def c(self):
            """
            Element c ftype=real(dbl) pytype=float
            
            
            Defined at src.f95 line 12
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _sim_lib.f90wrap_const_vec__array__c(self._handle)
            if array_handle in self._arrays:
                c = self._arrays[array_handle]
            else:
                c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _sim_lib.f90wrap_const_vec__array__c)
                self._arrays[array_handle] = c
            return c
        
        @c.setter
        def c(self, c):
            self.c[...] = c
        
        def __str__(self):
            ret = ['<const_vec>{\n']
            ret.append('    c : ')
            ret.append(repr(self.c))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def init_const_vec(self, n):
        """
        init_const_vec(self, n)
        
        
        Defined at src.f95 lines 16-20
        
        Parameters
        ----------
        dertype : Const_Vec
        n : int
        
        """
        _sim_lib.f90wrap_init_const_vec(dertype=self._handle, n=n)
    
    @staticmethod
    def destroy_const_vec(self):
        """
        destroy_const_vec(self)
        
        
        Defined at src.f95 lines 22-24
        
        Parameters
        ----------
        dertype : Const_Vec
        
        """
        _sim_lib.f90wrap_destroy_const_vec(dertype=self._handle)
    
    @property
    def dbl(self):
        """
        Element dbl ftype=integer pytype=int
        
        
        Defined at src.f95 line 8
        
        """
        return _sim_lib.f90wrap_kinds__get__dbl()
    
    @property
    def sgl(self):
        """
        Element sgl ftype=integer pytype=int
        
        
        Defined at src.f95 line 9
        
        """
        return _sim_lib.f90wrap_kinds__get__sgl()
    
    @property
    def lint(self):
        """
        Element lint ftype=integer pytype=int
        
        
        Defined at src.f95 line 10
        
        """
        return _sim_lib.f90wrap_kinds__get__lint()
    
    def __str__(self):
        ret = ['<kinds>{\n']
        ret.append('    dbl : ')
        ret.append(repr(self.dbl))
        ret.append(',\n    sgl : ')
        ret.append(repr(self.sgl))
        ret.append(',\n    lint : ')
        ret.append(repr(self.lint))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

kinds = Kinds()

class System(f90wrap.runtime.FortranModule):
    """
    Module system
    
    
    Defined at src.f95 lines 27-104
    
    """
    @staticmethod
    def init_system_arrays():
        """
        init_system_arrays()
        
        
        Defined at src.f95 lines 64-86
        
        
        """
        _sim_lib.f90wrap_init_system_arrays()
    
    @staticmethod
    def print_system_info():
        """
        print_system_info()
        
        
        Defined at src.f95 lines 88-95
        
        
        """
        _sim_lib.f90wrap_print_system_info()
    
    @staticmethod
    def switch_diameter(i, j):
        """
        switch_diameter(i, j)
        
        
        Defined at src.f95 lines 97-103
        
        Parameters
        ----------
        i : int
        j : int
        
        """
        _sim_lib.f90wrap_switch_diameter(i=i, j=j)
    
    @property
    def nthreads(self):
        """
        Element nthreads ftype=integer pytype=int
        
        
        Defined at src.f95 line 37
        
        """
        return _sim_lib.f90wrap_system__get__nthreads()
    
    @nthreads.setter
    def nthreads(self, nthreads):
        _sim_lib.f90wrap_system__set__nthreads(nthreads)
    
    @property
    def ndim(self):
        """
        Element ndim ftype=integer pytype=int
        
        
        Defined at src.f95 line 39
        
        """
        return _sim_lib.f90wrap_system__get__ndim()
    
    @ndim.setter
    def ndim(self, ndim):
        _sim_lib.f90wrap_system__set__ndim(ndim)
    
    @property
    def npart(self):
        """
        Element npart ftype=integer pytype=int
        
        
        Defined at src.f95 line 40
        
        """
        return _sim_lib.f90wrap_system__get__npart()
    
    @npart.setter
    def npart(self, npart):
        _sim_lib.f90wrap_system__set__npart(npart)
    
    @property
    def ntypes(self):
        """
        Element ntypes ftype=integer pytype=int
        
        
        Defined at src.f95 line 41
        
        """
        return _sim_lib.f90wrap_system__get__ntypes()
    
    @ntypes.setter
    def ntypes(self, ntypes):
        _sim_lib.f90wrap_system__set__ntypes(ntypes)
    
    @property
    def inv_ndim(self):
        """
        Element inv_ndim ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 42
        
        """
        return _sim_lib.f90wrap_system__get__inv_ndim()
    
    @inv_ndim.setter
    def inv_ndim(self, inv_ndim):
        _sim_lib.f90wrap_system__set__inv_ndim(inv_ndim)
    
    @property
    def real_npart(self):
        """
        Element real_npart ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 43
        
        """
        return _sim_lib.f90wrap_system__get__real_npart()
    
    @real_npart.setter
    def real_npart(self, real_npart):
        _sim_lib.f90wrap_system__set__real_npart(real_npart)
    
    @property
    def box_volume(self):
        """
        Element box_volume ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 45
        
        """
        return _sim_lib.f90wrap_system__get__box_volume()
    
    @box_volume.setter
    def box_volume(self, box_volume):
        _sim_lib.f90wrap_system__set__box_volume(box_volume)
    
    @property
    def boxl(self):
        """
        Element boxl ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 46
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_system__array__boxl(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            boxl = self._arrays[array_handle]
        else:
            boxl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_system__array__boxl)
            self._arrays[array_handle] = boxl
        return boxl
    
    @boxl.setter
    def boxl(self, boxl):
        self.boxl[...] = boxl
    
    @property
    def boxl2(self):
        """
        Element boxl2 ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 46
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_system__array__boxl2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            boxl2 = self._arrays[array_handle]
        else:
            boxl2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_system__array__boxl2)
            self._arrays[array_handle] = boxl2
        return boxl2
    
    @boxl2.setter
    def boxl2(self, boxl2):
        self.boxl2[...] = boxl2
    
    @property
    def box_strain(self):
        """
        Element box_strain ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 47
        
        """
        return _sim_lib.f90wrap_system__get__box_strain()
    
    @box_strain.setter
    def box_strain(self, box_strain):
        _sim_lib.f90wrap_system__set__box_strain(box_strain)
    
    @property
    def box_delrx(self):
        """
        Element box_delrx ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 47
        
        """
        return _sim_lib.f90wrap_system__get__box_delrx()
    
    @box_delrx.setter
    def box_delrx(self, box_delrx):
        _sim_lib.f90wrap_system__set__box_delrx(box_delrx)
    
    @property
    def box_rho(self):
        """
        Element box_rho ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 48
        
        """
        return _sim_lib.f90wrap_system__get__box_rho()
    
    @box_rho.setter
    def box_rho(self, box_rho):
        _sim_lib.f90wrap_system__set__box_rho(box_rho)
    
    @property
    def ptypes(self):
        """
        Element ptypes ftype=integer pytype=int
        
        
        Defined at src.f95 line 50
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_system__array__ptypes(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ptypes = self._arrays[array_handle]
        else:
            ptypes = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_system__array__ptypes)
            self._arrays[array_handle] = ptypes
        return ptypes
    
    @ptypes.setter
    def ptypes(self, ptypes):
        self.ptypes[...] = ptypes
    
    @property
    def mtypes(self):
        """
        Element mtypes ftype=integer pytype=int
        
        
        Defined at src.f95 line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_system__array__mtypes(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            mtypes = self._arrays[array_handle]
        else:
            mtypes = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_system__array__mtypes)
            self._arrays[array_handle] = mtypes
        return mtypes
    
    @mtypes.setter
    def mtypes(self, mtypes):
        self.mtypes[...] = mtypes
    
    @property
    def mass(self):
        """
        Element mass ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 52
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_system__array__mass(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            mass = self._arrays[array_handle]
        else:
            mass = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_system__array__mass)
            self._arrays[array_handle] = mass
        return mass
    
    @mass.setter
    def mass(self, mass):
        self.mass[...] = mass
    
    @property
    def d(self):
        """
        Element d ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 53
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_system__array__d(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            d = self._arrays[array_handle]
        else:
            d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_system__array__d)
            self._arrays[array_handle] = d
        return d
    
    @d.setter
    def d(self, d):
        self.d[...] = d
    
    @property
    def r(self):
        """
        Element r ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 54
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_system__array__r(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            r = self._arrays[array_handle]
        else:
            r = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_system__array__r)
            self._arrays[array_handle] = r
        return r
    
    @r.setter
    def r(self, r):
        self.r[...] = r
    
    @property
    def vel(self):
        """
        Element vel ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 55
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_system__array__vel(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            vel = self._arrays[array_handle]
        else:
            vel = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_system__array__vel)
            self._arrays[array_handle] = vel
        return vel
    
    @vel.setter
    def vel(self, vel):
        self.vel[...] = vel
    
    @property
    def f(self):
        """
        Element f ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 56
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_system__array__f(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            f = self._arrays[array_handle]
        else:
            f = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_system__array__f)
            self._arrays[array_handle] = f
        return f
    
    @f.setter
    def f(self, f):
        self.f[...] = f
    
    @property
    def s(self):
        """
        Element s ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 57
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_system__array__s(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            s = self._arrays[array_handle]
        else:
            s = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_system__array__s)
            self._arrays[array_handle] = s
        return s
    
    @s.setter
    def s(self, s):
        self.s[...] = s
    
    @property
    def typical_grad(self):
        """
        Element typical_grad ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 58
        
        """
        return _sim_lib.f90wrap_system__get__typical_grad()
    
    @typical_grad.setter
    def typical_grad(self, typical_grad):
        _sim_lib.f90wrap_system__set__typical_grad(typical_grad)
    
    @property
    def grad_scale(self):
        """
        Element grad_scale ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 58
        
        """
        return _sim_lib.f90wrap_system__get__grad_scale()
    
    @grad_scale.setter
    def grad_scale(self, grad_scale):
        _sim_lib.f90wrap_system__set__grad_scale(grad_scale)
    
    @property
    def thermo_temp(self):
        """
        Element thermo_temp ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 60
        
        """
        return _sim_lib.f90wrap_system__get__thermo_temp()
    
    @thermo_temp.setter
    def thermo_temp(self, thermo_temp):
        _sim_lib.f90wrap_system__set__thermo_temp(thermo_temp)
    
    @property
    def thermo_pot(self):
        """
        Element thermo_pot ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 61
        
        """
        return _sim_lib.f90wrap_system__get__thermo_pot()
    
    @thermo_pot.setter
    def thermo_pot(self, thermo_pot):
        _sim_lib.f90wrap_system__set__thermo_pot(thermo_pot)
    
    @property
    def thermo_pre(self):
        """
        Element thermo_pre ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 61
        
        """
        return _sim_lib.f90wrap_system__get__thermo_pre()
    
    @thermo_pre.setter
    def thermo_pre(self, thermo_pre):
        _sim_lib.f90wrap_system__set__thermo_pre(thermo_pre)
    
    @property
    def thermo_sigma(self):
        """
        Element thermo_sigma ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 61
        
        """
        return _sim_lib.f90wrap_system__get__thermo_sigma()
    
    @thermo_sigma.setter
    def thermo_sigma(self, thermo_sigma):
        _sim_lib.f90wrap_system__set__thermo_sigma(thermo_sigma)
    
    def __str__(self):
        ret = ['<system>{\n']
        ret.append('    nthreads : ')
        ret.append(repr(self.nthreads))
        ret.append(',\n    ndim : ')
        ret.append(repr(self.ndim))
        ret.append(',\n    npart : ')
        ret.append(repr(self.npart))
        ret.append(',\n    ntypes : ')
        ret.append(repr(self.ntypes))
        ret.append(',\n    inv_ndim : ')
        ret.append(repr(self.inv_ndim))
        ret.append(',\n    real_npart : ')
        ret.append(repr(self.real_npart))
        ret.append(',\n    box_volume : ')
        ret.append(repr(self.box_volume))
        ret.append(',\n    boxl : ')
        ret.append(repr(self.boxl))
        ret.append(',\n    boxl2 : ')
        ret.append(repr(self.boxl2))
        ret.append(',\n    box_strain : ')
        ret.append(repr(self.box_strain))
        ret.append(',\n    box_delrx : ')
        ret.append(repr(self.box_delrx))
        ret.append(',\n    box_rho : ')
        ret.append(repr(self.box_rho))
        ret.append(',\n    ptypes : ')
        ret.append(repr(self.ptypes))
        ret.append(',\n    mtypes : ')
        ret.append(repr(self.mtypes))
        ret.append(',\n    mass : ')
        ret.append(repr(self.mass))
        ret.append(',\n    d : ')
        ret.append(repr(self.d))
        ret.append(',\n    r : ')
        ret.append(repr(self.r))
        ret.append(',\n    vel : ')
        ret.append(repr(self.vel))
        ret.append(',\n    f : ')
        ret.append(repr(self.f))
        ret.append(',\n    s : ')
        ret.append(repr(self.s))
        ret.append(',\n    typical_grad : ')
        ret.append(repr(self.typical_grad))
        ret.append(',\n    grad_scale : ')
        ret.append(repr(self.grad_scale))
        ret.append(',\n    thermo_temp : ')
        ret.append(repr(self.thermo_temp))
        ret.append(',\n    thermo_pot : ')
        ret.append(repr(self.thermo_pot))
        ret.append(',\n    thermo_pre : ')
        ret.append(repr(self.thermo_pre))
        ret.append(',\n    thermo_sigma : ')
        ret.append(repr(self.thermo_sigma))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

system = System()

class Cell(f90wrap.runtime.FortranModule):
    """
    Module cell
    
    
    Defined at src.f95 lines 106-430
    
    """
    @f90wrap.runtime.register_class("sim_lib.cell_info")
    class cell_info(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=cell_info)
        
        
        Defined at src.f95 lines 115-122
        
        """
        def __init__(self, handle=None):
            """
            self = Cell_Info()
            
            
            Defined at src.f95 lines 115-122
            
            
            Returns
            -------
            this : Cell_Info
            	Object to be constructed
            
            
            Automatically generated constructor for cell_info
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _sim_lib.f90wrap_cell_info_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Cell_Info
            
            
            Defined at src.f95 lines 115-122
            
            Parameters
            ----------
            this : Cell_Info
            	Object to be destructed
            
            
            Automatically generated destructor for cell_info
            """
            if self._alloc:
                _sim_lib.f90wrap_cell_info_finalise(this=self._handle)
        
        @property
        def nt(self):
            """
            Element nt ftype=integer pytype=int
            
            
            Defined at src.f95 line 116
            
            """
            return _sim_lib.f90wrap_cell_info__get__nt(self._handle)
        
        @nt.setter
        def nt(self, nt):
            _sim_lib.f90wrap_cell_info__set__nt(self._handle, nt)
        
        @property
        def nc(self):
            """
            Element nc ftype=integer pytype=int
            
            
            Defined at src.f95 line 117
            
            """
            return _sim_lib.f90wrap_cell_info__get__nc(self._handle)
        
        @nc.setter
        def nc(self, nc):
            _sim_lib.f90wrap_cell_info__set__nc(self._handle, nc)
        
        @property
        def n_xyz(self):
            """
            Element n_xyz ftype=integer pytype=int
            
            
            Defined at src.f95 line 118
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _sim_lib.f90wrap_cell_info__array__n_xyz(self._handle)
            if array_handle in self._arrays:
                n_xyz = self._arrays[array_handle]
            else:
                n_xyz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _sim_lib.f90wrap_cell_info__array__n_xyz)
                self._arrays[array_handle] = n_xyz
            return n_xyz
        
        @n_xyz.setter
        def n_xyz(self, n_xyz):
            self.n_xyz[...] = n_xyz
        
        @property
        def diam(self):
            """
            Element diam ftype=real(kind=dbl) pytype=float
            
            
            Defined at src.f95 line 119
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _sim_lib.f90wrap_cell_info__array__diam(self._handle)
            if array_handle in self._arrays:
                diam = self._arrays[array_handle]
            else:
                diam = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _sim_lib.f90wrap_cell_info__array__diam)
                self._arrays[array_handle] = diam
            return diam
        
        @diam.setter
        def diam(self, diam):
            self.diam[...] = diam
        
        @property
        def ll_a(self):
            """
            Element ll_a ftype=integer pytype=int
            
            
            Defined at src.f95 line 120
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _sim_lib.f90wrap_cell_info__array__ll_a(self._handle)
            if array_handle in self._arrays:
                ll_a = self._arrays[array_handle]
            else:
                ll_a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _sim_lib.f90wrap_cell_info__array__ll_a)
                self._arrays[array_handle] = ll_a
            return ll_a
        
        @ll_a.setter
        def ll_a(self, ll_a):
            self.ll_a[...] = ll_a
        
        @property
        def hoc_a(self):
            """
            Element hoc_a ftype=integer pytype=int
            
            
            Defined at src.f95 line 120
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _sim_lib.f90wrap_cell_info__array__hoc_a(self._handle)
            if array_handle in self._arrays:
                hoc_a = self._arrays[array_handle]
            else:
                hoc_a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _sim_lib.f90wrap_cell_info__array__hoc_a)
                self._arrays[array_handle] = hoc_a
            return hoc_a
        
        @hoc_a.setter
        def hoc_a(self, hoc_a):
            self.hoc_a[...] = hoc_a
        
        @property
        def list_a(self):
            """
            Element list_a ftype=integer pytype=int
            
            
            Defined at src.f95 line 121
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _sim_lib.f90wrap_cell_info__array__list_a(self._handle)
            if array_handle in self._arrays:
                list_a = self._arrays[array_handle]
            else:
                list_a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _sim_lib.f90wrap_cell_info__array__list_a)
                self._arrays[array_handle] = list_a
            return list_a
        
        @list_a.setter
        def list_a(self, list_a):
            self.list_a[...] = list_a
        
        @property
        def maps(self):
            """
            Element maps ftype=integer pytype=int
            
            
            Defined at src.f95 line 122
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _sim_lib.f90wrap_cell_info__array__maps(self._handle)
            if array_handle in self._arrays:
                maps = self._arrays[array_handle]
            else:
                maps = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _sim_lib.f90wrap_cell_info__array__maps)
                self._arrays[array_handle] = maps
            return maps
        
        @maps.setter
        def maps(self, maps):
            self.maps[...] = maps
        
        def __str__(self):
            ret = ['<cell_info>{\n']
            ret.append('    nt : ')
            ret.append(repr(self.nt))
            ret.append(',\n    nc : ')
            ret.append(repr(self.nc))
            ret.append(',\n    n_xyz : ')
            ret.append(repr(self.n_xyz))
            ret.append(',\n    diam : ')
            ret.append(repr(self.diam))
            ret.append(',\n    ll_a : ')
            ret.append(repr(self.ll_a))
            ret.append(',\n    hoc_a : ')
            ret.append(repr(self.hoc_a))
            ret.append(',\n    list_a : ')
            ret.append(repr(self.list_a))
            ret.append(',\n    maps : ')
            ret.append(repr(self.maps))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def init_cell():
        """
        init_cell()
        
        
        Defined at src.f95 lines 129-174
        
        
        """
        _sim_lib.f90wrap_init_cell()
    
    @staticmethod
    def get_2d_icell(idx, idy, n_xyz):
        """
        icell = get_2d_icell(idx, idy, n_xyz)
        
        
        Defined at src.f95 lines 176-186
        
        Parameters
        ----------
        idx : int
        idy : int
        n_xyz : int array
        
        Returns
        -------
        icell : int
        
        """
        icell = _sim_lib.f90wrap_get_2d_icell(idx=idx, idy=idy, n_xyz=n_xyz)
        return icell
    
    @staticmethod
    def construct_2dcell_map():
        """
        construct_2dcell_map()
        
        
        Defined at src.f95 lines 188-209
        
        
        """
        _sim_lib.f90wrap_construct_2dcell_map()
    
    @staticmethod
    def update_shear_2dcell_map():
        """
        update_shear_2dcell_map()
        
        
        Defined at src.f95 lines 211-242
        
        
        """
        _sim_lib.f90wrap_update_shear_2dcell_map()
    
    @staticmethod
    def get_3d_icell(idx, idy, idz, n_xyz):
        """
        icell = get_3d_icell(idx, idy, idz, n_xyz)
        
        
        Defined at src.f95 lines 244-255
        
        Parameters
        ----------
        idx : int
        idy : int
        idz : int
        n_xyz : int array
        
        Returns
        -------
        icell : int
        
        """
        icell = _sim_lib.f90wrap_get_3d_icell(idx=idx, idy=idy, idz=idz, n_xyz=n_xyz)
        return icell
    
    @staticmethod
    def construct_3dcell_map():
        """
        construct_3dcell_map()
        
        
        Defined at src.f95 lines 257-348
        
        
        """
        _sim_lib.f90wrap_construct_3dcell_map()
    
    @staticmethod
    def update_shear_3dcell_map():
        """
        update_shear_3dcell_map()
        
        
        Defined at src.f95 lines 350-417
        
        
        """
        _sim_lib.f90wrap_update_shear_3dcell_map()
    
    @staticmethod
    def get_particle_cell(pos):
        """
        icel = get_particle_cell(pos)
        
        
        Defined at src.f95 lines 419-429
        
        Parameters
        ----------
        pos : float array
        
        Returns
        -------
        icel : int
        
        """
        icel = _sim_lib.f90wrap_get_particle_cell(pos=pos)
        return icel
    
    @property
    def frc_cell_rc_min(self):
        """
        Element frc_cell_rc_min ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 126
        
        """
        return _sim_lib.f90wrap_cell__get__frc_cell_rc_min()
    
    @frc_cell_rc_min.setter
    def frc_cell_rc_min(self, frc_cell_rc_min):
        _sim_lib.f90wrap_cell__set__frc_cell_rc_min(frc_cell_rc_min)
    
    def __str__(self):
        ret = ['<cell>{\n']
        ret.append('    frc_cell_rc_min : ')
        ret.append(repr(self.frc_cell_rc_min))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

cell = Cell()

class Sim_Tools(f90wrap.runtime.FortranModule):
    """
    Module sim_tools
    
    
    Defined at src.f95 lines 432-483
    
    """
    @staticmethod
    def apply_lepbc(pos, box2, box, drx_shear):
        """
        apply_lepbc(pos, box2, box, drx_shear)
        
        
        Defined at src.f95 lines 444-465
        
        Parameters
        ----------
        pos : float array
        box2 : float array
        box : float array
        drx_shear : float
        
        """
        _sim_lib.f90wrap_apply_lepbc(pos=pos, box2=box2, box=box, drx_shear=drx_shear)
    
    @staticmethod
    def put_back_in_box(pos, box):
        """
        put_back_in_box(pos, box)
        
        
        Defined at src.f95 lines 467-482
        
        Parameters
        ----------
        pos : float array
        box : float array
        
        """
        _sim_lib.f90wrap_put_back_in_box(pos=pos, box=box)
    
    _dt_array_initialisers = []
    

sim_tools = Sim_Tools()

class Potentials(f90wrap.runtime.FortranModule):
    """
    Module potentials
    
    
    Defined at src.f95 lines 485-632
    
    """
    @staticmethod
    def init_poly_force_field(type_bn, arg):
        """
        init_poly_force_field(type_bn, arg)
        
        
        Defined at src.f95 lines 523-554
        
        Parameters
        ----------
        type_bn : str
        arg : float array
        
        """
        _sim_lib.f90wrap_init_poly_force_field(type_bn=type_bn, arg=arg)
    
    @staticmethod
    def dij_swap_ipl(di, dj, arg):
        """
        dij2 = dij_swap_ipl(di, dj, arg)
        
        
        Defined at src.f95 lines 559-569
        
        Parameters
        ----------
        di : float
        dj : float
        arg : Const_Vec
        
        Returns
        -------
        dij2 : float
        
        """
        dij2 = _sim_lib.f90wrap_dij_swap_ipl(di=di, dj=dj, arg=arg._handle)
        return dij2
    
    @staticmethod
    def rcut_swap_ipl(dij2, arg):
        """
        rc2 = rcut_swap_ipl(dij2, arg)
        
        
        Defined at src.f95 lines 571-580
        
        Parameters
        ----------
        dij2 : float
        arg : Const_Vec
        
        Returns
        -------
        rc2 : float
        
        """
        rc2 = _sim_lib.f90wrap_rcut_swap_ipl(dij2=dij2, arg=arg._handle)
        return rc2
    
    @staticmethod
    def pot_swap_ipl(r2, dij2, arg):
        """
        pot = pot_swap_ipl(r2, dij2, arg)
        
        
        Defined at src.f95 lines 582-594
        
        Parameters
        ----------
        r2 : float
        dij2 : float
        arg : Const_Vec
        
        Returns
        -------
        pot : float
        
        """
        pot = _sim_lib.f90wrap_pot_swap_ipl(r2=r2, dij2=dij2, arg=arg._handle)
        return pot
    
    @staticmethod
    def phir_swap_ipl(r2, dij2, arg):
        """
        phir = phir_swap_ipl(r2, dij2, arg)
        
        
        Defined at src.f95 lines 596-609
        
        Parameters
        ----------
        r2 : float
        dij2 : float
        arg : Const_Vec
        
        Returns
        -------
        phir : float
        
        """
        phir = _sim_lib.f90wrap_phir_swap_ipl(r2=r2, dij2=dij2, arg=arg._handle)
        return phir
    
    @staticmethod
    def full_swap_ipl(r2, dij2, arg, arg_out):
        """
        full_swap_ipl(r2, dij2, arg, arg_out)
        
        
        Defined at src.f95 lines 611-628
        
        Parameters
        ----------
        r2 : float
        dij2 : float
        arg : Const_Vec
        arg_out : Const_Vec
        
        """
        _sim_lib.f90wrap_full_swap_ipl(r2=r2, dij2=dij2, arg=arg._handle, \
            arg_out=arg_out._handle)
    
    _dt_array_initialisers = []
    

potentials = Potentials()

class Pairwise(f90wrap.runtime.FortranModule):
    """
    Module pairwise
    
    
    Defined at src.f95 lines 634-850
    
    """
    @staticmethod
    def init_pairwise_arrays():
        """
        init_pairwise_arrays()
        
        
        Defined at src.f95 lines 654-667
        
        
        """
        _sim_lib.f90wrap_init_pairwise_arrays()
    
    @staticmethod
    def create_chain_list():
        """
        create_chain_list()
        
        
        Defined at src.f95 lines 669-686
        
        
        """
        _sim_lib.f90wrap_create_chain_list()
    
    @staticmethod
    def update_chain_list(id, old_cell, new_cell):
        """
        update_chain_list(id, old_cell, new_cell)
        
        
        Defined at src.f95 lines 688-714
        
        Parameters
        ----------
        id : int
        old_cell : int
        new_cell : int
        
        """
        _sim_lib.f90wrap_update_chain_list(id=id, old_cell=old_cell, new_cell=new_cell)
    
    @staticmethod
    def compute_poly_everything():
        """
        compute_poly_everything()
        
        
        Defined at src.f95 lines 716-789
        
        
        """
        _sim_lib.f90wrap_compute_poly_everything()
    
    @staticmethod
    def compute_hessian(hessian):
        """
        compute_hessian(hessian)
        
        
        Defined at src.f95 lines 791-839
        
        Parameters
        ----------
        hessian : float array
        
        """
        _sim_lib.f90wrap_compute_hessian(hessian=hessian)
    
    @staticmethod
    def kron_delta(a, b):
        """
        d = kron_delta(a, b)
        
        
        Defined at src.f95 lines 841-849
        
        Parameters
        ----------
        a : int
        b : int
        
        Returns
        -------
        d : int
        
        """
        d = _sim_lib.f90wrap_kron_delta(a=a, b=b)
        return d
    
    @property
    def pair_ncontacts(self):
        """
        Element pair_ncontacts ftype=integer pytype=int
        
        
        Defined at src.f95 line 645
        
        """
        return _sim_lib.f90wrap_pairwise__get__pair_ncontacts()
    
    @property
    def pair_count(self):
        """
        Element pair_count ftype=integer pytype=int
        
        
        Defined at src.f95 line 646
        
        """
        return _sim_lib.f90wrap_pairwise__get__pair_count()
    
    @pair_count.setter
    def pair_count(self, pair_count):
        _sim_lib.f90wrap_pairwise__set__pair_count(pair_count)
    
    @property
    def pair_lookupbond(self):
        """
        Element pair_lookupbond ftype=integer pytype=int
        
        
        Defined at src.f95 line 647
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_pairwise__array__pair_lookupbond(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pair_lookupbond = self._arrays[array_handle]
        else:
            pair_lookupbond = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_pairwise__array__pair_lookupbond)
            self._arrays[array_handle] = pair_lookupbond
        return pair_lookupbond
    
    @pair_lookupbond.setter
    def pair_lookupbond(self, pair_lookupbond):
        self.pair_lookupbond[...] = pair_lookupbond
    
    @property
    def pair_first(self):
        """
        Element pair_first ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 648
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_pairwise__array__pair_first(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pair_first = self._arrays[array_handle]
        else:
            pair_first = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_pairwise__array__pair_first)
            self._arrays[array_handle] = pair_first
        return pair_first
    
    @pair_first.setter
    def pair_first(self, pair_first):
        self.pair_first[...] = pair_first
    
    @property
    def pair_second(self):
        """
        Element pair_second ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 648
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_pairwise__array__pair_second(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pair_second = self._arrays[array_handle]
        else:
            pair_second = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_pairwise__array__pair_second)
            self._arrays[array_handle] = pair_second
        return pair_second
    
    @pair_second.setter
    def pair_second(self, pair_second):
        self.pair_second[...] = pair_second
    
    @property
    def pair_third(self):
        """
        Element pair_third ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 648
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_pairwise__array__pair_third(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pair_third = self._arrays[array_handle]
        else:
            pair_third = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_pairwise__array__pair_third)
            self._arrays[array_handle] = pair_third
        return pair_third
    
    @pair_third.setter
    def pair_third(self, pair_third):
        self.pair_third[...] = pair_third
    
    @property
    def pair_fourth(self):
        """
        Element pair_fourth ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 648
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_pairwise__array__pair_fourth(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pair_fourth = self._arrays[array_handle]
        else:
            pair_fourth = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_pairwise__array__pair_fourth)
            self._arrays[array_handle] = pair_fourth
        return pair_fourth
    
    @pair_fourth.setter
    def pair_fourth(self, pair_fourth):
        self.pair_fourth[...] = pair_fourth
    
    @property
    def pair_rij(self):
        """
        Element pair_rij ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 649
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_pairwise__array__pair_rij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pair_rij = self._arrays[array_handle]
        else:
            pair_rij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_pairwise__array__pair_rij)
            self._arrays[array_handle] = pair_rij
        return pair_rij
    
    @pair_rij.setter
    def pair_rij(self, pair_rij):
        self.pair_rij[...] = pair_rij
    
    @property
    def typical_contact_force(self):
        """
        Element typical_contact_force ftype=real(kind=dbl) pytype=float
        
        
        Defined at src.f95 line 650
        
        """
        return _sim_lib.f90wrap_pairwise__get__typical_contact_force()
    
    @typical_contact_force.setter
    def typical_contact_force(self, typical_contact_force):
        _sim_lib.f90wrap_pairwise__set__typical_contact_force(typical_contact_force)
    
    @property
    def pair_list_ni(self):
        """
        Element pair_list_ni ftype=integer pytype=int
        
        
        Defined at src.f95 line 651
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_pairwise__array__pair_list_ni(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pair_list_ni = self._arrays[array_handle]
        else:
            pair_list_ni = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_pairwise__array__pair_list_ni)
            self._arrays[array_handle] = pair_list_ni
        return pair_list_ni
    
    @pair_list_ni.setter
    def pair_list_ni(self, pair_list_ni):
        self.pair_list_ni[...] = pair_list_ni
    
    @property
    def pair_list_ij(self):
        """
        Element pair_list_ij ftype=integer pytype=int
        
        
        Defined at src.f95 line 652
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _sim_lib.f90wrap_pairwise__array__pair_list_ij(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            pair_list_ij = self._arrays[array_handle]
        else:
            pair_list_ij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _sim_lib.f90wrap_pairwise__array__pair_list_ij)
            self._arrays[array_handle] = pair_list_ij
        return pair_list_ij
    
    @pair_list_ij.setter
    def pair_list_ij(self, pair_list_ij):
        self.pair_list_ij[...] = pair_list_ij
    
    def __str__(self):
        ret = ['<pairwise>{\n']
        ret.append('    pair_ncontacts : ')
        ret.append(repr(self.pair_ncontacts))
        ret.append(',\n    pair_count : ')
        ret.append(repr(self.pair_count))
        ret.append(',\n    pair_lookupbond : ')
        ret.append(repr(self.pair_lookupbond))
        ret.append(',\n    pair_first : ')
        ret.append(repr(self.pair_first))
        ret.append(',\n    pair_second : ')
        ret.append(repr(self.pair_second))
        ret.append(',\n    pair_third : ')
        ret.append(repr(self.pair_third))
        ret.append(',\n    pair_fourth : ')
        ret.append(repr(self.pair_fourth))
        ret.append(',\n    pair_rij : ')
        ret.append(repr(self.pair_rij))
        ret.append(',\n    typical_contact_force : ')
        ret.append(repr(self.typical_contact_force))
        ret.append(',\n    pair_list_ni : ')
        ret.append(repr(self.pair_list_ni))
        ret.append(',\n    pair_list_ij : ')
        ret.append(repr(self.pair_list_ij))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

pairwise = Pairwise()

