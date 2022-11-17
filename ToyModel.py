import numpy as np
import numba as nb
import scipy.optimize as optimize

from EconModel import EconModelClass
from consav.grids import nonlinspace
from consav import linear_interp, linear_interp_1d
from consav import quadrature

# set gender indication as globals
woman = 1
man = 2

class HouseholdModelClass(EconModelClass):
    
    def settings(self):
        """ fundamental settings """

        # a. namespaces
        self.namespaces = []
        
        # b. other attributes
        self.other_attrs = []
        
        # c. savefolder
        self.savefolder = 'saved'
        
        # d. cpp
        self.cpp_filename = 'cppfuncs/solve.cpp'
        self.cpp_options = {'compiler':'vs'}
        
    def setup(self):
        par = self.par
        
        par.R = 1.03
        par.beta = 1.0/par.R # Discount factor
        
        # divorce
        par.div_A_share = 0.5 # divorce share of wealth to wife
        par.div_H_share =  1.0 # how much of the home capital survives a divorce
        
        # wages
        par.wage_const = 0.0
        par.wage_K = 0.05

        # Utility: 
        par.rho = 2.0 

        # home production
        par.home_power_1 = -0.1
        par.home_power_2 = 0.5
        par.home_load_1 = 1.0
        par.home_load_2 = 1.0   

        # state accumulation
        par.depre_K = 0.1
        par.depre_H = 0.1

        par.accum_K = 1.0
        par.accum_H = 1.0 
        
        # state variables
        par.T = 3
        par.max_time = 1.0
        
        # wealth
        par.num_A = 30
        par.max_A = 15.0

        # human capital
        par.num_K = 10
        par.max_K = 10.0

        # home capital
        par.num_H = 15
        par.max_H = 10.0        
        
        # bargaining power
        par.num_power = 21

        # love/match quality
        par.num_love = 41
        par.max_love = 1.0

        par.sigma_love = 0.1
        par.num_shock_love = 5

        # simulation
        par.seed = 9210
        par.simT = par.T
        par.simN = 50_000

        # cpp
        par.threads = 8

        # pre-computation
        par.num_pre_h = 50
        par.num_pre_C = 100
        
    def allocate(self):
        par = self.par
        sol = self.sol
        sim = self.sim

        # setup grids
        self.setup_grids()
        
        # singles
        shape_single = (par.T,par.num_H,par.num_K,par.num_A)
        sol.V_w_single = np.nan + np.ones(shape_single)
        sol.V_m_single = np.nan + np.ones(shape_single)
        
        sol.c_w_single = np.nan + np.ones(shape_single)
        sol.l_w_single = np.nan + np.ones(shape_single)
        sol.h_w_single = np.nan + np.ones(shape_single)
        sol.m_w_single = np.nan + np.ones(shape_single)
        
        sol.c_m_single = np.nan + np.ones(shape_single)
        sol.l_m_single = np.nan + np.ones(shape_single)
        sol.h_m_single = np.nan + np.ones(shape_single)
        sol.m_m_single = np.nan + np.ones(shape_single)

        # pre-computation
        shape_pre_single = (par.num_H,par.num_pre_h,par.num_pre_C)
        sol.pre_cons_w_single = np.nan + np.ones(shape_pre_single)
        sol.pre_market_w_single = np.nan + np.ones(shape_pre_single)
        sol.pre_cons_m_single = np.nan + np.ones(shape_pre_single)
        sol.pre_market_m_single = np.nan + np.ones(shape_pre_single)
                            

        # couples
        # shape_couple = (par.T,par.num_power,par.num_love,par.num_A)
        # sol.Vw_couple = np.nan + np.ones(shape_couple)
        # sol.Vm_couple = np.nan + np.ones(shape_couple)
        
        # sol.Cw_priv_couple = np.nan + np.ones(shape_couple)
        # sol.Cm_priv_couple = np.nan + np.ones(shape_couple)
        # sol.C_pub_couple = np.nan + np.ones(shape_couple)
        # sol.C_tot_couple = np.nan + np.ones(shape_couple)

        # sol.Vw_remain_couple = np.nan + np.ones(shape_couple)
        # sol.Vm_remain_couple = np.nan + np.ones(shape_couple)
        
        # sol.Cw_priv_remain_couple = np.nan + np.ones(shape_couple)
        # sol.Cm_priv_remain_couple = np.nan + np.ones(shape_couple)
        # sol.C_pub_remain_couple = np.nan + np.ones(shape_couple)
        # sol.C_tot_remain_couple = np.nan + np.ones(shape_couple)

        # sol.power_idx = np.zeros(shape_couple,dtype=np.int_)
        # sol.power = np.zeros(shape_couple)

        
    def setup_grids(self):
        par = self.par
        
        # wealth. Single grids are such to avoid interpolation
        par.grid_A = nonlinspace(1.0e-6,par.max_A,par.num_A,1.1)

        par.grid_Aw = par.div_A_share * par.grid_A
        par.grid_Am = (1.0 - par.div_A_share) * par.grid_A

        # market human capital
        par.grid_K = nonlinspace(0.0,par.max_K,par.num_K,1.1)

        # home capital
        par.grid_H = nonlinspace(0.0,par.max_H,par.num_H,1.1)

        par.grid_Hw = par.div_H_share * par.grid_H
        par.grid_Hm = par.div_H_share * par.grid_H

        # power. non-linear grid with more mass in both tails.
        odd_num = np.mod(par.num_power,2)
        first_part = nonlinspace(0.0,0.5,(par.num_power+odd_num)//2,1.3)
        last_part = np.flip(1.0 - nonlinspace(0.0,0.5,(par.num_power-odd_num)//2 + 1,1.3))[1:]
        par.grid_power = np.append(first_part,last_part)

        # love grid and shock
        if par.num_love>1:
            par.grid_love = np.linspace(-par.max_love,par.max_love,par.num_love)
        else:
            par.grid_love = np.array([0.0])

        if par.sigma_love<=1.0e-6:
            par.num_shock_love = 1
            par.grid_shock_love,par.grid_weight_love = np.array([0.0]),np.array([1.0])

        else:
            par.grid_shock_love,par.grid_weight_love = quadrature.normal_gauss_hermite(par.sigma_love,par.num_shock_love)

        # pre-computation
        par.grid_pre_Ctot = np.linspace(par.grid_A[0],par.grid_A[-1],par.num_pre_C)
        par.grid_pre_hours = np.linspace(0.0001,par.max_time,par.num_pre_h)

    def solve(self):
        sol = self.sol
        par = self.par 

        # setup grids
        self.setup_grids()

        # solve model
        self.cpp.solve(sol,par)

        

    # def simulate(self):
    #     sol = self.sol
    #     sim = self.sim
    #     par = self.par

    #     if par.do_cpp:
    #         self.cpp.simulate(sim,sol,par)

    #     else:
    #         simulate(sim,sol,par)

    #     # total consumption
    #     sim.Cw_tot = sim.Cw_priv + sim.Cw_pub
    #     sim.Cm_tot = sim.Cm_priv + sim.Cm_pub
    #     sim.C_tot = sim.Cw_priv + sim.Cm_priv + sim.Cw_pub
