
// functions for solving model for singles.
#ifndef MAIN
#define COUPLE
#include "myheader.cpp"
#endif

namespace couple {
    typedef struct {
        
        double A;
        double Hw;
        double Hm;
        int iP;
        par_struct *par;
        sol_struct *sol;

        double cons_w;
        double cons_m;
        double market;

    } solver_last_struct;

    double objfunc_last(unsigned n, const double *x, double *grad, void *solver_data_in){
        // fixed states/choices that does not matter here
        int iL = 0;

        // unpack
        solver_last_struct *solver_data = (solver_last_struct *) solver_data_in;
        double A = solver_data->A;
        double Hw = solver_data->Hw;
        double Hm = solver_data->Hm;
        int iP = solver_data->iP;
        par_struct* par = solver_data->par;
        sol_struct* sol = solver_data->sol;

        double leisure_w = x[0];
        double leisure_m = x[1];

        // implied variables
        double hours_w = par->max_time - leisure_w;
        double hours_m = par->max_time - leisure_m;

        double Qtot = utils::time_input(hours_w,hours_m,Hw,Hm,par);
        double Ctot = par->R*A;
        
        // interpolate pre-computed consumption allocation (TODO: speed up by interpolating simultaneously)
        int idx = index::index3(iP,0,0,par->num_power,par->num_pre_Q,par->num_pre_C);
        double cons_w = tools::interp_2d(par->grid_pre_Qtot,par->grid_pre_Ctot,par->num_pre_Q,par->num_pre_C,
                                        &sol->pre_cons_w[idx],Qtot,Ctot);
        
        double cons_m = tools::interp_2d(par->grid_pre_Qtot,par->grid_pre_Ctot,par->num_pre_Q,par->num_pre_C,
                                        &sol->pre_cons_m[idx],Qtot,Ctot);
        
        double market = Ctot - cons_w - cons_m;

        // store in solver data for later use
        solver_data->cons_w = cons_w;
        solver_data->cons_m = cons_m;
        solver_data->market = market;

        // value of choice
        double home_prod = utils::home_prod_Qtot(Qtot,market,par);
        return  - utils::util_couple(cons_w,cons_m,home_prod,leisure_w,leisure_m,iP,iL,par);

    }


    void solve_couple_last(sol_struct *sol,par_struct *par){
        #pragma omp parallel num_threads(1) //num_threads(par->threads)
        {
            // setup numerical solver
            solver_last_struct* solver_data = new solver_last_struct;
                    
            const int dim = 2;
            double lb[dim],ub[dim],x[dim];
            
            auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LN_BOBYQA NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
            nlopt_set_ftol_abs(opt,1.e-8);
            double minf=0.0;

            // fixed states
            int t = par->T-1;
            int iL = 0;
            int iK_w = 0;
            int iK_m = 0;

            // solve for values of reminaing a couple
            #pragma omp for
            for (int iP=0; iP<par->num_power; iP++){
                for (int iH_w=0; iH_w<par->num_H; iH_w++){
                    double Hw = par->grid_H[iH_w];
                    for (int iH_m=0; iH_m<par->num_H; iH_m++){
                        double Hm = par->grid_H[iH_m];
                        for (int iA=0; iA<par->num_A; iA++){
                            double A = par->grid_A[iA];
                            
                            int idx = index::index8(t,iP,iL,iH_w,iH_m,iK_w,iK_m,iA,
                                                    par->T,par->num_power,par->num_love,par->num_H,par->num_H,par->num_K,par->num_K,par->num_A);

                            // solve for leisure_w and leisure_m 
                            solver_data->A = A;
                            solver_data->Hw = Hw;
                            solver_data->Hm = Hm;
                            solver_data->iP = iP;
                            solver_data->par = par;
                            solver_data->sol = sol;
                            nlopt_set_min_objective(opt, objfunc_last, solver_data);

                            // bounds on leisure time
                            lb[0] = 0.000000001;
                            ub[0] = par->max_time;
                            lb[1] = 0.000000001;
                            ub[1] = par->max_time;
                            nlopt_set_lower_bounds(opt, lb);
                            nlopt_set_upper_bounds(opt, ub);

                            // optimize over leisure time
                            x[0] = 0.5*ub[0];
                            x[1] = 0.5*ub[0];
                            nlopt_optimize(opt, x, &minf);

                            // store results 
                            sol->leisure_w_remain[idx] = x[0];
                            sol->leisure_m_remain[idx] = x[1];
                            sol->hours_w_remain[idx] = par->max_time - sol->leisure_w_remain[idx];
                            sol->hours_m_remain[idx] = par->max_time - sol->leisure_m_remain[idx];
                            sol->labor_w_remain[idx] = 0.0; // not working in the last period 
                            sol->labor_m_remain[idx] = 0.0; // not working in the last period 

                            sol->cons_w_remain[idx] = solver_data->cons_w;
                            sol->cons_m_remain[idx] = solver_data->cons_m;
                            sol->market_remain[idx] = solver_data->market;

                            // value functions
                            double home_prod = utils::home_prod(sol->hours_w_remain[idx], sol->hours_m_remain[idx], Hw, Hm,sol->market_remain[idx],par);
                            double love = par->grid_love[iL];
                            sol->Vw_remain[idx] = utils::util(sol->cons_w_remain[idx],sol->leisure_w_remain[idx],home_prod, woman,par) + love;
                            sol->Vm_remain[idx] = utils::util(sol->cons_m_remain[idx],sol->leisure_m_remain[idx],home_prod, man,par) + love;

                            // copy solution into other states
                            for (int iL_now=0; iL_now<par->num_love; iL_now++){
                                double love_now = par->grid_love[iL_now];

                                for (int iK_w_now=0; iK_w_now<par->num_K; iK_w_now++){
                                    for (int iK_m_now=0; iK_m_now<par->num_K; iK_m_now++){
                                        int idx_now = index::index8(t,iP,iL_now,iH_w,iH_m,iK_w_now,iK_m_now,iA,
                                                                    par->T,par->num_power,par->num_love,par->num_H,par->num_H,par->num_K,par->num_K,par->num_A);

                                        sol->leisure_w_remain[idx_now] = sol->leisure_w_remain[idx];
                                        sol->leisure_m_remain[idx_now] = sol->leisure_m_remain[idx];
                                        sol->hours_w_remain[idx_now] = sol->hours_w_remain[idx];
                                        sol->hours_m_remain[idx_now] = sol->hours_m_remain[idx];
                                        sol->labor_w_remain[idx_now] = sol->labor_w_remain[idx];
                                        sol->labor_m_remain[idx_now] = sol->labor_m_remain[idx];

                                        sol->cons_w_remain[idx_now] = sol->cons_w_remain[idx];
                                        sol->cons_m_remain[idx_now] = sol->cons_m_remain[idx];
                                        sol->market_remain[idx_now] = sol->market_remain[idx];

                                        sol->Vw_remain[idx_now] = sol->Vw_remain[idx] - love + love_now;
                                        sol->Vm_remain[idx_now] = sol->Vm_remain[idx] - love + love_now;

                                    } // human capital, man
                                } // human capital, woman
                            } // love

                        } // wealth
                    } // home capital, man
                } // home capital, woman
            } // power

        } // pragma
    }

    
    void solve_couple(int t,sol_struct *sol,par_struct *par){
        
        
        if (t==par->T-1){
            solve_couple_last(sol,par);

        } else {
            
        }

        // check participation

    }


    
}