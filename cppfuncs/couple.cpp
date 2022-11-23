
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

    void cons_from_C_and_Q(double Ctot, double Qtot, int iP,sol_struct* sol, par_struct* par, double* cons_w, double* cons_m, double* market){
        // interpolate pre-computed consumption allocation (TODO: speed up by interpolating simultaneously)
        int idx = index::index3(iP,0,0,par->num_power,par->num_pre_Q,par->num_pre_C);
        cons_w[0] = tools::interp_2d(par->grid_pre_Qtot,par->grid_pre_Ctot,par->num_pre_Q,par->num_pre_C,
                                        &sol->pre_cons_w[idx],Qtot,Ctot);
        
        cons_m[0] = tools::interp_2d(par->grid_pre_Qtot,par->grid_pre_Ctot,par->num_pre_Q,par->num_pre_C,
                                        &sol->pre_cons_m[idx],Qtot,Ctot);
        
        market[0] = Ctot - cons_w[0] - cons_m[0];
    }

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
        double Ctot = par->R*A; // last period -> consume everything and not working.

        // interpolate consumption allocation from pre-computed intra-temporal problem
        cons_from_C_and_Q(Ctot, Qtot,iP,sol,par,&solver_data->cons_w,&solver_data->cons_m,&solver_data->market);

        // value of choice
        double home_prod = utils::home_prod_Qtot(Qtot,solver_data->market,par);
        return  - utils::util_couple(solver_data->cons_w,solver_data->cons_m,home_prod,leisure_w,leisure_m,iP,iL,par);

    }


    void solve_couple_last(sol_struct *sol,par_struct *par){
        #pragma omp parallel num_threads(par->threads)
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

    double calc_marital_surplus(double V_remain,double V_single,par_struct* par){
        return V_remain - V_single;
    }

    
    void check_participation_constraints(int t,sol_struct* sol, par_struct* par){
        // loop through all states except power
        #pragma omp parallel num_threads(par->threads)
        {

            // allocate memory to store relevant objects for the participation constraint check
            int shape_tmp = par->num_power;
            double* tmp_Vw = new double[shape_tmp];
            double* tmp_Vm = new double[shape_tmp];
            double* tmp_cons_w = new double[shape_tmp];
            double* tmp_cons_m = new double[shape_tmp];
            double* tmp_market = new double[shape_tmp];
            double* tmp_leisure_w = new double[shape_tmp];
            double* tmp_leisure_m = new double[shape_tmp];
            double* tmp_hours_w = new double[shape_tmp];
            double* tmp_hours_m = new double[shape_tmp];

            int num_list = 5;
            double** list_couple_w = new double*[num_list]; 
            double** list_couple_m = new double*[num_list]; 
            double** list_raw_w = new double*[num_list]; 
            double** list_raw_m = new double*[num_list]; 
            double* list_single_w = new double[num_list]; 
            double* list_single_m = new double[num_list]; 

            double* Sw = new double[par->num_power];
            double* Sm = new double[par->num_power];

            index::index_couple_struct* idx_couple = new index::index_couple_struct;   

            // solve for values of reminaing a couple
            #pragma omp for
            for (int iL=0; iL<par->num_love; iL++){
                for (int iH_w=0; iH_w<par->num_H; iH_w++){
                    for (int iH_m=0; iH_m<par->num_H; iH_m++){
                        for (int iK_w=0; iK_w<par->num_K; iK_w++){
                            for (int iK_m=0; iK_m<par->num_K; iK_m++){
                                for (int iA=0; iA<par->num_A; iA++){

                                    // indices
                                    int idx_single_w = index::index4(t,iH_w,iK_w,iA,par->T,par->num_H,par->num_K,par->num_A);
                                    int idx_single_m = index::index4(t,iH_m,iK_m,iA,par->T,par->num_H,par->num_K,par->num_A);

                                    idx_couple->t = t;
                                    idx_couple->iL = iL;
                                    idx_couple->iH_w = iH_w;
                                    idx_couple->iH_m = iH_m;
                                    idx_couple->iK_w = iK_w;
                                    idx_couple->iK_m = iK_m;
                                    idx_couple->iA = iA;
                                    idx_couple->par = par; 

                                    // setup temporary arrays with the one dimension being power
                                    for (int iP=0; iP<par->num_power; iP++){
                                        int idx_tmp = idx_couple->idx(iP);

                                        tmp_Vw[iP] = sol->Vw_remain[idx_tmp];
                                        tmp_Vm[iP] = sol->Vm_remain[idx_tmp];
                                        tmp_cons_w[iP] = sol->cons_w_remain[idx_tmp];
                                        tmp_cons_m[iP] = sol->cons_m_remain[idx_tmp];
                                        tmp_market[iP] = sol->market_remain[idx_tmp];
                                        tmp_leisure_w[iP] = sol->leisure_w_remain[idx_tmp];
                                        tmp_leisure_m[iP] = sol->leisure_m_remain[idx_tmp];
                                        tmp_hours_w[iP] = sol->hours_w_remain[idx_tmp];
                                        tmp_hours_m[iP] = sol->hours_m_remain[idx_tmp];

                                        // marital surplus
                                        Sw[iP] = calc_marital_surplus(tmp_Vw[iP],sol->Vw_single[idx_single_w],par);
                                        Sm[iP] = calc_marital_surplus(tmp_Vm[iP],sol->Vm_single[idx_single_m],par);
                                    }
                            
                                    // setup relevant lists
                                    int i = 0;
                                    list_couple_w[i] = sol->Vw_couple; i++;
                                    list_couple_w[i] = sol->cons_w_couple; i++;
                                    list_couple_w[i] = sol->market_couple; i++;
                                    list_couple_w[i] = sol->leisure_w_couple; i++;
                                    list_couple_w[i] = sol->hours_w_couple; i++;
                                    i = 0;
                                    list_couple_m[i] = sol->Vm_couple; i++;
                                    list_couple_m[i] = sol->cons_m_couple; i++;
                                    list_couple_m[i] = sol->market_couple; i++;
                                    list_couple_m[i] = sol->leisure_m_couple; i++;
                                    list_couple_m[i] = sol->hours_m_couple; i++;

                                    i = 0;
                                    list_raw_w[i] = tmp_Vw; i++;
                                    list_raw_w[i] = tmp_cons_w; i++;
                                    list_raw_w[i] = tmp_market; i++;
                                    list_raw_w[i] = tmp_leisure_w; i++;
                                    list_raw_w[i] = tmp_hours_w; i++;
                                    i = 0;
                                    list_raw_m[i] = tmp_Vm; i++;
                                    list_raw_m[i] = tmp_cons_m; i++;
                                    list_raw_m[i] = tmp_market; i++;
                                    list_raw_m[i] = tmp_leisure_m; i++;
                                    list_raw_m[i] = tmp_hours_m; i++;

                                    // TODO: reduction in value from transitioning? prefereably directly in the grids. 
                                    i = 0;
                                    list_single_w[i] = sol->Vw_single[idx_single_w]; i++;
                                    list_single_w[i] = sol->cons_w_single[idx_single_w]; i++;
                                    list_single_w[i] = sol->market_w_single[idx_single_w]; i++;
                                    list_single_w[i] = sol->leisure_w_single[idx_single_w]; i++;
                                    list_single_w[i] = sol->hours_w_single[idx_single_w]; i++;
                                    i = 0;
                                    list_single_m[i] = sol->Vm_single[idx_single_m]; i++;
                                    list_single_m[i] = sol->cons_m_single[idx_single_m]; i++;
                                    list_single_m[i] = sol->market_m_single[idx_single_m]; i++;
                                    list_single_m[i] = sol->leisure_m_single[idx_single_m]; i++;
                                    list_single_m[i] = sol->hours_m_single[idx_single_m]; i++;

                                    // update solution
                                    if (par->do_bargaining){
                                        tools::update_bargaining(sol->power_idx,sol->power,Sw,Sm,
                                                                            idx_couple,
                                                                            list_couple_w,list_couple_m,
                                                                            list_raw_w,list_raw_m,
                                                                            list_single_w,list_single_m,
                                                                            num_list, par);
                                    } else {
                                        for (int iP=0; iP<par->num_power; iP++){
                                            int idx_c = idx_couple->idx(iP);

                                            if (Sw[iP]>0.0 & Sm[iP] > 0.0){
                                                for (int i=0; i< num_list; i++){
                                                    list_couple_w[i][idx_c] = list_raw_w[i][iP];
                                                    list_couple_m[i][idx_c] = list_raw_m[i][iP];
                                                }

                                                sol->power_idx[idx_c] = par->grid_power[iP];
                                                sol->power[idx_c] = iP;   

                                            } else { // divorce
                                                for (int i=0; i< num_list; i++){
                                                    list_couple_w[i][idx_c] = list_single_w[i];
                                                    list_couple_m[i][idx_c] = list_single_m[i];
                                                }

                                                sol->power_idx[idx_c] = -1.0;
                                                sol->power[idx_c] = -1;

                                            }

                                        }
                                    }

                                } // wealth
                            } // human capital, man
                        } // human capital, woman
                    } // home capital, man
                } // home capital, woman
            } // love

            // delete pointers
            delete[] list_couple_w;
            delete[] list_couple_m;
            delete[] list_raw_w;
            delete[] list_raw_m;
            delete list_single_w;
            delete list_single_m;

            delete Sw;
            delete Sm;

            delete tmp_cons_w;
            delete tmp_cons_m;
            delete tmp_market;
            delete tmp_leisure_w;
            delete tmp_leisure_m;
            delete tmp_hours_w;
            delete tmp_hours_m;

        } // pragma

    }

    typedef struct {
        
        double A;
        double Hw;
        double Hm;
        double Kw;
        double Km;
        double love;
        int iP;
        int t;
        par_struct *par;
        sol_struct *sol;

        double cons_w;
        double cons_m;
        double market;

        double leisure_w, hours_w, labor_w;
        double leisure_m, hours_m, labor_m;

        double resources;
        double Qtot;

        double Vw, Vm;

    } solver_struct;

    double time_constraint_w(unsigned n, const double *x, double *grad, void *solver_data_in){
        solver_struct *solver_data = (solver_struct *) solver_data_in;
        
        double leisure = x[0];
        double hours = x[1];

        return -(solver_data->par->max_time - leisure - hours);
    }
    double time_constraint_m(unsigned n, const double *x, double *grad, void *solver_data_in){
        solver_struct *solver_data = (solver_struct *) solver_data_in;
        
        double leisure = x[2];
        double hours = x[3];

        return -(solver_data->par->max_time - leisure - hours);
    }

    double objfunc_Ctot(unsigned n, const double *x, double *grad, void *solver_data_in){
        // unpack
        solver_struct *solver_data = (solver_struct *) solver_data_in;
        double A = solver_data->A;
        double Hw = solver_data->Hw;
        double Hm = solver_data->Hm;
        double Kw = solver_data->Kw;
        double Km = solver_data->Km;
        double love = solver_data->love;
        int iP = solver_data->iP;
        int t = solver_data->t;
        par_struct* par = solver_data->par;
        sol_struct* sol = solver_data->sol;

        double Qtot = solver_data->Qtot;
        double resources = solver_data->resources;

        double leisure_w = solver_data->leisure_w;
        double hours_w = solver_data->hours_w;
        double labor_w = solver_data->labor_w;
        double leisure_m = solver_data->leisure_m;
        double hours_m = solver_data->hours_m;
        double labor_m = solver_data->labor_m;
        
        double Ctot = x[0];

        // instantaneous utility
        cons_from_C_and_Q(Ctot, Qtot,iP,sol,par,&solver_data->cons_w,&solver_data->cons_m,&solver_data->market);
        double cons_w = solver_data->cons_w;
        double cons_m = solver_data->cons_m;
        double market = solver_data->market;

        double home_prod = utils::home_prod(hours_w,hours_m,Hw,Hm,market,par);
        double Vw = utils::util(cons_w,leisure_w,home_prod,woman,par) + love;
        double Vm = utils::util(cons_m,leisure_m,home_prod,man,par) + love;

        // continuation value [be be sped up by re-using everything across genders]
        // next period states (power is fixed/end-of-period and potentially updated in subsequent barganing)
        int idx_next = index::index8(t,iP,0,0,0,0,0,0,
                                    par->T,par->num_power,par->num_love,par->num_H,par->num_H,par->num_K,par->num_K,par->num_A);

        double Hw_next = utils::trans_H(Hw, hours_w, par);
        double Hm_next = utils::trans_H(Hm, hours_m, par);
        double Kw_next = utils::trans_K(Kw, labor_w, par);
        double Km_next = utils::trans_K(Km, labor_m, par);
        double A_next = resources - Ctot;

        double EVw_plus = 0.0;
        double EVm_plus = 0.0;
        for (int i_love = 0; i_love<par->num_shock_love; i_love++){
            double love_next = love + par->grid_shock_love[i_love];
            double weight = par->grid_weight_love[i_love];

            EVw_plus += weight * tools::interp_6d(
                            par->grid_love,par->grid_H,par->grid_H,par->grid_K,par->grid_K,par->grid_A,
                            par->num_love,par->num_H,par->num_H,par->num_K,par->num_K,par->num_A,
                            &sol->Vw_couple[idx_next], 
                            love_next,Hw_next,Hm_next,Kw_next,Km_next,A_next);

            EVm_plus += weight * tools::interp_6d(
                            par->grid_love,par->grid_H,par->grid_H,par->grid_K,par->grid_K,par->grid_A,
                            par->num_love,par->num_H,par->num_H,par->num_K,par->num_K,par->num_A,
                            &sol->Vm_couple[idx_next], 
                            love_next,Hw_next,Hm_next,Kw_next,Km_next,A_next);
        }

        // individual values
        Vw += par->beta*EVw_plus;
        Vm += par->beta*EVm_plus;

        // store on solver_data to be used to store solution
        solver_data->Vw = Vw;
        solver_data->Vm = Vm;

        // return negative weighted values
        double power = par->grid_power[iP];
        return -(power*Vw + (1.0-power)*Vm);

    }

    double objfunc(unsigned n, const double *x, double *grad, void *solver_data_in){
        // unpack
        solver_struct *solver_data = (solver_struct *) solver_data_in;
        double A = solver_data->A;
        double Hw = solver_data->Hw;
        double Hm = solver_data->Hm;
        double Kw = solver_data->Kw;
        double Km = solver_data->Km;
        par_struct* par = solver_data->par;

        double leisure_w = x[0];
        double hours_w = x[1];
        double leisure_m = x[2];
        double hours_m = x[3];

        // add implied variables to solver_data
        double labor_w = utils::labor_implied(leisure_w,hours_w,par);
        double labor_m = utils::labor_implied(leisure_m,hours_m,par);

        solver_data->Qtot = utils::time_input(hours_w,hours_m,Hw,Hm,par);
        solver_data->resources = utils::resources_couple(labor_w, labor_m, Hw, Hm, Kw, Km, A,par);

        solver_data->leisure_w = leisure_w;
        solver_data->hours_w = hours_w;
        solver_data->labor_w = labor_w;
        solver_data->leisure_m = leisure_m;
        solver_data->hours_m = hours_m;
        solver_data->labor_m = labor_m;

        // solve for optimal total consumption, Ctot
        const int dim = 1;
        double lb[dim],ub[dim],Ctot[dim];
        auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LN_BOBYQA NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
        nlopt_set_ftol_abs(opt,1.e-8);
        double minf=0.0;

        nlopt_set_min_objective(opt, objfunc_Ctot, solver_data);

        // bounds
        lb[0] = 0.000000001;
        ub[0] = solver_data->resources;
        nlopt_set_lower_bounds(opt, lb);
        nlopt_set_upper_bounds(opt, ub);

        Ctot[0] = 0.5*ub[0];
        nlopt_optimize(opt, Ctot, &minf);

        // penalty on the inequality constraints
        labor_w = par->max_time - leisure_w - hours_w;
        labor_m = par->max_time - leisure_m - hours_m;
        double penalty = 0.0;
        if (labor_w<0.0){
            penalty += labor_w*1000.0;
        }
        if (labor_m<0.0){
            penalty += labor_m*1000.0;
        }

        return minf + penalty;

    }

    void solve_couple_period(int t,sol_struct* sol, par_struct* par){
        // loop through all states except power
        #pragma omp parallel num_threads(par->threads)
        {
            // setup numerical solver to solve for time allocation
            solver_struct* solver_data = new solver_struct;
                    
            const int dim = 4;
            double lb[dim],ub[dim],x[dim];
            
            auto opt = nlopt_create(NLOPT_LN_COBYLA, dim); // NLOPT_LN_BOBYQA NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
            nlopt_set_ftol_abs(opt,1.e-8);
            double minf=0.0;

            // solve for values of reminaing a couple
            #pragma omp for
            for (int iP=0; iP<par->num_power; iP++){
                for (int iL=0; iL<par->num_love; iL++){
                    double love = par->grid_love[iL];
                    for (int iH_w=0; iH_w<par->num_H; iH_w++){
                        double Hw = par->grid_H[iH_w];
                        for (int iH_m=0; iH_m<par->num_H; iH_m++){
                            double Hm = par->grid_H[iH_m];
                            for (int iK_w=0; iK_w<par->num_K; iK_w++){
                                double Kw = par->grid_K[iK_w];
                                for (int iK_m=0; iK_m<par->num_K; iK_m++){
                                    double Km = par->grid_K[iK_m];
                                    for (int iA=0; iA<par->num_A; iA++){
                                        double A = par->grid_A[iA];

                                        int idx = index::index8(t,iP,iL,iH_w,iH_m,iK_w,iK_m,iA,
                                                    par->T,par->num_power,par->num_love,par->num_H,par->num_H,par->num_K,par->num_K,par->num_A);


                                        // solve for (leisure_w, hours_w, leisure_m, hours_m)
                                        solver_data->A = A;
                                        solver_data->Hw = Hw;
                                        solver_data->Hm = Hm;
                                        solver_data->Kw = Kw;
                                        solver_data->Km = Km;
                                        solver_data->love = love;
                                        solver_data->iP = iP;
                                        solver_data->t = t;
                                        solver_data->par = par;
                                        solver_data->sol = sol;
                                        nlopt_set_min_objective(opt, objfunc, solver_data);

                                        // bounds and cosntraints
                                        for (int i=0; i<dim; i++){
                                            lb[i] = 0.000000001;
                                            ub[i] = par->max_time;
                                        }
                                        nlopt_set_lower_bounds(opt, lb);
                                        nlopt_set_upper_bounds(opt, ub);

                                        nlopt_add_inequality_constraint(opt, time_constraint_w, solver_data, 1e-8);
                                        nlopt_add_inequality_constraint(opt, time_constraint_m, solver_data, 1e-8);
                                        printf("%d,%d,%d,%d,%d,%d,%d\n",iP,iL,iH_w,iH_m,iK_w,iK_m,iA);
                                        // optimize over leisure and home production time
                                        for (int i=0; i<dim; i++){
                                            x[i] = 0.3*ub[i];
                                        }
                                        nlopt_optimize(opt, x, &minf);

                                        // store results
                                        double leisure_w = x[0];
                                        double hours_w = x[1];
                                        double leisure_m = x[2];
                                        double hours_m = x[3];

                                        sol->leisure_w_remain[idx] = leisure_w;
                                        sol->hours_w_remain[idx] = hours_w;
                                        sol->leisure_m_remain[idx] = leisure_m;
                                        sol->hours_m_remain[idx] = hours_m;

                                        sol->cons_w_remain[idx] = solver_data->cons_w;
                                        sol->cons_m_remain[idx] = solver_data->cons_m;
                                        sol->market_remain[idx] = solver_data->market;

                                        sol->Vw_remain[idx] = solver_data->Vw;
                                        sol->Vm_remain[idx] = solver_data->Vm;

                                    } // wealth
                                } // human capital, man
                            } // human capital, woman
                        } // home capital, man
                    } // home capital, woman
                } // love
            } // power


        } // pragma

    }
    
    void solve_couple(int t,sol_struct *sol,par_struct *par){
        
        
        if (t==par->T-1){
            solve_couple_last(sol,par);
        } else {
            solve_couple_period(t,sol,par);
        }

        // check participation constraints
        check_participation_constraints(t,sol,par);

    }


    
}