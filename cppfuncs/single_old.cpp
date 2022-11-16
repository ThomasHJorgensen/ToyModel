#ifndef MAIN
#define SINGLE
#include "myheader.cpp"
#endif

namespace single {
    typedef struct {
        
        double A;
        double H;
        int gender;
        par_struct *par;

    } solver_single_last_struct;

    double objfunc_single_last(unsigned n, const double *x, double *grad, void *solver_data_in){
        double love = 0.0;

        // unpack
        solver_single_last_struct *solver_data = (solver_single_last_struct *) solver_data_in;
        
        double hours = x[0];
        double market = x[1];
        int gender = solver_data->gender;
        double A = solver_data->A;
        double H = solver_data->H;
        par_struct *par = solver_data->par;

        // consumption and leisure
        double cons = A - market;
        double leisure = par->max_time - hours;

        // home production
        double hours_w = hours;
        double hours_m = hours;
        double Hw = H;
        double Hm = H;
        if (gender==man){
            hours_w = 0.0;
        } else {
            hours_m = 0.0;
        }
        double home_prod = utils::home_prod(hours_w,hours_m,Hw,Hm,market,par);

        return - utils::util(cons,leisure,home_prod,gender,par);

    }

    typedef struct {
        
        double A;
        double K;
        double H;
        double* V_next;
        int gender;
        par_struct *par;

    } solver_single_struct;

    double value_of_choice(double cons,double leisure,double hours,double market,double A,double K,double H,int gender,double* V_next,par_struct *par){

        // flow-utility
        double hours_w = hours;
        double hours_m = hours;
        double Hw = H;
        double Hm = H;
        if (gender==man){
            hours_w = 0.0;
        } else {
            hours_m = 0.0;
        }
        double home_prod = utils::home_prod(hours_w,hours_m,Hw,Hm,market,par);
        double Util = utils::util(cons,leisure,home_prod,gender,par);
        
        // continuation value
        // grids
        double *grid_A = par->grid_Aw;
        double *grid_H = par->grid_Hw;
        double *grid_K = par->grid_K;
        if (gender==man){
            grid_A = par->grid_Am;
            grid_H = par->grid_Hm;
        }

        // states
        double labor = par->max_time - hours - leisure;
        double income = utils::wage_func(H,par)*labor;

        double H_next = (1.0-par->depre_H)*H + par->accum_H * hours;
        double K_next = (1.0-par->depre_H)*K + par->accum_K * labor;
        double A_next = par->R*( A + income - cons - market);
        
        double V_next_interp = tools::interp_3d(grid_H,grid_K,grid_A,par->num_H,par->num_K,par->num_A,V_next,H_next,K_next,A_next);
        
        // return discounted sum
        return Util + par->beta*V_next_interp;
    }

    double objfunc_single(unsigned n, const double *x, double *grad, void *solver_data_in){

        // unpack
        solver_single_struct *solver_data = (solver_single_struct *) solver_data_in;
        
        double cons = x[0];
        double leisure = x[1];
        double hours = x[2];
        double market = x[3];
        int gender = solver_data->gender;
        double A = solver_data->A;
        double K = solver_data->K;
        double H = solver_data->H;
        par_struct *par = solver_data->par;

        printf("1:%2.2f 2:%2.2f 3:%2.2f 4:%2.2f (constr2:%2.2f)\n",x[0],x[1],x[2],x[3] , par->max_time - hours - leisure);
printf("val:%2.4f\n",value_of_choice(cons,leisure,hours,market,A,K,H,gender,solver_data->V_next,par));

        return - value_of_choice(cons,leisure,hours,market,A,K,H,gender,solver_data->V_next,par);

    }

    typedef struct {
        
        double A;
        double K;
        double H;
        
        par_struct *par;

    } constr_single_struct;

    double constraint1(unsigned n, const double *x, double *grad, void *constr_data_in){
        // savings must be non-negative (constraint is that it is positive by default)
        // unpack
        constr_single_struct *constr_data = (constr_single_struct *) constr_data_in;

        double cons = x[0];
        double leisure = x[1];
        double hours = x[2];
        double market = x[3];
        double A = constr_data->A;
        double K = constr_data->K;
        double H = constr_data->H;
        par_struct *par = constr_data->par;

        double labor = par->max_time - hours - leisure;
        double income = utils::wage_func(H,par)*labor;

        double savings = A + income - cons - market;

        return - savings;
        // return savings;

    }
    
    double constraint2(unsigned n, const double *x, double *grad, void *constr_data_in){
        // hours must not be more than available
        // unpack
        constr_single_struct *constr_data = (constr_single_struct *) constr_data_in;

        double cons = x[0];
        double leisure = x[1];
        double hours = x[2];
        double market = x[3];
        double A = constr_data->A;
        double K = constr_data->K;
        double H = constr_data->H;
        par_struct *par = constr_data->par;

        double labor = par->max_time - hours - leisure;

        return - labor;
        // return labor;

    }

    void solve_single(int t,sol_struct *sol,par_struct *par){
        
        // terminal period
        if (t == (par->T-1)){
            
            #pragma omp parallel num_threads(par->threads)
            {

                // 1. allocate objects for solver
                solver_single_last_struct* solver_data = new solver_single_last_struct;
                
                int dim = 2;
                double lb[2],ub[2],x[2];
                
                auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); //NLOPT_LN_BOBYQA NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
                double minf=0.0;

                // 2. loop over states
                #pragma omp for
                for (int iH=0; iH<par->num_H;iH++){
                    int iK = 0; // human capital does not matter here. So re-use this value below
                        for (int iA=0; iA<par->num_A;iA++){
                            int idx = index::index4(t,iH,iK,iA,par->T,par->num_H,par->num_K,par->num_A);

                            // state variables
                            double Aw = par->grid_Aw[iA];
                            double Am = par->grid_Am[iA];

                            double Hw = par->grid_Hw[iH];
                            double Hm = par->grid_Hm[iH];

                            // solve for women and men
                            // WOMEN
                            // settings
                            solver_data->A = Aw;
                            solver_data->H = Hw;
                            solver_data->gender = woman;
                            solver_data->par = par;
                            nlopt_set_min_objective(opt, objfunc_single_last, solver_data);
                                
                            // bounds on home production hours and market purchases
                            lb[0] = 1.0e-8;
                            ub[0] = par->max_time;
                            lb[1] = 1.0e-8;
                            ub[1] = solver_data->A;
                            nlopt_set_lower_bounds(opt, lb);
                            nlopt_set_upper_bounds(opt, ub);

                            // optimize
                            if (iA==0){
                                x[0] = ub[0]/1.1;
                                x[1] = ub[1]/1.1;
                            }
                            
                            nlopt_optimize(opt, x, &minf);

                            // store results
                            sol->h_w_single[idx] = x[0];
                            sol->m_w_single[idx] = x[1];

                            sol->c_w_single[idx] = solver_data->A - sol->m_w_single[idx];
                            sol->l_w_single[idx] = par->max_time - sol->h_w_single[idx];
                            
                            sol->V_w_single[idx] = -minf;

                            // insert for all values of human capital
                            for (int iK_now=1; iK_now<par->num_K;iK_now++){
                                int idx_now = index::index4(t,iH,iK_now,iA,par->T,par->num_H,par->num_K,par->num_A);
                                sol->c_w_single[idx_now] = sol->c_w_single[idx];
                                sol->l_w_single[idx_now] = sol->l_w_single[idx];
                                sol->h_w_single[idx_now] = sol->h_w_single[idx];
                                sol->m_w_single[idx_now] = sol->m_w_single[idx];
                                sol->V_w_single[idx_now] = sol->V_w_single[idx];
                            }


                            // MEN
                            // settings
                            solver_data->A = Am;
                            solver_data->H = Hm;
                            solver_data->gender = man;
                            solver_data->par = par;
                            nlopt_set_min_objective(opt, objfunc_single_last, solver_data);
                                
                            // bounds on home production hours and market purchases
                            lb[0] = 1.0e-8;
                            ub[0] = par->max_time;
                            lb[1] = 1.0e-8;
                            ub[1] = solver_data->A;
                            nlopt_set_lower_bounds(opt, lb);
                            nlopt_set_upper_bounds(opt, ub);

                            // optimize
                            if (iA==0){
                                x[0] = ub[0]/1.1;
                                x[1] = ub[1]/1.1;
                            }
                            
                            nlopt_optimize(opt, x, &minf);

                            // store results
                            sol->h_m_single[idx] = x[0];
                            sol->m_m_single[idx] = x[1];

                            sol->c_m_single[idx] = solver_data->A - sol->m_m_single[idx];
                            sol->l_m_single[idx] = par->max_time - sol->h_m_single[idx];
                            
                            sol->V_m_single[idx] = -minf;

                            // insert for all values of human capital
                            for (int iK_now=1; iK_now<par->num_K;iK_now++){
                                int idx_now = index::index4(t,iH,iK_now,iA,par->T,par->num_H,par->num_K,par->num_A);
                                sol->c_m_single[idx_now] = sol->c_m_single[idx];
                                sol->l_m_single[idx_now] = sol->l_m_single[idx];
                                sol->h_m_single[idx_now] = sol->h_m_single[idx];
                                sol->m_m_single[idx_now] = sol->m_m_single[idx];
                                sol->V_m_single[idx_now] = sol->V_m_single[idx];
                            }

                        }
                    
                }


                nlopt_destroy(opt);

            } // pragma
        } else {
                printf("t:%d\n",t);
//             #pragma omp parallel num_threads(par->threads)
//             {

//                 // 1. allocate objects for solver
//                 solver_single_struct* solver_data = new solver_single_struct;
//                 constr_single_struct* constr_data = new constr_single_struct;
                
//                 int dim = 4;
//                 double lb[4],ub[4],x[4];
                
//                 auto opt = nlopt_create(NLOPT_LN_COBYLA, dim); //NLOPT_LN_BOBYQA NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
//                 double minf=0.0;

//                 // 2. loop over states
//                 #pragma omp for
//                 for (int iH=0; iH<par->num_H;iH++){
//                     for (int iK=1; iK<par->num_K;iK++){
//                         for (int iA=0; iA<par->num_A;iA++){
//                             int idx = index::index4(t,iH,iK,iA,par->T,par->num_H,par->num_K,par->num_A);

//                             // state variables
//                             double Aw = par->grid_Aw[iA];
//                             double Am = par->grid_Am[iA];

//                             double Kw = par->grid_K[iK];
//                             double Km = par->grid_K[iK];

//                             double Hw = par->grid_Hw[iH];
//                             double Hm = par->grid_Hm[iH];

//                             // solve for women and men
//                             // WOMEN
//                             // settings
//                             solver_data->A = Aw;
//                             solver_data->K = Kw;
//                             solver_data->H = Hw;
//                             solver_data->V_next = &sol->V_w_single[index::index4(t+1,0,0,0,par->T,par->num_H,par->num_K,par->num_A)];
                            
//                             solver_data->gender = woman;
//                             solver_data->par = par;
//                             nlopt_set_min_objective(opt, objfunc_single, solver_data);
//                             nlopt_set_xtol_rel(opt, 1.0e-5);
                                
//                             // bounds on consumption, leisure, home production hours and market purchases
//                             lb[0] = 1.0e-8;
//                             ub[0] = 10000000000000000000000.0;
//                             lb[1] = 1.0e-8;
//                             ub[1] = par->max_time;
//                             lb[2] = 1.0e-8;
//                             ub[2] = par->max_time;
//                             lb[3] = 1.0e-8;
//                             ub[3] = 10000000000000000000000.0;
//                             nlopt_set_lower_bounds(opt, lb);
//                             nlopt_set_upper_bounds(opt, ub);

//                             // budget and time constraints
//                             constr_data->A = Aw;
//                             constr_data->K = Kw;
//                             constr_data->H = Hw;
//                             constr_data->par = par;
//                             nlopt_add_inequality_constraint(opt, constraint1, constr_data,1.0e-6);
//                             nlopt_add_inequality_constraint(opt, constraint2, constr_data,1.0e-6);

//                             // optimize
//                             // if (iA==0){
//                             //     x[0] = 0.1;
//                             //     x[1] = ub[1]/2.0;
//                             //     x[2] = ub[0]/2.0;
//                             //     x[3] = ub[1]/1.1;
//                             // }
                            
//                             if (iA==0){
//                                 x[0] = 0.1;
//                                 x[1] = 0.001;
//                                 x[2] = 0.001;
//                                 x[3] = 0.001;
//                             }
//                             // printf("t:%2.4f\n",x[0]);
//                             if (iH==6 & iK==3 & iA==26){
//                                 nlopt_optimize(opt, x, &minf);
//                             }
//                             // store results
//                             sol->c_w_single[idx] = x[0];
//                             sol->l_w_single[idx] = x[1];
//                             sol->h_w_single[idx] = x[2];
//                             sol->m_w_single[idx] = x[3];

//                             sol->V_w_single[idx] = -minf;

//                             // printf("iH:%d,iK:%d,iA:%d,cons:%2.3f\n",iH,iK,iA, x[0]);

// // HER: lav objective function + constraint
                            


//                         }
//                     }
                    
//                 }


//                 nlopt_destroy(opt);

//             } // pragma
            
        }   
        
    }


}