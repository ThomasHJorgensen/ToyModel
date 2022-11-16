#define MAIN
#include "myheader.h"

// Pre-compute: SINGLES
typedef struct {
        
    double C_tot;
    double hours;
    double leisure;
    double H;
    int gender;
    par_struct *par;

} solver_single_pre_struct;

double objfunc_pre_single(unsigned n, const double *x, double *grad, void *solver_data_in){
    double love = 0.0;

    // unpack
    solver_single_pre_struct *solver_data = (solver_single_pre_struct *) solver_data_in;
    
    double cons_share = x[0];
    double C_tot = solver_data->C_tot;
    double hours = solver_data->hours;
    double leisure = solver_data->leisure;
    double H = solver_data->H;
    int gender = solver_data->gender;
    par_struct *par = solver_data->par;

    // consumption and leisure
    double cons = cons_share * C_tot;
    double market = C_tot - cons;

    // home production
    double home_prod = utils::home_prod_single(hours,market,H,gender,par);

    return - utils::util(cons,leisure,home_prod,gender,par);

}
void solve_intraperiod_single(double* cons,double* market, double C_tot,double leisure, double hours, double H,int gender,par_struct *par){
    // setup numerical solver
    solver_single_pre_struct* solver_data = new solver_single_pre_struct;
            
    int dim = 1;
    double lb[1],ub[1],x[1];
    
    auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LN_BOBYQA NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
    double minf=0.0;

    solver_data->C_tot = C_tot;
    solver_data->hours = hours;
    solver_data->leisure = leisure;
    solver_data->H = H;
    solver_data->gender = gender;
    solver_data->par = par;
    nlopt_set_min_objective(opt, objfunc_pre_single, solver_data);
        
    // bounds on share of total spending on consumption
    lb[0] = 0.000000001;
    ub[0] = 0.999999999;
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    // optimize
    x[0] = 0.5;
    nlopt_optimize(opt, x, &minf);
    nlopt_destroy(opt);

    // store results in pointers
    cons[0] = x[0] * C_tot;
    market[0] = C_tot - cons[0];
}


void precompute_single(sol_struct* sol, par_struct* par){
    double leisure = 1.0; // independent of leisure
    #pragma omp parallel num_threads(par->threads)
    {
        #pragma omp for
        for (int iH=0; iH<par->num_H;iH++){
            for (int i_h=0; i_h < par->num_pre_h; i_h++){
                for (int i_c=0; i_c < par->num_pre_C; i_c++){
                    int idx = index::index3(iH,i_h,i_c,par->num_H,par->num_pre_h,par->num_pre_C);

                    double C_tot = par->grid_pre_Ctot[i_c];
                    double hours = par->grid_pre_hours[i_h];

                    // women
                    double H = par->grid_Hw[iH];
                    solve_intraperiod_single(&sol->pre_cons_w_single[idx], &sol->pre_market_w_single[idx] , C_tot,leisure,hours,H,woman,par);

                    // men
                    H = par->grid_Hm[iH];
                    solve_intraperiod_single(&sol->pre_cons_m_single[idx], &sol->pre_market_m_single[idx] , C_tot,leisure,hours,H,man,par);
                    
                }
            }
        }
    }
}

/////////////
// 5. MAIN //
/////////////

EXPORT void solve(sol_struct *sol, par_struct *par){

    // pre-compute intra-temporal optimal allocation
    precompute_single(sol,par);

    // loop backwards
    for (int t = par->T-1; t >= 0; t--){

        single::solve_single(t,sol,par);
        // couple::solve_couple(t,sol,par);

    }
}


EXPORT void simulate(sim_struct *sim, sol_struct *sol, par_struct *par){
    
    sim::model(sim,sol,par);

}
