
// functions for solving model for singles.
#ifndef MAIN
#define COUPLE
#include "myheader.cpp"
#endif

namespace couple {
    
    void solve_couple(int t,sol_struct *sol,par_struct *par){
        
        #pragma omp parallel num_threads(par->threads)
        {
            // allocate memory to store relevant objects for the participation constraint check
            int shape_tmp = par->num_power;
            double* tmp_Vw = new double[shape_tmp];
            double* tmp_Vm = new double[shape_tmp];
            double* tmp_Cw_priv = new double[shape_tmp];
            double* tmp_Cm_priv = new double[shape_tmp];
            double* tmp_C_pub = new double[shape_tmp];
            double* tmp_marg_V = new double[shape_tmp];

            int num = 6;
            double** list_couple = new double*[num]; 
            double** list_single = new double*[num]; 
            double** list_raw = new double*[num]; 

            double* Sw = new double[par->num_power];
            double* Sm = new double[par->num_power];

            index_couple_struct* idx_couple = new index_couple_struct;

            // a. solve for values of reminaing a couple
            #pragma omp for
            for (int iP=0; iP<par->num_power; iP++){

                // continuation values
                int idx_next = index::index4(t+1,iP,0,0,par->T,par->num_power,par->num_love,par->num_A);
                if (t==(par->T-1)){ // does not matter in last period-> fix at some valid index
                    idx_next = 0;
                }
                double *Vw_next = &sol->Vw_couple[idx_next];
                double *Vm_next = &sol->Vm_couple[idx_next];
                double *marg_V_next = &sol->marg_V_couple[idx_next]; 

                for (int iL=0; iL<par->num_love; iL++){

                    // solve for all values in grid_A.
                    if (par->do_egm){
                        solve_remain_Agrid_egm(t,iP,iL,Vw_next,Vm_next,marg_V_next,sol,par);

                    } else {
                        solve_remain_Agrid_vfi(t,iP,iL,Vw_next,Vm_next,sol,par);

                    }
 
                } // love
            } // power
        }
    }


    
}