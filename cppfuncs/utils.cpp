// functions related to utility and environment.
#ifndef MAIN
#define UTILS
#include "myheader.cpp"
#endif

namespace utils {
    double labor_implied(double leisure, double hours, par_struct* par){
        double labor = par->max_time - leisure - hours;
        if (labor<0){
            labor = 0.0;
        }
        return labor;
    }
    double trans_K(double K, double labor, par_struct* par){
        return (1.0-par->depre_K)*K + par->accum_K*labor;
    }
    double trans_H(double H, double hours, par_struct* par){
        return (1.0-par->depre_H)*H + par->accum_H*hours;
    }

    double wage_func(double K, par_struct* par){
        return std::exp(par->wage_const + par->wage_K * K);
    }

    double home_prod_load(double H,int gender,par_struct *par){
        return par->home_load_1 + par->home_load_2*H;

    }
    double time_input(double hours_w, double hours_m, double Hw, double Hm,par_struct *par){
        double input_w = 0.0;
        if(hours_w>0.0) {
            input_w = home_prod_load(Hw,woman,par)*pow(hours_w,par->home_power_2);
        }
        double input_m = 0.0;
        if(hours_m>0.0) {
            input_m = home_prod_load(Hm,man,par)*pow(hours_m,par->home_power_2);
        }

        double Qtot = input_w + input_m;
        return Qtot;
    }
    double home_prod_Qtot(double Qtot,double market, par_struct *par){
        
        double input_time = pow(Qtot , par->home_power_1/par->home_power_2);
        double input_market = pow(market,par->home_power_1);

        double Q = pow(input_time + input_market , 1.0/par->home_power_1);
        
        return Q;
    }

    double home_prod(double hours_w, double hours_m, double Hw, double Hm,double market, par_struct *par){
        
        double Qtot = time_input(hours_w,hours_m,Hw,Hm,par);
        return home_prod_Qtot(Qtot,market,par);

    }

    double home_prod_single(double hours, double market, double H, int gender, par_struct* par){
        double hours_w = hours;
        double hours_m = hours;
        double Hw = H;
        double Hm = H;
        if (gender==man){
            hours_w = 0.0;
        } else {
            hours_m = 0.0;
        }
        return home_prod(hours_w,hours_m,Hw,Hm,market,par);
    }

    double util(double cons,double leisure,double home_prod, int gender,par_struct *par){
        
        double curv = 1.0 - par->rho;
        double curv_inv = 1.0/curv; 

        return pow(cons,curv) * curv_inv + pow(leisure,curv) * curv_inv + pow(home_prod,curv) * curv_inv;
    }

    double util_couple(double cons_w, double cons_m, double home_prod,double leisure_w, double leisure_m, int iP,int iL,par_struct* par){
        double power = par->grid_power[iP];
        double love = par->grid_love[iP];

        double Uw = util(cons_w,leisure_w,home_prod,woman,par);
        double Um = util(cons_m,leisure_m,home_prod,man,par);

        return power*Uw + (1.0-power)*Um + love;
    }
}