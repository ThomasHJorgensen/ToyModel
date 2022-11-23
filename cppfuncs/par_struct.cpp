typedef struct par_struct
{
 double R;
 double beta;
 double div_A_share;
 double div_H_share;
 double wage_const;
 double wage_K;
 double rho;
 double home_power_1;
 double home_power_2;
 double home_load_1;
 double home_load_2;
 double depre_K;
 double depre_H;
 double accum_K;
 double accum_H;
 int T;
 double max_time;
 int num_A;
 double max_A;
 int num_K;
 double max_K;
 int num_H;
 double max_H;
 int num_power;
 int num_love;
 double max_love;
 double sigma_love;
 int num_shock_love;
 bool do_bargaining;
 int seed;
 int simT;
 int simN;
 int threads;
 int num_pre_h;
 int num_pre_C;
 int num_pre_Q;
 double max_pre_Q;
 double* grid_A;
 double* grid_Aw;
 double* grid_Am;
 double* grid_K;
 double* grid_H;
 double* grid_Hw;
 double* grid_Hm;
 double* grid_power;
 double* grid_love;
 double* grid_shock_love;
 double* grid_weight_love;
 double* grid_pre_Ctot;
 double* grid_pre_hours;
 double* grid_pre_Qtot;
} par_struct;

double get_double_par_struct(par_struct* x, char* name){

 if( strcmp(name,"R") == 0 ){ return x->R; }
 else if( strcmp(name,"beta") == 0 ){ return x->beta; }
 else if( strcmp(name,"div_A_share") == 0 ){ return x->div_A_share; }
 else if( strcmp(name,"div_H_share") == 0 ){ return x->div_H_share; }
 else if( strcmp(name,"wage_const") == 0 ){ return x->wage_const; }
 else if( strcmp(name,"wage_K") == 0 ){ return x->wage_K; }
 else if( strcmp(name,"rho") == 0 ){ return x->rho; }
 else if( strcmp(name,"home_power_1") == 0 ){ return x->home_power_1; }
 else if( strcmp(name,"home_power_2") == 0 ){ return x->home_power_2; }
 else if( strcmp(name,"home_load_1") == 0 ){ return x->home_load_1; }
 else if( strcmp(name,"home_load_2") == 0 ){ return x->home_load_2; }
 else if( strcmp(name,"depre_K") == 0 ){ return x->depre_K; }
 else if( strcmp(name,"depre_H") == 0 ){ return x->depre_H; }
 else if( strcmp(name,"accum_K") == 0 ){ return x->accum_K; }
 else if( strcmp(name,"accum_H") == 0 ){ return x->accum_H; }
 else if( strcmp(name,"max_time") == 0 ){ return x->max_time; }
 else if( strcmp(name,"max_A") == 0 ){ return x->max_A; }
 else if( strcmp(name,"max_K") == 0 ){ return x->max_K; }
 else if( strcmp(name,"max_H") == 0 ){ return x->max_H; }
 else if( strcmp(name,"max_love") == 0 ){ return x->max_love; }
 else if( strcmp(name,"sigma_love") == 0 ){ return x->sigma_love; }
 else if( strcmp(name,"max_pre_Q") == 0 ){ return x->max_pre_Q; }
 else {return NAN;}

}


int get_int_par_struct(par_struct* x, char* name){

 if( strcmp(name,"T") == 0 ){ return x->T; }
 else if( strcmp(name,"num_A") == 0 ){ return x->num_A; }
 else if( strcmp(name,"num_K") == 0 ){ return x->num_K; }
 else if( strcmp(name,"num_H") == 0 ){ return x->num_H; }
 else if( strcmp(name,"num_power") == 0 ){ return x->num_power; }
 else if( strcmp(name,"num_love") == 0 ){ return x->num_love; }
 else if( strcmp(name,"num_shock_love") == 0 ){ return x->num_shock_love; }
 else if( strcmp(name,"seed") == 0 ){ return x->seed; }
 else if( strcmp(name,"simT") == 0 ){ return x->simT; }
 else if( strcmp(name,"simN") == 0 ){ return x->simN; }
 else if( strcmp(name,"threads") == 0 ){ return x->threads; }
 else if( strcmp(name,"num_pre_h") == 0 ){ return x->num_pre_h; }
 else if( strcmp(name,"num_pre_C") == 0 ){ return x->num_pre_C; }
 else if( strcmp(name,"num_pre_Q") == 0 ){ return x->num_pre_Q; }
 else {return -9999;}

}


bool get_bool_par_struct(par_struct* x, char* name){

 if( strcmp(name,"do_bargaining") == 0 ){ return x->do_bargaining; }
 else {return false;}

}


double* get_double_p_par_struct(par_struct* x, char* name){

 if( strcmp(name,"grid_A") == 0 ){ return x->grid_A; }
 else if( strcmp(name,"grid_Aw") == 0 ){ return x->grid_Aw; }
 else if( strcmp(name,"grid_Am") == 0 ){ return x->grid_Am; }
 else if( strcmp(name,"grid_K") == 0 ){ return x->grid_K; }
 else if( strcmp(name,"grid_H") == 0 ){ return x->grid_H; }
 else if( strcmp(name,"grid_Hw") == 0 ){ return x->grid_Hw; }
 else if( strcmp(name,"grid_Hm") == 0 ){ return x->grid_Hm; }
 else if( strcmp(name,"grid_power") == 0 ){ return x->grid_power; }
 else if( strcmp(name,"grid_love") == 0 ){ return x->grid_love; }
 else if( strcmp(name,"grid_shock_love") == 0 ){ return x->grid_shock_love; }
 else if( strcmp(name,"grid_weight_love") == 0 ){ return x->grid_weight_love; }
 else if( strcmp(name,"grid_pre_Ctot") == 0 ){ return x->grid_pre_Ctot; }
 else if( strcmp(name,"grid_pre_hours") == 0 ){ return x->grid_pre_hours; }
 else if( strcmp(name,"grid_pre_Qtot") == 0 ){ return x->grid_pre_Qtot; }
 else {return NULL;}

}


