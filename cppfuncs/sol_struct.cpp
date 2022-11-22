typedef struct sol_struct
{
 double* Vw_single;
 double* Vm_single;
 double* cons_w_single;
 double* leisure_w_single;
 double* hours_w_single;
 double* market_w_single;
 double* cons_m_single;
 double* leisure_m_single;
 double* hours_m_single;
 double* market_m_single;
 double* pre_cons_w_single;
 double* pre_market_w_single;
 double* pre_cons_m_single;
 double* pre_market_m_single;
 double* Vw_remain;
 double* Vm_remain;
 double* cons_w_remain;
 double* cons_m_remain;
 double* market_remain;
 double* leisure_w_remain;
 double* leisure_m_remain;
 double* hours_w_remain;
 double* hours_m_remain;
 double* labor_w_remain;
 double* labor_m_remain;
 int* power_idx;
 double* power;
 double* pre_cons_w;
 double* pre_cons_m;
 double* pre_market;
} sol_struct;

double* get_double_p_sol_struct(sol_struct* x, char* name){

 if( strcmp(name,"Vw_single") == 0 ){ return x->Vw_single; }
 else if( strcmp(name,"Vm_single") == 0 ){ return x->Vm_single; }
 else if( strcmp(name,"cons_w_single") == 0 ){ return x->cons_w_single; }
 else if( strcmp(name,"leisure_w_single") == 0 ){ return x->leisure_w_single; }
 else if( strcmp(name,"hours_w_single") == 0 ){ return x->hours_w_single; }
 else if( strcmp(name,"market_w_single") == 0 ){ return x->market_w_single; }
 else if( strcmp(name,"cons_m_single") == 0 ){ return x->cons_m_single; }
 else if( strcmp(name,"leisure_m_single") == 0 ){ return x->leisure_m_single; }
 else if( strcmp(name,"hours_m_single") == 0 ){ return x->hours_m_single; }
 else if( strcmp(name,"market_m_single") == 0 ){ return x->market_m_single; }
 else if( strcmp(name,"pre_cons_w_single") == 0 ){ return x->pre_cons_w_single; }
 else if( strcmp(name,"pre_market_w_single") == 0 ){ return x->pre_market_w_single; }
 else if( strcmp(name,"pre_cons_m_single") == 0 ){ return x->pre_cons_m_single; }
 else if( strcmp(name,"pre_market_m_single") == 0 ){ return x->pre_market_m_single; }
 else if( strcmp(name,"Vw_remain") == 0 ){ return x->Vw_remain; }
 else if( strcmp(name,"Vm_remain") == 0 ){ return x->Vm_remain; }
 else if( strcmp(name,"cons_w_remain") == 0 ){ return x->cons_w_remain; }
 else if( strcmp(name,"cons_m_remain") == 0 ){ return x->cons_m_remain; }
 else if( strcmp(name,"market_remain") == 0 ){ return x->market_remain; }
 else if( strcmp(name,"leisure_w_remain") == 0 ){ return x->leisure_w_remain; }
 else if( strcmp(name,"leisure_m_remain") == 0 ){ return x->leisure_m_remain; }
 else if( strcmp(name,"hours_w_remain") == 0 ){ return x->hours_w_remain; }
 else if( strcmp(name,"hours_m_remain") == 0 ){ return x->hours_m_remain; }
 else if( strcmp(name,"labor_w_remain") == 0 ){ return x->labor_w_remain; }
 else if( strcmp(name,"labor_m_remain") == 0 ){ return x->labor_m_remain; }
 else if( strcmp(name,"power") == 0 ){ return x->power; }
 else if( strcmp(name,"pre_cons_w") == 0 ){ return x->pre_cons_w; }
 else if( strcmp(name,"pre_cons_m") == 0 ){ return x->pre_cons_m; }
 else if( strcmp(name,"pre_market") == 0 ){ return x->pre_market; }
 else {return NULL;}

}


int* get_int_p_sol_struct(sol_struct* x, char* name){

 if( strcmp(name,"power_idx") == 0 ){ return x->power_idx; }
 else {return NULL;}

}


