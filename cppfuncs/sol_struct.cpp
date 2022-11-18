typedef struct sol_struct
{
 double* V_w_single;
 double* V_m_single;
 double* c_w_single;
 double* l_w_single;
 double* h_w_single;
 double* m_w_single;
 double* c_m_single;
 double* l_m_single;
 double* h_m_single;
 double* m_m_single;
 double* pre_cons_w_single;
 double* pre_market_w_single;
 double* pre_cons_m_single;
 double* pre_market_m_single;
 double* pre_cons_w;
 double* pre_cons_m;
 double* pre_market;
} sol_struct;

double* get_double_p_sol_struct(sol_struct* x, char* name){

 if( strcmp(name,"V_w_single") == 0 ){ return x->V_w_single; }
 else if( strcmp(name,"V_m_single") == 0 ){ return x->V_m_single; }
 else if( strcmp(name,"c_w_single") == 0 ){ return x->c_w_single; }
 else if( strcmp(name,"l_w_single") == 0 ){ return x->l_w_single; }
 else if( strcmp(name,"h_w_single") == 0 ){ return x->h_w_single; }
 else if( strcmp(name,"m_w_single") == 0 ){ return x->m_w_single; }
 else if( strcmp(name,"c_m_single") == 0 ){ return x->c_m_single; }
 else if( strcmp(name,"l_m_single") == 0 ){ return x->l_m_single; }
 else if( strcmp(name,"h_m_single") == 0 ){ return x->h_m_single; }
 else if( strcmp(name,"m_m_single") == 0 ){ return x->m_m_single; }
 else if( strcmp(name,"pre_cons_w_single") == 0 ){ return x->pre_cons_w_single; }
 else if( strcmp(name,"pre_market_w_single") == 0 ){ return x->pre_market_w_single; }
 else if( strcmp(name,"pre_cons_m_single") == 0 ){ return x->pre_cons_m_single; }
 else if( strcmp(name,"pre_market_m_single") == 0 ){ return x->pre_market_m_single; }
 else if( strcmp(name,"pre_cons_w") == 0 ){ return x->pre_cons_w; }
 else if( strcmp(name,"pre_cons_m") == 0 ){ return x->pre_cons_m; }
 else if( strcmp(name,"pre_market") == 0 ){ return x->pre_market; }
 else {return NULL;}

}


