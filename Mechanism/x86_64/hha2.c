/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__hha2
#define _nrn_initial _nrn_initial__hha2
#define nrn_cur _nrn_cur__hha2
#define _nrn_current _nrn_current__hha2
#define nrn_jacob _nrn_jacob__hha2
#define nrn_state _nrn_state__hha2
#define _net_receive _net_receive__hha2 
#define rates rates__hha2 
#define states states__hha2 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define ar2 _p[0]
#define gnabar _p[1]
#define gkbar _p[2]
#define gl _p[3]
#define el _p[4]
#define ina _p[5]
#define ik _p[6]
#define il _p[7]
#define m _p[8]
#define h _p[9]
#define n _p[10]
#define s _p[11]
#define ena _p[12]
#define ek _p[13]
#define Dm _p[14]
#define Dh _p[15]
#define Dn _p[16]
#define Ds _p[17]
#define _g _p[18]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define _ion_ek	*_ppvar[3]._pval
#define _ion_ik	*_ppvar[4]._pval
#define _ion_dikdv	*_ppvar[5]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_alpv(void);
 static void _hoc_alpr(void);
 static void _hoc_betr(void);
 static void _hoc_rates(void);
 static void _hoc_varss(void);
 static void _hoc_vartau(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_hha2", _hoc_setdata,
 "alpv_hha2", _hoc_alpv,
 "alpr_hha2", _hoc_alpr,
 "betr_hha2", _hoc_betr,
 "rates_hha2", _hoc_rates,
 "varss_hha2", _hoc_varss,
 "vartau_hha2", _hoc_vartau,
 0, 0
};
#define alpv alpv_hha2
#define alpr alpr_hha2
#define betr betr_hha2
#define varss varss_hha2
#define vartau vartau_hha2
 extern double alpv( double );
 extern double alpr( double );
 extern double betr( double );
 extern double varss( double , double );
 extern double vartau( double , double );
 /* declare global and static user variables */
#define a0r a0r_hha2
 double a0r = 0.0003;
#define b0r b0r_hha2
 double b0r = 0.0003;
#define gmr gmr_hha2
 double gmr = 0.2;
#define inf inf_hha2
 double inf[4];
#define taumin taumin_hha2
 double taumin = 10;
#define tau tau_hha2
 double tau[4];
#define vvh vvh_hha2
 double vvh = -58;
#define vhalfr vhalfr_hha2
 double vhalfr = -60;
#define vvs vvs_hha2
 double vvs = 2;
#define zetas zetas_hha2
 double zetas = 12;
#define zetar zetar_hha2
 double zetar = 12;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "a0r_hha2", "/ms",
 "b0r_hha2", "/ms",
 "taumin_hha2", "ms",
 "vvs_hha2", "mV",
 "vhalfr_hha2", "mV",
 "vvh_hha2", "mV",
 "tau_hha2", "ms",
 "gnabar_hha2", "mho/cm2",
 "gkbar_hha2", "mho/cm2",
 "gl_hha2", "mho/cm2",
 "el_hha2", "mV",
 "ina_hha2", "mA/cm2",
 "ik_hha2", "mA/cm2",
 "il_hha2", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double n0 = 0;
 static double s0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "a0r_hha2", &a0r_hha2,
 "b0r_hha2", &b0r_hha2,
 "zetar_hha2", &zetar_hha2,
 "zetas_hha2", &zetas_hha2,
 "gmr_hha2", &gmr_hha2,
 "taumin_hha2", &taumin_hha2,
 "vvs_hha2", &vvs_hha2,
 "vhalfr_hha2", &vhalfr_hha2,
 "vvh_hha2", &vvh_hha2,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 "inf_hha2", inf_hha2, 4,
 "tau_hha2", tau_hha2, 4,
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[6]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"hha2",
 "ar2_hha2",
 "gnabar_hha2",
 "gkbar_hha2",
 "gl_hha2",
 "el_hha2",
 0,
 "ina_hha2",
 "ik_hha2",
 "il_hha2",
 0,
 "m_hha2",
 "h_hha2",
 "n_hha2",
 "s_hha2",
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 19, _prop);
 	/*initialize range parameters*/
 	ar2 = 1;
 	gnabar = 0;
 	gkbar = 0;
 	gl = 0;
 	el = -70;
 	_prop->param = _p;
 	_prop->param_size = 19;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _hha2_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	ion_reg("k", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 19, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 hha2 /Users/danielchapman/PythonDev/Final/Mechanism/x86_64/hha2.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "HH channel that includes both a sodium and a delayed rectifier channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double, double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[4], _dlist1[4];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v , ar2 ) ;
   Dm = ( inf [ 0 ] - m ) / tau [ 0 ] ;
   Dh = ( inf [ 1 ] - h ) / tau [ 1 ] ;
   Dn = ( inf [ 2 ] - n ) / tau [ 2 ] ;
   Ds = ( inf [ 3 ] - s ) / tau [ 3 ] ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v , ar2 ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau[0] )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau[1] )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau[2] )) ;
 Ds = Ds  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau[3] )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargscomma_ v , ar2 ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau[0])))*(- ( ( ( inf[0] ) ) / tau[0] ) / ( ( ( ( - 1.0 ) ) ) / tau[0] ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau[1])))*(- ( ( ( inf[1] ) ) / tau[1] ) / ( ( ( ( - 1.0 ) ) ) / tau[1] ) - h) ;
    n = n + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau[2])))*(- ( ( ( inf[2] ) ) / tau[2] ) / ( ( ( ( - 1.0 ) ) ) / tau[2] ) - n) ;
    s = s + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau[3])))*(- ( ( ( inf[3] ) ) / tau[3] ) / ( ( ( ( - 1.0 ) ) ) / tau[3] ) - s) ;
   }
  return 0;
}
 
static int  rates (  double _lv , double _la2 ) {
   double _ltmp , _lc ;
 {int  _li ;for ( _li = 0 ; _li <= 2 ; _li ++ ) {
     tau [ _li ] = vartau ( _threadargscomma_ _lv , ((double) _li ) ) ;
     inf [ _li ] = varss ( _threadargscomma_ _lv , ((double) _li ) ) ;
     } }
   tau [ 3 ] = betr ( _threadargscomma_ _lv ) / ( a0r * ( 1.0 + alpr ( _threadargscomma_ _lv ) ) ) ;
   if ( tau [ 3 ] < taumin ) {
     tau [ 3 ] = taumin ;
     }
   _lc = alpv ( _threadargscomma_ _lv ) ;
   inf [ 3 ] = _lc + _la2 * ( 1.0 - _lc ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double varss (  double _lv , double _li ) {
   double _lvarss;
 if ( _li  == 0.0 ) {
     _lvarss = 1.0 / ( 1.0 + exp ( ( _lv + 44.0 ) / ( - 3.0 ) ) ) ;
     }
   else if ( _li  == 1.0 ) {
     _lvarss = 1.0 / ( 1.0 + exp ( ( _lv + 49.0 ) / ( 3.5 ) ) ) ;
     }
   else if ( _li  == 2.0 ) {
     _lvarss = 1.0 / ( 1.0 + exp ( ( _lv + 46.3 ) / ( - 3.0 ) ) ) ;
     }
   
return _lvarss;
 }
 
static void _hoc_varss(void) {
  double _r;
   _r =  varss (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double alpv (  double _lv ) {
   double _lalpv;
 _lalpv = 1.0 / ( 1.0 + exp ( ( _lv - vvh ) / vvs ) ) ;
   
return _lalpv;
 }
 
static void _hoc_alpv(void) {
  double _r;
   _r =  alpv (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alpr (  double _lv ) {
   double _lalpr;
  _lalpr = exp ( 1.e-3 * zetar * ( _lv - vhalfr ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
    
return _lalpr;
 }
 
static void _hoc_alpr(void) {
  double _r;
   _r =  alpr (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double betr (  double _lv ) {
   double _lbetr;
  _lbetr = exp ( 1.e-3 * zetar * gmr * ( _lv - vhalfr ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
    
return _lbetr;
 }
 
static void _hoc_betr(void) {
  double _r;
   _r =  betr (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double vartau (  double _lv , double _li ) {
   double _lvartau;
 double _ltmp ;
 if ( _li  == 0.0 ) {
     _lvartau = 0.05 ;
     }
   else if ( _li  == 1.0 ) {
     _lvartau = 1.0 ;
     }
   else if ( _li  == 2.0 ) {
     _lvartau = 3.5 ;
     }
   
return _lvartau;
 }
 
static void _hoc_vartau(void) {
  double _r;
   _r =  vartau (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 4;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
     _ode_spec1 ();
   }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 4; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
  n = n0;
  s = s0;
 {
   rates ( _threadargscomma_ v , ar2 ) ;
   m = inf [ 0 ] ;
   h = inf [ 1 ] ;
   n = inf [ 2 ] ;
   s = inf [ 3 ] ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
  ek = _ion_ek;
 initmodel();
  }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ina = gnabar * m * m * h * s * ( v - ena ) ;
   ik = gkbar * n * n * ( v - ek ) ;
   il = gl * ( v - el ) ;
   }
 _current += ina;
 _current += ik;
 _current += il;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
  ek = _ion_ek;
 { error =  states();
 if(error){fprintf(stderr,"at line 79 in file hha2.mod:\n	SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 }  }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(n) - _p;  _dlist1[2] = &(Dn) - _p;
 _slist1[3] = &(s) - _p;  _dlist1[3] = &(Ds) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/danielchapman/PythonDev/Final/Mechanism/hha2.mod";
static const char* nmodl_file_text = 
  "TITLE HH channel that includes both a sodium and a delayed rectifier channel\n"
  ": and accounts for sodium conductance attenuation\n"
  ": Bartlett Mel-modified Hodgkin - Huxley conductances (after Ojvind et al.)\n"
  ": Terrence Brannon-added attenuation \n"
  ": Yiota Poirazi-modified Kdr and Na threshold and time constants\n"
  ": to make it more stable, 2000, poirazi@LNC.usc.edu\n"
  ": Used in all BUT somatic and axon sections. The spike threshold is about -50 mV\n"
  ":\n"
  ": Modified to use CVode --Carl Gold 08/12/03\n"
  ":  Updated by Maria Markaki  12/05/03\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX hha2\n"
  "	USEION na READ ena WRITE ina \n"
  "	USEION k READ ek WRITE ik\n"
  "	NONSPECIFIC_CURRENT il\n"
  "	RANGE gnabar, gkbar, gl, el,ik,ina\n"
  "	RANGE ar2, vhalfs\n"
  "	GLOBAL inf, tau, taumin\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "}\n"
  "\n"
  "\n"
  "PARAMETER {   : parameters that can be entered when function is called in cell-setup\n"
  "        a0r = 0.0003 (/ms)\n"
  "        b0r = 0.0003 (/ms)\n"
  "       : a0r = 0.0003 (ms)\n"
  "        :b0r = 0.0003 (ms)\n"
  "        zetar = 12    \n"
  "	zetas = 12   \n"
  "        gmr = 0.2   \n"
  "	ar2 = 1.0               :initialized parameter for location-dependent\n"
  "                                :Na-conductance attenuation, \"s\", (ar=1 -> zero attenuation)\n"
  ":	taumin = 3   (ms)       :min activation time for \"s\" attenuation system\n"
  "	taumin = 10   (ms)       : smax=10		(ms)\n"
  "        vvs  = 2     (mV)       :slope for \"s\" attenuation system\n"
  "        vhalfr = -60 (mV)       :half potential for \"s\" attenuation system\n"
  "        vvh=-58		(mV) \n"
  "        gnabar = 0   (mho/cm2)  :initialized conductances\n"
  "	gkbar = 0    (mho/cm2)  :actual values set in hoc code\n"
  "	gl = 0       (mho/cm2)  : keep this to zero if using the NEURON pas mechanism\n"
  "	el = -70.0   (mV)       :steady state \n"
  "\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "	celsius      (degC)\n"
  "	v            (mV)\n"
  "	ena          (mV)       :Na reversal potential (also reset in\n"
  "	ek           (mV)       :K reversal potential  cell-setup.hoc)\n"
  "	ina          (mA/cm2)\n"
  "	ik           (mA/cm2)\n"
  "	il           (mA/cm2)\n"
  "	inf[4]\n"
  "	tau[4]	     (ms)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	m \n"
  "	h \n"
  "	n \n"
  "	s\n"
  "}\n"
  "\n"
  "INITIAL {                       : initialize the following parameter using states()\n"
  "	rates(v,ar2)\n"
  "	m = inf[0]\n"
  "	h = inf[1]\n"
  "	n = inf[2]\n"
  "	s = inf[3]\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "	ina = gnabar*m*m*h*s*(v - ena) :Sodium current\n"
  "	ik = gkbar*n*n*(v - ek)        :Potassium current\n"
  "	il = gl*(v - el)               :leak current\n"
  "}\n"
  "\n"
  "\n"
  "DERIVATIVE states {\n"
  "	rates(v,ar2)\n"
  "	m' = (inf[0]-m)/tau[0]\n"
  "	h' = (inf[1]-h)/tau[1]\n"
  "	n' = (inf[2]-n)/tau[2]\n"
  "	s' = (inf[3]-s)/tau[3]\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v(mV),a2) {\n"
  "	LOCAL tmp, c\n"
  "	FROM i=0 TO 2 {\n"
  "		tau[i] = vartau(v,i)\n"
  "		inf[i] = varss(v,i)\n"
  "	}\n"
  "	tau[3] = betr(v)/(a0r*(1+alpr(v))) \n"
  "	if (tau[3]<taumin) {tau[3]=taumin} :s activation tau\n"
  "	c = alpv(v)\n"
  "	inf[3] = c+a2*(1-c) \n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION varss(v(mV), i) { :steady state values\n"
  "	if (i==0) {\n"
  "		varss = 1 / (1 + exp((v + 44)/(-3(mV)))) :Na activation\n"
  "	}\n"
  "	else if (i==1) {\n"
  "		varss = 1 / (1 + exp((v + 49)/(3.5(mV))))  :Na inactivation\n"
  "	}\n"
  "	else if (i==2) {	\n"
  "		varss = 1 / (1 + exp((v + 46.3)/(-3(mV)))) :K activation\n"
  "\n"
  "	}\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION alpv(v(mV)) {\n"
  "         alpv = 1/(1+exp((v-vvh)/vvs))\n"
  "}\n"
  "\n"
  "FUNCTION alpr(v(mV)) {       :used in \"s\" activation system tau\n"
  "UNITSOFF\n"
  "  alpr = exp(1.e-3*zetar*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) \n"
  "UNITSON\n"
  "}\n"
  "\n"
  "FUNCTION betr(v(mV)) {       :used in \"s\" activation system tau\n"
  "UNITSOFF\n"
  "  betr = exp(1.e-3*zetar*gmr*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) \n"
  "UNITSON\n"
  "}\n"
  "\n"
  "FUNCTION vartau(v(mV), i) (ms){ :estimate tau values\n"
  "	LOCAL tmp\n"
  "	if (i==0) {\n"
  "	   vartau = 0.05(ms)      :Na activation tau\n"
  "	}\n"
  "	else if (i==1) {\n"
  "           vartau = 1(ms)       :Na inactivation tau\n"
  "        }\n"
  "	else if (i==2) {\n"
  "            vartau = 3.5(ms)      :K activation tau\n"
  "       	}: else {\n"
  "	  :   tmp =1* betr(v)/(a0r+b0r*alpr(v)) \n"
  "	   :  if (tmp<taumin) {tmp=taumin}\n"
  "	 :    vartau = tmp      :s activation tau\n"
  "       :}\n"
  "}	\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  ;
#endif
