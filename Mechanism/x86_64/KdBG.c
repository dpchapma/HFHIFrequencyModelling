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
 
#define nrn_init _nrn_init__kdBG
#define _nrn_initial _nrn_initial__kdBG
#define nrn_cur _nrn_cur__kdBG
#define _nrn_current _nrn_current__kdBG
#define nrn_jacob _nrn_jacob__kdBG
#define nrn_state _nrn_state__kdBG
#define _net_receive _net_receive__kdBG 
#define rates rates__kdBG 
#define states states__kdBG 
 
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
#define gbar _p[0]
#define ik _p[1]
#define xs _p[2]
#define ys _p[3]
#define q10 _p[4]
#define T _p[5]
#define Dxs _p[6]
#define Dys _p[7]
#define _g _p[8]
#define _ion_ik	*_ppvar[0]._pval
#define _ion_dikdv	*_ppvar[1]._pval
 
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
 static void _hoc_rates(void);
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
 "setdata_kdBG", _hoc_setdata,
 "rates_kdBG", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
#define Ky Ky_kdBG
 double Ky = 0.0002;
#define gammay gammay_kdBG
 double gammay = 0;
#define taoy taoy_kdBG
 double taoy = 0;
#define taox taox_kdBG
 double taox = 1;
#define vhalfy vhalfy_kdBG
 double vhalfy = -73;
#define vhalfx vhalfx_kdBG
 double vhalfx = -63;
#define xinf xinf_kdBG
 double xinf = 0;
#define xtau xtau_kdBG
 double xtau = 0;
#define yinf yinf_kdBG
 double yinf = 0;
#define ytau ytau_kdBG
 double ytau = 0;
#define zettay zettay_kdBG
 double zettay = -2.5;
#define zettax zettax_kdBG
 double zettax = 3;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Ky_kdBG", "1/ms",
 "gammay_kdBG", "1",
 "zettax_kdBG", "1",
 "zettay_kdBG", "1",
 "vhalfx_kdBG", "mV",
 "vhalfy_kdBG", "mV",
 "taox_kdBG", "ms",
 "taoy_kdBG", "ms",
 "xtau_kdBG", "ms",
 "ytau_kdBG", "ms",
 "xinf_kdBG", "1",
 "yinf_kdBG", "1",
 "gbar_kdBG", "S/cm2",
 "ik_kdBG", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double v = 0;
 static double xs0 = 0;
 static double ys0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Ky_kdBG", &Ky_kdBG,
 "gammay_kdBG", &gammay_kdBG,
 "zettax_kdBG", &zettax_kdBG,
 "zettay_kdBG", &zettay_kdBG,
 "vhalfx_kdBG", &vhalfx_kdBG,
 "vhalfy_kdBG", &vhalfy_kdBG,
 "taox_kdBG", &taox_kdBG,
 "taoy_kdBG", &taoy_kdBG,
 "xtau_kdBG", &xtau_kdBG,
 "ytau_kdBG", &ytau_kdBG,
 "xinf_kdBG", &xinf_kdBG,
 "yinf_kdBG", &yinf_kdBG,
 0,0
};
 static DoubVec hoc_vdoub[] = {
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
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"kdBG",
 "gbar_kdBG",
 0,
 "ik_kdBG",
 0,
 "xs_kdBG",
 "ys_kdBG",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 9, _prop);
 	/*initialize range parameters*/
 	gbar = 0.001;
 	_prop->param = _p;
 	_prop->param_size = 9;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[1]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
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

 void _KdBG_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 9, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 kdBG /Users/danielchapman/PythonDev/Final/Mechanism/x86_64/KdBG.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96485.3;
 static double R = 8.3145;
static int _reset;
static char *modelname = "Kd current";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates();
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargs_ ) ;
   Dxs = ( xinf - xs ) / xtau ;
   Dys = ( yinf - ys ) / ytau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargs_ ) ;
 Dxs = Dxs  / (1. - dt*( ( ( ( - 1.0 ) ) ) / xtau )) ;
 Dys = Dys  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ytau )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargs_ ) ;
    xs = xs + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / xtau)))*(- ( ( ( xinf ) ) / xtau ) / ( ( ( ( - 1.0 ) ) ) / xtau ) - xs) ;
    ys = ys + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ytau)))*(- ( ( ( yinf ) ) / ytau ) / ( ( ( ( - 1.0 ) ) ) / ytau ) - ys) ;
   }
  return 0;
}
 
static int  rates (  ) {
   double _la , _lb ;
 _la = q10 * exp ( ( 1.0e-3 ) * zettax * ( v - vhalfx ) * FARADAY / ( R * T ) ) ;
   _lb = q10 * exp ( ( 1.0e-3 ) * - zettax * ( v - vhalfx ) * FARADAY / ( R * T ) ) ;
   xinf = _la / ( _la + _lb ) ;
   xtau = taox ;
   _la = q10 * Ky * exp ( ( 1.0e-3 ) * zettay * gammay * ( v - vhalfy ) * FARADAY / ( R * T ) ) ;
   _lb = q10 * Ky * exp ( ( 1.0e-3 ) * - zettay * ( 1.0 - gammay ) * ( v - vhalfy ) * FARADAY / ( R * T ) ) ;
   yinf = _la / ( _la + _lb ) ;
   ytau = 1.0 / ( _la + _lb ) + taoy ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
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
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  xs = xs0;
  ys = ys0;
 {
   T = celsius + 273.15 ;
   q10 = pow( 1.0 , ( ( celsius - 35.0 ) / 10.0 ) ) ;
   rates ( _threadargs_ ) ;
   xs = xinf ;
   ys = yinf ;
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
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ik = gbar * pow( xs , 4.0 ) * pow( ys , 4.0 ) * ( v + 95.0 ) ;
   }
 _current += ik;

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
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
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
 { error =  states();
 if(error){fprintf(stderr,"at line 56 in file KdBG.mod:\n	SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(xs) - _p;  _dlist1[0] = &(Dxs) - _p;
 _slist1[1] = &(ys) - _p;  _dlist1[1] = &(Dys) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/danielchapman/PythonDev/Final/Mechanism/KdBG.mod";
static const char* nmodl_file_text = 
  "TITLE Kd current\n"
  "\n"
  "COMMENT Equations from \n"
  "		  Lyle J Borg-Graham Interpretation of Data and Mechanisms for Hippocampal Pyramidal Cell Models A Chapter in \"Cerebral Cortex, Volumne 13: Cortical Models\" Edited by P.S.Ulinski, E.G.Jones and A.Peters,New York:plenum Press,1998\n"
  "		  \n"
  "		  The Krasnow Institute\n"
  "		  George Mason University\n"
  "\n"
  "Copyright	  Maciej Lazarewicz, 2001\n"
  "		  (mlazarew@gmu.edu)\n"
  "		  All rights reserved.\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX kdBG\n"
  "	USEION k WRITE ik\n"
  "	RANGE  gbar,ik\n"
  "	GLOBAL xtau, ytau, xinf, yinf\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(S)	= (siemens)\n"
  "	(mA)	= (milliamp)\n"
  "	(mV)	= (millivolt)\n"
  "	FARADAY	= (faraday) (coulombs)\n"
  "	R	= (k-mole)  (joule/degC)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gbar	=   1.0e-3	(S/cm2)\n"
  "	Ky	=   2.0e-4	(1/ms)\n"
  "	gammay	=   0.0		(1)\n"
  "	zettax	=   3.0		(1)\n"
  "	zettay	=  -2.5		(1)\n"
  "	vhalfx	= -63.0		(mV)\n"
  "	vhalfy	= -73.0		(mV)\n"
  "	taox	=   1.0		(ms)\n"
  "	taoy	=   0.0		(ms)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v       (mV)\n"
  "	ik     	(mA/cm2)\n"
  "	celsius			(degC)\n"
  "	xtau    (ms)\n"
  "	ytau    (ms)\n"
  "	xinf	(1)\n"
  "	yinf	(1)\n"
  "	q10	(1)\n"
  "	T     	(K)\n"
  "}\n"
  "\n"
  "STATE { xs ys }\n"
  "\n"
  "BREAKPOINT { \n"
  "	SOLVE states METHOD cnexp\n"
  "	ik= gbar * xs^4 * ys^4 * ( v + 95.0 ) \n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "	rates()\n"
  "	xs'= (xinf- xs)/ xtau	\n"
  "	ys'= (yinf- ys)/ ytau\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	T  = celsius + 273.15\n"
  "	q10= 1.0^( (celsius-35.0) / 10.0(K) )\n"
  "	rates()\n"
  "	xs= xinf\n"
  "	ys= yinf\n"
  "}\n"
  "\n"
  "PROCEDURE rates() { LOCAL a, b  \n"
  "	a = q10*exp( (1.0e-3)*  zettax*(v-vhalfx)*FARADAY/(R*T) )\n"
  "	b = q10*exp( (1.0e-3)* -zettax*(v-vhalfx)*FARADAY/(R*T) )\n"
  "	xinf = a / ( a + b )\n"
  "	xtau = taox\n"
  "\n"
  "	a = q10*Ky*exp( (1.0e-3)*  zettay*     gammay *(v-vhalfy)*FARADAY/(R*T) )\n"
  "	b = q10*Ky*exp( (1.0e-3)* -zettay*(1.0-gammay)*(v-vhalfy)*FARADAY/(R*T) )\n"
  "	yinf = a   / ( a + b )\n"
  "	ytau = 1.0 / ( a + b ) + taoy\n"
  "}\n"
  ;
#endif
