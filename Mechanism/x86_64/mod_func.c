#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _KdBG_reg(void);
extern void _Kir_reg(void);
extern void _cad_reg(void);
extern void _cadL_reg(void);
extern void _cadN_reg(void);
extern void _cal_reg(void);
extern void _can_reg(void);
extern void _car_reg(void);
extern void _cat_reg(void);
extern void _dists_reg(void);
extern void _eff_reg(void);
extern void _exc_reg(void);
extern void _h_reg(void);
extern void _hha2_reg(void);
extern void _ic_reg(void);
extern void _id_reg(void);
extern void _inh_reg(void);
extern void _kad_reg(void);
extern void _kap_reg(void);
extern void _kca_reg(void);
extern void _kdr_reg(void);
extern void _km_reg(void);
extern void _na3_reg(void);
extern void _nmdaSyn_reg(void);
extern void _syns_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," KdBG.mod");
    fprintf(stderr," Kir.mod");
    fprintf(stderr," cad.mod");
    fprintf(stderr," cadL.mod");
    fprintf(stderr," cadN.mod");
    fprintf(stderr," cal.mod");
    fprintf(stderr," can.mod");
    fprintf(stderr," car.mod");
    fprintf(stderr," cat.mod");
    fprintf(stderr," dists.mod");
    fprintf(stderr," eff.mod");
    fprintf(stderr," exc.mod");
    fprintf(stderr," h.mod");
    fprintf(stderr," hha2.mod");
    fprintf(stderr," ic.mod");
    fprintf(stderr," id.mod");
    fprintf(stderr," inh.mod");
    fprintf(stderr," kad.mod");
    fprintf(stderr," kap.mod");
    fprintf(stderr," kca.mod");
    fprintf(stderr," kdr.mod");
    fprintf(stderr," km.mod");
    fprintf(stderr," na3.mod");
    fprintf(stderr," nmdaSyn.mod");
    fprintf(stderr," syns.mod");
    fprintf(stderr, "\n");
  }
  _KdBG_reg();
  _Kir_reg();
  _cad_reg();
  _cadL_reg();
  _cadN_reg();
  _cal_reg();
  _can_reg();
  _car_reg();
  _cat_reg();
  _dists_reg();
  _eff_reg();
  _exc_reg();
  _h_reg();
  _hha2_reg();
  _ic_reg();
  _id_reg();
  _inh_reg();
  _kad_reg();
  _kap_reg();
  _kca_reg();
  _kdr_reg();
  _km_reg();
  _na3_reg();
  _nmdaSyn_reg();
  _syns_reg();
}
