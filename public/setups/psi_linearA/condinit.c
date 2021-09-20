#include "fargo3d.h"

void CondInit() {
  int id_gas = 0;
  int feedback = YES;
  //We first create the gaseous fluid and store it in the array Fluids
  Fluids[id_gas] = CreateFluid("gas",GAS);

  //We now select the fluid
  SelectFluid(id_gas);
  //and fill its fields
  int i,j,k;
  real *rho  = Density->field_cpu;
  real *cs   = Energy->field_cpu;
  real *vx   = Vx->field_cpu;
  real *vy   = Vy->field_cpu;
  real *vz   = Vz->field_cpu;
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
        rho[l] = 1.0 + 0.0001*(2.2733602246716968e-05*cos(60*Ymed(j) + 60*Zmed(k)) - (3.167583397097529e-05)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.3334116724847603 + 0.0001*(0.7973654573327341*cos(60*Ymed(j) + 60*Zmed(k)) - (0.6762546707509746)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.01101308694128686 + 0.0001*(0.391109550601909*cos(60*Ymin(j) + 60*Zmed(k)) - (0.33281392786638453)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 0.0001*(0.5983087535871898*cos(60*Ymed(j) + 60*Zmin(k)) - (0.18673418560371335)*sin(60*Ymed(j) + 60*Zmin(k)));
        cs[l]  = 20.0;      // H*Omega*Omega/eta
      }
    }
  }
  
  //We repeat the process for the dust fluids
  char dust_name[MAXNAMELENGTH];
  sprintf(dust_name,"dust%d",1); //We assign different names to the dust fluids
  Fluids[1]  = CreateFluid(dust_name, DUST);
  SelectFluid(1);
  rho  = Density->field_cpu;
  cs   = Energy->field_cpu;
  vx   = Vx->field_cpu;
  vy   = Vy->field_cpu;
  vz   = Vz->field_cpu;
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
        rho[l] = 0.48050615 + 0.0001*(1.0*cos(60*Ymed(j) + 60*Zmed(k)) - (0.0)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.33342575145459785 + 0.0001*(0.248245714629571*cos(60*Ymed(j) + 60*Zmed(k)) - (0.9488811518333182)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.008904315770667167 + 0.0001*(0.6672374531003724*cos(60*Ymin(j) + 60*Zmed(k)) - (0.09589793559411208)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 0.0001*(0.4418396661678128*cos(60*Ymed(j) + 60*Zmin(k)) - (0.8864799193275177)*sin(60*Ymed(j) + 60*Zmin(k)));
        cs[l]  = 0.0;
      }
    }
  }

  sprintf(dust_name,"dust%d",2); //We assign different names to the dust fluids
  Fluids[2]  = CreateFluid(dust_name, DUST);
  SelectFluid(2);
  rho  = Density->field_cpu;
  cs   = Energy->field_cpu;
  vx   = Vx->field_cpu;
  vy   = Vy->field_cpu;
  vz   = Vz->field_cpu;
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
        rho[l] = 1.51949385 + 0.0001*(0.6974534998820221*cos(60*Ymed(j) + 60*Zmed(k)) - (0.3264728640701121)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.3332525520737937 + 0.0001*(0.7339281633300665*cos(60*Ymed(j) + 60*Zmed(k)) - (0.22013495554548623)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = -0.010063657336049383 + 0.0001*(0.08159456954220812*cos(60*Ymin(j) + 60*Zmed(k)) - (0.15989560107504752)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 0.0001*(0.3401001849547053*cos(60*Ymed(j) + 60*Zmin(k)) - (0.46519315370205094)*sin(60*Ymed(j) + 60*Zmin(k)));
        cs[l]  = 0.0;
      }
    }
  }

  /*We now fill the collision matrix (Feedback from dust included)
  Note: ColRate() moves the collision matrix to the device.
  If feedback=NO, gas does not feel the drag force.*/
  ColRate(INVSTOKES1, id_gas, 1, feedback);
  ColRate(INVSTOKES2, id_gas, 2, feedback);
}
