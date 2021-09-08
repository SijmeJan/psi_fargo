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
        rho[l] = 1.0 + 0.0*exp(-((Ymed(j)-0.05)*(Ymed(j)-0.05) + (Zmed(k)-0.05)*(Zmed(k)-0.05))/(0.025*0.025)) + 0.0001*(7.4637e-06*cos(30.0*Ymed(j) + 30.0*Zmed(k)) - 7.0677e-06*sin(30.0*Ymed(j) + 30.0*Zmed(k)));
        vx[l]  = -0.2504684572142411 + 0.0001*(0.0445570113*cos(30.0*Ymed(j) + 30.0*Zmed(k)) - 0.0197224299*sin(30.0*Ymed(j) + 30.0*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.03747657713928795 + 0.0001*(-0.0563787907*cos(30.0*Ymin(j) + 30.0*Zmed(k)) - 0.0120535455*sin(30.0*Ymin(j) + 30.0*Zmed(k)));
        vz[l]  = 0.0 + 0.0001*(0.0563784989*cos(30.0*Ymed(j) + 30.0*Zmin(k)) - -0.0120536242*sin(30.0*Ymed(j) + 30.0*Zmin(k)));
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
        rho[l] = 3.0 + 0.0001*(1.0*cos(30.0*Ymed(j) + 30.0*Zmed(k)) - 0.0*sin(30.0*Ymed(j) + 30.0*Zmed(k)));
        vx[l]  = -0.249843847595253 + 0.0001*(0.0435211557*cos(30.0*Ymed(j) + 30.0*Zmed(k)) - 0.0213517453*sin(30.0*Ymed(j) + 30.0*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = -0.012492192379762646 + 0.0001*(-0.0466198076*cos(30.0*Ymin(j) + 30.0*Zmed(k)) - 0.0124333223*sin(30.0*Ymin(j) + 30.0*Zmed(k)));
        vz[l]  = 0.0 + 0.0001*(0.0546507401*cos(30.0*Ymed(j) + 30.0*Zmin(k)) - -0.0077776652*sin(30.0*Ymed(j) + 30.0*Zmin(k)));
        cs[l]  = 0.0;
      }
    }
  }

  /*We now fill the collision matrix (Feedback from dust included)
  Note: ColRate() moves the collision matrix to the device.
  If feedback=NO, gas does not feel the drag force.*/
  ColRate(INVSTOKES1, id_gas, 1, feedback);
}
