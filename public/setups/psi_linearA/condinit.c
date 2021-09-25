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
        rho[l] = 1.0 + 1e-05*(-5.842314036904915e-05*cos(60*Ymed(j) + 60*Zmed(k)) - (-2.5669732914177673e-05)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.33362329689650616 + 1e-05*(-0.31548326185418685*cos(60*Ymed(j) + 60*Zmed(k)) - (-0.8962485443130421)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.016412768082719808 + 1e-05*(1.0*cos(60*Ymin(j) + 60*Zmed(k)) - (0.0)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(-0.999999486521875*cos(60*Ymed(j) + 60*Zmin(k)) - (1.121571940598809e-07)*sin(60*Ymed(j) + 60*Zmin(k)));
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
        rho[l] = 0.07386066065566742 + 1e-05*(0.0028207913290286967*cos(60*Ymed(j) + 60*Zmed(k)) - (0.02393727593324721)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.333633646991661 + 1e-05*(-0.31662312683662974*cos(60*Ymed(j) + 60*Zmed(k)) - (-0.895947891794636)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.015522952845211426 + 1e-05*(0.9990236908972593*cos(60*Ymin(j) + 60*Zmed(k)) - (-0.0029634057074862316)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(-0.9998689397661567*cos(60*Ymed(j) + 60*Zmin(k)) - (0.0005747605210755163)*sin(60*Ymed(j) + 60*Zmin(k)));
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
        rho[l] = 0.09849477397807484 + 1e-05*(0.006012766552709017*cos(60*Ymed(j) + 60*Zmed(k)) - (0.06262409625439008)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.3336408810992883 + 1e-05*(-0.3174206920859448*cos(60*Ymed(j) + 60*Zmed(k)) - (-0.8957417956299747)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.014830393657574532 + 1e-05*(0.9982584989627448*cos(60*Ymin(j) + 60*Zmed(k)) - (-0.005169125530298975)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(-0.9997670962239121*cos(60*Ymed(j) + 60*Zmin(k)) - (0.0009232966443552898)*sin(60*Ymed(j) + 60*Zmin(k)));
        cs[l]  = 0.0;
      }
    }
  }

  sprintf(dust_name,"dust%d",3); //We assign different names to the dust fluids
  Fluids[3]  = CreateFluid(dust_name, DUST);
  SelectFluid(3);
  rho  = Density->field_cpu;
  cs   = Energy->field_cpu;
  vx   = Vx->field_cpu;
  vy   = Vy->field_cpu;
  vz   = Vz->field_cpu;
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
        rho[l] = 0.13134489205584537 + 1e-05*(0.007805780264315092*cos(60*Ymed(j) + 60*Zmed(k)) - (0.18116467681228884)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.3336519696667888 + 1e-05*(-0.3186456250768081*cos(60*Ymed(j) + 60*Zmed(k)) - (-0.8954370002452439)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.013598770703308923 + 1e-05*(0.9968900359011932*cos(60*Ymin(j) + 60*Zmed(k)) - (-0.008875151466176748)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(-0.9995860528922722*cos(60*Ymed(j) + 60*Zmin(k)) - (0.0013298336824458096)*sin(60*Ymed(j) + 60*Zmin(k)));
        cs[l]  = 0.0;
      }
    }
  }

  sprintf(dust_name,"dust%d",4); //We assign different names to the dust fluids
  Fluids[4]  = CreateFluid(dust_name, DUST);
  SelectFluid(4);
  rho  = Density->field_cpu;
  cs   = Energy->field_cpu;
  vx   = Vx->field_cpu;
  vy   = Vy->field_cpu;
  vz   = Vz->field_cpu;
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
        rho[l] = 0.17515122856164805 + 1e-05*(-0.09687295992498843*cos(60*Ymed(j) + 60*Zmed(k)) - (0.6744581744579063)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.3336660726730801 + 1e-05*(-0.32021360738230725*cos(60*Ymed(j) + 60*Zmed(k)) - (-0.8950953317000779)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.011408482967754903 + 1e-05*(0.9944484678485207*cos(60*Ymin(j) + 60*Zmed(k)) - (-0.014787247437926587)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(-0.9992657629817074*cos(60*Ymed(j) + 60*Zmin(k)) - (0.0013791675779289182)*sin(60*Ymed(j) + 60*Zmin(k)));
        cs[l]  = 0.0;
      }
    }
  }

  sprintf(dust_name,"dust%d",5); //We assign different names to the dust fluids
  Fluids[5]  = CreateFluid(dust_name, DUST);
  SelectFluid(5);
  rho  = Density->field_cpu;
  cs   = Energy->field_cpu;
  vx   = Vx->field_cpu;
  vy   = Vy->field_cpu;
  vz   = Vz->field_cpu;
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
        rho[l] = 0.23356791715669456 + 1e-05*(-3.051201831063354*cos(60*Ymed(j) + 60*Zmed(k)) - (-0.3852356756704245)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.3336733943338212 + 1e-05*(-0.32107594500225456*cos(60*Ymed(j) + 60*Zmed(k)) - (-0.8951402028991791)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.007513555628983113 + 1e-05*(0.9901620402123752*cos(60*Ymin(j) + 60*Zmed(k)) - (-0.023189861528739086)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(-0.9986984121306737*cos(60*Ymed(j) + 60*Zmin(k)) - (-0.0006586067071176905)*sin(60*Ymed(j) + 60*Zmin(k)));
        cs[l]  = 0.0;
      }
    }
  }

  sprintf(dust_name,"dust%d",6); //We assign different names to the dust fluids
  Fluids[6]  = CreateFluid(dust_name, DUST);
  SelectFluid(6);
  rho  = Density->field_cpu;
  cs   = Energy->field_cpu;
  vx   = Vx->field_cpu;
  vy   = Vy->field_cpu;
  vz   = Vz->field_cpu;
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
        rho[l] = 0.31146782339420004 + 1e-05*(-0.8567136350564517*cos(60*Ymed(j) + 60*Zmed(k)) - (-1.445049562146815)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.33363028683222773 + 1e-05*(-0.3165126863189757*cos(60*Ymed(j) + 60*Zmed(k)) - (-0.8972852184254327)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.0005895262905937792 + 1e-05*(0.9830599032592415*cos(60*Ymin(j) + 60*Zmed(k)) - (-0.03164376011040021)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(-0.997568205530217*cos(60*Ymed(j) + 60*Zmin(k)) - (-0.01097377703527285)*sin(60*Ymed(j) + 60*Zmin(k)));
        cs[l]  = 0.0;
      }
    }
  }

  sprintf(dust_name,"dust%d",7); //We assign different names to the dust fluids
  Fluids[7]  = CreateFluid(dust_name, DUST);
  SelectFluid(7);
  rho  = Density->field_cpu;
  cs   = Energy->field_cpu;
  vx   = Vx->field_cpu;
  vy   = Vy->field_cpu;
  vz   = Vz->field_cpu;
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
        rho[l] = 0.41534901792542683 + 1e-05*(-0.5201226982690978*cos(60*Ymed(j) + 60*Zmed(k)) - (-1.4750196963580133)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.3333765206396656 + 1e-05*(-0.2888338966259539*cos(60*Ymed(j) + 60*Zmed(k)) - (-0.9065928460028239)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = -0.011703974533066923 + 1e-05*(0.9730007850209015*cos(60*Ymin(j) + 60*Zmed(k)) - (-0.02700781821012077)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(-0.9933642853468078*cos(60*Ymed(j) + 60*Zmin(k)) - (-0.05016903664565876)*sin(60*Ymed(j) + 60*Zmin(k)));
        cs[l]  = 0.0;
      }
    }
  }

  sprintf(dust_name,"dust%d",8); //We assign different names to the dust fluids
  Fluids[8]  = CreateFluid(dust_name, DUST);
  SelectFluid(8);
  rho  = Density->field_cpu;
  cs   = Energy->field_cpu;
  vx   = Vx->field_cpu;
  vy   = Vy->field_cpu;
  vz   = Vz->field_cpu;
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
        rho[l] = 0.5538768172315454 + 1e-05*(-0.07461030342107641*cos(60*Ymed(j) + 60*Zmed(k)) - (-1.7292458898754814)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.332369637062138 + 1e-05*(-0.17642916458512045*cos(60*Ymed(j) + 60*Zmed(k)) - (-0.9241899581148348)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = -0.03343564515544567 + 1e-05*(0.9587533631758136*cos(60*Ymin(j) + 60*Zmed(k)) - (0.04129398889985225)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(-0.959328378808431*cos(60*Ymed(j) + 60*Zmin(k)) - (-0.17900364978160715)*sin(60*Ymed(j) + 60*Zmin(k)));
        cs[l]  = 0.0;
      }
    }
  }

  /*We now fill the collision matrix (Feedback from dust included)
  Note: ColRate() moves the collision matrix to the device.
  If feedback=NO, gas does not feel the drag force.*/
  ColRate(INVSTOKES1, id_gas, 1, feedback);
  ColRate(INVSTOKES2, id_gas, 2, feedback);
  ColRate(INVSTOKES3, id_gas, 3, feedback);
  ColRate(INVSTOKES4, id_gas, 4, feedback);
  ColRate(INVSTOKES5, id_gas, 5, feedback);
  ColRate(INVSTOKES6, id_gas, 6, feedback);
  ColRate(INVSTOKES7, id_gas, 7, feedback);
  ColRate(INVSTOKES8, id_gas, 8, feedback);
}
