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
        rho[l] = 1.0 + 1e-05*(2.2733602246716968e-05*cos(60*Ymed(j) + 60*Zmed(k)) - (3.167583397097529e-05)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.33362329689650616 + 1e-05*(0.7973654573327341*cos(60*Ymed(j) + 60*Zmed(k)) - (0.6762546707509746)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.016412768082719808 + 1e-05*(0.391109550601909*cos(60*Ymin(j) + 60*Zmed(k)) - (0.33281392786638453)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(0.5983087535871898*cos(60*Ymed(j) + 60*Zmin(k)) - (0.18673418560371335)*sin(60*Ymed(j) + 60*Zmin(k)));
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
        rho[l] = 0.02711010553657533 + 1e-05*(1.0*cos(60*Ymed(j) + 60*Zmed(k)) - (0.0)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.3336318884354722 + 1e-05*(0.248245714629571*cos(60*Ymed(j) + 60*Zmed(k)) - (0.9488811518333182)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.01568161604812898 + 1e-05*(0.6672374531003724*cos(60*Ymin(j) + 60*Zmed(k)) - (0.09589793559411208)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(0.4418396661678128*cos(60*Ymed(j) + 60*Zmin(k)) - (0.8864799193275177)*sin(60*Ymed(j) + 60*Zmin(k)));
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
        rho[l] = 0.07190145629123526 + 1e-05*(0.6974534998820221*cos(60*Ymed(j) + 60*Zmed(k)) - (0.3264728640701121)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.33363555233587905 + 1e-05*(0.7339281633300665*cos(60*Ymed(j) + 60*Zmed(k)) - (0.22013495554548623)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.015347065934369976 + 1e-05*(0.08159456954220812*cos(60*Ymin(j) + 60*Zmed(k)) - (0.15989560107504752)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(0.3401001849547053*cos(60*Ymed(j) + 60*Zmin(k)) - (0.46519315370205094)*sin(60*Ymed(j) + 60*Zmin(k)));
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
        rho[l] = 0.13858961669912775 + 1e-05*(0.26642102829077097*cos(60*Ymed(j) + 60*Zmed(k)) - (0.815776403424807)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.33364479974663425 + 1e-05*(0.1932943892894945*cos(60*Ymed(j) + 60*Zmed(k)) - (0.12946907617720027)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.014423093985888579 + 1e-05*(0.09166475154493592*cos(60*Ymin(j) + 60*Zmed(k)) - (0.5985680136649132)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(0.8547419043740013*cos(60*Ymed(j) + 60*Zmin(k)) - (0.6016212416937131)*sin(60*Ymed(j) + 60*Zmin(k)));
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
        rho[l] = 0.23756584144318396 + 1e-05*(0.9319883611359835*cos(60*Ymed(j) + 60*Zmed(k)) - (0.7247813610920201)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.333662752489953 + 1e-05*(0.8605513173932924*cos(60*Ymed(j) + 60*Zmed(k)) - (0.9293378015753163)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.012038523312311686 + 1e-05*(0.546186009082353*cos(60*Ymin(j) + 60*Zmed(k)) - (0.9376729587677569)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(0.4949879400788243*cos(60*Ymed(j) + 60*Zmin(k)) - (0.2737731824899875)*sin(60*Ymed(j) + 60*Zmin(k)));
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
        rho[l] = 0.3624254092490066 + 1e-05*(0.4517787074747607*cos(60*Ymed(j) + 60*Zmed(k)) - (0.6650389233995303)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.3336708335132149 + 1e-05*(0.33089093046705464*cos(60*Ymed(j) + 60*Zmed(k)) - (0.9034540068082391)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = 0.00623194514488969 + 1e-05*(0.2570741752765343*cos(60*Ymin(j) + 60*Zmed(k)) - (0.33982833761031983)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(0.25885339864292733*cos(60*Ymed(j) + 60*Zmin(k)) - (0.355446479944286)*sin(60*Ymed(j) + 60*Zmin(k)));
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
        rho[l] = 0.4647967723375903 + 1e-05*(0.005022333717131788*cos(60*Ymed(j) + 60*Zmed(k)) - (0.6286045440996787)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.33352338155378347 + 1e-05*(0.2823827074251183*cos(60*Ymed(j) + 60*Zmed(k)) - (0.06808768948794575)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = -0.005958401552861245 + 1e-05*(0.6168289772563805*cos(60*Ymin(j) + 60*Zmed(k)) - (0.17632632028120343)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(0.3043883872195896*cos(60*Ymed(j) + 60*Zmin(k)) - (0.44088681087611803)*sin(60*Ymed(j) + 60*Zmin(k)));
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
        rho[l] = 0.45019862483388096 + 1e-05*(0.1502023410627008*cos(60*Ymed(j) + 60*Zmed(k)) - (0.21792886308543502)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.33283228213151184 + 1e-05*(0.4743331153335445*cos(60*Ymed(j) + 60*Zmed(k)) - (0.47636885508119187)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = -0.025266675824230492 + 1e-05*(0.25523235381950027*cos(60*Ymin(j) + 60*Zmed(k)) - (0.29756526814804807)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(0.27906711981376664*cos(60*Ymed(j) + 60*Zmin(k)) - (0.26057921249129756)*sin(60*Ymed(j) + 60*Zmin(k)));
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
        rho[l] = 0.24741217360939965 + 1e-05*(0.48276159279931574*cos(60*Ymed(j) + 60*Zmed(k)) - (0.2119790363515106)*sin(60*Ymed(j) + 60*Zmed(k)));
        vx[l]  = -0.33161032889397407 + 1e-05*(0.49563059667304066*cos(60*Ymed(j) + 60*Zmed(k)) - (0.24626132583073757)*sin(60*Ymed(j) + 60*Zmed(k))) - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);
        vy[l]  = -0.044114058087172156 + 1e-05*(0.8384826524669448*cos(60*Ymin(j) + 60*Zmed(k)) - (0.18013059009503507)*sin(60*Ymin(j) + 60*Zmed(k)));
        vz[l]  = 0.0 + 1e-05*(0.8621562915092364*cos(60*Ymed(j) + 60*Zmin(k)) - (0.17829944484518745)*sin(60*Ymed(j) + 60*Zmin(k)));
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
