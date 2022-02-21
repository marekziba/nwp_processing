//
//  thermo.h
//  
//
//  Created by Marek Zięba on 27/12/2019.
//

#ifndef thermo_h
#define thermo_h

#include <stdio.h>

#define KAPPA 0.28562982892500527
#define INT_NLEVS 2000
#define HMAX 20000

#define G 9.8076
#define HV 2501000.0
#define RSD 287.0
#define EPS 0.622

typedef struct tph {
    double T;
    double p;
    double h;
} tph;

double Tv();

double satPressureBuck(double T);

double mixr(double Td, double p);

double relative_humidity(double T, double Td,double p);

tph lcl_temp(double T, double Td, double p);

double dT(double T, double p);

double* profile_malr(double* p, double* height, double T, int nlevs);

double* profile_dalr(double* p, double T, double start_pres, int nlevs);

double cape(double* T, double* Td, double* p, double* h, double T0, double Td0, double p0, int nlevs);

#endif /* thermo_h */
