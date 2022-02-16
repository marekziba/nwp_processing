//
//  thermo.c
//  
//
//  Created by Marek Zięba on 27/12/2019.
//

#include "thermo.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

//  double satPressureBuck(double T)
//
//  Calculates vapour pressure of saturated air with given temperature
//
//  Arguments:
//      T (double) - temperature of air parcel (in K)
//
//  Return value:
//      saturation vapour pressure for the parcel (in Pa)
double satPressureBuck(double T){
    // the Buck formula operates on temperature in deg C,
    // but for convenience the function gets temperature in
    // Kelvins, therefore we need to convert it
    T -= 273.15;
    
    double exponent = (18.678 - (T/234.5))*(T/(257.14+T));
    double saturation_pressure = 0.61121 * exp(exponent);
    
    // the result is in kPa, but since we want to stick to
    // standard SI units we need to multiply by 10^3 to get
    // the value in Pa
    return saturation_pressure*1000;
    //return exp(20.386 - 5132.0/T) * 133.322;
}

//  double mixr(double T, double Td)
//
//  Calculates mixing ratio/specific humidity for air parcel
//  of given temperature, dewpoint and pressure
//
//  Arguments:
//      T (double) - temperature of air parcel (in K)
//      Td (double) - dewpoint temperature of air parcel (in K)
//      p (double) - air pressure (in Pa)
//
//  Return value:
//      mixing ratio of this parcel (dimensionless [kg/kg])
double mixr(double Td, double p){
    double pg = satPressureBuck(Td);
    return (0.622*pg)/(p-pg);
}

//  double relative_humidity(double T, double Td, double p)
//
//  Calculate relative humidity from temperature, dewpoint and pressure
//
//  Arguments:
//      T (double)  -   air temperature (in K)
//      Td (double) -   dewpoint (in K)
//      p (double)  -   air pressure (in Pa)
//
//  Return value:
//      relative humidity, fractional (dimensionless)
double relative_humidity(double T, double Td,double p){
    //  simply a ratio of mixing ratio of given parcel of air to saturation
    //  mixing ratio of this parcel
    return mixr(Td,p)/mixr(T,p);
}

//  double lcl_height(double T, double Td, double p)
//
//  Numerically compute the pressure at Lifted Condensation Level
//  Function uses iterative approach to approximate the value of
//  LCL pressure described in Stackpole, J. D. 1967
//
//  Arguments:
//      T (double)  -   temperature of lifted parcel (in K)
//      Td (double) -   dewpoint of lifted parcel (in K)
//      p (double)  -   pressure at level from which parcel is lifted
//                      (in Pa)
//
//  Return value:
//      air pressure at LCL level (in Pa)
tph lcl_temp(double T, double Td, double p){
    
    int n_iter = 0;
    double T0 = Td;
    double mix_ratio = mixr(Td,p);
    
    // initial conditions
    double T1 = 1 / ( 1/Td - 6.561*0.0001*(1+0.25*mix_ratio) * log(T0/T) );
    
    //  numerically approximate the solution for the equation until the difference
    //  between two subsequent solutions is lower than the given threshold, in this
    //  case |T1-T0| < 0.001
    while(fabs(T1-T0) > 0.00001 && n_iter < 500){
        T0 = T1;
        T1 = 1 / ( 1/Td - 6.561*0.0001*(1+0.25*mix_ratio) * log(T0/T) );
        n_iter++;
    }
    
    //  calculate LCL pressure and height from LCL temperature using equations proposed in
    //  [source needed]
    double cpm = (1 - mix_ratio)*(719+287.04) + mix_ratio*(1418+461);
    double rm = (1-mix_ratio)*287.04 + mix_ratio*461;
    double lcl_pres = p*pow((T1/T),(cpm/rm));
    
    double lcl_height = (cpm/9.8076)*(T-T1);
    
    //printf("Pressure at LCL: %lf\nTemperate at LCL: %lf\nLCL height: %lf\n",lcl_pres,T1,lcl_height);

    tph res; res.T = T1; res.p = lcl_pres; res.h = lcl_height;
    return res;
}

double dT(double T, double p){
    
    //  Calculate the moist adiabatic lapse rate for air parcel of given
    //  temperature and pressure (in K/m)
    //
    //  Arguments:
    //      T (double)  -   temperature of lifted parcel (in K)
    //      p (double)  -   pressure level from which the parcel is
    //                      lifted, assuming it's same as parcel
    //                      pressure (in Pa)
    //  Return value:
    //      temperature lapse rate for given parcel (in K/m)
    
    //  currently constants are defined here as variables as a temporary solution
    //  future plan is to #define all of those constans in header file
    double g, Rsd, Hv, r, Cpd, eps, e, rs;
    
    g = 9.8076;
    Hv = 2501000.0;
    Rsd = 287.0;
    eps = 0.622;
    
    e = satPressureBuck(T);
    r = (eps*e)/(p-e);
    rs = mixr(T,p);
    
    //printf("r = %lf, rs = %lf\n",r,rs);
    Cpd = 1003.5;
    
    double dt = g*( Rsd*(T*T) + Hv*rs*T )/( Cpd*Rsd*T*T + Hv*Hv*rs*eps );
    
    return dt;
}

double* profile_malr(double* p, double* height, double T, int nlevs){
    //  Calculate the vertical temperature profile of moist adiabatically lifted
    //  parcel of given start temperature. The parcel temperature is calculated
    //  at each level for which data are provided in p and height arrays.
    //
    //  Arguments:
    //      p (array: double)           vertical profile of pressure (in Pa)
    //      height (array: double)      height of every pressure level (in m)
    //      T (double)                  starting (LCL) temperature of lifted parcel
    //      nlevs (int)                 number of pressure levels
    
    // allocate memory for output data array
    double* profile = (double*)malloc(nlevs*sizeof(double));
    
    for(int j = 0; j < nlevs; j++){
        printf("%lf - %lf\n",p[j]/100,height[j]);
    }
    
    //  assuming output profile includes LCL, we need to include
    //  parcel starting temperature in the profile
    profile[0] = T;
    // printf("p0 = %lf, T0 = %lf\n",p[0],profile[0]-273.15);
    // printf("----------------\n");
    int n; double h_step, temp, pres;
    
    for(int i=1; i<nlevs; i++){
        
        // n = abs((p[i]-p[i-1])/2500);
        h_step = height[i]-height[i-1];
        printf("n = %d, h_step = %lf\n", n, h_step);
        
        temp = profile[i-1]; pres = p[i-1];
        // printf("n = %d\n",n);
        // for(int i=0;i<n;i++){
        //     temp -= dT(temp,pres-2500)*h_step;
        // }

        temp -= dT(temp,pres)*h_step;
        
        // profile[i] = profile[i-1] - (dT(profile[i-1],p[i]) * (height[i] - height[i-1]));
        
        profile[i] = temp;
        
        // printf("%lf - %lf\n",p[i],profile[i]-273.15);
    }
    
    return profile;
}


 double* profile_dalr(double* p, double T, double start_pres, int nlevs){
     //  Calculate the vertical temperature profile of dry adiabatically lifted
     //  parcel of given start temperature and pressure. The parcel temperature is
     //  calculated at each level for which data are provided in p.
     //
     //  Arguments:
     //     p (array: double):
     //         vertical profile of pressure (in Pa)
     //     T (double):
     //         starting temperature of lifted parcel
     //     start_pres (double):
     //         starting pressure of lifted parcel
     //     nlevs (int):
     //         number of pressure levels
     
     // allocate memory for output data array
     double* profile = (double*)malloc(nlevs*sizeof(double));
     
     for(int i = 0; i < nlevs; i++){
         profile[i] = T * pow((p[i]/start_pres),KAPPA);
     }
     return profile;
 }

double cape(double* T, double* p, double* h, double T0, double Td0, double P0, int nlevs){
    
    tph lcl = lcl_temp(T0,Td0,P0);
    printf("Pressure at LCL: %lf\nTemperate at LCL: %lf\nLCL height: %lf\n",lcl.p,lcl.T-273.15,lcl.h);
    int i=0;
    while(h[i] < lcl.h && i < nlevs){
        i++;
    }
    
    printf("test = %lf\n",p[i]/100);
    
    if(i == nlevs){
        return 0.0;
    }
    
    //  insert LCL temperature, pressure and height into given profiles
    
    double* temperature = (double*) malloc((nlevs+1)*sizeof(double));
    double* pressure = (double*) malloc((nlevs+1)*sizeof(double));
    double* height = (double*) malloc((nlevs+1)*sizeof(double));
    
    memcpy(temperature,T,i*sizeof(double));
    memcpy(temperature+i,&lcl.T,sizeof(double));
    memcpy(temperature+i+1,T+i,(nlevs-i)*sizeof(double));
    
    memcpy(pressure,p,i*sizeof(double));
    memcpy(pressure+i,&lcl.p,sizeof(double));
    memcpy(pressure+i+1,p+i,(nlevs-i)*sizeof(double));
    
    memcpy(height,h,i*sizeof(double));
    memcpy(height+i,&lcl.h,sizeof(double));
    memcpy(height+i+1,h+i,(nlevs-i)*sizeof(double));
    
    // double* dry_profile = profile_dalr(pressure,T0,P0,i);
    // double* moist_profile = profile_malr(pressure+i,height+i,lcl.T,nlevs-i+1);
    
    // double* profile = (double*) malloc(nlevs*sizeof(double));
    
    // memcpy(profile,dry_profile,i*sizeof(double));
    // memcpy(profile+i,moist_profile,(nlevs-i)*sizeof(double));

    // printf("\n++++++ FINAL PROFILE ++++++\n");

    // for(int i = 0; i < nlevs; i++){
    //     printf("%lf\t%lf\n",pressure[i],profile[i] - 273.15);
    // }

    // ++++++ INTERPOLATE PROFILE ++++++
    // double dh = 0.1;
    double* interpolated_temperature = (double*) malloc(INT_NLEVS*sizeof(double));
    double* interpolated_pressure = (double*) malloc(INT_NLEVS*sizeof(double));
    double* interpolated_height = (double*) malloc(INT_NLEVS*sizeof(double));
    
    double p0 = pressure[0], p1 = pressure[1]; 
    double t0 = temperature[0], t1 = temperature[1];
    double h0 = height[0], h1 = height[1];

    double hdiff = (HMAX - h0)/INT_NLEVS;
    double initial_height = h0;
    printf("Initial height = %lf\n", initial_height);

    int hindex = 0;

    int interpolated_profile_size = 0;

    for(int i = 0; i < INT_NLEVS; i++){
        double h = initial_height + i * hdiff;
        if(i == 0) printf("%lf\n", h);

        //TODO: if this if fires, then interpolation and assignment to arrays is not performed - fix this!
        if(h > h1){
            if(hindex + 1 < nlevs){
                hindex++;

                h0 = h1;
                h1 = height[hindex+1];

                p0 = p1;
                p1 = pressure[hindex+1];

                t0 = t1;
                t1 = temperature[hindex+1];
            } else {
                break;
            }
        }

        interpolated_height[i] = h;
        
        if(h == h1){
            // handle situation
            interpolated_temperature[i] = temperature[hindex+1];
            interpolated_pressure[i] = pressure[hindex+1];
        }
        else if(h == h0){
            interpolated_temperature[i] = temperature[hindex];
            interpolated_pressure[i] = pressure[hindex];
        }
        else {
            double w0 = h - h0;
            double w1 = h1 - h;
            

            w0 = w0/(h1 - h0);
            w1 = w1/(h1 - h0);

            interpolated_pressure[i] = p0 * w1 + p1 * w0;
            // printf("interpolated_pressure[%d] = %lf * %lf + %lf * %lf = %lf\n", i, p0, w1, p1, w0, interpolated_pressure[i]);
            interpolated_temperature[i] = t0 * w1 + t1 * w0;
            // printf("interpolated_temperature[%d] = %lf * %lf + %lf * %lf = %lf\n", i, t0, w1, t1, w0, interpolated_temperature[i]);
        }

        interpolated_profile_size++;
    }

    int iprofile_lcl_index = 0;
    while(interpolated_height[iprofile_lcl_index] < lcl.h && iprofile_lcl_index < interpolated_profile_size){
        iprofile_lcl_index++;
    }

    printf("+++++++++++++++++\niprofile_lcl_index = %d\n+++++++++++++++++\n", iprofile_lcl_index);

    double* full_interpolated_temperature = (double*) malloc((interpolated_profile_size+1)*sizeof(double));
    double* full_interpolated_pressure = (double*) malloc((interpolated_profile_size+1)*sizeof(double));
    double* full_interpolated_height = (double*) malloc((interpolated_profile_size+1)*sizeof(double));

    memcpy(full_interpolated_temperature,interpolated_temperature,iprofile_lcl_index*sizeof(double));
    memcpy(full_interpolated_temperature+iprofile_lcl_index,&lcl.T,sizeof(double));
    memcpy(full_interpolated_temperature+iprofile_lcl_index+1,interpolated_temperature+iprofile_lcl_index,(interpolated_profile_size-iprofile_lcl_index)*sizeof(double));
    
    memcpy(full_interpolated_pressure,interpolated_pressure,iprofile_lcl_index*sizeof(double));
    memcpy(full_interpolated_pressure+iprofile_lcl_index,&lcl.p,sizeof(double));
    memcpy(full_interpolated_pressure+iprofile_lcl_index+1,interpolated_pressure+iprofile_lcl_index,(interpolated_profile_size-iprofile_lcl_index)*sizeof(double));
    
    memcpy(full_interpolated_height,interpolated_height,iprofile_lcl_index*sizeof(double));
    memcpy(full_interpolated_height+iprofile_lcl_index,&lcl.h,sizeof(double));
    memcpy(full_interpolated_height+iprofile_lcl_index+1,interpolated_height+iprofile_lcl_index,(interpolated_profile_size-iprofile_lcl_index)*sizeof(double));

    double* dry_profile = profile_dalr(full_interpolated_pressure,T0,P0,iprofile_lcl_index);
    double* moist_profile = profile_malr(full_interpolated_pressure+iprofile_lcl_index,full_interpolated_height+iprofile_lcl_index,lcl.T,interpolated_profile_size-iprofile_lcl_index+1);
    
    double* profile = (double*) malloc((interpolated_profile_size+1)*sizeof(double));
    
    memcpy(profile,dry_profile,iprofile_lcl_index*sizeof(double));
    memcpy(profile+iprofile_lcl_index,moist_profile,(interpolated_profile_size+1-iprofile_lcl_index)*sizeof(double));

    printf("\n++++++ INTERPOLATED PROFILE ++++++\n");

    for(int i = 0; i < interpolated_profile_size+1; i++){
        printf("%lf\t%lf\t%lf\t%lf\n",full_interpolated_height[i], full_interpolated_pressure[i], full_interpolated_temperature[i] - 273.15, profile[i] - 273.15);
        
    }

    return 4200.0;
}
 

void printarray(int* array, int length){
    for(int i=0;i<length;i++){
        printf("%d ",array[i]);
    }
    printf("\n");
}

int main(int argc, const char* argv[]){
    /*
     printf("Saturation pressure: %lf \n",satPressureBuck(temp+273.15));
     printf("Mixing ratio: %lf \n", mixr(20+273.15,pressure)*1000);
     printf("Relative humidity: %lf%% \n \n",relative_humidity(30+273.15,20+273.15,1000)*100);
     double huj = lcl_height(34.9+273.15,22.9+273.15,100400);
     huj+=1;
     double pressure[4] = {85000,80000,75000,70000};
     double height[4] = {1551,2059,2595,3160};
     double* profile = profile_malr(pressure,height,20+273.15,4);
     for(int i=0; i<4; i++){
     printf("%lf, ",profile[i]);
     }
     printf("\n");
     */
    //int array[10] = {1,2,3,4,5,6,7,8,9,10};
    //printarray(array+2,8);
    double Temp[6] = {25.0+273.15,12.0+273.15,2.0+273.15,-16.0+273.15,-45.0+273.15,-55.0+273.15};
    /*
     for(int i=0;i<6;i++){
     Temp[i] = Temp[i]+273.15;
     }
     */
    double pres[6] = {100000.0,85000.0,70000.0,50000.0,30000.0,20000.0};
    double height[6] = {0.0,1500.0,3000.0,5600.0,9500.0,11600.0};
    double T0 = 25.0 + 273.15, Td0 = 15.0 + 273.15, p0 = 100000;
    double luj = cape(Temp,pres,height,T0,Td0,p0,6); luj += 1;
}
