#ifndef FUNCTION_H
#define FUNCTION_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define Euler 2.71828182845904523536
#define G 6.67430e-11
#define M_Earth 5.9722e24
#define omega_Earth 7.292115e-5

typedef struct {
    double t, t_oe, sqrt_a_s, e_s, i_0, OMEGA_0, omega_satelite, M_0, delta_n, i_dot, OMEGA_dot, C_uc, C_us, C_rc, C_rs, C_ic, C_is;
    char tps[50];
} stateLiteData;

typedef struct {
    double T_k, n_0, n, M_k, E_k, f_k, phi_k, sigma_u, sigma_r, sigma_i, u_k, r_k, i_k, X_orbitalPlane, Y_orbitalPlane, OMEGA_k, X_ECEF, Y_ECEF, Z_ECEF;
} caculatingData;

int Skip (FILE *fileName) {
    char dustbin = 'a';
    while( dustbin != '#' ){
        dustbin = fgetc(fileName);
        }
        dustbin = fgetc(fileName); // Skip the \n
        dustbin = 1;
    return 1;
}

int SkipPro (FILE *fileName, char endTag, int mode) // mode = 0: don't pass \n, mode = 1: pass \n 
{
    char dustbin = 'a';
    while( dustbin != endTag ){
        dustbin = fgetc(fileName);
        }
    if (mode == 0){
        return 1;
    }else if (mode == 1){
        dustbin = fgetc(fileName); // Skip the \n
        dustbin = 1;
        return 1;
    }else{
        printf("Error mode in SkipPro function.\n");
        return -1;
    }
}

int readInData ( FILE *dataFile, stateLiteData *data_P ) {
    char dustbin = 'a';
    int intBuffer = 0;
    double *doubleMemberPointer = (double *)data_P; // Trans the sateliteData pointer to double pointer

    if (dataFile == NULL) {
        printf("Error reading dataFile.\n");
        return -1;
    }
    Skip(dataFile);
    while( dustbin != '#' ){
        dustbin = fgetc(dataFile);
        SkipPro(dataFile, '=', 0);
        fscanf(dataFile, "%lf\n", doubleMemberPointer + intBuffer);
        intBuffer++;
        dustbin = fgetc(dataFile); // get the first char of the next line
        }
    intBuffer = 0;
    return 1;
}

double julianCenturies(double T_send, double T_oe) {
    double T_k = 0.000; // julianCenturies
    T_k = T_send - T_oe;
    if ( T_k > 302400.000 )
        T_k -= 604800.000;
    else if ( T_k < -302400.000 ) 
        T_k += 604800.000;
    return T_k;
}

double meanMotion(double sqrt_a, double delta_n) {
    double n = 0.000;
    double n_0 = 0.000;
    n_0 = sqrt( (G * M_Earth) / pow(sqrt_a,6) );
    n = n_0 + delta_n;
    return n;
}

double meanAnomaly(double M_0, double T_k, double n) {
    double M_k = 0.000;
    M_k = M_0 + T_k * n;
    return M_k;
}

double solveKepler(double M, double e) {
    #define MAX_ITER 100
    #define EPSILON 1e-12

    double E = M; // 初始猜值
    for (int i = 0; i < MAX_ITER; i++) {
        double f = E - e * sin(E) - M;
        double f_prime = 1 - e * cos(E);
        double delta = -f / f_prime;
        E += delta;

        if (fabs(delta) < EPSILON) {
            break;
        }
    }
    return E;
}

double trueAnomaly(double e_s, double E_k) {
    double f_k = 0.000; // true anomaly
    f_k = atan2(sqrt(1 - pow(e_s, 2)) * sin(E_k), cos(E_k) - e_s);
    return f_k;
}

double argumentOfLatitude(double f_k, double omega_satelite) {
    double phi_k = 0.000;
    phi_k = f_k + omega_satelite;
    return phi_k;
}

void perturbationCorrection(stateLiteData *data_p, caculatingData *res_p) {
    res_p->sigma_u = data_p->C_uc * cos(2*res_p->phi_k) + data_p->C_us * sin(2*res_p->phi_k);
    res_p->sigma_r = data_p->C_rc * cos(2*res_p->phi_k) + data_p->C_rs * sin(2*res_p->phi_k);
    res_p->sigma_i = data_p->C_ic * cos(2*res_p->phi_k) + data_p->C_is * sin(2*res_p->phi_k);
    res_p->u_k = res_p->phi_k + res_p->sigma_u;
    res_p->r_k = pow( data_p->sqrt_a_s, 2) * ( 1 - data_p->e_s * cos(res_p->E_k) ) + res_p->sigma_r;
    // printf("%lf\n",pow( data_p->sqrt_a_s, 2));
    res_p->i_k = data_p->i_0 + res_p->sigma_i + data_p->i_dot * ( data_p->t - data_p->t_oe);

}

void satelliteCoordinateComputation(stateLiteData *data_p, caculatingData *res_p) {
    //Compute the satellite's position in the orbital plane
    res_p->X_orbitalPlane = res_p->r_k * cos(res_p->u_k);
    res_p->Y_orbitalPlane = res_p->r_k * sin(res_p->u_k);
    //Calculate the longitude of the ascending node (Ω)
    res_p->OMEGA_k = data_p->OMEGA_0 + ( data_p->OMEGA_dot - omega_Earth)*( data_p->t - data_p->t_oe ) - omega_Earth * data_p->t_oe;
    //Compute the satellite's position in the instantaneous Earth-fixed frame
    res_p->X_ECEF = res_p->X_orbitalPlane * cos(res_p->OMEGA_k) - res_p->Y_orbitalPlane * cos(res_p->i_k) * sin(res_p->OMEGA_k);
    res_p->Y_ECEF = res_p->X_orbitalPlane * sin(res_p->OMEGA_k) + res_p->Y_orbitalPlane * cos(res_p->i_k) * cos(res_p->OMEGA_k);
    res_p->Z_ECEF = res_p->Y_orbitalPlane * sin(res_p->i_k);
}

void writeOut (FILE *result, stateLiteData *data_p, caculatingData *cul_p) {
    fprintf(result, "■■■■■■-Satellite Parameter Calculation Program-■■■■■■\n");
    fprintf(result,"Dev by autumnal_leaf in 2025\n");
    fprintf(result,"Version: 1.0.0 - SNAPSHOT\n");
    fprintf(result, "■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■\n");
    fprintf(result,"【STEP 1】 计算规划时间\n");
    fprintf(result,"T_k = T - T_oe =\t%.9lf\n",cul_p->T_k);
    fprintf(result, "-----------------------------------------------------\n");
    fprintf(result,"【STEP 2】 计算卫星运行的平均角速度\n");
    fprintf(result,"n_0 =\t%.9le\n",cul_p->n_0);
    fprintf(result,"n =\t%.9le\n",cul_p->n);
    fprintf(result, "-----------------------------------------------------\n");
    fprintf(result,"【STEP 3】 观测时刻的平近点角\n");
    fprintf(result,"M_k =\t%.9le\n",cul_p->M_k);
    fprintf(result, "-----------------------------------------------------\n");
    fprintf(result,"【STEP 4】 迭代法计算偏近点角\n");
    fprintf(result,"E_k =\t%.9le\n",cul_p->E_k);
    fprintf(result, "-----------------------------------------------------\n");
    fprintf(result,"【STEP 5】 计算真近点角\n");
    fprintf(result,"f_k =\t%.9le\n",cul_p->f_k);
    fprintf(result, "-----------------------------------------------------\n");
    fprintf(result,"【STEP 6】 计算升交距角\n");
    fprintf(result,"Φ_k =\t%.9le\n",cul_p->phi_k);
    fprintf(result, "-----------------------------------------------------\n");
    fprintf(result,"【STEP 7】 计算摄动改正\n");
    fprintf(result,"σu =\t%.9le\n",cul_p->sigma_u);
    fprintf(result,"σr =\t%.9le\n",cul_p->sigma_r);
    fprintf(result,"σi =\t%.9le\n",cul_p->sigma_i);
    fprintf(result, "-----------------------------------------------------\n");
    fprintf(result,"【STEP 8】 进行摄动改正\n");
    fprintf(result,"u_k =\t%.9le\n",cul_p->u_k);
    fprintf(result,"r_k =\t%.9le\n",cul_p->r_k);
    fprintf(result,"i_k =\t%.9le\n",cul_p->i_k);
    fprintf(result, "-----------------------------------------------------\n");
    fprintf(result,"【STEP 9】 计算卫星的轨道平面坐标\n");
    fprintf(result,"x_k =\t%.9le\n",cul_p->X_orbitalPlane);
    fprintf(result,"y_k =\t%.9le\n",cul_p->Y_orbitalPlane);
    fprintf(result, "-----------------------------------------------------\n");
    fprintf(result,"【STEP 10】 计算升交点经度\n");
    fprintf(result,"Ω_k =\t%.9le\n",cul_p->OMEGA_k);
    fprintf(result, "-----------------------------------------------------\n");
    fprintf(result,"【STEP 11】 计算卫星在瞬时地球坐标系中的位置\n");
    fprintf(result,"X_k =\t%.9le\n",cul_p->X_ECEF);
    fprintf(result,"Y_k =\t%.9le\n",cul_p->Y_ECEF);
    fprintf(result,"Z_k =\t%.9le\n",cul_p->Z_ECEF);
    fprintf(result, "■■■■■■■■■■■■■■■■■■■■■■   DONE  ■■■■■■■■■■■■■■■■■■■■■■■\n");
}
#endif //FUNCTION_H