#include "function.h"

void main() {

    FILE *dataFilePointer = fopen("files\\data.txt", "r");
    if (dataFilePointer == NULL) {
        printf("Error opening dataFile.\n");
        fclose(dataFilePointer);
        return;
    }
    FILE *resultFilePointer = fopen("files\\result.txt", "w");
    if (resultFilePointer == NULL) {
        printf("Error opening resultFile.\n");
        fclose(dataFilePointer);
        fclose(resultFilePointer);
        return;
    }
    
    stateLiteData data_0, *data_P;
    data_P = &data_0;

    caculatingData cul_0, *cul_p;
    cul_p = &cul_0;

    readInData( dataFilePointer, data_P );
    
    cul_0.T_k = julianCenturies( data_0.t, data_0.t_oe);
    cul_0.n = meanMotion( data_0.sqrt_a_s, data_0.delta_n);
    cul_0.n_0 = sqrt( (G * M_Earth) / pow(data_0.sqrt_a_s,6) );
    cul_0.M_k = meanAnomaly(data_0.M_0, cul_0.T_k, cul_0.n);
    cul_0.E_k = solveKepler(cul_0.M_k, data_0.e_s);
    cul_0.f_k = trueAnomaly(data_0.e_s, cul_0.E_k);
    cul_0.phi_k = argumentOfLatitude(cul_0.f_k, data_0.omega_satelite);
    perturbationCorrection( data_P, cul_p);
    satelliteCoordinateComputation( data_P,cul_p);

    writeOut(resultFilePointer, data_P, cul_p);
    //printf("Hello, World!\n");
    
    /* for(int i = 0; i < 17; i++){
        double *doubleMemberPointer = (double *)data_P;
        printf("%lf\n", *(doubleMemberPointer + i) );
    }
    */

    fclose(dataFilePointer);
    fclose(resultFilePointer);
}