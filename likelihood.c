/* 
*Author: Maike Jung
*Date: 15.11.2016

*Purpose: Calculate the likelihood for the random events generated with generateEvents.c to belong to a certain mass spectrum

SN - Model: Lawrence-Livermore
    time spectrum is convoluted with the first hit distribution, to account for not knowing the absolute arrival times

UNITS: mass: eV
       energy: MeV
       distance: Mpc
       time: s

BINNING: see spectrum.h
*/

#include "spectrum.h"

void calcLLH(double mass, double distance, double events, bool triggEff, bool energyRes, int filenumber){

    /*read in trigger efficiency*/
    int i;
    double triggerEnergy[RESE+1];
    double triggerEfficiency[RESE+1];
    /* initialize with 1 in case triggerEfficiency is not used*/
    for(i = 0; i < RESE+1 ; triggerEfficiency[i++] = 1);
    if(triggEff){
        FILE *myFile;
        if (RESE==600){
            myFile = fopen("trigger_efficiency_100keV_steps.txt", "r");
        }
        else if (RESE==6000){
            myFile = fopen("trigger_efficiency_10keV_steps.txt", "r");
        }
        else if (RESE==60000){
            myFile = fopen("trigger_efficiency_1keV_steps.txt", "r");
        }
        else {
            printf("Invalid grid size for the energy resolution.");
        }
        for (i = 0; i < RESE+1; i++) {
            fscanf(myFile, "%lf %lf", &triggerEnergy[i], &triggerEfficiency[i]);
        }
        fclose(myFile);
    }
   
    /*load events & store energy and time in arrays*/
    char filename[sizeof "1eV_ideal/eventsGenerated_1.34eV_10.5Mpc_1000Events_ideal_1111.txt"];
    if (triggEff && energyRes){
        sprintf(filename, "%.1feV_real_1/events_%.2feV_%.1fMpc_%.0fEvents_real_%d.txt",mass,  mass, distance, events, filenumber);
    }
    else {
        sprintf(filename, "%.1feV_ideal_test2/events_%.2feV_%.1fMpc_%.0fEvents_ideal_%d.txt", mass, mass, distance, events, filenumber);
    }
    FILE *f = fopen(filename, "r");
    int eventEnergy[(int) events];
    int eventTime[(int) events];
    for(i = 0; i < events; eventEnergy[i++] = 1);
    for(i = 0; i < events; eventTime[i++] = 1);
    for (i = 0; i < events; i++){
        fscanf(f, "%d %d", &eventEnergy[i], &eventTime[i]);
    }

    // calculate the likelihood
    double llh;
    double testMass;
    // store current value of the minimumLLH and the corresponding mass
    double minLLH = INFINITY;
    double massOfMinLLH;
    // go over all the spectra around a certain range of the input mass & calculate the likelihood for each spectrum
    double *testSpectrum= (double*) malloc((RESE-1) * REST * sizeof(double));
    // first go over broad range - there are no negative entries in the spectrum!!!!!
    for (testMass = mass - 0.3; testMass <= mass + 0.3; testMass+=0.1){
        llh = 0.0;
        generateDist(testMass, distance, events, testSpectrum, triggerEfficiency, energyRes);
        for (i = 0; i < events; i++){
            if (testSpectrum[eventTime[i]*(RESE-1)+eventEnergy[i]] < pow(10,-200)){
                llh += -10000000;   
            }
            else llh += log(testSpectrum[eventTime[i]*(RESE-1)+eventEnergy[i]]);
            //printf("tset %d %f \n", i, llh);
        }
        llh*=-1;
        printf("mass %f, llh %f\n",testMass, llh);
        if (llh < minLLH) {
            minLLH = llh;
            massOfMinLLH = testMass;
        }
    }
    double currentMinimum = massOfMinLLH;

    // check around the minimum mass in finer steps
    for (testMass = currentMinimum - 0.05; testMass <= currentMinimum + 0.05; testMass += 0.01){
        llh = 0.0;
        generateDist(testMass, distance, events, testSpectrum, triggerEfficiency, energyRes);
        for (i = 0; i < events; i++){
            if (testSpectrum[eventTime[i]*(RESE-1)+eventEnergy[i]] < pow(10,-200)){
                llh += -10000000;   
            }
            else llh += log(testSpectrum[eventTime[i]*(RESE-1)+eventEnergy[i]]);
            //printf("tset %d %f \n", i, llh);
        }
        llh*=-1;
        printf("mass %f, llh %f \n",testMass, llh);
        if (llh < minLLH) {
            minLLH = llh;
            massOfMinLLH = testMass;
        }
    }

    printf("mass found: %f eV\n", massOfMinLLH);

    fclose(f);
    free(testSpectrum);

    /*write to file*/
    char filename2[sizeof "Results_Likelihood/test_masses_1.55eV_1Mpc_real.txt"];
    if (triggEff && energyRes){
        sprintf(filename2, "masses_%.2feV_%.1fMpc_real_1.txt", mass, distance);
    }
    else {
        sprintf(filename2, "Results_Likelihood/masses_%.2feV_%.1fMpc_ideal_test2.txt", mass, distance);
    }
    FILE *g = fopen(filename2, "a+");
    fprintf(g, "%f \n", massOfMinLLH);
    fclose(g);
}

int main(void){
	/*set parameters*/
    /*flag for trigger efficiency*/
    bool triggEff = true;
    bool energyRes = true;
    double mass = 1.0;
    double distance = 1.0;
    double events = 160;
    int filenumber;

    /*calculate uncertainty for certain configuration*/
    for (filenumber=1001; filenumber<2001; filenumber++){ 
        printf("evaluating file %d \n", filenumber);
        calcLLH(mass, distance, events, triggEff, energyRes, filenumber);
    }

    printf("DONE\n");
    return 0;
}
