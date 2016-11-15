/* 
*Author: Maike Jung
*Date: 2.11.2016

*Purpose: Draw random events from a certain mass spectrum

SN - Model: Lawrence-Livermore
    time spectrum is convoluted with the first hit distribution, to account for not knowing the absolute arrival times

UNITS: mass: eV
       energy: MeV
       distance: Mpc
       time: s

BINNING defined in .h file
*/
#include "spectrum.h"
#include <time.h>

void createEvents(double mass, double distance, double events, double stepSize, bool triggEff, bool energyRes, bool checkEntireSpectrum, int filenumber){

    /*read in trigger efficiency*/
    int i;
    double triggerEnergy[RESE+1];
    double triggerEfficiency[RESE+1];
    /* initialize with 1 in case triggerEfficiency is not used*/
    for(i = 0; i < RESE+1 ; triggerEfficiency[i++] = 1.0);
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

	double *myTrue= (double*) malloc((RESE-1) * REST * sizeof(double));
    /*create dist*/
	generateDist(mass, distance, events, myTrue, triggerEfficiency, energyRes);

    /*generate random events E=0.1-60.0 0.1
                             t=0.0-10.0 0.01*/
   
    /*storing time & energy in file*/
    char filename[sizeof "1.5eV_ideal/eventsGenerated_1.45eV_10.5Mpc_1000Events_real_1111.txt"];
    if (triggEff && energyRes){
        sprintf(filename, "%.1feV_real/events_%.2feV_%.1fMpc_%.0fEvents_real_%d.txt",mass, mass, distance, events, filenumber);
    }
    else {
        sprintf(filename, "%.1feV_ideal/events_%.2feV_%.1fMpc_%.0fEvents_ideal_%d.txt", mass, mass, distance, events, filenumber);
    }
    FILE *f = fopen(filename, "a+");
    if (f == NULL){
        printf("Error opening file!\n");
        exit(1);
    }

    int eventsGenerated = 0;
    int eventsTest = 0;
    srand( (unsigned)time( NULL ) );
    int randE, randT;
    double randCheck;
    // find the maximum in the list and use this for the random number gen
    /*double max = 0.0;
    for(i=0; i<((RESE-1)*REST);i++){
        if(myTrue[i]>max) max = myTrue[i];
    }
    printf("max val %f \n", max);*/
    while(eventsGenerated < 160){
        randE = rand() % (RESE+1);
        randT = rand() % (REST+1);
        randCheck = rand()*19.0/RAND_MAX;
        if (myTrue[randT*(RESE-1)+randE] >= randCheck){
            //printf("myTrue rand %f %f \n", myTrue[randE*(RESE-1)+randT], randCheck);
            fprintf(f, "%d %d\n", randE, randT);
            eventsGenerated += 1;
        }
        eventsTest+=1;
    }
    //printf("%d eventgen %d eventest \n", eventsGenerated, eventsTest);

    fclose(f);
    free(myTrue);
}

int main(void){
	/*set parameters*/
    /*flag for trigger efficiency*/
    bool triggEff = false;
    bool energyRes = false;
    bool checkEntireSpectrum = false;
    double mass;
    double distance = 1.0;
    double events = 160;
    double stepSize = 0.001;
    int filenumber;
    
    /*calculate uncertainty for certain configuration*/
    for (filenumber=101; filenumber<1001; filenumber++){ 
        printf("creating file %d \n", filenumber);
        createEvents(1.5, distance, events, stepSize, triggEff, energyRes, checkEntireSpectrum, filenumber);
        createEvents(1.0, distance, events, stepSize, triggEff, energyRes, checkEntireSpectrum, filenumber);
    }
    for (filenumber=1; filenumber<1001; filenumber++){ 
        printf("creating file %d \n", filenumber);
        createEvents(0.5, distance, events, stepSize, triggEff, energyRes, checkEntireSpectrum, filenumber);
    }

    printf("DONE\n");
    return 0;
}
