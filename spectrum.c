#include "spectrum.h"

/* time shift due to neutrino mass */
double time_shift(double t, double E, double mass, double dist){
    double tDelta = dist*51.4635*(mass/E)*(mass/E);
    double time = t - tDelta;
    if (time <= 0){
        return 0.0;
    }   
    return LL_time_spectrum(time);
}

/* arrival time probability for a certain mass/distance - normalized */
double LLSpectrumTotal (double t, double E, double mass, double dist){
    return time_shift(t, E, mass, dist)*LL_energy_spectrum(E);
}

/* calculate the probability to get the first hit after a certain amount of time */
void ProbFirstHitDist (double mass, double dist, double events, double *result){
    /*arrival time distribution of all the hits (for a certain mass) - project the E,t spectrum
    on the t axis - t in 0.01 steps from 0 to 10 seconds*/
    double totalArrivalTimeDist[REST];
    int i;
    double sum;
    double y, j;
    for (i = 0; i < REST; i++){
        sum = 0.0;     
        /*integrate over the energy part for every part of the spectrum*/
        for (j = 0.01; j < 60.0; j += 0.01) {
            y = LLSpectrumTotal((i*10.0/REST), j, mass, dist);
            sum += y * 0.01; 
        } 
        totalArrivalTimeDist[i] = sum*(10.0/REST);
    }
    /*calculate the cumulative sum of the array*/
    double cumulative[REST];
    int k, l, m;
    double cum;
    for (k = 0; k <REST; k++){
        cum = 0.0;
        for (l = 0; l <= k; l++){
            cum += totalArrivalTimeDist[l];
        }
        cumulative[k] = cum;    
    }
    /*calculate the 1st Hit Distribution - and normalize it*/
    double count = 0.0;
    for (m = 0; m < REST; m++){
        result[m] = totalArrivalTimeDist[m]*events*pow((1 - cumulative[m]), events-1);
        count += result[m];
    }
    for (m = 0; m < REST; m++){
        result[m] = result[m]/count;
    }
}

/*calculate the correlation - new spectrum between -3 and 10s*/
/*this is stored in an array so newSpec[0] corresponds to a time of -3s 
and newSpec[1.3*REST-1] to 10s*/
void correlation(double mass, double dist, double events, double *newSpec){
    double hitDist[REST];
    ProbFirstHitDist(mass, dist, events, hitDist);
    int i, j;
    double pNew;
    /*perform the convolution*/
    for (i = 0; i < REST*1.3; i++){
        pNew = 0.0;
        for (j = 0; j < REST; j++){
            if ((i-0.3*REST + j) < REST && (i-0.3*REST + j) > 0){
                pNew += hitDist[j] * LL_time_spectrum( (j+i-0.3*REST)/(0.1*REST) );
            }
        newSpec[i] = pNew;
        }
    }
}

/*generate the proper distribution*/
void generateDist(double mass, double dist, double events, double *distribution, double *triggerEfficiency, bool energyRes){
	double timeArray[(int) (1.3*REST)];
	correlation(mass, dist, events, timeArray);
	double timeShift;
	double time;
	int t, e, f, g;
	int arrayIndex;
    
    /*energy resolution taken into account*/
    if (energyRes){
        for (t=0; t<REST;t++){
            double energySpectrum[RESE];
            energySpectrum[0] = 0.0;
            for (e=1; e<RESE; e++){  
			    timeShift = dist*(mass/e)*(mass/e)*51.4635*100;
			    time = t/(0.1*REST) - timeShift;
			    arrayIndex = (int) (time*(0.1*REST) + 0.3*REST);       
			    if (arrayIndex <= 0){
				    arrayIndex = 0;
			    }
                energySpectrum[e] = LL_energy_spectrum(e/(RESE/60.0))*timeArray[arrayIndex]*triggerEfficiency[e];
            }
            for (f=1; f<RESE; f+=1){
                double pNew = 0.0;
                for (g=-RESE; g<RESE+1; g+=5){
                    if (f-g >=0 && f-g <=RESE){
                        pNew += GAUSS(g/(RESE/60.0),f/(RESE/60.0))*energySpectrum[f-g];
                    }
                distribution[t*(RESE-1) +f] = pNew/(RESE/60.0);
                }
            }
        }
    }
    /*without energy resolution*/
    else {
        for (t=0; t<REST;t++){
                for (e=1; e<RESE; e++){  
			        timeShift = dist*(mass/e)*(mass/e)*51.4635*(RESE/60.0)*(RESE/60.0);
			        time = t/(0.1*REST) - timeShift;
			        arrayIndex = (int) (time*(0.1*REST) + (0.3*REST));
			        if (arrayIndex <= 0){
				        arrayIndex = 0;
			        }
                    distribution[t*(RESE-1) +e] = LL_energy_spectrum(e/(RESE/60.0))*timeArray[arrayIndex]*triggerEfficiency[e];
                }
            }
    }
    int k;
    double normalize = 0;
    for (k=0; k<(RESE-1)*REST;k++){
        normalize += distribution[k]*(1/(REST*0.1))*(1/(RESE/60.0));
    }
    // todo include events for binned likelihood
    for (k=0; k<(RESE-1)*REST;k++){
        distribution[k] *= 1.0/normalize;
    }
}
