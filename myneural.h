#ifndef __MY_NEURAL__
#define __MY_NEURAL__
#include <stdint.h>
#include <math.h>
#include "kbhit.h"

#define CHEBY_ORDER 8
#define CTRLF_NORMPARAMS    0x01
#define CTRLF_SQUASHED       0x02
#define CTRLF_NORMALIZE     0x04

typedef unsigned long int t_neuronum;


typedef struct
{
    unsigned long int N;
    double *T; // cada valor del polinomios de chebyshev 
    double *d; // cada valor de cada derivada 
    double *p; // cada parámetro
    double lrate;
    unsigned char ControlFlags;
    
} ts_fparams;

typedef ts_fparams* tps_fparams;

typedef struct {
        t_neuronum ninputs; // Número de entradas
        double b; // bias;
        double *w; //coeficientes de ponderación (pesos)
        double **x; // entradas
        double v; 
        double o;   // salida
        double delta;                
        ts_fparams fpars;
        unsigned char ControlFlags;
        
        double (*act_function)(double, tps_fparams); // funcion de activacion
        double (*der_act_func)(double, tps_fparams); // derivada 
        void (*der_par_func)(double, tps_fparams); // derivada 
        
} ts_neurona;

typedef ts_neurona* tps_neurona;

// Capa 
typedef struct {
    t_neuronum nneuronas;
    tps_neurona n; // Puntero Neuronas
} ts_capa;

typedef ts_capa* tps_capa;

typedef struct {
    t_neuronum ncapas;
    tps_capa l;   
    double lrate;
    double reg;
} ts_red;

typedef ts_red* tps_red;

double ReLUSoftplus(double, tps_fparams );
double dReLUSoftplus(double, tps_fparams );
void ShowNeuralInfo(tps_red rnn);
#endif 
