#ifndef __MY_NEURAL__
#define __MY_NEURAL__
#include <stdint.h>
#include <math.h>
#include "kbhit.h"

typedef unsigned long int t_neuronum;
typedef struct 
{
    double (*act_function)(double); // funcion de activacion
    double (*der_act_func)(double); // derivada 
} ts_neurona_func;

typedef struct {
        t_neuronum ninputs; // Número de entradas
        double b; // bias;
        double *w; //coeficientes de ponderación (pesos)
        double **x; // entradas
        double v; 
        double o;   // salida
        double delta;
        double (*act_function)(double); // funcion de activacion
        double (*der_act_func)(double); // derivada 
        
        //ts_neurona_func v;
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

double ReLUSoftplus(double);
double dReLUSoftplus(double);
void ShowNeuralInfo(tps_red rnn);
#endif 