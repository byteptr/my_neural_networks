#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fcntl.h>   /* File control */
#include <curses.h>
#include "myneural.h"
#include "caracteres_ansi.h"

double ReLUSoftplus(double x)
{
   return log(1.0+exp(x));   
}
double dReLUSoftplus(double x)
{
    return (1.0/(1.0+exp(-x)));
}
double dtanh(double x)
{
   return (1-tanh(x)*tanh(x));
}
double logistic(double x)
{
    return (1.0/(1.0+exp(-x)));
}
double dlogistic(double x)
{
    double y; 
    y = logistic(x);
    return y*(1-y);
}


void CreateDefaulNeural(tps_red rnn, 
                        t_neuronum Ninputs, // numero de entradas
                        t_neuronum Noutputs, // numero de salida (igual a neuronas en salida)
                        t_neuronum NNIlayer, // numero de neuronas en la entrada
                        t_neuronum NHlayers, // numero de capas ocultas
                        t_neuronum NHiddenNeuron) // numero de neuronas en las capas ocultas
{
    t_neuronum  i,j,k;
    
    rnn->l = (tps_capa)malloc(sizeof(ts_capa)*(NHlayers+2));
    if(rnn->l == NULL)
    {
        printf("Out of memory");
        exit(0);
    }
    rnn->ncapas = NHlayers+2;
    rnn->lrate = 1.0;
    
    rnn->l[0].n = (tps_neurona)malloc(sizeof(ts_neurona)*NNIlayer);
    rnn->l[0].nneuronas = NNIlayer;
    
    srand (time(NULL)); // Inicia generador de numeros pseudaleatorios 
    
    // Para todas las neuronas de la capa de entrada
    for(i = 0; i < NNIlayer; i++)
    {
        rnn->l[0].n[i].ninputs = Ninputs;
        rnn->l[0].n[i].w = (double *)malloc(sizeof(double)*Ninputs);
        rnn->l[0].n[i].x = (double **)malloc(sizeof(double*)*Ninputs);
        rnn->l[0].n[i].act_function = tanh;
        rnn->l[0].n[i].der_act_func = dtanh;
        rnn->l[0].n[i].b = 0.0;
        rnn->l[0].n[i].o = 0.0;
        rnn->l[0].n[i].v = 0.0;
        for(j = 0; j < Ninputs; j++)
        {
            rnn->l[0].n[i].w[j] = ((double)rand()/RAND_MAX-0.5)*2.0; 
            rnn->l[0].n[i].x[j] = NULL;
        }
        
        
    }
    
    // Capas ocultas
    for(j = 1; j < (NHlayers+1); j++)
    {
        rnn->l[j].nneuronas = NHiddenNeuron; // cuantas neuronas en esta capa 
        rnn->l[j].n = (tps_neurona)malloc(sizeof(ts_neurona)*NHiddenNeuron);
        
        // Para todas las neuronas de la capa intermedia
        for(i = 0; i < NHiddenNeuron; i++)
        {
            rnn->l[j].n[i].ninputs = rnn->l[j-1].nneuronas;
            
            rnn->l[j].n[i].act_function = tanh;
            rnn->l[j].n[i].der_act_func = dtanh;
            rnn->l[j].n[i].b = 0.0;
            rnn->l[j].n[i].o = 0.0;
            rnn->l[j].n[i].v = 0.0;
            rnn->l[j].n[i].w = (double *)malloc(sizeof(double)*rnn->l[j-1].nneuronas);            
            rnn->l[j].n[i].x = (double **)malloc(sizeof(double*)*rnn->l[j-1].nneuronas);
            
            //conectar las entradas con las salidas de las neuronas anteriores
            // e inicializar pesos 
            for(k = 0; k < rnn->l[j-1].nneuronas; k++)
            {
                rnn->l[j].n[i].w[k] = ((double)rand()/RAND_MAX-0.5)*2.0;
                rnn->l[j].n[i].x[k] = &(rnn->l[j-1].n[k].o);
            }
        }        
    }
   //Capa de salida 
   rnn->l[NHlayers+1].nneuronas = Noutputs;
   rnn->l[NHlayers+1].n = (tps_neurona)malloc(sizeof(ts_neurona)*Noutputs);
   
   for(i = 0; i < Noutputs; i++)
   {
       rnn->l[NHlayers+1].n[i].ninputs = rnn->l[NHlayers].nneuronas;       
       rnn->l[NHlayers+1].n[i].act_function = tanh;
       rnn->l[NHlayers+1].n[i].der_act_func = dtanh;
       rnn->l[NHlayers+1].n[i].w = (double *)malloc(sizeof(double)*rnn->l[NHlayers].nneuronas);            
       rnn->l[NHlayers+1].n[i].x = (double **)malloc(sizeof(double*)*rnn->l[NHlayers].nneuronas);       
       rnn->l[NHlayers+1].n[i].b = 0;
       rnn->l[NHlayers+1].n[i].o = 0;
       rnn->l[NHlayers+1].n[i].v = 0;
       for(k = 0; k < rnn->l[NHlayers].nneuronas; k++)
       {
           rnn->l[NHlayers+1].n[i].w[k] = ((double)rand()/RAND_MAX-0.5)*2.0;
           rnn->l[NHlayers+1].n[i].x[k] = &(rnn->l[NHlayers].n[k].o);
       }       
   }
}

void EvaluateNeural(tps_red rnn)
{
    t_neuronum i,j,k; 
    
    for(i = 0; i < rnn->ncapas; i++)
    {
        for(j = 0; j < rnn->l[i].nneuronas; j++)
        {
            rnn->l[i].n[j].v = 0.0;
            for(k = 0; k < rnn->l[i].n[j].ninputs; k++)
            {
                if(rnn->l[i].n[j].x[k] != NULL)
                    rnn->l[i].n[j].v  += rnn->l[i].n[j].w[k]*((double)*(rnn->l[i].n[j].x[k]));                
            }
            rnn->l[i].n[j].v += rnn->l[i].n[j].b;
            rnn->l[i].n[j].o = rnn->l[i].n[j].act_function(rnn->l[i].n[j].v);                 
        }
    }
}

// d es el vector de salidas a probar para una entrada dada
void BackPropagIteration(tps_red rnn, double *d, double *errsq)
{
    double err, dw,scale;
    t_neuronum i,j,k;
    
    //primera fase;
    EvaluateNeural(rnn);
    // primero evaluar la capa de salida

    err = 0.0;
    if(errsq != NULL)
        *errsq = 0;    
    // calculo de las deltas, capa de salida
    for(j = 0; j < rnn->l[rnn->ncapas-1].nneuronas; j++)
    {
        err = d[j]-rnn->l[rnn->ncapas-1].n[j].o;
        
        if(errsq != NULL)
            *errsq += err*err;
        rnn->l[rnn->ncapas-1].n[j].delta = err;
        
        rnn->l[rnn->ncapas-1].n[j].delta *= 
        rnn->l[rnn->ncapas-1].n[j].der_act_func(rnn->l[rnn->ncapas-1].n[j].v);
                
    }
    // calculo de las deltas, capas intrmedias   
    for(i = rnn->ncapas-1; i>0; i--) // direccionar con [i-1] unsigned long int. Macarra
    {
        for(j = 0; j < rnn->l[i-1].nneuronas; j++)
        {
           
            rnn->l[i-1].n[j].delta = 0.0;
            // neuronas de la capa posterior
            for(k = 0; k < rnn->l[i].nneuronas; k++) 
            {
                rnn->l[i-1].n[j].delta += rnn->l[i].n[k].delta * rnn->l[i].n[k].w[j]; // peso j a neurona j
            }                
                        
            rnn->l[i-1].n[j].delta *= rnn->l[i-1].n[j].der_act_func(rnn->l[i-1].n[j].v); 
        }
        
    }

    // ahora es cuando se actualizan los putos pesos, joder 
    
    for(j = 0; j < rnn->l[rnn->ncapas-1].nneuronas; j++)
    {
        scale = 0.0;
        for(k = 0; k < rnn->l[rnn->ncapas-1].n[j].ninputs; k++)
            scale += abs((double)*(rnn->l[rnn->ncapas-1].n[j].x[k]));
        scale *= scale;
        scale += 1.0;                

        for(k = 0; k < rnn->l[rnn->ncapas-1].n[j].ninputs; k++)
        {
            dw = (rnn->lrate)*rnn->l[rnn->ncapas-1].n[j].delta*
                 ((double)*(rnn->l[rnn->ncapas-1].n[j].x[k]));            
            rnn->l[rnn->ncapas-1].n[j].w[k] += dw/scale;
        }
       
        dw = (rnn->lrate)*rnn->l[rnn->ncapas-1].n[j].delta;  
        rnn->l[rnn->ncapas-1].n[j].b += dw/scale;        
    }    
    
    for(i = rnn->ncapas-1; i>0; i--) // direccionar con [i-1] unsigned long int. Macarra
    {
        for(j = 0; j < rnn->l[i-1].nneuronas; j++)
        {
            scale = 0.0;
            for(k = 0; k < rnn->l[i-1].n[j].ninputs; k++)
                scale += abs((double)*(rnn->l[i-1].n[j].x[k]));
            scale *= scale;
            scale += 1.0;

            //para todas las entradas de esa puta neurona 
            for(k = 0; k < rnn->l[i-1].n[j].ninputs; k++)
            {
                dw = (rnn->lrate)* rnn->l[i-1].n[j].delta * 
                     ((double)*(rnn->l[i-1].n[j].x[k]))/scale;
                rnn->l[i-1].n[j].w[k] += dw;
            }
            rnn->l[i-1].n[j].b += (rnn->lrate)* rnn->l[i-1].n[j].delta/scale;            
        }
    }
    
}

void ShowNeuralInfo(tps_red rnn)
{
    t_neuronum i,j,k;
    printf("\n");
    printf("Imprimiendo estructura de la red neuronal convolutiva\n");
    printf("Neurona(i,j): Neurona j-esima de capa i-esima (capa 0 = entrada)\n");
    printf("&x[n]-> Indica donde esta conectada la entrada de la neurona j-esima\n");
    printf("x[n] = Indica valor de la entrada de la neurona j-esima\n");
    printf("Salida: Dir-> , Valor = Indica punto de conexion con salida de neurona j-esima y su valor\n");
    
    if(rnn->l != NULL)
    {
        printf("Numero de capas: %d\n", rnn->ncapas);
        for(i = 0; i < rnn->ncapas; i++)
        {
            printf("Capa: %d, Neuronas: %d\n", i, rnn->l[i].nneuronas);
            
            for(j = 0; j < rnn->l[i].nneuronas; j++)
            {
                //rnn->l[i].n[j]
                printf("Neurona: (%d,%d), Entradas %d\n",i,j,rnn->l[i].n[j].ninputs);
                printf("Bias: %.15le\n", rnn->l[i].n[j].b);
                for(k = 0; k < rnn->l[i].n[j].ninputs; k++)
                {
                    if(rnn->l[i].n[j].x != NULL)
                    {
                        if(rnn->l[i].n[j].x[k] != NULL)
                        {
                            printf("&x[%d] -> 0x%04X",k,rnn->l[i].n[j].x[k]);
                            printf(", x[%d] = %.15le, ",k,*(rnn->l[i].n[j].x[k]));
                        } else 
                        {
                            printf("&x[%d] -> NULL, ",k);
                        }
                    } else 
                    {
                        printf("x[%d] -> NULL, ",k);
                    }
                    printf("w[%d] = %.15le\n",k,rnn->l[i].n[j].w[k]);
                }
                printf("Salida: Dir -> 0x%04X, Valor = %.15le\n",&(rnn->l[i].n[j].o),rnn->l[i].n[j].o);
                printf("____________________________________________\n");
                while(!kbhit()); //ungetch();
            }
            while(!kbhit()); //ungetch();
        }
        
    } else 
    {
        printf("No hay red\n");
    }
}
void DumpNeuralGraph(tps_red rnn, unsigned int NGraphs, FILE *f)
{
    //FILE *f;    
    t_neuronum i,j,k,r;
    unsigned int GraphCounter;
    
    //f = fopen("neur.dot","w");
    
    fprintf(f,"digraph G {\n");
    
    fprintf(f,"  rankdir=LR;\n");
    //fprintf(f,"size=\"8,5\"\n");
    //GraphCounter = 0;
    for(GraphCounter = 0; GraphCounter < NGraphs; GraphCounter++)
    {
        if(NGraphs > 1)
        {
            fprintf(f,"\tsubgraph cluster_%d {\n",GraphCounter);
            fprintf(f,"\tlabel = \"subnet %d of %d\"\n",GraphCounter+1,NGraphs);
        }
        for(k = 0; k < rnn->l[0].n[0].ninputs; k++)
        {   
            {
                fprintf(f,"\t\tx%d_%d [shape = doublecircle; label = \"x%d,%d\";];\n",GraphCounter,k,GraphCounter,k);            
            }       
        }   
    
        fprintf(f,"/*_____________________________*/\n");    
    
        for(i = 0; i < rnn->ncapas; i++)
        {
            for(j = 0; j < rnn->l[i].nneuronas; j++)
            {
                fprintf(f,"\t\tn%d_%d_%d [shape = circle; label = \"n%d,%d,%d\"];\n",GraphCounter,i,j,GraphCounter,i,j);
            }
            fprintf(f,"/******************************/\n");
        }
    
        for(i = 0; i < rnn->ncapas; i++)
        {
            if(i == 0)
            {
                for(k = 0; k < rnn->l[0].n[0].ninputs; k++)
                {
                    // cuidado aqui con la j
                    for(j = 0; j < rnn->l[i+1].n[0].ninputs; j++ )
                    {
                        fprintf(f,"\t\tx%d_%d->n%d_%d_%d;\n",GraphCounter,k,GraphCounter,i,j);
                    }
                }
            } else  
            {
                for(j = 0; j < rnn->l[i].nneuronas; j++)
                {
                    for(k = 0; k < rnn->l[i].n[j].ninputs; k++)
                    {
                        fprintf(f,"\t\tn%d_%d_%d->n%d_%d_%d;\n",GraphCounter,i-1,k,GraphCounter, i,j);
                    }
                }
            }
        }
        if(NGraphs > 1)
        {
            fprintf(f,"\t}\n");
            rnn++; // incrementamos el puntero
        }
    }
    
    fprintf(f,"}\n");    
}

DumpNeuralNet(tps_red rnn)
{
    t_neuronum i,j,k;
    FILE *f;    
    f = fopen("neur.rnn","w");
    fprintf(f,"%d\n",rnn->ncapas-2); // cuantas capas ocultas
    fprintf(f,"%d\n",rnn->l[0].nneuronas); // cuantas neuronas en la entrada
    fprintf(f,"%d\n",rnn->l[0].n[0].ninputs); // cuantas entradas
    fprintf(f,"%d\n",rnn->l[rnn->ncapas-1].nneuronas); // cuantas salidas
    
    for(i = 0; i<rnn->ncapas; i++)
    {
        for(j = 0; j<rnn->l[i].nneuronas; j++)
        {
            for(k = 0; k < rnn->l[i].n[j].ninputs; k++)
            {
                fprintf(f,"%.15le\n",rnn->l[i].n[j].w[k]);
            }
        }
    }
    fclose(f);
    
}
int main(void)
{
#define ORDEN 3
#define NETS 3
    double random_value;
    unsigned long int k,j,r,kk, epoch_counter;
    double input[ORDEN];
    double err[NETS],errm[NETS], errp[NETS], errt, phi;
    //double itestv[4][2] = {{0.0,0.0},{0.0,1.0},{1.0,0.0},{1.0,1.0}};
    double otestv[NETS];// deberia ser doble, nets y nsalidas
    FILE *f,*fres, *neural_graph, *pipe[NETS]; 
    // ficheros 
    neural_graph    = fopen("neur.dot","w"); // grafico de la red
    f               = fopen("neur.txt","w"); // pesos y estructura
    fres            = fopen("resul.txt","w"); // resultados
    
    ts_red neural[NETS]; 
    
    for(k = 0; k < NETS; k++)
    {
        if(k == 0)
        CreateDefaulNeural(&neural[k], 
                            ORDEN, // numero de entradas
                            1, // numero de salidas (igual a neuronas en salida)
                            3, // numero de neuronas en la entrada
                            1, // numero de capas ocultas
                            11); // numero de neuronas en las capas ocultas    
        else 
        CreateDefaulNeural(&neural[k], 
                            ORDEN+1, // numero de entradas
                            1, // numero de salidas (igual a neuronas en salida)
                            3, // numero de neuronas en la entrada
                            1, // numero de capas ocultas
                            3); // numero de neuronas en las capas ocultas                            
    }
        
    for(k = 0; k < NETS; k++)
    {
        for(j = 0; j < neural[k].l[0].nneuronas; j++)
        {
            for(kk = 0; kk < ORDEN; kk++)
                neural[k].l[0].n[j].x[kk] = &input[kk]; 
        }
        if(k > 0)
        {
            for(j = 0; j < neural[k].l[0].nneuronas; j++)
                neural[k].l[0].n[j].x[ORDEN] = &errp[k-1];
                /* neural[k-1].l[neural[k-1].ncapas-1].n[0].o-*/
        }
    }
    
    
    j = 1;
    
    for(k = 0; k < NETS; k++)
    {
        errm[k] = 0.0;
        errp[k] = 0.0;
        err[k]  = 0.0;
    }
    clrscr();
    neural[0].lrate = 0.5;
    DumpNeuralGraph(&neural[0],NETS, neural_graph);    
    fclose(neural_graph);
    system("dot neur.dot -Tps -o neur.ps");
    //for (k = 0; k < 200; k++)    
#if 1    
k = 0;
epoch_counter = 0;

    while((!kbhit())&&(epoch_counter < 50000))
    {
        
        phi = ((rand()%500)/500.0);
        input[0] = phi;
        
        for(kk = 1; kk < ORDEN; kk++)
            input[kk] = input[kk-1]*input[0];

        otestv[0] = (1+cos(2*M_PI*phi))/2;//+ 0.01*((rand()/((double)RAND_MAX))-0.5)*2; 
        
        errt = otestv[0];
        
        for(kk = 0; kk < NETS; kk++)
        {
            for(r = 0; r < 8; r++)
            {
                EvaluateNeural(&neural[kk]);
                
                if(kk == 0)
                {
                    BackPropagIteration(&neural[kk],(double *)&otestv[0], &err[kk]);
                } else 
                {
                    BackPropagIteration(&neural[kk],(double *)&errp[kk-1], &err[kk]);
                }
                
                if(kk == 0)
                {
                    errp[kk] = (otestv[0]-neural[kk].l[neural[kk].ncapas-1].n[0].o);
                } else 
                {
                    errp[kk] = (errp[kk-1]-neural[kk].l[neural[kk].ncapas-1].n[0].o);
                }
                                
            }
            
            //EvaluateNeural(&neural[kk]);
            errm[kk] += err[kk];
            
            if(kk == 0)
                errt = otestv[0]-neural[kk].l[neural[kk].ncapas-1].n[0].o;
            else 
                errt += pow(-1,kk)*neural[kk].l[neural[kk].ncapas-1].n[0].o;
        }

         if(j == 0)
         {
             k++;
             printf("      f: %.15le\n",kk,otestv[0]);
             for(kk = 0; kk < NETS; kk++)
                printf(" Out[%d]: %.15le\n",kk,neural[kk].l[neural[kk].ncapas-1].n[0].o);
             {
                 double p;
                 p = 0.0;
                for(kk = 0; kk < NETS; kk++)
                {   if(kk > 0)               
                        p += neural[kk].l[neural[kk].ncapas-1].n[0].o/10.0;
                    else 
                        p += neural[kk].l[neural[kk].ncapas-1].n[0].o;
                }
                printf("SumOuts: %.15le\n",p);
                printf("DifErrs: %.15le\n",p-otestv[0]);
             }
             for(kk = 0; kk < NETS; kk++)                         
                printf("Errp[%d]: %.15le\n",kk, errp[kk]);             
             printf("Error t: %.15le\n",errt);
             for(kk = 0; kk < NETS; kk++)
                printf("Errm[%d]: %.15le  \t epoch %d\n",kk, sqrt(errm[kk]/500), epoch_counter);
             
             printf("\n");
             
             if(k >= 100)
             {                
                k = 0;
                epoch_counter++;
             }
             
             for(kk = 0; kk < NETS; kk++)
                errm[kk]=0.0;
             
         }
        j++; j%=500;
    }
    fclose(f);
    //DumpNeuralNet(&neural[0]);
    
    //if(epoch_counter < 10000)
      //  ShowNeuralInfo(&neural[0]);
    
    for(j = 0; j < 1000; j++)
    {
        
        input[0] = j/1000.0;
        for(kk = 1; kk < ORDEN; kk++)
            input[kk] = input[kk-1]*input[0];        
        
        for(kk = 1; kk < NETS; kk++)
        {
            EvaluateNeural(&neural[kk]);
                if(kk == 0)
                {
                    errp[kk] = otestv[0]-neural[kk].l[neural[kk].ncapas-1].n[0].o;
                } else 
                {
                    errp[kk] = (errp[kk-1]-neural[kk].l[neural[kk].ncapas-1].n[0].o);
                }
        }
        
        for(k = 0; k < NETS; k++)
        {
            fprintf(fres,"%.15le",errp[k]);     
            if(k < (NETS-1))
                fprintf(fres,"\t");     
        }
        fprintf(fres,"\n");
    }
    fclose(fres);
#endif     
    printf("\n");
    return 0;
}