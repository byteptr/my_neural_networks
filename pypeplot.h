#ifndef __PYPEPLOT_C__

#define __PYPEPLOT_C__ 
#define PYPEPLOT_VERSION "0.1b"
#define PYPEPLOT_FIFOPATH "/tmp/pypeplot.fifo.0"

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>

typedef enum {PYPE_SUNDEFINED=0, PYPE_SOK=1, PYPE_SERROR = -1} te_PypeOpStatus;
typedef struct s_pypeplot{
    unsigned long int PointsWindow;     
    unsigned long int TotalPoints;
    unsigned long int counter; 
    FILE * stream; 
    te_PypeOpStatus OpStatus;     
} ts_pypeplot;
typedef ts_pypeplot * tps_pypeplot;

void PypePlotOpenStream(tps_pypeplot s);
void PypePlotCloseStream(tps_pypeplot s);
void PypePlotReset(tps_pypeplot s);
void PypePlotNumPoints(tps_pypeplot s);
void PypePlotNewPoint(tps_pypeplot s, double y);
void PypePlotQuit(tps_pypeplot s);

#endif 
