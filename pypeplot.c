#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include "pypeplot.h"

void PypePlotReset(tps_pypeplot s)
{
    s->PointsWindow = 0;
    s->TotalPoints = 0;
    s->counter = 0; 
    fprintf (s->stream, "reset\n");
    if (ferror (s->stream))    
        s->OpStatus = PYPE_SERROR;
    else     
        s->OpStatus = PYPE_SOK;
}

void PypePlotNumPoints(tps_pypeplot s)
{    
    fprintf(s->stream, "npoints\n");
    fprintf(s->stream, "%d\n",s->PointsWindow);
    if (ferror (s->stream))    
        s->OpStatus = PYPE_SERROR;
    else     
        s->OpStatus = PYPE_SOK;
}

void PypePlotNewPoint(tps_pypeplot s, double y)
{
    fprintf(s->stream, "%g\n",y);
    if (ferror (s->stream))    
        s->OpStatus = PYPE_SERROR;
    else     
        s->OpStatus = PYPE_SOK;
}


void PypePlotOpenStream(tps_pypeplot s)
{
    s->stream = fopen((const char *)PYPEPLOT_FIFOPATH, "w+"); // ,0600);
    if (!s->stream )
        s->OpStatus = PYPE_SERROR;
    else     
        s->OpStatus = PYPE_SOK;
}

void PypePlotCloseStream(tps_pypeplot s)
{
    
    if (fclose(s->stream) != 0)
        s->OpStatus = PYPE_SERROR;
    else     
        s->OpStatus = PYPE_SOK;
}

void PypePlotQuit(tps_pypeplot s)
{
    fprintf (s->stream, "quit\n");
    if (ferror (s->stream))    
        s->OpStatus = PYPE_SERROR;
    else     
        s->OpStatus = PYPE_SOK;    
}
