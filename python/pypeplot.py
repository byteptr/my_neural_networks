#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 20:01:33 2018
V 1.0
@author: byteptr 
"""
import os 
import sys
import ast
import numpy as np
import matplotlib 
from matplotlib import pyplot as plt

from enum import Enum
from pathlib import Path

path = "/tmp/pypeplot.fifo.0"
ftest = Path(path); 
if ftest.exists(): 
    print("fichero detectado, borrar")
    os.remove(path) 
else:
    print("fichero NO detectado");


PRINT_DEBUG = False

CMD_QUIT = "quit\n"
CMD_NPOINTS = "npoints\n"
#CMD_ROLL_MODE = "rollmode\n"
#CMD_CLEAR = "clear\n"
CMD_RESET = "reset\n"
CMD_NONE = ''; 

cmdlist = [CMD_NPOINTS, CMD_QUIT, CMD_RESET, CMD_NONE]

class status(Enum):
    INIT=0 
    WAIT = 1
    GETNPOINTS = 2
    GETPOINTS = 3
    QUIT = 4

s = {'status': status.INIT, 'recv' : '', 'N': 0,'x': [], 'y': [],'n': 0}    

def TreatMachine(r):
    r0 = r.copy();
    # Intenta leer el estado, si hubiese algun error de diccionario se reinicia la maquina

    #este es el case propiamente dicho
        
    #estado INICIAL
    if r['status'] == status.INIT :
        r = {'status': status.WAIT, 'recv' : '','N': 0,'x': [],'y': [],'n': 0}
        
    #estado ESPERA
    elif r['status'] == status.WAIT :
        if r['recv'] == CMD_NPOINTS :
            r['status'] = status.GETNPOINTS
        elif r['recv'] == CMD_QUIT :
            r['status'] = status.QUIT
    #estado OBTENER NUMERO DE PUNTOS
    
    elif r['status'] == status.GETNPOINTS :
        if not r['recv'] in cmdlist:
            r['status'] = status.GETPOINTS;
            if PRINT_DEBUG:
                print('RECIBIDO: ',r['recv'])
            r['N'] = ast.literal_eval(r['recv']);
            r['x'] = np.linspace(0,r['N']-1,r['N'],dtype=float);            
            r['y'] = np.zeros((r['N'],)) 
            plt.ion()
        else :
            if r['recv'] == CMD_QUIT :
                r['status'] = status.QUIT
            elif r['recv'] == CMD_RESET :
                r = {'status': status.WAIT, 'recv' : '','N': 0,'x': [],'y': [],'n': 0}         
        return r    
    #estado OBTENER PUNTOS PARA GRAFICAR 
    elif r['status'] == status.GETPOINTS :
            
        r['status'] = status.GETPOINTS
        if not r['recv'] in cmdlist:

            num = ast.literal_eval(r['recv']);
                    
            if r['n'] < r['N']:
                #r['y'] = np.append(r['y'],num); # YEEEEEEEEEEEEEEEEEEEEEEEAAAAAHHHH MADAFAKAAAAAA
                r['y'][r['n']] = num;
                r['n'] += 1
            else :
                r['y'] = np.roll(r['y'],-1);
                r['y'][r['N']-1] = num
                r['x'] = np.roll(r['x'],-1);
                r['x'][r['N']-1] = r['x'][r['N']-2]+1;

            if(r['n'] == 1):                
                r['fig'], r['axes'] = plt.subplots(nrows=1, ncols=1);
                r['axes'].grid(True);
                r['axes'].set_xlim((0,r['n']))                
                r['line'], = r['axes'].plot(r['x'][0:r['n']],r['y'][0:r['n']]);                                
                r['ymin'] = np.min(r['y'])
                r['ymax'] = np.max(r['y'])                
            elif(r['n'] > 1):
                
                r['change'] = True; 
                
                if r['ymin'] > np.min(r['y'][0:r['n']]):
                    r['ymin'] = np.min(r['y'][0:r['n']])
                    r['change'] = True                    
                if r['ymax'] < np.max(r['y'][0:r['n']]):
                    r['ymax'] = np.max(r['y'][0:r['n']])
                    r['change'] = True                
                if(r['change']):
                    r['axes'].set_ylim(r['ymin'],r['ymax'])
                if(r['n'] < r['N']):
                    r['axes'].set_xlim((0,r['n']))                    
                else :
                    r['axes'].set_xlim((r['x'][0],r['x'][r['N']-1]))     
                    
                r['line'].set_ydata(r['y'][0:r['n']])
                r['line'].set_xdata(r['x'][0:r['n']])
                r['fig'].canvas.draw();                                                 
        else :
            if r['recv'] == CMD_QUIT :
                r['status'] = status.QUIT
            elif r['recv'] == CMD_RESET :
                r = {'status': status.WAIT, 'recv' : '','N': 0,'x': [],'y': [],'n': 0} 
    #estado de SALIDA
    elif r['status'] == status.QUIT :
        pass
    else :
        if PRINT_DEBUG:
            print ("SE VIENE AQUI r:", r);
        #r['status'] = status.INIT

    if PRINT_DEBUG:    
        print('Old status:',r0, '\nNew status:',r,'\n');
    
    return r


os.mkfifo(path)
fifo = open(path,"r")
s = TreatMachine(s)

while not s['status'] == status.QUIT:    
    for s['recv'] in fifo:
        s = TreatMachine(s)
input("Esto ya ha terminado. Pulsa una tecla para continuar...");        

if ftest.exists(): 
    fifo.close()
    os.remove(path)
