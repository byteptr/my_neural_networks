#ifndef __CARACTERES__ANSI__
#define __CARACTERES__ANSI__

/* Necesito la librería estándar para usar printf, pero de todas formas esto queda relegado al usuario */
/* #include <stdio.h> */

/* Defino los colores que puedo utilizar */
#define NEGRO    0
#define ROJO     1
#define VERDE    2
#define AMARILLO 3
#define AZUL     4
#define MAGENTA  5
#define CYAN     6
#define BLANCO   7

/* Defino algunos de los atributos posibles */
#define NEGRITA   1
#define SUBRAYADO 4
#define PARPADEO  5
#define INVERSO   7

/* Las macros */
/* gotoxy(x,y): Toma dos argumentos, siendo 'x' el nro. de columnas [ 1 ... 80 ] e 'y' el nro. de filas [ 1 ... 25 ] */
/* clrscr(): Borra la pantalla */
/* color(x,y): Toma dos argumentos, siendo 'x' el color de primer plano e 'y' el color de fondo */
/* atrr(x): Toma un único argumento que representa el atributo a establecer */
#define gotoxy(x,y) printf ( "\x1b[%d;%dH", (y), (x) )
#define clrscr()    printf ( "\x1b[2J" )
#define clrline()   printf ( "\x1b[2K" )
#define color(x,y)  printf ( "\x1b[%d;%dm", (x) + 30, (y) + 40  )
#define attr(x)     printf ( "\x1b[%dm", (x) )

#endif  