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
/*
cursorup(n) CUU       Move cursor up n lines                 ^[[<n>A
cursordn(n) CUD       Move cursor down n lines               ^[[<n>B
cursorrt(n) CUF       Move cursor right n lines              ^[[<n>C
cursorlf(n) CUB       Move cursor left n lines               ^[[<n>D
 */
#define cursorup(x) printf("\x1b[%dA",(x))
#define cursordn(x) printf("\x1b[%dB",(x))
#define cursorrt(x) printf("\x1b[%dC",(x))
#define cursorlt(x) printf("\x1b[%dD",(x))
//nextline NEL          Move to next line                      ^[E
#define cursornwline()  printf("\x1bE")

#define clrscr()    printf ( "\x1b[2J" )
#define clrline()   printf ( "\x1b[2K" )
#define colorbg(x,y)  printf ( "\x1b[%d;%dm", (x) + 30, (y) + 40  )
#define color(x)  printf ( "\x1b[0;%dm", (x) + 30 )
#define attr(x)     printf ( "\x1b[%dm", (x) )
#define restorescreen() printf("\033[m")
#endif  
