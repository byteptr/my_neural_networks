#include <stdio.h>
#include <termios.h>
#include <unistd.h>
#include <fcntl.h>
#include "kbhit.h"

int kbhit(void)
{
  struct termios oldt, newt;
  int ch;
  int oldf;

  tcgetattr(STDIN_FILENO, &oldt);
  newt = oldt;
  newt.c_lflag &= ~(ICANON | ECHO);
  tcsetattr(STDIN_FILENO, TCSANOW, &newt);
  oldf = fcntl(STDIN_FILENO, F_GETFL, 0);
  fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);

  ch = getchar();

  tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
  fcntl(STDIN_FILENO, F_SETFL, oldf);
  
    //fseek(stdin,0,SEEK_END);
    
  if(ch != EOF)
  {
    //ungetc(ch, stdin); // opcional
    return 1;
  }

  return 0;
}
