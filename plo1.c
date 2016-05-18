
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  
  FILE *gp1 = popen("gnuplot -persist","w");
  char * cmds[] = { "set term pdf", "set yrange [0:]","plot \"plot.txt\" with lines" };

  for (int i=0;i<3;i++)
    fprintf(gp1,"%s\n",cmds[i]);

  fclose(gp1);  

  return 0;
}
