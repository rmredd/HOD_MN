#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define pi 3.14159265358979323846

double nu;
double func(double);
double func_gam(double);
double func_gauss(double);

int main(int argc, char **argv)
{
  int i,n;
  double x,f,qromo(),midinf(),fgam,midpnt(),fgauss;
  
  nu=2;
  if(argc>1)
    nu=atof(argv[1]);

  fgam = qromo(func_gam,0.0,0.1,midpnt) + qromo(func_gam,0.1,1.0E+30,midinf);

  for(i=1;i<=200;++i)
    {
      x=i/10.0;
      f=qromo(func,x,1.0E+30,midinf);
      fgauss = 2*qromo(func_gauss,0.0,x,midpnt);
      f*=1/(pow(2.0,nu/2.0)*fgam);
      printf("%f %e %e\n",x,1-f,fgauss);
    }
}

double func(double x)
{
  return(pow(x,nu/2-1)*exp(-x/2));
}
double func_gam(double x)
{
  return(pow(x,nu/2-1)*exp(-x));
}
double func_gauss(double x)
{
  return(exp(-x*x/2)/sqrt(2*pi));
}
