#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/* This calculates and tabulates both the linear and non-linear
 * matter correlation function.
 *
 * If the parameter BOX_SIZE is set, then the lower limits of the
 * Fourier transform in 1/BOX_SIZE, else it starts at k=0.
 *
 * The integral (Smith etal 2003 MNRAS.341.1311S  Eq [4]) to transform:
 *
 *   Int_0^\infty Delta(k) sin(rk)/rk/k dk
 * 
 * This integral sometimes gives qromo problems when calculating xi(r) at large r.
 * Therefore, we cut off the integral at k = 10^3, which has negligible effect on the
 * correlation function at the relevent scales.
 */

double r_g4;
double xi_int(double xk);
double tabulated_power_spectrum(double xk);
void discrete_fourier_transform(void);

/* Calculates and tabulates the non-linear matter correlation function.
 * Since this is only really needed for the scale-dependence of the bias,
 * which is basically 1 at scales r>~8 Mpc/h, I won't calculate this much
 * past that value. 
 */
double xi_interp(double r)
{
  static int flag=0,prev_cosmo=0;
  static double *x,*y,*y2;
  int n=30,i,j;
  double a,xi_int(),rhi=95,rlo=0.1,dlogr,klo,s1,s2,tolerance=1.0e-6;

  if(!flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      // discrete_fourier_transform();

      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=1;
      dlogr = (log(rhi)-log(rlo))/(n-1);

      for(i=1;i<=n;++i)
	{
	  klo = 0;
	  if(BOX_SIZE>0)klo = TWOPI/BOX_SIZE;
	  r_g4 = x[i] = exp((i-1)*dlogr)*rlo;
	  
	  j=1;
	  s1 = qromo(xi_int,klo,j/r_g4,midpnt);
	  s2 = s1;
	  klo = j/r_g4;
	  while(mabs(s1)>tolerance*mabs(s2)) {
	    j+=16;
	    s1 = qromo(xi_int,klo,j/r_g4,midpnt);
	    s2 += s1;
	    klo = j/r_g4;
	  }
	  y[i]=s2;
	  //fprintf(stderr,"%f %f\n",r_g4,s2);
	  
	}
      check_for_smoothness(x,y,n,0);
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      prev_cosmo=RESET_COSMOLOGY;
    }

  splint(x,y,y2,n,r,&a);
  return(a);
}


/* This is the integrand of Smith et al Eq. [4]
 */
double xi_int(double xk)
{
  double xk1,xk2,psp;
    
  if(xk==0)return(0);

  /* power spectrum at xk
   */
  psp=nonlinear_power_spectrum(xk);
  //psp = tabulated_power_spectrum(xk)*(xk*xk*xk)/TWOPI/PI;

  /* Integrand of Fourier transform
   */
  xk1=r_g4*xk;
  psp*=sin(xk1)/xk1/xk;
  return psp;
}


double xi_linear_interp(double r)
{
  static int flag=0,prev_cosmo=0;
  static double *x,*y,*y2;
  int n=100,i;
  double a,rlo=0.1,rhi=150,dlogr,klo;
  double xi_linear_int();

  if(!flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=1;

      dlogr = (log(rhi)-log(rlo))/(n-1);
      klo = 0;
      if(BOX_SIZE>0)klo = 1/BOX_SIZE;
      for(i=1;i<=n;++i)
	{
	  r_g4 = x[i] = exp((i-1)*dlogr)*rlo;
	  y[i] = qromo(xi_linear_int,klo,1.0/r_g4,midpnt)+
	    qromo(xi_linear_int,1.0/r_g4,1.0E+3,midpnt);
	}
      check_for_smoothness(x,y,n,0);
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      prev_cosmo=RESET_COSMOLOGY;
    }

  splint(x,y,y2,n,r,&a);
  return(a);
}

double xi_linear_int(xk)
double xk;
{
  double psp;
  double xk1,xk2;
    
  /* power spectrum at xk
   */
  psp=linear_power_spectrum(xk);

  /* Now take Fourier transform
   */
  xk1=r_g4*xk;
  psp*=sin(xk1)/xk1/xk;

  return psp;
}

/* Sometimes one or two of the correlation function values will
 * totally crap out, [a side-effect of qromo] so here we check for that and if
 * we find it then we interpolate the values in log-space.
 */
void check_for_smoothness(double *x, double *y, int n, double r)
{
  int i,flag,flag2=0;
  double m,b,new;

  for(i=2;i<n;++i)
    {
      flag=0;
      if(y[i]>0)flag2=1;
      if(y[i]<0 && !flag2)continue;
      if(x[i]<r)continue;
      if(fabs(y[i]/y[i-1])>3.0 &&  fabs(y[i]/y[i+1])>3.0)flag=1; 
      if(fabs(y[i]/y[i-1])<0.2 && fabs(y[i]/y[i+1])<0.2)flag=1; 
      if(y[i]<0 && (y[i-1]>=0 && y[i+1]>=0))flag=1;
      if(y[i+1]<0)flag=0;
      if(!flag)continue;
      
      m=(log(y[i+1])-log(y[i-1]))/(log(x[i+1])-log(x[i-1]));
      b=log(y[i+1])-m*log(x[i+1]);
      new=m*log(x[i])+b;
      
      fprintf(stderr,"SMOOTHING: %e %e %e r= %f new=%e\n",y[i-1],y[i],y[i+1],x[i],exp(new));
      y[i]=exp(new);
    }


}

double tabulated_power_spectrum(double xk)
{
  static double *x,*y,*y2;
  static int flag=1,n;
  int i;
  double t,x0;
  FILE *fp;
  char a[1000];
  float x1,x2;

  if(flag)
    {
      flag=0;

      fp=openfile("tmp.particles");
      n=filesize(fp);

      x=dvector(1,n);
      y=dvector(1,n);
      y2=dvector(1,n);
      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f",&x1,&x2);
	  x[i]=x1;
	  y[i]=x2;
	  x[i]=log(x[i]);
	  y[i]=log(y[i]);
	}
      fclose(fp);
      spline(x,y,n,1.0E+30,1.0E+30,y2);
    }
  xk=log(xk);
  splint(x,y,y2,n,xk,&t);
  return(exp(t));

}

void discrete_fourier_transform(void)
{
  static double *x,*y,*y2;
  static int flag=1,n;
  int i,j;
  double t,x0;
  FILE *fp;
  char a[1000];
  float x1,x2;

  double rlo = 0.1, rhi = 50, dlogr, r, xi, xkr;
  int nr=50;

  if(flag)
    {
      flag=0;

      fp=openfile("pk.smith");
      n=filesize(fp);

      x=dvector(1,n);
      y=dvector(1,n);
      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f",&x1,&x2);
	  x[i]=x1;
	  y[i]=x2;
	}
      fclose(fp);
    }

  dlogr = (log(rhi)-log(rlo))/(nr-1);

  n = 64;
  nr = n;
  for(i=1;i<=nr;++i)
    {
      r = exp(dlogr*(i-1))*rlo;
      r = TWOPI/x[i];
      xi = 0;
      for(j=1;j<=n;++j)
	{
	  xkr = r*x[j];
	  xi += y[j]/TWOPI/PI*sin(xkr)/xkr*x[j]*x[j];
	}
      xi/=n;
      printf("DFT %f %f\n",r,xi);
    }
}

