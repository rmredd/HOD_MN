  if(Work.iquad && Work.iquad_covar)
    {
      calc_rquad(Work.r_quad,rquad,Work.n_quad);
      for(i=0;i<Work.n_quad;++i)
	for(j=0;j<Work.n_quad;++j)
	  {
	    if(Work.r_quad[i]<Work.rhlo || Work.r_quad[j]<Work.rhlo)continue;
	    e=esys_quad(Work.r_quad[i])*rquad[i]*esys_quad(Work.r_quad[j]);
	    chi2c+=(rquad[i]-Work.data_h[i])*(rquad[j]-Work.data_h[j])*
	      (Work.covar_q[i][j]);
	    if(!ThisTask && OUTPUTZ)
	      printf("CHIQUAD%d %d %d %d %f %e %e %e %e %e %e\n",Work.imodel,iter,i,j,Work.r_quad[j],
		     rquad[i],Work.data_h[i],rquad[j],Work.data_h[j],Work.covar_h[i][j],e);
	  }
    }
