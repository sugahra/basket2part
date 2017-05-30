#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include <oxdraw.h>
#import<solvenle>
#import<maximize>
#include"myconst_indep.h"

static decl s_vD,s_mX;


fn_Logit_j(const vgamma, const Func, const Score, const Hess)
{
  decl vexp,vlikelihood;
  vlikelihood=zeros(c_iN,1);
  vexp=exp(s_mX*vgamma);
  vlikelihood=s_vD.*(vexp./(ones(c_iN,1)+vexp))
    +(1-s_vD).*(ones(c_iN,1)./(ones(c_iN,1)+vexp));    
  Func[0]=sumc(log(vlikelihood));
  return 1;
}

fn_Reg(const vY,const mX, const fl_onlyest)
{
  decl iN,vbeta,dsigma2,mvar_beta,vsd_beta,vp_beta,mest;
  iN=sizer(vY);
  vbeta=invert(mX'*mX)*mX'*vY;
  dsigma2=sumsqrc(vY-mX*vbeta)/iN;
  if(fl_onlyest==1)
    {
      return vbeta|dsigma2;
    }
  else
    {
      mvar_beta=dsigma2*invert(mX'*mX);
      vsd_beta=sqrt(diagonal(mvar_beta))';
      vp_beta=tailn(fabs(vbeta./vsd_beta));
      mest=zeros(c_iP+1,3);
      mest[0:c_iP-1][0]=vbeta;
      mest[0:c_iP-1][1]=vsd_beta;
      mest[0:c_iP-1][2]=vp_beta;
      mest[c_iP][0]=dsigma2;
      return mest;
    }
}

main()
{
  decl time,file,mdata,mZ,mD,mY,vinsure_drug, vinsure_dental,
    mbeta,vsigma2,cj,maxi_gamma,dopt,mgamma,msd_gamma,msd_beta,
    mp_beta,mp_gamma,mhess,mvar,mest_reg,dloglikelihood,dAIC,dBIC,
    mloglikelihood,dexp,dlogphi,vgamma,vY,iJ,cij,ci,mX;
  time=timer();
  //Load data
  file=fopen("data/mdata.dat","r");
  fscan(file,"%#m",c_iN,c_ivar_file,&mdata);
  fclose(file);

  mZ=mdata[][0:c_iJ-1];
  s_mX=ones(c_iN,1)~mdata[][c_iJ:c_iJ+c_iP-2];
  vinsure_dental=mdata[][c_iJ+c_iP-3];
  vinsure_drug=mdata[][c_iJ+c_iP-2];

  mD=(mZ.>0);
  mY=log(mZ).*mD;
  mgamma=mbeta=msd_gamma=msd_beta=mp_beta=mp_gamma=zeros(c_iP,c_iJ);
  vsigma2=zeros(c_iP,1);

  //Estimation
  for(cj=0;cj<c_iJ;cj++)
    {
      s_vD=mD[][cj];
      //Logit
      ///Initial values
      vgamma=0.1*ones(c_iP,1);
      maxi_gamma=MaxBFGS(fn_Logit_j, &vgamma, &dopt, 0,TRUE);
      Num2Derivative(fn_Logit_j,vgamma,&mhess);
      mvar=invert(-mhess);
      mgamma[][cj]=vgamma;
      msd_gamma[][cj]=sqrt(diagonal(mvar))';
      mp_gamma[][cj]=tailn(fabs(mgamma[][cj]./msd_gamma[][cj]));

      //Regression
      iJ=sumc(s_vD);
      vY=zeros(iJ,1);
      mX=zeros(iJ,c_iP);
      cij=0;
      for(ci=0;ci<c_iN;ci++)
	{
	  if(s_vD[ci]==1)
	    {
	      vY[cij]=mY[ci][cj];
	      mX[cij][]=s_mX[ci][];
	      cij++;
	    }
	}
      mest_reg=fn_Reg(vY,mX,0);
      //Last argument==0 if we need inference results
      mbeta[][cj]=mest_reg[0:c_iP-1][0];
      vsigma2[cj]=mest_reg[c_iP][0];
      msd_beta[][cj]=mest_reg[0:c_iP-1][1];
      mp_beta[][cj]=mest_reg[0:c_iP-1][2];
    }
  println("time lapsed=",timespan(time));

  //Information criterion

  mloglikelihood=zeros(c_iN,c_iJ);
  for(ci=0;ci<c_iN;ci++)
    {
      for(cj=0;cj<c_iJ;cj++)
	{
	  dexp=exp(s_mX[ci][]*mgamma[][cj]);
	  if(mD[ci][cj]==1)
	    {
	      dlogphi=-0.5*log(M_2PI*vsigma2[cj])
		-sumsqrc(mY[ci][cj]-s_mX[ci][]*mbeta[][cj])/(2*vsigma2[cj]);
	      mloglikelihood[ci][cj]=log(dexp)-log(1+dexp)+dlogphi;
	    }
	  else
	    {
	      mloglikelihood[ci][cj]=-log(1+dexp);
	    }
	}
    }
  dloglikelihood=sumc(sumr(mloglikelihood));
  dAIC=-2*dloglikelihood+2*c_iq_indep;
  dBIC=-2*dloglikelihood+c_iq_indep*log(c_iN);

  //Report estimation results
  for(cj=0;cj<c_iJ;cj++)
    {
      println("Gamma",cj);
      println(mgamma[][cj]~ msd_gamma[][cj]~mp_gamma[][cj]);
      println("Beta",cj, " Sample size=",sumc(mD[][cj]));
      println(mbeta[][cj]~ msd_beta[][cj]~mp_beta[][cj]);
      println("sigma^{2}",cj);
      println(vsigma2[cj]);
    }
  println("LogL, AIC, BIC");
  println(dloglikelihood~dAIC~dBIC);

  format(2000);
  file=fopen("results/mle_indep.dat","w");
  for(cj=0;cj<c_iJ;cj++)
    {
      fprintln(file,"Gamma",cj);
      fprintln(file,"%15.5f",mgamma[][cj]~ msd_gamma[][cj]~mp_gamma[][cj]);
      fprintln(file,"Beta",cj, " Sample size=",sumc(mD[][cj]));
      fprintln(file,"%15.5f",mbeta[][cj]~ msd_beta[][cj]~mp_beta[][cj]);
      fprintln(file,"sigma^{2}",cj);
      fprintln(file,"%15.5f",vsigma2[cj]);
    }
  fprintln(file,"LogL, AIC, BIC");
  fprintln(file,"%15.5f",dloglikelihood~dAIC~dBIC);
  fclose(file);
  format(2000);
  file=fopen("results/gamma_indep.dat","w");
  fprintln(file,"%15.5f",mgamma);
  fclose(file);  
}
