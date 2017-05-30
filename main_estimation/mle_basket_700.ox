#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include <oxdraw.h>
#import<solvenle>
#import<maximize>
#include"myconst_basket_700.h"

static decl s_mX,s_ml,  //For multinomial logit for m
  s_mc_m, s_mX_m, //For k, case 2 and 3
  s_vmulti_m; //For k, only case 3
fn_Likelihood_k(const vY, const vX, const vbeta, const mSigma)
{
  decl im,dphi,mbeta;
  im=sizec(mSigma);
  if(im==1)
    {
      dphi=(1/sqrt(mSigma))
	*densn((1/sqrt(mSigma))*(vY-vbeta'*vX));
    }
  if(im!=1)
    {
      mbeta=reshape(vbeta,im,c_iP)';
      dphi=(1/sqrt(determinant(mSigma)))
	*prodc(densn( invert(choleski(mSigma))* (vY-mbeta'*vX)));
    }
  return dphi;
}

fn_Mlogit_m(const vgamma_l, const Func, const Score, const Hess)
{
  decl mexp,vlikelihood,mgamma;
  mgamma=reshape(vgamma_l,c_iM,c_iP)';
  mexp=exp(s_mX*mgamma);
  mexp=ones(c_iN,1)~mexp;
  vlikelihood=sumr(mexp.*s_ml)./sumr(mexp);
  Func[0]=sumc(log(vlikelihood));
  return 1;
}

fn_K_case2(const vgamma_km2, const Func, const Score, const Hess)
{
  decl mexp,vlikelihood,mgamma,iK_m,iN_m;
  iN_m=sizer(s_mX_m);
  iK_m=sizec(s_mc_m);
  mgamma=reshape(vgamma_km2,iK_m-1,c_iP)';
  mexp=ones(iN_m,1)~exp(s_mX_m*mgamma);
  vlikelihood=sumr(mexp.*s_mc_m)./sumr(mexp);
  Func[0]=sumc(log(vlikelihood));
  return 1;
}

fn_K_case3(const vgamma_km3, const Func, const Score, const Hess)
{
  decl mgamma,mexp,mlikelihood,ci,ck,iK_m,iN_m;
  iN_m=sizer(s_mX_m);
  iK_m=sizec(s_mc_m);
  mgamma=reshape(vgamma_km3,iK_m,c_iP)';
  mexp=exp(s_mX_m*mgamma);
  mlikelihood=zeros(iN_m,iK_m);
  for(ci=0;ci<iN_m;ci++)
    {
      if(s_vmulti_m[ci]==0)
	{
	  for(ck=0;ck<iK_m;ck++)
	    {
	      if(s_mc_m[ci][ck]==0)
		{
		  mlikelihood[ci][ck]=1/(1+mexp[ci][ck]);
		}
	      if(s_mc_m[ci][ck]==1)
		{
		  mlikelihood[ci][ck]=mexp[ci][ck]/sumr(mexp[ci][]);
		}
	    }
	}
      if(s_vmulti_m[ci]==1)
	{
	  for(ck=0;ck<iK_m;ck++)
	    {
	      if(s_mc_m[ci][ck]==0)
		{
		  mlikelihood[ci][ck]=1/(1+mexp[ci][ck]);
		}
	      if(s_mc_m[ci][ck]==1)
		{
		  mlikelihood[ci][ck]=mexp[ci][ck]/(1+mexp[ci][ck]);
		}
	    }
	}
    }
  Func[0]=sumc(sumr(log(mlikelihood)));
  return 1;
}

fn_Reg(const vY,const mX, const fl_onlyest)
{
  decl iN,vbeta,dsigma2,mvar_beta,vsd_beta,vp_beta,mest;
  iN=sizer(vY);
  println(iN);
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


fn_MReg(const mY, const mX, const fl_onlyest)
{
  //If fl_onlyest==1, return only estimates.
  //If fl_onlyest==0, return estimates, s.e. and p-values.
  decl iN_k,iJ_k,mXX,mXY,ci,mRSS,mbeta,mSigma,vbeta,vsigma_vec,
    mest,mvar_beta,vsd_beta,vp_beta;
  iN_k=sizer(mY);
  iJ_k=sizec(mY);
  mXX=zeros(c_iP,c_iP);
  mXY=zeros(c_iP,iJ_k);
  for(ci=0;ci<iN_k;ci++)
    {
      mXX+=mX[ci][]'*mX[ci][];
      mXY+=mX[ci][]'*mY[ci][];
    }
  mbeta=invert(mXX)*mXY;
  mRSS=zeros(iJ_k,iJ_k);
  for(ci=0;ci<iN_k;ci++)
    {
      mRSS+=(mY[ci][]-mX[ci][]*mbeta)'*(mY[ci][]-mX[ci][]*mbeta);
    }
  mSigma=mRSS/iN_k;
  vbeta=vec(mbeta);
  vsigma_vec=vec(mSigma);
  if(fl_onlyest==1)
    {
      return vbeta|vsigma_vec;
    }
  else
    {
      mvar_beta=mSigma**invert(mXX);
      vsd_beta=sqrt(diagonal(mvar_beta))';
      vp_beta=tailn(fabs(vbeta./vsd_beta));
      mest=zeros(iJ_k*c_iP+iJ_k*iJ_k,3);
      mest[0:iJ_k*c_iP-1][0]=vbeta;
      mest[0:iJ_k*c_iP-1][1]=vsd_beta;
      mest[0:iJ_k*c_iP-1][2]=vp_beta;
      mest[iJ_k*c_iP: iJ_k*c_iP+iJ_k*iJ_k-1][0]=vsigma_vec;
      return mest;
    }
}

main()
{
  decl time,file,mdata,mZ,vcv,mD,mY,
    vinsure_drug,vinsure_dental,
    mcluster_i_all,mcluster_i,mcluster_index,vncluster,
    ci,incluster,cki,ik,
    mcluster_k,mcluster_k_all,vN_k,vJ_k,amY_k,amX_k,ck,
    iJ_k,iN_k,mY_k,mX_k,cik,
    vl,icluster,
    vN_m,vK_m,im,avc_m,
    cm,vc_m,ckm,mkm_to_k,
    iK_m,amX_m,amc_m,avmulti_m,iN_m,mX_m,mc_m,cim,vmulti_m,
    vfl_m,mfl_k,vfl_m_count,
    vgamma_l,avgamma_c,aestgamma_c,avbeta,avsd_beta,avp_beta,amSigma,
    maxi_m,dopt,mhess,mvar,vsd_gamma_l,vp_gamma_l,
    vgamma_cm,maxi_k,vsd_gamma_cm,vp_gamma_cm,
    mest_reg,vsigma_vec,
    vlikelihood,mgamma_l,mexp_l,vci_m,dLambda_l,mgamma_cm,ikm,
    vexp_cm,dphi,vbeta,mbeta,mSigma,
    dloglikelihood,vci_k,dlikelihood_k,dlikelihood_m,dlikelihood_mk,
    iq_sc,dAIC,dBIC,vsd_beta,vp_beta,dq_secondterm,dq_thirdterm,
    mest_gamma,cj,avnonzero_k,vnonzero,cjk;

  time=timer();
  //Load data
  ///Data of X, Y and splitting for cross validation
  file=fopen("data/mdata.dat","r");
  fscan(file,"%#m",c_iN,c_ivar_file,&mdata);
  fclose(file);

  mZ=mdata[][0:c_iJ-1];
  s_mX=ones(c_iN,1)~mdata[][c_iJ:c_iJ+c_iP-2];
  vinsure_dental=mdata[][c_iJ+c_iP-3];
  vinsure_drug=mdata[][c_iJ+c_iP-2];
  vcv=mdata[][c_iJ+c_iP-1];
  //vcv_raw=mdata[][c_iJ+c_iP-1];
  //vcv=floor((vcv_raw-1)/10)+1;

  mD=(mZ.>0);
  mY=log(mZ).*mD;

  //Data of cluster belongings for each individual
  //mcluster_i_all: for each i
  //raw 1: i=0,...,N
  //raw 2 to 2+iK_max_i-1 : Belonging clusters, 9999 for cells > m
  file=fopen(sprint("data/",c_itau,"mcluster_i.dat"),"r");
  fscan(file,"%#m",c_iN, 1+c_iK_max_i,&mcluster_i_all);
  fclose(file);

  //Number of clusters in which an individual belongs
  mcluster_i=mcluster_i_all[][1:1+c_iK_max_i-1]; 

  // Dummies for cluster belongings of individuals
  mcluster_index=zeros(c_iN,c_iK); 
  //Number ofclusters to which an individual belongs. 

  vncluster=zeros(c_iN,1);

  for(ci=0;ci<c_iN;ci++)
    {
      incluster=0;
      for(cki=0;cki<c_iK_max_i;cki++)
	{
	  ik=mcluster_i[ci][cki];
	  if(ik!=9999)
	    {
	      mcluster_index[ci][ik]=1;
	      incluster++;
	    }
	}
      vncluster[ci]=incluster;
    }

  //Load data about clusters, on k 
  file=fopen(sprint("data/",c_itau,"mcluster_k.dat"),"r");
  fscan(file,"%#m",c_iK,2+c_iKm_max,&mcluster_k_all);
  fclose(file);
  //mcluster_k.dat:  For each k, 
  //raw 1: k: cluster ID, start from 0
  //raw 2: Nk: number of individuals in k th cluster
  //raw 3 - : index of non-zero element (correspond to j=1,2,...,J)
  //1st column, 3rd- raws are filled by 0s, which collesponds to  y=(0,0,,,0)
  //Manual sort is required!
  vN_k=mcluster_k_all[][1];
  mcluster_k=mcluster_k_all[][2:2+c_iKm_max-1];
  //Corresponding m for each cluster k
  vJ_k=zeros(c_iK,1);
  for(ck=1;ck<c_iK;ck++)
    {
      cm=0;
      while (cm<c_iKm_max)
	{
	  if(mcluster_k[ck][cm]!=9999)
	    {
	      vJ_k[ck]++;
	    }
	  cm++;
	}
    }
  //Non-zero elements for each k
  avnonzero_k=new array[c_iK];
  ///For (0,...,0)
  avnonzero_k[0]=0;
  for(ck=1;ck<c_iK;ck++)
    {
      iJ_k=vJ_k[ck];
      vnonzero=zeros(iJ_k,1);
      cj=0;
      cjk=0;
      while (cj<c_iKm_max)
	{
	  if(mcluster_k[ck][cj]!=9999)
	    {
	      vnonzero[cjk]=mcluster_k[ck][cj]-9;
	      //Label for j starts from 9 in our data
	      cjk++;
	    }
	  cj++;
	}
      avnonzero_k[ck]=vnonzero;
    }


  //For k=1,...,K
  amY_k=new array[c_iK];
  amX_k=new array[c_iK];

  ///For (0,...,0)
  amY_k[0]=0;
  amX_k[0]=0; // Not required
  for(ck=1;ck<c_iK;ck++)
    {
      iN_k=vN_k[ck];
      iJ_k=vJ_k[ck];
      mY_k=zeros(iN_k,iJ_k);
      mX_k=zeros(iN_k,c_iP);
      vnonzero=avnonzero_k[ck];
      cik=0;
      for(ci=0;ci<c_iN;ci++)
	{
	  if(mcluster_index[ci][ck]==1)
	    {
	      cj=0;
	      cjk=0;
	      while (cjk<iJ_k)
		{
		  if(cj==vnonzero[cjk])
		    {
		      mY_k[cik][cjk]=mY[ci][cj];
		      cjk++;
		    }
		  cj++;
		}
	      mX_k[cik][]=s_mX[ci][];
	      cik++;
	    }	      
	}
      amY_k[ck]=mY_k;
      amX_k[ck]=mX_k;
    }
  //For which m individuals belong to
  vl=zeros(c_iN,1);
  s_ml=zeros(c_iN,c_iM+1);
  for(ci=0;ci<c_iN;ci++)
    {
      icluster=mcluster_i[ci][0];
      im=vl[ci]=vJ_k[icluster];
      s_ml[ci][im]=1;
    }

  //M specific variables
  vN_m=sumc(s_ml)';
  vK_m=zeros(c_iM+1,1);

  for(ck=1;ck<c_iK;ck++)
    {
      im=vJ_k[ck];
      vK_m[im]++;
    }

  //avc_m: Which clusters have m non-zero elements
  avc_m=new array[c_iM+1]; 
  avc_m[0]=0;  //(0,...,0)
  for(cm=1;cm<c_iM+1;cm++)
    {
      iK_m=vK_m[cm];
      vc_m=zeros(iK_m,1);
      ckm=0;
      for(ck=1;ck<c_iK;ck++)
	{
	  if(vJ_k[ck]==cm)
	    {
	      vc_m[ckm]=ck;
	      ckm++;
	    }
	}
      avc_m[cm]=vc_m;
    }

  //matrix to transform k_m to k
  mkm_to_k=zeros(c_iM+1,c_iK); 
  mkm_to_k[0][0]=0; // (0,...,0) is k=0
  for(cm=1;cm<c_iM+1;cm++)
    {
      vc_m=avc_m[cm];
      iK_m=sizer(vc_m);
      for(ckm=0;ckm<iK_m;ckm++)
	{
	  mkm_to_k[cm][ckm]=vc_m[ckm];
	}
    }

  //Vectors and matrices for each m=1,...,M and i=1,...,N_m 
  amX_m=new array[c_iM+1];  //X
  amc_m=new array[c_iM+1];  //Clusters
  avmulti_m=new array[c_iM+1]; //Multiple cluster belongings

  //For (0,...,0)
  avmulti_m[0]=0;

  for(cm=1;cm<c_iM+1;cm++)
    {
      iN_m=vN_m[cm];
      iK_m=vK_m[cm];
      mX_m=zeros(iN_m,c_iP);
      mc_m=zeros(iN_m,iK_m);
      vmulti_m=zeros(iN_m,1);
      cim=0;
      for(ci=0;ci<c_iN;ci++)
	{
	  if(vl[ci]==cm)
	    {
	      mX_m[cim][]=s_mX[ci][];
	      if(vncluster[ci]>=2)
		{
		  vmulti_m[cim]=1;
		}
	      for(ckm=0;ckm<iK_m;ckm++)
		{
		  ck=mkm_to_k[cm][ckm]; 
		  if(mcluster_index[ci][ck]==1)
		    {
		      mc_m[cim][ckm]=1;
		    }
		}
	      cim++;
	    }
	}
      amc_m[cm]=mc_m;
      amX_m[cm]=mX_m;
      avmulti_m[cm]=vmulti_m;
    }
      
  //fl_m: cases  
  /// ==1: Km==1, 
  /// ==2: no individual belong to multiple clusters
  /// ==3: with individuals who belong to multiple clusters
  //mfl_k: for each k,
  //raw 1: #Individuals who belong only to the k-th cluster
  //raw 2: #Individuals who belong to the k-th cluster and to the other cluster
  //vfl_m_count: for each m, count number of individuals 
  // who belong to multiple clusters

  mfl_k=zeros(c_iK,2);
  for(ci=0;ci<c_iN;ci++)
    {
      incluster=vncluster[ci];
      if(incluster==1)
	{
	  ik=mcluster_i[ci][0];
	  mfl_k[ik][0]++;
	}
      if(incluster>1)
	{
	  for(cki=0;cki<incluster;cki++)
	    {
	      ik=mcluster_i[ci][cki];
	      mfl_k[ik][1]++;
	    }
	}
    }
  vfl_m_count=zeros(c_iM+1,1);

  for(ck=1;ck<c_iK;ck++)
    {
      im=vJ_k[ck];
      if(mfl_k[ck][1]>1)
	{
	  vfl_m_count[im]++;
	}
    }

  vfl_m=ones(c_iM+1,1).*(vK_m.==1)
    +(vK_m.>1).*(2*(vfl_m_count.==0)+3*(vfl_m_count.>0));

  //Set dimensions of paramaters
  avgamma_c=new array[c_iM+1];
  aestgamma_c=new array[c_iK];
  avbeta=new array[c_iK];
  avsd_beta=new array[c_iK];
  avp_beta=new array[c_iK];
  amSigma=new array[c_iK];

  aestgamma_c[0]=0;
  avbeta[0]=0;
  avsd_beta[0]=0;
  avp_beta[0]=0;
  amSigma[0]=0;

  //Estimation
  //1st step: Multinomial Logit for m
  vgamma_l=0.1*ones(c_iP*c_iM,1);
  maxi_m=MaxBFGS(fn_Mlogit_m, &vgamma_l, &dopt,0,1);
  Num2Derivative(fn_Mlogit_m,vgamma_l,&mhess);
  mvar=invert(-mhess);
  vsd_gamma_l=sqrt(diagonal(mvar))';
  vp_gamma_l=tailn(fabs(vgamma_l./vsd_gamma_l));

  for(cm=1;cm<c_iM+1;cm++)
    {
      println("m=",cm);
      iK_m=vK_m[cm];
      iN_m=vN_m[cm];
      s_mc_m=amc_m[cm];
      s_mX_m=amX_m[cm];
      s_vmulti_m=avmulti_m[cm];
      //2nd step: Choice of a cluster given m
      ///Case 1: No estimation
      ///Case 2: 
      if(vfl_m[cm]==2)
	{
	  //For multinomial choice
	  vgamma_cm=0.1*ones(c_iP*(iK_m-1),1);
	  maxi_k=MaxBFGS(fn_K_case2, &vgamma_cm, &dopt,0,1);
	  Num2Derivative(fn_K_case2,vgamma_cm,&mhess);
	  mvar=invert(-mhess);
	  vsd_gamma_cm=sqrt(diagonal(mvar))';
	  vp_gamma_cm=tailn(fabs(vgamma_cm./vsd_gamma_cm));
	  avgamma_c[cm]=vgamma_cm;
	  aestgamma_c[cm]=vgamma_cm~vsd_gamma_cm~vp_gamma_cm;
	}
      ///Case 3: 
      if(vfl_m[cm]==3)
	{
	  //For multivarialte choice
	  vgamma_cm=0.1*ones(c_iP*iK_m,1);
	  maxi_k=MaxBFGS(fn_K_case3, &vgamma_cm, &dopt,0,1);
	  Num2Derivative(fn_K_case3,vgamma_cm,&mhess);
	  mvar=invert(-mhess);
	  vsd_gamma_cm=sqrt(diagonal(mvar))';
	  vp_gamma_cm=tailn(fabs(vgamma_cm./vsd_gamma_cm));
	  avgamma_c[cm]=vgamma_cm;
	  aestgamma_c[cm]=vgamma_cm~vsd_gamma_cm~vp_gamma_cm;
	}
      println("gamma");
      println(vgamma_cm~vsd_gamma_cm~vp_gamma_cm);
      //3rd step: Regression
      for(ckm=0;ckm<iK_m;ckm++)
	{
	  ck=mkm_to_k[cm][ckm];
	  println("k=",ck);
	  mY_k=amY_k[ck];
	  mX_k=amX_k[ck];
	  if(cm==1)
	    {
	      //For m=1, regression for single dependent variable
	      mest_reg=fn_Reg(mY_k,mX_k,0); 
	      avbeta[ck]=mest_reg[0:c_iP-1][0];
	      amSigma[ck]=mest_reg[c_iP][0];
	      avsd_beta[ck]=mest_reg[0:c_iP-1][1];
	      avp_beta[ck]=mest_reg[0:c_iP-1][2];
	    }
	  if(cm>1)
	    {
	      //Multivariate Regression
	      mest_reg=fn_MReg(mY_k,mX_k,0); 
	      //Last argment is 0 when we want s.e and p-values
	      avbeta[ck]=mest_reg[0:cm*c_iP-1][0];
	      avsd_beta[ck]=mest_reg[0:cm*c_iP-1][1];
	      avp_beta[ck]=mest_reg[0:cm*c_iP-1][2];
	      vsigma_vec=mest_reg[cm*c_iP:cm*c_iP+cm*cm-1][0];
	      amSigma[ck]=reshape(vsigma_vec,cm,cm);
	    }
	  println("beta");
	  println(mest_reg);

	}
    }
  //Report estimation results
  println("Results");
  for(cm=1;cm<c_iM+1;cm++)
    {
      println("m=",cm);
      println("Gamma_l",cm);
      println(vgamma_l[(cm-1)*c_iP:cm*c_iP-1]
	      ~vsd_gamma_l[(cm-1)*c_iP:cm*c_iP-1]
	      ~vp_gamma_l[(cm-1)*c_iP:cm*c_iP-1]);
      mest_gamma=aestgamma_c[cm];
      iK_m=vK_m[cm];
      if(vfl_m[cm]==2)
	{
	  for(ckm=0;ckm<iK_m-1;ckm++)
	    {
	      ck=mkm_to_k[cm][ckm];
	      println("Gamma_c, m=" , cm);
	      println("k=", ck);
	      println(mest_gamma[ckm*c_iP:(ckm+1)*c_iP-1][]);
	    }
	}      
      if(vfl_m[cm]==3)
	{
	  for(ckm=0;ckm<iK_m;ckm++)
	    {
	      ck=mkm_to_k[cm][ckm];
	      println("Gamma_c, m=" , cm);
	      println("k=", ck);
	      println(mest_gamma[ckm*c_iP:(ckm+1)*c_iP-1][]);
	    }
	}      
    }

  for(ck=1;ck<c_iK;ck++)
    {
      iJ_k=vJ_k[ck];
      vbeta=avbeta[ck];
      vsd_beta=avsd_beta[ck];
      vp_beta=avp_beta[ck];
      mSigma=amSigma[ck];
      for(cj=0;cj<iJ_k;cj++)
	{
	  println("Beta",ck~cj);
	  println(vbeta[cj*c_iP:(cj+1)*c_iP-1]
		  ~vsd_beta[cj*c_iP:(cj+1)*c_iP-1]
		  ~vp_beta[cj*c_iP:(cj+1)*c_iP-1]);
	}
      println("Sigma",ck);
      println(mSigma);
    }

  //Save output

  format(2000);
  file=fopen(sprint("results/",c_itau,"gamma_l.dat"),"w");
  fprintln(file,"%15.5f",vgamma_l~vsd_gamma_l~vp_gamma_l);
  fclose(file);

  for(cm=1;cm<c_iM+1;cm++)
    {
      mest_gamma=aestgamma_c[cm];
      file=fopen(sprint("results/",c_itau,"gamma_c_",cm, ".dat"),"w");
      fprintln(file,"%15.5f",mest_gamma);
      fclose(file);

    }
  for(ck=1;ck<c_iK;ck++)
    {
      iJ_k=vJ_k[ck];
      vbeta=avbeta[ck];
      vsd_beta=avsd_beta[ck];
      vp_beta=avp_beta[ck];
      mSigma=amSigma[ck];
      file=fopen(sprint("results/",c_itau,"beta_",ck, ".dat"),"w");
      fprintln(file,"%15.5f",vbeta~vsd_beta~vp_beta);
      fclose(file);
      file=fopen(sprint("results/",c_itau,"Sigma_",ck, ".dat"),"w");
      fprintln(file,"%15.5f",mSigma);
      fclose(file);
    }

  //Readable results
  file=fopen(sprint("results/",c_itau,"esimates.dat"),"w");
  for(cm=1;cm<c_iM+1;cm++)
    {
      fprintln(file,"Gamma_l");
      fprint(file,"m=");
      fprintln(file,"%5.0",cm);   
      fprintln(file,"%15.5f",vgamma_l[(cm-1)*c_iP:cm*c_iP-1]
	      ~vsd_gamma_l[(cm-1)*c_iP:cm*c_iP-1]
	      ~vp_gamma_l[(cm-1)*c_iP:cm*c_iP-1]);
      mest_gamma=aestgamma_c[cm];
      iK_m=vK_m[cm];
      if(vfl_m[cm]==2)
	{
	  for(ckm=0;ckm<iK_m-1;ckm++)
	    {
	      fprint(file,"Gamma_ck,k=");
	      ck=mkm_to_k[cm][ckm];
	      fprintln(file,"%5.0f",ck);
	      fprintln(file,"%15.5f",mest_gamma[ckm*c_iP:(ckm+1)*c_iP-1][]);
	    }
	}      
      if(vfl_m[cm]==3)
	{
	  for(ckm=0;ckm<iK_m;ckm++)
	    {
	      fprint(file,"Gamma_ck,k=");
	      ck=mkm_to_k[cm][ckm];
	      fprintln(file,"%5.0f",ck);
	      fprintln(file,"%15.5f",mest_gamma[ckm*c_iP:(ckm+1)*c_iP-1][]);
	    }
	}      
    }

  for(ck=1;ck<c_iK;ck++)
    {
      iJ_k=vJ_k[ck];
      vbeta=avbeta[ck];
      vsd_beta=avsd_beta[ck];
      vp_beta=avp_beta[ck];
      mSigma=amSigma[ck];
      for(cj=0;cj<iJ_k;cj++)
	{
	  fprint(file,"Beta");
	  fprintln(file,"%5.0f",ck~cj);
	  fprintln(file,"%15.5f",vbeta[cj*c_iP:(cj+1)*c_iP-1]
		  ~vsd_beta[cj*c_iP:(cj+1)*c_iP-1]
		  ~vp_beta[cj*c_iP:(cj+1)*c_iP-1]);
	}
      fprint(file,"Sigma");
      fprintln(file,"%5.0f",ck);
      fprintln(file,"%15.5f",mSigma);
    }
  fclose(file);

  //Likelihood function
  vlikelihood=zeros(c_iN,1);
  mgamma_l=reshape(vgamma_l,c_iM,c_iP)';
  mexp_l=ones(c_iN,1)~exp(s_mX*mgamma_l);
  //vci_m: # individuals with m before i
  vci_m=zeros(c_iM+1,1);
  vci_k=zeros(c_iK,1);

  for(ci=0;ci<c_iN;ci++)
    {
      im=vl[ci];
      cim=vci_m[im];
      iK_m=vK_m[im];
      dLambda_l=mexp_l[ci][im]/sumr(mexp_l[ci][]);
      vlikelihood[ci]=dLambda_l;
      if(im!=0)
	{
	  mc_m=amc_m[im];
	  vc_m=mc_m[cim][]';
	  dlikelihood_mk=0;
	  if(vfl_m[im]==1)
	    {
	      for(ckm=0;ckm<iK_m;ckm++)
		{
		  if(vc_m[ckm]==1)
		    {
		      ikm=ckm;
		    }
		}
	      ik=mkm_to_k[im][ikm];
	      cik=vci_k[ik];
	      mY_k=amY_k[ik];
	      mX_k=amX_k[ik];
	      vbeta=avbeta[ik];
	      mSigma=amSigma[ik];
	      dlikelihood_mk=fn_Likelihood_k(mY_k[cik][]',mX_k[cik][]',vbeta,mSigma);
	      vci_k[ik]++;
	    }
	  if(vfl_m[im]==2)
	    {
	      vgamma_cm=avgamma_c[im];
	      mgamma_cm=reshape(vgamma_cm,iK_m-1, c_iP)';
	      vexp_cm=(1~exp(s_mX[ci][]*mgamma_cm))';
	      for(ckm=0;ckm<iK_m;ckm++)
		{
		  if(vc_m[ckm]==1)
		    {
		      ikm=ckm;
		    }
		}
	      dlikelihood_m=vexp_cm[ikm]/sumc(vexp_cm);
	      ik=mkm_to_k[im][ikm];
	      cik=vci_k[ik];
	      mY_k=amY_k[ik];
	      mX_k=amX_k[ik];
	      vbeta=avbeta[ik];
	      mSigma=amSigma[ik];
	      dlikelihood_k=fn_Likelihood_k(mY_k[cik][]',mX_k[cik][]',vbeta,mSigma);
	      dlikelihood_mk=dlikelihood_m*dlikelihood_k;
	      vci_k[ik]++;
	    }
	  if(vfl_m[im]==3)
	    {
	      vgamma_cm=avgamma_c[im];
	      mgamma_cm=reshape(vgamma_cm,iK_m, c_iP)';
	      vexp_cm=exp(s_mX[ci][]*mgamma_cm)';
	      vmulti_m=avmulti_m[im];
	      if(vmulti_m[cim]==0)
		{
		  for(ckm=0;ckm<iK_m;ckm++)
		    {
		      if(vc_m[ckm]==1)
			{
			  ikm=ckm;
			}
		    }
		  dlikelihood_m=vexp_cm[ikm]/sumc(vexp_cm);
		  ik=mkm_to_k[im][ikm];
		  cik=vci_k[ik];
		  mY_k=amY_k[ik];
		  mX_k=amX_k[ik];
		  vbeta=avbeta[ik];
		  mSigma=amSigma[ik];
		  dlikelihood_k=fn_Likelihood_k(mY_k[cik][]',mX_k[cik][]',vbeta,mSigma);
		  dlikelihood_mk=dlikelihood_m*dlikelihood_k;
		  vci_k[ik]++;
		}
	      if(vmulti_m[cim]==1)
		{
		  for(ckm=0;ckm<iK_m;ckm++)
		    {
		      if(vc_m[ckm]==1)
			{
			  dlikelihood_m=vexp_cm[ckm]/(1+vexp_cm[ckm]);
			  ik=mkm_to_k[im][ckm];
			  cik=vci_k[ik];
			  mY_k=amY_k[ik];
			  mX_k=amX_k[ik];
			  vbeta=avbeta[ik];
			  mSigma=amSigma[ik];
			  dlikelihood_k=fn_Likelihood_k(mY_k[cik][]',mX_k[cik][]',vbeta,mSigma);
			  dlikelihood_mk+=dlikelihood_m*dlikelihood_k;
			  vci_k[ik]++;
			}
		    }
		}
	    }
	  vlikelihood[ci]*=dlikelihood_mk;
	}
      vci_m[im]++;
    }
  dloglikelihood=sumc(log(vlikelihood));

  //Information criterion
  dq_secondterm=0;
  dq_thirdterm=0;
  for(cm=1;cm<c_iM+1;cm++)
    {
      iK_m=vK_m[cm];
      if(vfl_m==2)
	{
	  dq_secondterm+=iK_m-1;
	}
      if(vfl_m==3)
	{
	  dq_secondterm+=iK_m;
	}
      for(ckm=0;ckm<iK_m;ckm++)
	{
	  dq_thirdterm+=c_iP*cm+cm*(cm+1)/2;
	}
    }

  iq_sc=c_iM*c_iP+dq_secondterm+dq_thirdterm; 
  dAIC=-2*dloglikelihood+2*iq_sc;
  dBIC=-2*dloglikelihood+iq_sc*log(c_iN);

  //Report AIC
  println("LogL, AIC, BIC");
  println(dloglikelihood~dAIC~dBIC);
  file=fopen(sprint("results/",c_itau,"IC_interdep.dat"),"w");
  fprintln(file,"%15.5f",dloglikelihood~dAIC~dBIC);
  fclose(file);

}
