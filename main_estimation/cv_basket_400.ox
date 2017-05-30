#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include <oxdraw.h>
#import<solvenle>
#import<maximize>
#include"myconst_basket_400.h"

static decl s_mX,s_ml,  //For multinomial logit for m
  s_mc_m, s_mX_m, //For k, case 2 and 3
  s_vmulti_m; //For k, only case 3



fn_Mlogit_m(const vgamma_l, const Func, const Score, const Hess)
{
  decl mexp,vlikelihood,mgamma,iN;
  iN=sizer(s_mX);
  mgamma=reshape(vgamma_l,c_iM,c_iP)';
  mexp=exp(s_mX*mgamma);
  mexp=ones(iN,1)~mexp;
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
  decl time,file,mdata,mZ_all,vcv,mD_all,mY_all,mX_all,
    vinsure_drug_all,vinsure_dental_all,
    mcluster_i_file,mcluster_i_all,mcluster_index_all,vncluster_all,incluster_all,
    ci,cki,cj,cik,cm,ckm,cim,cjk,ck,
    im,ik,iJ_k,iN_k,icluster,iK_m,iN_m,
    incluster,mcluster_k,mcluster_k_all,vN_k_all,vJ_k,amY_k,amX_k,
    vl,mY_k,mX_k,
    vN_m,vK_m,avc_m,
    vc_m,mkm_to_k,
    amX_m,amc_m,avmulti_m,mX_m,mc_m,vmulti_m,
    vfl_m,mfl_k,vfl_m_count,avnonzero_k,vnonzero,
    avgamma_c_init,
    mest_gamma_l,vgamma_l_init,mest_gamma_c,
    vgamma_l,avgamma_c,aestgamma_c,avbeta,avsd_beta,avp_beta,amSigma,
    maxi_m,dopt,mhess,mvar,vsd_gamma_l,vp_gamma_l,
    vgamma_cm,maxi_k,vsd_gamma_cm,vp_gamma_cm,
    mest_reg,vsigma_vec,
    vbeta,mbeta,mSigma,
    //CV
    cv,vN_k_test,vN_k_train,
    mD_cv_mean,mD_cv_mse,mZ_cv_mean,mZ_cv_mse,vZ_cv_mean_agg,vZ_cv_mse_agg,
    mD_cv_mean_drug,mD_cv_meanchange_drug,mD_cv_mse_drug,
    mD_cv_mean_dental,mD_cv_meanchange_dental,mD_cv_mse_dental,
    mZ_cv_mean_drug,mZ_cv_meanchange_drug,mZ_cv_mse_drug,
    mZ_cv_mean_dental,mZ_cv_meanchange_dental,mZ_cv_mse_dental,
    vD_star,mD_cv,vD_drug_star,vD_drug_cv,vD_dental_star,vD_dental_cv,mZ_cv,
    mY_true_one,mZ_true,mD_true,
    mD_true_drug,mD_true_dental,
    mZ_true_drug,mZ_true_dental,vZ_true_agg,
    mY_mean_one,mY_mse_one,mZ_mean,mZ_mse,
    mZ_meanchange_drug,mZ_mse_drug,
    mZ_meanchange_dental,mZ_mse_dental,
    mZ_cv_ndrug,mZ_cv_ndental,mY_all_adjusted,
    mZ_mean_dental,mZ_mean_drug,mindex_nan,
    mZ_train,mY_train,mD_train,mX_train,
    mZ_test,mD_test,mX_test,ci_train,ci_test,
    iN_ndrug_test,iN_ndental_test,
    mZ_ndrug_test,mD_ndrug_test,mX_ndrug_test,
    mZ_ndental_test,mD_ndental_test,mX_ndental_test,
    ci_ndrug,ci_ndental,vZ_agg_test,vD_agg_test,mCP_v,
    maxi_gamma,vexp,vexp_ndrug,vexp_ndental,vtheta_reg,
    vN_cv_test,iv,
    indental,indrug,vN_ndrug,vN_ndental,iN_train,iN_test,
    mY_test,vgamma,iJ,mCP_v_ndrug,mCP_v_ndental,
    vY,mX,cij,vsigma2_ndrug,
    vsigma2_ndental,vcv_raw,mY_v,
    inD_test_one,inD_test_zero,
    iinonnan,cnonnan,vY_mean_one_adj,vY_mse_one_adj,
    vY_mean_one_mean,vY_rmse_one,
    //Only for this file
    vl_all,ml_all,
    mcluster_i_test,mcluster_i_train,vl_train,
    mcluster_index_train,vncluster_train,
    vl_test,ml_test,ml_train,
    mcluster_index_test,vncluster_test,
    amX_m_train,vN_m_train,amY_k_train,amX_k_train,amc_m_train,avmulti_m_train,
    mL_star,vL_cv,vC_star,iC_cv_km,vY_cv_i,
    vC_star_unity,inC_unity,
    mnonzero_cv,vnonzero_count_cv,mY_cv_multi,
    mgamma_c,vest_reg,mgamma_l,cji,iji,mY_cv,ifl_unity,vY_true_denom,
    mY_cv_drug,mD_cv_drug,mZ_cv_drug,mY_cv_multi_drug,
    mY_cv_dental,mD_cv_dental,mZ_cv_dental,mY_cv_multi_dental,
    vN11,vN01,vN10,mZ_cv_mean_11,mZ_cv_mse_11,mZ_cv_mean_01,mZ_cv_mse_01,mZ_cv_mean_10,mZ_cv_mse_10,
    mZ_true_11,mZ_true_10,mN11,mN10,mN01,mN00,cr,amD_cv_r,mD_cv_r,amZ_cv_r,mZ_cv_r;
  time=timer();
  //Load data
  ///Data of X, Y and splitting for cross validation
  file=fopen("data/mdata.dat","r");
  fscan(file,"%#m",c_iN,c_ivar_file,&mdata);
  fclose(file);
  mZ_all=mdata[][0:c_iJ-1];
  mX_all=ones(c_iN,1)~mdata[][c_iJ:c_iJ+c_iP-2];
  vinsure_dental_all=mdata[][c_iJ+c_iP-3];
  vinsure_drug_all=mdata[][c_iJ+c_iP-2];
  indental=c_iP-2;
  indrug=c_iP-1;

  vcv_raw=mdata[][c_iJ+c_iP-1];
  vcv=floor((vcv_raw-1)/2)+1;
  //vcv=mdata[][c_iJ+c_iP-1];

  mD_all=(mZ_all.>0);
  mY_all=mD_all.*log(mZ_all);
  mY_all_adjusted=mY_all;
  for(ci=0;ci<c_iN;ci++)
    {
      for(cj=0;cj<c_iJ;cj++)
	{
	  if(mD_all[ci][cj]==0)
	    {
	      mY_all_adjusted[ci][cj]=0;
	    }
	}
    }

  //Data of cluster belongings for each individual
  //mcluster_i_all: for each i
  //raw 1: i=0,...,N
  //raw 2 to 2+iK_max_i-1 : Belonging clusters, 9999 for cells > m
  file=fopen(sprint("data/",c_itau,"mcluster_i.dat"),"r");
  fscan(file,"%#m",c_iN, 1+c_iK_max_i,&mcluster_i_file);
  fclose(file);
  //Number of clusters in which an individual belongs
  mcluster_i_all=mcluster_i_file[][1:1+c_iK_max_i-1]; 

  // Dummies for cluster belongings of individuals
  mcluster_index_all=zeros(c_iN,c_iK); 
  //Number ofclusters to which an individual belongs. 

  vncluster_all=zeros(c_iN,1);

  for(ci=0;ci<c_iN;ci++)
    {
      incluster_all=0;
      for(cki=0;cki<c_iK_max_i;cki++)
	{
	  ik=mcluster_i_all[ci][cki];
	  if(ik!=9999)
	    {
	      mcluster_index_all[ci][ik]=1;
	      incluster_all++;
	    }
	}
      vncluster_all[ci]=incluster_all;
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

  vN_k_all=mcluster_k_all[][1];
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

  //For which m individuals belong to
  vl_all=zeros(c_iN,1);
  ml_all=zeros(c_iN,c_iM+1);
  for(ci=0;ci<c_iN;ci++)
    {
      icluster=mcluster_i_all[ci][0];
      im=vl_all[ci]=vJ_k[icluster];
      ml_all[ci][im]=1;
    }
  //M specific variables
  vK_m=zeros(c_iM+1,1);

  for(ck=1;ck<c_iK;ck++)
    {
      im=vJ_k[ck];
      vK_m[im]++;
    }

  //avc_m: Which clusters have m non-zero elements
  avc_m=new array[c_iM+1]; 
  avc_m[0]=0;  
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
  mkm_to_k[0][0]=0; 
  for(cm=1;cm<c_iM+1;cm++)
    {
      vc_m=avc_m[cm];
      iK_m=sizer(vc_m);
      for(ckm=0;ckm<iK_m;ckm++)
  	{
  	  mkm_to_k[cm][ckm]=vc_m[ckm];
  	}
    }

  //fl_m: cases  
  // ==1: Km==1, 
  // ==2: no individual belong to multiple clusters
  // ==3: with individuals who belong to multiple clusters
  // mfl_k: for each k,
  // raw 1: #Individuals who belong only to the k-th cluster
  // raw 2: #Individuals who belong to the k-th cluster and to the other cluster
  // vfl_m_count: for each m, count number of individuals 
  // who belong to multiple clusters

  mfl_k=zeros(c_iK,2);
  for(ci=0;ci<c_iN;ci++)
    {
      incluster=vncluster_all[ci];
      if(incluster==1)
  	{
  	  ik=mcluster_i_all[ci][0];
  	  mfl_k[ik][0]++;
  	}
      if(incluster>1)
  	{
  	  for(cki=0;cki<incluster;cki++)
  	    {
  	      ik=mcluster_i_all[ci][cki];
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

  // Set dimensions of paramaters
  avgamma_c_init=new array[c_iM+1];
  //Load output
  file=fopen(sprint("results/",c_itau,"gamma_l.dat"),"r");
  fscan(file,"%#m",c_iM*c_iP,3,&mest_gamma_l);
  fclose(file);
  vgamma_l_init=mest_gamma_l[][0];

  for(cm=1;cm<c_iM+1;cm++)
    {
      iK_m=vK_m[cm];
      if(vfl_m[cm]==2)
  	{
  	  file=fopen(sprint("results/",c_itau,"gamma_c_",cm, ".dat"),"r");
  	  fscan(file,"%#m",c_iP*(iK_m-1),3,&mest_gamma_c);
  	  fclose(file);
  	}
      if(vfl_m[cm]==3)
  	{
  	  file=fopen(sprint("results/",c_itau,"gamma_c_",cm, ".dat"),"r");
  	  fscan(file,"%#m",c_iP*iK_m,3,&mest_gamma_c);
  	  fclose(file);
  	}
      avgamma_c_init[cm]=mest_gamma_c[][0];
    }
  mD_cv_mean=zeros(c_iV_2p,c_iJ); 
  mD_cv_mse=zeros(c_iV_2p,c_iJ); 
  mY_mean_one=zeros(c_iV_2p,c_iJ); 
  mY_mse_one=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mean=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mse=zeros(c_iV_2p,c_iJ); 
  vZ_cv_mean_agg=zeros(c_iV_2p,1);
  vZ_cv_mse_agg=zeros(c_iV_2p,1); 

  mD_cv_mean_drug=zeros(c_iV_2p,c_iJ); 
  mD_cv_mean_dental=zeros(c_iV_2p,c_iJ); 
  mD_cv_meanchange_drug=zeros(c_iV_2p,c_iJ); 
  mD_cv_meanchange_dental=zeros(c_iV_2p,c_iJ); 
  mD_cv_mse_drug=zeros(c_iV_2p,c_iJ); 
  mD_cv_mse_dental=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mean_drug=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mean_dental=zeros(c_iV_2p,c_iJ); 
  mZ_cv_meanchange_drug=zeros(c_iV_2p,c_iJ); 
  mZ_cv_meanchange_dental=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mse_drug=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mse_dental=zeros(c_iV_2p,c_iJ); 

  mZ_cv_mean_11=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mean_01=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mean_10=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mse_11=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mse_01=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mse_10=zeros(c_iV_2p,c_iJ); 

  //Summary for true data
  mY_true_one=mZ_true=mD_true=mD_true_drug=mD_true_dental=
    mZ_true_drug=mZ_true_dental=zeros(c_iV_2p,c_iJ);
  vZ_true_agg=zeros(c_iV_2p,1);
  mZ_true_11=mZ_true_10=zeros(c_iV_2p,c_iJ);
  mN11=mN10=mN01=mN00=zeros(c_iV_2p,c_iJ);

  //Index for test data, which takes unity when all test data show D_j=0
  mindex_nan=zeros(c_iV_2p,c_iJ);

  vN_ndrug=vN_ndental=zeros(c_iV_2p,1);
  //Size of test data
  vN_cv_test=zeros(c_iV_2p,1);
  for(ci=0;ci<c_iN;ci++)
    {
      iv=vcv[ci]-1;
      vN_cv_test[iv]++;
    } 

  //Loop start

  for(cv=0;cv<c_iV_2p;cv++)
    {
      println("cv=",cv);
      println("time lapsed=",timespan(time));
      //Separating training and test data
      iN_test=vN_cv_test[cv];
      iN_train=c_iN-iN_test;
      vN_k_test=vN_k_train=zeros(c_iK,1);
      mY_train=mD_train=zeros(iN_train,c_iJ);
      mX_train=zeros(iN_train,c_iP);
      mY_test=mZ_test=mD_test=zeros(iN_test,c_iJ);
      mX_test=zeros(iN_test,c_iP);

      mcluster_i_test=zeros(iN_test,c_iK);
      mcluster_index_test=zeros(iN_test,c_iK);
      vncluster_test=zeros(iN_test,1);
      mcluster_i_train=zeros(iN_train,c_iK);
      mcluster_index_train=zeros(iN_train,c_iK);
      vncluster_train=zeros(iN_train,1);

      vl_test=zeros(iN_test,1);
      ml_test=zeros(iN_test,c_iM+1);
      vl_train=zeros(iN_train,1);
      ml_train=zeros(iN_train,c_iM+1);

      ci_train=ci_test=iN_ndrug_test=iN_ndental_test=0;
      for(ci=0;ci<c_iN;ci++)
  	{
  	  if(vcv[ci]==cv+1)
  	    {
  	      mX_test[ci_test][]=mX_all[ci][];
  	      mZ_test[ci_test][]=mZ_all[ci][];
  	      mY_test[ci_test][]=mY_all_adjusted[ci][];
  	      mD_test[ci_test][]=mD_all[ci][];
  	      mcluster_i_test[ci_test][]=mcluster_i_all[ci][];
  	      mcluster_index_test[ci_test][]=mcluster_index_all[ci][];
  	      vncluster_test[ci_test]=vncluster_all[ci];
  	      vl_test[ci_test]=vl_all[ci];
  	      ml_test[ci_test][]=ml_all[ci][];
  	      if(vinsure_drug_all[ci]==1)
  		{
  		  iN_ndrug_test++;
  		}
  	      if(vinsure_dental_all[ci]==1)
  		{
  		  iN_ndental_test++;
  		}
  	      ci_test++;
  	    }
  	  else
  	    {
  	      mX_train[ci_train][]=mX_all[ci][];
  	      mY_train[ci_train][]=mY_all_adjusted[ci][];
  	      mD_train[ci_train][]=mD_all[ci][];
  	      mcluster_i_train[ci_train][]=mcluster_i_all[ci][];
  	      mcluster_index_train[ci_train][]=mcluster_index_all[ci][];
  	      vncluster_train[ci_train]=vncluster_all[ci];
  	      vl_train[ci_train]=vl_all[ci];
  	      ml_train[ci_train][]=ml_all[ci][];
  	      ci_train++;
  	      for(ck=0;ck<c_iK;ck++)
  		{
  		  if(mcluster_index_all[ci][ck]==1)
  		    {
  		      vN_k_train[ck]++;
  		    }
  		}
  	    }
  	}
      vY_true_denom=sumc(mD_test);
      vY_true_denom=vY_true_denom+(vY_true_denom.==0);
      mY_true_one[cv][]=sumc(mD_test.*mY_test) ./ vY_true_denom;
      mD_true[cv][]= meanc(mD_test);
      mZ_true[cv][]= meanc(mZ_test);
      vZ_true_agg[cv]=meanc(sumr(mZ_test)); //Aggregate medical expenditure

      //On counter-factual
      vN_ndrug[cv]=iN_ndrug_test;
      vN_ndental[cv]=iN_ndental_test;
      mZ_ndrug_test=mD_ndrug_test=zeros(iN_ndrug_test,c_iJ);
      mX_ndrug_test=zeros(iN_ndrug_test,c_iP);
      mZ_ndental_test=mD_ndental_test=zeros(iN_ndental_test,c_iJ);
      mX_ndental_test=zeros(iN_ndental_test,c_iP);
      ci_ndrug=ci_ndental=0;
      for(ci=0;ci<c_iN;ci++)
  	{
  	  if(vcv[ci]==cv+1 && vinsure_drug_all[ci]==1)
  	    {
  	      mZ_ndrug_test[ci_ndrug][]=mZ_all[ci][];
  	      mD_ndrug_test[ci_ndrug][]=mD_all[ci][];
  	      mX_ndrug_test[ci_ndrug][]=mX_all[ci][];
  	      ci_ndrug++;
  	    }
  	  if(vcv[ci]==cv+1 && vinsure_dental_all[ci]==1)
  	    {
  	      mZ_ndental_test[ci_ndental][]=mZ_all[ci][];
  	      mD_ndental_test[ci_ndental][]=mD_all[ci][];
  	      mX_ndental_test[ci_ndental][]=mX_all[ci][];
  	      ci_ndental++;
  	    }	  
  	}
      mX_ndrug_test[][indrug]=ones(iN_ndrug_test,1);
      mX_ndental_test[][indental]=ones(iN_ndental_test,1);


      // K-specific variables for training data
      // For k=1,...,K
      amY_k_train=new array[c_iK];
      amX_k_train=new array[c_iK];
      
      // For (0,...,0)
      amY_k_train[0]=0;
      amX_k_train[0]=0; 
      for(ck=1;ck<c_iK;ck++)
  	{
  	  iN_k=vN_k_train[ck];
  	  iJ_k=vJ_k[ck];
  	  mY_k=zeros(iN_k,iJ_k);
  	  mX_k=zeros(iN_k,c_iP);
  	  vnonzero=avnonzero_k[ck];
  	  cik=0;
  	  for(ci=0;ci<iN_train;ci++)
  	    {
  	      if(mcluster_index_train[ci][ck]==1)
  		{
  		  cj=0;
  		  cjk=0;
  		  while (cjk<iJ_k)
  		    {
  		      if(cj==vnonzero[cjk])
  			{
  			  mY_k[cik][cjk]=mY_train[ci][cj];
  			  cjk++;
  			}
  		      cj++;
  		    }
  		  mX_k[cik][]=mX_train[ci][];
  		  cik++;
  		}	      
  	    }
  	  amY_k_train[ck]=mY_k;
  	  amX_k_train[ck]=mX_k;
  	}

      // M-specific variables for training data
      vN_m_train=sumc(ml_train)';
      // Vectors and matrices for each m=1,...,M and i=1,...,N_m 
      amX_m_train=new array[c_iM+1];  
      amc_m_train=new array[c_iM+1];  
      avmulti_m_train=new array[c_iM+1]; 
      // For (0,...,0)
      avmulti_m_train[0]=0;
      for(cm=1;cm<c_iM+1;cm++)
  	{
  	  iN_m=vN_m_train[cm];
  	  iK_m=vK_m[cm];
  	  mX_m=zeros(iN_m,c_iP);
  	  mc_m=zeros(iN_m,iK_m);
  	  vmulti_m=zeros(iN_m,1);
  	  cim=0;
  	  for(ci=0;ci<iN_train;ci++)
  	    {
  	      if(vl_train[ci]==cm)
  		{
  		  mX_m[cim][]=mX_train[ci][];
  		  if(vncluster_train[ci]>=2)
  		    {
  		      vmulti_m[cim]=1;
  		    }
  		  for(ckm=0;ckm<iK_m;ckm++)
  		    {
  		      ck=mkm_to_k[cm][ckm]; 
  		      if(mcluster_index_train[ci][ck]==1)
  			{
  			  mc_m[cim][ckm]=1;
  			}
  		    }
  		  cim++;
  		}
  	    }
  	  amX_m_train[cm]=mX_m;
  	  amc_m_train[cm]=mc_m;
  	  avmulti_m_train[cm]=vmulti_m;
  	}


      // On counter-factual
      vN_ndrug[cv]=iN_ndrug_test;
      vN_ndental[cv]=iN_ndental_test;
      mZ_ndrug_test=mD_ndrug_test=zeros(iN_ndrug_test,c_iJ);
      mX_ndrug_test=zeros(iN_ndrug_test,c_iP);
      mZ_ndental_test=mD_ndental_test=zeros(iN_ndental_test,c_iJ);
      mX_ndental_test=zeros(iN_ndental_test,c_iP);
      ci_ndrug=ci_ndental=0;
      for(ci=0;ci<c_iN;ci++)
  	{
  	  if(vcv[ci]==cv+1 && vinsure_drug_all[ci]==1)
  	    {
  	      mZ_ndrug_test[ci_ndrug][]=mZ_all[ci][];
  	      mD_ndrug_test[ci_ndrug][]=mD_all[ci][];
  	      mX_ndrug_test[ci_ndrug][]=mX_all[ci][];
  	      ci_ndrug++;
  	    }
  	  if(vcv[ci]==cv+1 && vinsure_dental_all[ci]==1)
  	    {
  	      mZ_ndental_test[ci_ndental][]=mZ_all[ci][];
  	      mD_ndental_test[ci_ndental][]=mD_all[ci][];
  	      mX_ndental_test[ci_ndental][]=mX_all[ci][];
  	      ci_ndental++;
  	    }	  
  	}
      mX_ndrug_test[][indrug]=ones(iN_ndrug_test,1);
      mX_ndental_test[][indental]=ones(iN_ndental_test,1);

      mD_true_dental[cv][]= meanc(mD_ndental_test);
      mZ_true_dental[cv][]= meanc(mZ_ndental_test);
      mD_true_drug[cv][]= meanc(mD_ndrug_test);
      mZ_true_drug[cv][]= meanc(mZ_ndrug_test);

      //Estimation
      //Set dimensions of paramaters
      avgamma_c=new array[c_iM+1];
      avbeta=new array[c_iK];
      amSigma=new array[c_iK];
      avbeta[0]=0;
      amSigma[0]=0;

      //1st step: Multinomial Logit for m
      s_ml=ml_train;
      s_mX=mX_train;
      vgamma_l=vgamma_l_init;
      //      vgamma_l=0.1*ones(c_iP,1);
      maxi_m=MaxBFGS(fn_Mlogit_m, &vgamma_l, &dopt,0,1);

      for(cm=1;cm<c_iM+1;cm++)
	{
	  iK_m=vK_m[cm];
	  s_mX_m=amX_m_train[cm];
	  s_vmulti_m=avmulti_m_train[cm];
	  s_mc_m=amc_m_train[cm];

	  //2nd step: Choice of a cluster given m
	  ///Case 1: No estimation
	  ///Case 2: 
	  if(vfl_m[cm]==2)
	    {
	      //For multinomial choice
	      vgamma_cm=avgamma_c_init[cm];
	      //	      vgamma_cm=0.1*ones(c_iP*(iK_m-1),1);
	      maxi_k=MaxBFGS(fn_K_case2, &vgamma_cm, &dopt,0,1);
	      mgamma_c=reshape(vgamma_cm,iK_m-1,c_iP)';
	      avgamma_c[cm]=mgamma_c;
	    }
	  ///Case 3: 
	  if(vfl_m[cm]==3)
	    {
	      //For multivarialte choice
	      vgamma_cm=avgamma_c_init[cm];
	      //	      vgamma_cm=0.1*ones(c_iP*iK_m,1);
	      maxi_k=MaxBFGS(fn_K_case3, &vgamma_cm, &dopt,0,1);
	      mgamma_c=reshape(vgamma_cm,iK_m,c_iP)';
	      avgamma_c[cm]=mgamma_c;
	    }

	  //3rd step: Regression
	  for(ckm=0;ckm<iK_m;ckm++)
	    {
	      ck=mkm_to_k[cm][ckm];
	      println("k=",ck);
	      mY_k=amY_k_train[ck];
	      mX_k=amX_k_train[ck];

	      if(cm==1)
		{
		  //For m=1, regression for single dependent variable
		  vest_reg=fn_Reg(mY_k,mX_k,1); 
		  //Last argment is 1 when we only want estimates
		  avbeta[ck]=vest_reg[0:c_iP-1];
		  amSigma[ck]=vest_reg[c_iP];
		}
	      if(cm>1)
		{
		  //Multivariate Regression
		  mest_reg=fn_MReg(mY_k,mX_k,1); 
		  //Last argment is 1 when we only want estimates
		  avbeta[ck]=mest_reg[0:cm*c_iP-1];
		  vsigma_vec=mest_reg[cm*c_iP:cm*c_iP+cm*cm-1];
		  amSigma[ck]=reshape(vsigma_vec,cm,cm);
		}
	    }
	}
       // CV targets:
       // L: 
      mgamma_l=reshape(vgamma_l,c_iM,c_iP)';
      amD_cv_r=amZ_cv_r=new array[c_iR];
      for(cr=0;cr<c_iR;cr++)
	{
	  mL_star=mX_test*mgamma_l+ranextremevalue(iN_test,c_iM,0,1);
	  mL_star=zeros(iN_test,1)~mL_star;
	  vL_cv=maxcindex(mL_star')';
	  mZ_cv=mY_cv=mD_cv=zeros(iN_test,c_iJ);
	  mZ_cv_drug=zeros(iN_ndrug_test,c_iJ);
	  mZ_cv_dental=zeros(iN_ndental_test,c_iJ);
	  mD_cv_drug=zeros(iN_ndrug_test,c_iJ); 
	  mD_cv_dental=zeros(iN_ndental_test,c_iJ); 
	  mY_cv_drug=zeros(iN_ndrug_test,c_iJ); 
	  mY_cv_dental=zeros(iN_ndental_test,c_iJ); 

	  for(ci=0;ci<iN_test;ci++)
	    {
	      im=vL_cv[ci];
	      iK_m=vK_m[im];
	      if(im!=0)
		{
		  // 	      //K
		  mgamma_c=avgamma_c[im];
		  if(iK_m==1) 
		    {
		      ik=mkm_to_k[im][0];
		      // 		  //Regression
		      vbeta=avbeta[ik];
		      mSigma=amSigma[ik];
		      vnonzero=avnonzero_k[ik];
		      if(im==1)
			{
			  vY_cv_i=mX_test[ci][]*vbeta+sqrt(mSigma)*rann(1,1);
			  // 		      //Cross varidation
			  iji=vnonzero;
			  mY_cv[ci][iji]=vY_cv_i;
			  mD_cv[ci][iji]=1;
			  mZ_cv[ci][iji]=mD_cv[ci][iji]*exp(mY_cv[ci][iji]);
			}
		      else
			{
			  mbeta=reshape(vbeta,im,c_iP)';
			  vY_cv_i=(mX_test[ci][]*mbeta)'+choleski(mSigma)*rann(im,1);
			  // 		      //Cross varidation
			  for(cji=0;cji<im;cji++)
			    {
			      iji=vnonzero[cji];
			      mY_cv[ci][iji]=vY_cv_i[cji];
			      mD_cv[ci][iji]=1;
			      mZ_cv[ci][iji]=mD_cv[ci][iji]*exp(mY_cv[ci][iji]);
			    }
			}
		    }
		  if(vfl_m[im]==2)
		    {
		      // 		  //Logit
		      vC_star=mX_test[ci][]*mgamma_c+ranextremevalue(1,iK_m-1,0,1);
		      vC_star=0|vC_star';
		      iC_cv_km=maxcindex(vC_star);
		      ik=mkm_to_k[im][iC_cv_km];
		      // 		  //Regression
		      vbeta=avbeta[ik];
		      mSigma=amSigma[ik];
		      vnonzero=avnonzero_k[ik];
		      if(im==1)
			{
			  vY_cv_i=mX_test[ci][]*vbeta+sqrt(mSigma)*rann(1,1);
			  // 		      //Cross varidation
			  iji=vnonzero;
			  mY_cv[ci][iji]=vY_cv_i;
			  mD_cv[ci][iji]=1;
			  mZ_cv[ci][iji]=mD_cv[ci][iji]*exp(mY_cv[ci][iji]);
			}
		      else
			{
			  mbeta=reshape(vbeta,im,c_iP)';
			  vY_cv_i=(mX_test[ci][]*mbeta)'+choleski(mSigma)*rann(im,1);
			  // 		      //Cross varidation
			  for(cji=0;cji<im;cji++)
			    {
			      iji=vnonzero[cji];
			      mY_cv[ci][iji]=vY_cv_i[cji];
			      mD_cv[ci][iji]=1;
			      mZ_cv[ci][iji]=mD_cv[ci][iji]*exp(mY_cv[ci][iji]);
			    }
			}
		    }

		  if(vfl_m[im]==3)
		    {
		      vC_star=mX_test[ci][]*mgamma_c+ranextremevalue(1,iK_m,0,1);
		      vC_star_unity=(vC_star'.>=0);
		      inC_unity=sumc(vC_star_unity); //#k C^{*}_{ki}>=0
		      if(inC_unity<=1)
			{
			  // 		      Same procedure for inC_unity==0 and 1 
			  iC_cv_km=maxcindex(vC_star');
			  ik=mkm_to_k[im][iC_cv_km];
			  // 		      //Regression
			  vbeta=avbeta[ik];
			  mSigma=amSigma[ik];
			  vnonzero=avnonzero_k[ik];
			  if(im==1)
			    {
			      vY_cv_i=mX_test[ci][]*vbeta+sqrt(mSigma)*rann(1,1);
			      // // 			  Cross varidation
			      iji=vnonzero;
			      mY_cv[ci][iji]=vY_cv_i;
			      mD_cv[ci][iji]=1;
			      mZ_cv[ci][iji]=mD_cv[ci][iji]*exp(mY_cv[ci][iji]);
			    }
			  else
			    {
			      mbeta=reshape(vbeta,im,c_iP)';
			      vY_cv_i=(mX_test[ci][]*mbeta)'+choleski(mSigma)*rann(im,1);
			      // 			  //Cross varidation
			      for(cji=0;cji<im;cji++)
				{
				  iji=vnonzero[cji];
				  mD_cv[ci][iji]=1;
				  mY_cv[ci][iji]=vY_cv_i[cji];
				  mZ_cv[ci][iji]=mD_cv[ci][iji]*exp(mY_cv[ci][iji]);
				}
			    }
			}
		      if(inC_unity>1)
			{
			  mnonzero_cv=zeros(im,inC_unity); //Index for non-zero element
			  vnonzero_count_cv=zeros(c_iJ,1); //Count for non-zero element
			  mY_cv_multi=zeros(c_iJ,inC_unity); //Matrix of Y from multiple clusters
			  cki=0;
			  for(ck=0;ck<iK_m;ck++)
			    {
			      ifl_unity=vC_star_unity[ck];
			      if(ifl_unity==1)
				{
				  ik=mkm_to_k[im][ck];
				  mnonzero_cv[][cki]=avnonzero_k[ik];
				  vbeta=avbeta[ik];
				  mSigma=amSigma[ik];
				  vnonzero=avnonzero_k[ik];
				  if(im==1)
				    {
				      vY_cv_i=mX_test[ci][]*vbeta+sqrt(mSigma)*rann(1,1);
				      // 				  //Cross varidation
				      iji=vnonzero;
				      vnonzero_count_cv[iji]++;
				      mY_cv_multi[iji][cki]=vY_cv_i;
				      mD_cv[ci][iji]=1;
				    }
				  else
				    {
				      mbeta=reshape(vbeta,im,c_iP)';
				      vY_cv_i=(mX_test[ci][]*mbeta)'+choleski(mSigma)*rann(im,1);
				      // 				  //Cross varidation
				      vnonzero=avnonzero_k[ik];
				      for(cji=0;cji<im;cji++)
					{
					  iji=vnonzero[cji];
					  vnonzero_count_cv[iji]++;
					  mY_cv_multi[iji][cki]=vY_cv_i[cji];
					  mD_cv[ci][iji]=1;
					}
				    }
				  cki++;
				}
			    }
			  vnonzero_count_cv=vnonzero_count_cv+(vnonzero_count_cv.==0); //Avoinding 0 division
			  mY_cv[ci][]=(sumr(mY_cv_multi)./vnonzero_count_cv)';
			  mZ_cv[ci][]=mD_cv[ci][].*exp(mY_cv[ci][]);
			}
		    }
		}
	    }
	  amD_cv_r[cr]=mD_cv;
	  amZ_cv_r[cr]=mZ_cv;
	}
      
      mD_cv=mZ_cv=zeros(iN_test,c_iJ);
      for(cr=0;cr<c_iR;cr++)
	{
	  mD_cv_r=amD_cv_r[cr];
	  mZ_cv_r=amZ_cv_r[cr];
	  mD_cv+=mD_cv_r;
	  mZ_cv+=mZ_cv_r;
	}

      mD_cv=mD_cv/c_iR;
      mZ_cv=mZ_cv/c_iR;

      mD_cv_mean[cv][]=meanc(mD_cv);
      mD_cv_mse[cv][]=(1/iN_test)*sumsqrc(mD_test-mD_cv);
      mZ_cv_mean[cv][]=meanc(mZ_cv);
      mZ_cv_mse[cv][]=(1/iN_test)*sumsqrc(mZ_test-mZ_cv);
      vZ_cv_mean_agg[cv]=meanc(sumr(mZ_cv));
      vZ_cv_mse_agg[cv]=(1/iN_test)*sumsqrc(sumr(mZ_test)-sumr(mZ_cv));
    }


  //Report cross validation results
  println("D: True mean, Est mean, RMSE");
  println(meanc(mD_true)'~meanc(mD_cv_mean)'~sqrt(meanc(mD_cv_mse))');
  println("Z: True mean, Est mean, RMSE");
  println(meanc(mZ_true)'~meanc(mZ_cv_mean)'~sqrt(meanc(mZ_cv_mse))');
  println("Aggregate Z: True mean, Est mean, RMSE");
  println(meanc(vZ_true_agg)'~meanc(vZ_cv_mean_agg)'~sqrt(meanc(vZ_cv_mse_agg))');
 
  format(2000);
  file=fopen(sprint("results/",c_itau,"cv_interdep_r2p.dat"),"w");
  //Report cross validation results
  fprintln(file,"D: True mean, Est mean, RMSE");
  fprintln(file,"%15.5f",meanc(mD_true)'~meanc(mD_cv_mean)'~sqrt(meanc(mD_cv_mse))');
  fprintln(file,"Z: True mean, Est mean, RMSE");
  fprintln(file,"%15.5f",meanc(mZ_true)'~meanc(mZ_cv_mean)'~sqrt(meanc(mZ_cv_mse))');
  fprintln(file,"Aggregate Z: True mean, Est mean, RMSE");
  fprintln(file,"%15.5f",meanc(vZ_true_agg)'~meanc(vZ_cv_mean_agg)'~sqrt(meanc(vZ_cv_mse_agg))');
  fclose(file);
}