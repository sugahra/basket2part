#include<oxstd.h>
#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include <oxdraw.h>
#import<solvenle>
#import<maximize>
#include"myconst_indep_check.h"

static decl s_vD,s_mX,s_iN_train;
fn_Logit_j(const vgamma, const Func, const Score, const Hess)
{
  decl vexp,vlikelihood;
  vlikelihood=zeros(s_iN_train,1);
  vexp=exp(s_mX*vgamma);
  vlikelihood=s_vD.*(vexp./(ones(s_iN_train,1)+vexp))
    +(1-s_vD).*(ones(s_iN_train,1)./(ones(s_iN_train,1)+vexp));    
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
  decl time,file,mdata,mZ,vcv,mD,mY,vinsure_drug, vinsure_dental,
    mbeta,vsigma2,mgamma,
    mgamma_init,mY_all,mD_all,mX_all,cv,
    mZ_train,mY_train,mD_train,mX_train,
    mZ_test,mD_test,mX_test,ci_train,ci_test,ci,
    iN_ndrug_test,iN_ndental_test,
    mZ_ndrug_test,mD_ndrug_test,mX_ndrug_test,
    mZ_ndental_test,mD_ndental_test,mX_ndental_test,
    ci_ndrug,ci_ndental,vZ_agg_test,vD_agg_test,mCP_v,
    maxi_gamma,dopt,vexp,vexp_ndrug,vexp_ndental,vtheta_reg,
    vN_cv_test,iv,
    indental,indrug,vN_ndrug,vN_ndental,iN_train,iN_test,
    mY_test,mZ_all,cj,vgamma,iJ,mCP_v_ndrug,mCP_v_ndental,
    vY,mX,cij,vsigma2_ndrug,
    vsigma2_ndental,vcv_raw,mY_v,
    inD_test_one,inD_test_zero,
    mY_true_one,mZ_true,mD_true,
    mD_true_drug,mD_true_dental,
    mZ_true_drug,mZ_true_dental,vZ_true_agg,
    mCP_mse,mCP_mean,mCP_mean_zero,mCP_mean_one,
    mY_mean_one,mY_mse_one,mZ_mean,mZ_mse,vZ_mean_agg,vZ_mse_agg,
    mCP_mean_drug,mCP_meanchange_drug,mCP_mse_drug,mZ_meanchange_drug,mZ_mse_drug,
    mCP_mean_dental,mCP_meanchange_dental,mCP_mse_dental,mZ_meanchange_dental,mZ_mse_dental,
    mZ_v,mZ_v_ndrug,mZ_v_ndental,mY_all_adjusted,
    mZ_mean_dental,mZ_mean_drug,mindex_nan,
    iinonnan,cnonnan,vCP_mean_one_adj,vY_mean_one_adj,vY_mse_one_adj,
    vCP_mean_one_mean,vY_mean_one_mean,vY_rmse_one,
    mD_cv_mean,mD_cv_mse,mZ_cv_mean,mZ_cv_mse,vZ_cv_mean_agg,vZ_cv_mse_agg,
    mD_cv_mean_drug,mD_cv_meanchange_drug,mD_cv_mse_drug,
    mD_cv_mean_dental,mD_cv_meanchange_dental,mD_cv_mse_dental,
    mZ_cv_mean_drug,mZ_cv_meanchange_drug,mZ_cv_mse_drug,
    mZ_cv_mean_dental,mZ_cv_meanchange_dental,mZ_cv_mse_dental,
    vD_star,mD_cv,vD_drug_star,vD_drug_cv,vD_dental_star,vD_dental_cv,mZ_cv,vY_true_denom,
    vN11,vN01,vN10,mZ_cv_mean_11,mZ_cv_mse_11,mZ_cv_mean_01,mZ_cv_mse_01,mZ_cv_mean_10,mZ_cv_mse_10,
    mZ_true_11,mZ_true_10,mN11,mN10,mN01,mN00,mD_cv_r,mZ_v_r,cr,
    mD_cv_drug_r,mD_cv_dental_r,mD_cv_drug,mD_cv_dental,
    mZ_cv_drug_r,mZ_cv_dental_r,mZ_cv_drug,mZ_cv_dental;

  time=timer();

  //Load data
  file=fopen("data/mdata.dat","r");
  fscan(file,"%#m",c_iN,c_ivar_file,&mdata);
  fclose(file);
  mZ_all=mdata[][0:c_iJ-1];
  mX_all=ones(c_iN,1)~mdata[][c_iJ:c_iJ+c_iP-2];
  vinsure_dental=mdata[][c_iJ+c_iP-3];
  vinsure_drug=mdata[][c_iJ+c_iP-2];
  //Culumns for dental or drug insurance variables in X
  indental=c_iP-2;
  indrug=c_iP-1;

  //vcv=mdata[][c_iJ+c_iP-1];
  vcv_raw=mdata[][c_iJ+c_iP-1];
  vcv=floor((vcv_raw-1)/2)+1;

  mD_all=(mZ_all.>0);
  mY_all=log(mZ_all).*mD_all;
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
  //Target of prediction
  mCP_mean=zeros(c_iV_2p,c_iJ); //Mean CP in test data
  mCP_mean_zero=zeros(c_iV_2p,c_iJ); //Mean CP for those with D=1 in test data
  mCP_mean_one=zeros(c_iV_2p,c_iJ); //Mean CP for those with D=0 in test data
  mCP_mse=zeros(c_iV_2p,c_iJ); //Mean squared error of CP
  mY_mean_one=zeros(c_iV_2p,c_iJ); //Prediction mean for Y for those with D=1 in test data
  mY_mse_one=zeros(c_iV_2p,c_iJ); //Prediction error for Y for those with D=1 in test data
  mZ_mean=zeros(c_iV_2p,c_iJ); //Prediction mean for X
  mZ_mse=zeros(c_iV_2p,c_iJ); //Prediction error for Z
  vZ_mean_agg=zeros(c_iV_2p,1); //Prediction mean for aggregate expenditure
  vZ_mse_agg=zeros(c_iV_2p,1); //Prediction error for aggregate expenditure

  //Counter-factual
  mCP_mean_drug=zeros(c_iV_2p,c_iJ); //Mean CP for those changes drug insurance coverage 
  mCP_mean_dental=zeros(c_iV_2p,c_iJ); //Mean CP for those changes dental insurance coverage 
  mCP_meanchange_drug=zeros(c_iV_2p,c_iJ); //Mean CP change for those changes drug insurance coverage 
  mCP_meanchange_dental=zeros(c_iV_2p,c_iJ); //Mean CP change for those changes dental insurance coverage 
  mCP_mse_drug=zeros(c_iV_2p,c_iJ); //MSE of CP for those changes drug insurance coverage 
  mCP_mse_dental=zeros(c_iV_2p,c_iJ); //MSE of CP for those changes dental insurance coverage 
  mZ_mean_drug=zeros(c_iV_2p,c_iJ); //Mean of Z by changes drug insurance coverage 
  mZ_mean_dental=zeros(c_iV_2p,c_iJ); //Mean of Z by changes dental insurance coverage 
  mZ_meanchange_drug=zeros(c_iV_2p,c_iJ); //Mean change of Z by changes drug insurance coverage 
  mZ_meanchange_dental=zeros(c_iV_2p,c_iJ); //Mean change of Z by changes dental insurance coverage 
  mZ_mse_drug=zeros(c_iV_2p,c_iJ); //MSE for Z among those changes drug insurance coverage 
  mZ_mse_dental=zeros(c_iV_2p,c_iJ); //MSE for Z among those changes dental insurance coverage 

  // 
  mD_cv_mean=zeros(c_iV_2p,c_iJ); //Mean CP in test data
  mD_cv_mse=zeros(c_iV_2p,c_iJ); //Mean squared error of CP
  mZ_cv_mean=zeros(c_iV_2p,c_iJ); //Prediction mean for X
  mZ_cv_mse=zeros(c_iV_2p,c_iJ); //Prediction error for Z
  vZ_cv_mean_agg=zeros(c_iV_2p,1); //Prediction mean for aggregate expenditure
  vZ_cv_mse_agg=zeros(c_iV_2p,1); //Prediction error for aggregate expenditure
  mD_cv_mean_drug=zeros(c_iV_2p,c_iJ); //Mean D for those changes drug insurance coverage 
  mD_cv_mean_dental=zeros(c_iV_2p,c_iJ); //Mean D for those changes dental insurance coverage 
  mD_cv_meanchange_drug=zeros(c_iV_2p,c_iJ); //Mean D change for those changes drug insurance coverage 
  mD_cv_meanchange_dental=zeros(c_iV_2p,c_iJ); //Mean D change for those changes dental insurance coverage 
  mD_cv_mse_drug=zeros(c_iV_2p,c_iJ); //MSE of D for those changes drug insurance coverage 
  mD_cv_mse_dental=zeros(c_iV_2p,c_iJ); //MSE of D for those changes dental insurance coverage 
  mZ_cv_mean_drug=zeros(c_iV_2p,c_iJ); //Mean of Z by changes drug insurance coverage 
  mZ_cv_mean_dental=zeros(c_iV_2p,c_iJ); //Mean of Z by changes dental insurance coverage 
  mZ_cv_meanchange_drug=zeros(c_iV_2p,c_iJ); //Mean change of Z by changes drug insurance coverage 
  mZ_cv_meanchange_dental=zeros(c_iV_2p,c_iJ); //Mean change of Z by changes dental insurance coverage 
  mZ_cv_mse_drug=zeros(c_iV_2p,c_iJ); //MSE for Z among those changes drug insurance coverage 
  mZ_cv_mse_dental=zeros(c_iV_2p,c_iJ); //MSE for Z among those changes dental insurance coverage 

  mZ_cv_mean_11=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mean_01=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mean_10=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mse_11=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mse_01=zeros(c_iV_2p,c_iJ); 
  mZ_cv_mse_10=zeros(c_iV_2p,c_iJ); 

  mZ_true_11=mZ_true_10=zeros(c_iV_2p,c_iJ);
  mN11=mN10=mN01=mN00=zeros(c_iV_2p,c_iJ);
  //Summary for true data
  mY_true_one=mZ_true=mD_true=mD_true_drug=mD_true_dental=
    mZ_true_drug=mZ_true_dental=zeros(c_iV_2p,c_iJ);
  vZ_true_agg=zeros(c_iV_2p,1); 

  //Index for test data, which takes unity when all test data show D_j=0
  mindex_nan=zeros(c_iV_2p,c_iJ);

  vN_ndrug=vN_ndental=zeros(c_iV_2p,1);
  ///Size of test data
  vN_cv_test=zeros(c_iV_2p,1);
  for(ci=0;ci<c_iN;ci++)
    {
      iv=vcv[ci]-1;
      vN_cv_test[iv]++;
    } 

  //Set estimates above for the initial values of MLE iteration for gamma
  file=fopen("results/gamma_indep.dat","r");
  fscan(file,"%#m",c_iP,c_iJ,&mgamma_init);
  fclose(file);
  mbeta=zeros(c_iP,c_iJ);
  vsigma2=zeros(c_iP,1);

  //Loop start

  for(cv=0;cv<c_iV_2p;cv++)
    {
      println("cv=",cv);
      println("time lapsed=",timespan(time));
      //Separating training and test data
      iN_test=vN_cv_test[cv];
      iN_train=c_iN-iN_test;

      mY_train=mD_train=zeros(iN_train,c_iJ);
      mX_train=zeros(iN_train,c_iP);
      mY_test=mZ_test=mD_test=zeros(iN_test,c_iJ);
      mX_test=zeros(iN_test,c_iP);

      ci_train=ci_test=iN_ndrug_test=iN_ndental_test=0;
      for(ci=0;ci<c_iN;ci++)
	{
	  if(vcv[ci]==cv+1)
	    {
	      mX_test[ci_test][]=mX_all[ci][];
	      mZ_test[ci_test][]=mZ_all[ci][];
	      mY_test[ci_test][]=mY_all_adjusted[ci][];
	      mD_test[ci_test][]=mD_all[ci][];
	      if(vinsure_drug[ci]==1)
		{
		  iN_ndrug_test++;
		}
	      if(vinsure_dental[ci]==1)
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
	      ci_train++;
	    }
	}
      vY_true_denom=sumc(mD_test);
      vY_true_denom=vY_true_denom+(vY_true_denom.==0);
      mY_true_one[cv][]=sumc(mD_test.*mY_test) ./ vY_true_denom;
      mD_true[cv][]= meanc(mD_test);
      mZ_true[cv][]= meanc(mZ_test);
      vZ_true_agg[cv]=meanc(sumr(mZ_test)); //Aggregate medical expenditure
      mZ_v=mY_v=mZ_cv=mD_cv=zeros(iN_test,c_iJ);

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
	  if(vcv[ci]==cv+1 && vinsure_drug[ci]==1)
	    {
	      mZ_ndrug_test[ci_ndrug][]=mZ_all[ci][];
	      mD_ndrug_test[ci_ndrug][]=mD_all[ci][];
	      mX_ndrug_test[ci_ndrug][]=mX_all[ci][];
	      ci_ndrug++;
	    }
	  if(vcv[ci]==cv+1 && vinsure_dental[ci]==1)
	    {
	      mZ_ndental_test[ci_ndental][]=mZ_all[ci][];
	      mD_ndental_test[ci_ndental][]=mD_all[ci][];
	      mX_ndental_test[ci_ndental][]=mX_all[ci][];
	      ci_ndental++;
	    }	  
	}
      mX_ndrug_test[][indrug]=ones(iN_ndrug_test,1);
      mX_ndental_test[][indental]=ones(iN_ndental_test,1);
      mZ_v=mY_v=mZ_cv=mD_cv=zeros(iN_test,c_iJ);
      mD_cv_drug=mZ_cv_drug=zeros(iN_ndrug_test,c_iJ);
      mD_cv_dental=mZ_cv_dental=zeros(iN_ndental_test,c_iJ);

      mD_true_dental[cv][]= meanc(mD_ndental_test);
      mZ_true_dental[cv][]= meanc(mZ_ndental_test);
      mD_true_drug[cv][]= meanc(mD_ndrug_test);
      mZ_true_drug[cv][]= meanc(mZ_ndrug_test);

      //Estimation
      for(cj=0;cj<c_iJ;cj++)
	{
	  inD_test_one=sumc(mD_test[][cj]);
	  inD_test_zero=iN_test-inD_test_one;
	  println("cj=",cj);
	  s_vD=mD_train[][cj];
	  s_mX=mX_train;
	  s_iN_train=iN_train;
	  //Logit
	  vgamma=mgamma_init[][cj];
	  maxi_gamma=MaxBFGS(fn_Logit_j, &vgamma, &dopt,0,1);
	  vexp=exp(mX_test*vgamma);	  
	  //Regression
	  iJ=sumc(s_vD);
	  vY=zeros(iJ,1);
	  mX=zeros(iJ,c_iP);
	  cij=0;

	  for(ci_train=0;ci_train<iN_train;ci_train++)
	    {
	      if(s_vD[ci_train]==1)
		{
		  vY[cij]=mY_train[ci_train][cj];
		  mX[cij][]=s_mX[ci_train][];
		  cij++;
		}
	    }

	  vtheta_reg=fn_Reg(vY,mX,1);
	  //Last argument takes 1 when we only need estimates
	  mbeta[][cj]=vtheta_reg[0:c_iP-1];
	  vsigma2[cj]=vtheta_reg[c_iP]; 

	  
	  //CV targets
	  mD_cv_r=mZ_v_r=zeros(iN_test,c_iR);
	  mD_cv_drug_r=mZ_cv_drug_r=zeros(vN_ndrug[cv],c_iR);
	  mD_cv_dental_r=mZ_cv_dental_r=zeros(vN_ndental[cv],c_iR);
	  for(cr=0;cr<c_iR;cr++)
	    {
	      vD_star=mX_test*vgamma+ranextremevalue(iN_test,1,0,1);
	      mD_cv_r[][cr]=(vD_star.>=0);
	      mZ_v_r[][cr]=mD_cv_r[][cr].*exp(mX_test*mbeta[][cj]+sqrt(vsigma2[cj])*rann(iN_test,1));

	      vD_drug_star=mX_ndrug_test*vgamma+ranextremevalue(vN_ndrug[cv],1,0,1);
	      mD_cv_drug_r[][cr]=(vD_drug_star.>=0);
	      mZ_cv_drug_r[][cr]=mD_cv_drug_r[][cr].*exp(mX_ndrug_test*mbeta[][cj]+sqrt(vsigma2[cj])*rann(vN_ndrug[cv],1));

	      vD_dental_star=mX_ndental_test*vgamma+ranextremevalue(vN_ndental[cv],1,0,1);
	      mD_cv_dental_r[][cr]=(vD_dental_star.>=0);
	      mZ_cv_dental_r[][cr]=mD_cv_dental_r[][cr].*exp(mX_ndental_test*mbeta[][cj]+sqrt(vsigma2[cj])*rann(vN_ndental[cv],1));
	    }
	  mD_cv[][cj]=meanr(mD_cv_r);
	  mZ_cv[][cj]=meanr(mZ_v_r);
	  mD_cv_mean[cv][cj]=meanc(mD_cv[][cj]);
	  mD_cv_mse[cv][cj]=(1/iN_test)*sumsqrc(mD_test[][cj]-mD_cv[][cj]);
	  mZ_cv_mean[cv][cj]=meanc(mZ_cv[][cj]);
	  mZ_cv_mse[cv][cj]=(1/iN_test)*sumsqrc(mZ_test[][cj]-mZ_cv[][cj]);

	  mD_cv_drug[][cj]=meanr(mD_cv_drug_r);
	  mZ_cv_drug[][cj]=meanr(mZ_cv_drug_r);
	  mZ_cv_mean_drug[cv][cj]=meanc(mZ_cv_drug[][cj]);
	  mZ_cv_meanchange_drug[cv][cj]=meanc(mZ_cv_drug[][cj]-mZ_ndrug_test[][cj]);
	  mZ_cv_mse_drug[cv][cj]=(1/vN_ndrug[cv])
	    *sumsqrc(
		     (mZ_cv_drug[][cj]-mZ_ndrug_test[][cj])
		     -mZ_cv_meanchange_drug[cv][cj]);

	  mD_cv_dental[][cj]=meanr(mD_cv_dental_r);
	  mZ_cv_dental[][cj]=meanr(mZ_cv_dental_r);
	  mZ_cv_mean_dental[cv][cj]=meanc(mZ_cv_dental[][cj]);
	  mZ_cv_meanchange_dental[cv][cj]=meanc(mZ_cv_dental[][cj]-mZ_ndental_test[][cj]);
	  mZ_cv_mse_dental[cv][cj]=(1/vN_ndental[cv])
	    *sumsqrc(
		     (mZ_cv_dental[][cj]-mZ_ndental_test[][cj])
		     -mZ_cv_meanchange_dental[cv][cj]);
	}

      vZ_cv_mean_agg[cv]=meanc(sumr(mZ_cv));
      vZ_cv_mse_agg[cv]=(1/iN_test)*sumsqrc(sumr(mZ_test)-sumr(mZ_cv)) ;
    }

  //Report cross validation results
  println("D: True mean, Est mean, RMSE");
  println(meanc(mD_true)'~meanc(mD_cv_mean)'~sqrt(meanc(mD_cv_mse))');
  println("Z: True mean, Est mean, RMSE");
  println(meanc(mZ_true)'~meanc(mZ_cv_mean)'~sqrt(meanc(mZ_cv_mse))');
  println("Aggregate Z: True mean, Est mean, RMSE");
  println(meanc(vZ_true_agg)'~meanc(vZ_cv_mean_agg)'~sqrt(meanc(vZ_cv_mse_agg))');

  //Report counter-factual result
  println("Counterfactual");
  println("Drug insurance, D: True mean, Est mean,  RMSE");
  println(<0;0;0;0;0;0;0;0;1;0;0;0;0>~meanc(mD_true_drug)'~meanc(mD_cv_mean_drug)'~sqrt(meanc(mD_cv_mse_drug))');
  println("Drug insurance, Z: True mean, Est mean, Est mean-change, RMSE");
  println(<0;0;0;0;0;0;0;0;1;0;0;0;0>~meanc(mZ_true_drug)'~meanc(mZ_cv_mean_drug)'~meanc(mZ_cv_meanchange_drug)'~sqrt(meanc(mZ_cv_mse_drug))');
  println("Dental insurance, D: True mean, Est mean, Est mean-change, RMSE");
  println(<0;0;0;0;0;0;0;0;0;1;0;0;0>~meanc(mD_true_dental)'~meanc(mD_cv_mean_dental)'~sqrt(meanc(mD_cv_mse_dental))');
  println("Dental insurance, Z: True mean, Est mean, Est mean-change, RMSE");
  println(<0;0;0;0;0;0;0;0;0;1;0;0;0>~meanc(mZ_true_dental)'~meanc(mZ_cv_mean_dental)'~meanc(mZ_cv_meanchange_dental)'~sqrt(meanc(mZ_cv_mse_dental))');

  format(2000);
  file=fopen("results/cv_indep_check_r2p.dat","w");
  //Report cross validation results
  fprintln(file,"D: True mean, Est mean, RMSE");
  fprintln(file,"%15.5f",meanc(mD_true)'~meanc(mD_cv_mean)'~sqrt(meanc(mD_cv_mse))');
  fprintln(file,"Z: True mean, Est mean, RMSE");
  fprintln(file,"%15.5f",meanc(mZ_true)'~meanc(mZ_cv_mean)'~sqrt(meanc(mZ_cv_mse))');
  fprintln(file,"Aggregate Z: True mean, Est mean, RMSE");
  fprintln(file,"%15.5f",meanc(vZ_true_agg)'~meanc(vZ_cv_mean_agg)'~sqrt(meanc(vZ_cv_mse_agg))');

  //Report counter-factual result
  fprintln(file,"Counterfactual");
  fprintln(file,"Drug insurance, D: True mean, Est mean,  RMSE");
  fprintln(file,"%15.5f",<0;0;0;0;0;0;0;0;1;0;0;0;0>~meanc(mD_true_drug)'~meanc(mD_cv_mean_drug)'~sqrt(meanc(mD_cv_mse_drug))');
  fprintln(file,"Drug insurance, Z: True mean, Est mean, Est mean-change, RMSE");
  fprintln(file,"%15.5f",<0;0;0;0;0;0;0;0;1;0;0;0;0>~meanc(mZ_true_drug)'~meanc(mZ_cv_mean_drug)'~meanc(mZ_cv_meanchange_drug)'~sqrt(meanc(mZ_cv_mse_drug))');
  fprintln(file,"Dental insurance, D: True mean, Est mean, RMSE");
  fprintln(file,"%15.5f",<0;0;0;0;0;0;0;0;0;1;0;0;0>~meanc(mD_true_dental)'~meanc(mD_cv_mean_dental)'~sqrt(meanc(mD_cv_mse_dental))');
  fprintln(file,"Dental insurance, Z: True mean, Est mean, Est mean-change, RMSE");
  fprintln(file,"%15.5f",<0;0;0;0;0;0;0;0;0;1;0;0;0>~meanc(mZ_true_dental)'~meanc(mZ_cv_mean_dental)'~meanc(mZ_cv_meanchange_dental)'~sqrt(meanc(mZ_cv_mse_dental))');
  fclose(file);
}
