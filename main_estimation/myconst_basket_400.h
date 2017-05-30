//tau
const decl c_itau=400;
//Constants
const decl c_iJ=13;
const decl c_iP=29;
const decl c_ivar_file=c_iJ+c_iP; // # variables in data file

//Sample size 
const decl c_iN=23147;

//For cross validation
///#Subsamples in cross varidation
const decl c_iV_2p=50;
const decl c_iR=50;
//Number of parameters
const decl c_iq_indep=(2*c_iP+1)*c_iJ;

//For soft clustering
const decl c_iM=4; //Maximum number of non-zero element
const decl c_iK=23; //Total number of clusters
const decl c_iKm_max=4; //Maximum number of clusters (wrt. m)
const decl c_iK_max_i=8; //Maximum number of clusters (wrt. i)
