set more off
use ../raw_data/H155.dta

//ID
gen id=DUPERSID
gen hhid=DUID
gen memberid=PID

//Demographic characteristics
//Region
drop if REGION53==-1
gen northea=0
gen midwest=0
gen south=0
gen west=0
replace northea=1 if REGION53==1
replace midwest=1 if REGION53==2
replace south=1 if REGION53==3
replace west=1 if REGION53==4

///Base
drop if AGE53X==-1
gen age=AGE53X 

gen male=.
replace male=1 if SEX==1
replace male=0 if SEX==2

//Race
gen white=0
gen black=0
gen hispanic=0
gen asian=0
gen otrace=0
gen multirace=0
replace white=1 if RACEV1X==1 & RACETHX!=1
replace black=1 if RACEV1X==2 & RACETHX!=1
replace hispanic=1 if RACETHX==1
replace asian=1 if RACEV1X==4 & RACETHX!=1
replace otrace=1 if ( RACEV1X==3 | RACEV1X==5 ) & RACETHX!=1
replace multirace=1 if RACEV1X==6
drop if multirace==1

//Marital Status
gen marry=0
gen widow=0
gen divorce=0
gen single=0
gen otmarry=0
replace marry=1 if MARRY53X==1
replace widow=1 if MARRY53X==2
replace divorce=1 if MARRY53X==3 | MARRY53X==4
replace single=1 if MARRY53X==5 | MARRY53X==6
replace otmarry=1 if MARRY53X<0 | MARRY53X>=7
drop if otmarry==1

//Educationq
gen eduyears=.
replace eduyears=EDRECODE if EDRECODE>=0
drop if eduyears==.
gen edulhs=0
gen eduhs=0
gen edusomec=0
gen educolle=0
replace edulhs=1 if EDRECODE<=12 & EDRECODE>=0 
replace eduhs=1 if EDRECODE==13 
replace edusomec=1 if EDRECODE==14 
replace educolle=1 if EDRECODE>=15 

//Income
gen pincome=TTLP12X
gen fincome=FAMINC12
gen pov1=0
gen pov2=0
gen pov3=0
gen pov4=0
gen pov5=0

replace pov1=1 if POVCAT12==1
replace pov2=1 if POVCAT12==2
replace pov3=1 if POVCAT12==3
replace pov4=1 if POVCAT12==4
replace pov5=1 if POVCAT12==5

//Health

//BMI
drop if BMINDX53<0
gen bmi=BMINDX53
gen underwe=0
gen normalwe=0
gen overwe=0
gen obese=0
replace underwe=1 if bmi<18.5
replace normalwe=1 if bmi>=18.5 & bmi<25
replace overwe=1 if bmi>=25 & bmi<30
replace obese=1 if bmi>=30 

//Smoking:
///Different from Finkelstein et al (2009)
///Note that this is asked only in the 4/2 th round
drop if ADSMOK42<0
gen smoke=0
replace smoke=1 if ADSMOK42==1

//Insurance
gen uninsure=0
gen privatei=0
gen medicaid=0
gen medicare=0
gen otherins=0
replace uninsure=1 if UNINS12==1
replace medicaid=1 if MCDEV12==1
replace medicare=1 if MCREV12==1
replace privatei=1 if PRVEV12==1
replace otherins=1 if OPAEV12==1 | OPBEV12==1

drop if PMEDIN31==-1 | DENTIN31==-1 | ///
  PMEDIN42==-1 | DENTIN42==-1 | ///
  PMEDIN53==-1 | DENTIN53==-1 
gen dentali=0
gen drugi=0
replace dentali=1 if  DENTIN31==1 & DENTIN42==1 & DENTIN53==1 
replace drugi=1 if PMEDIN31==1 & PMEDIN42==1 & PMEDIN53==1 

//Medical purchase
gen yall=TOTEXP12
gen yinpat=IPTEXP12
gen yemer=ERFEXP12
gen yphis=OPVEXP12+OPSEXP12+OBDEXP12
gen ychiro=AMCEXP12
gen yambu=AMNEXP12
gen yopt=AMEEXP12
gen yassist=AMAEXP12
gen ytherapy=AMTEXP12
gen ydrug=RXEXP12
gen ydental=DVTEXP12
gen yhome=HHAEXP12+HHNEXP12
gen yvision=VISEXP12
gen yequip=OTHEXP12

//Wear seat belt
gen seatbelt=SEATBE53

//Sampling weight
gen pweight=PERWT12F

keep if age>=18
keep id hhid memberid ///
  bmi underwe normalwe overwe obese ///
  yall yinpat yemer yphis ychiro yambu yopt yassist ///
  ytherapy ydrug ydental yhome yvision yequip ///
  age male white black hispanic asian otrace ///
  marry widow divorce single ///
  eduyears edulhs eduhs edusomec educolle ///
  pincome fincome pov1 pov2 pov3 pov4 pov5 ///
  smoke ///
  northea midwest south west ///
  uninsure privatei medicaid medicare otherins ///
  dentali drugi seatbelt pweight

order id hhid memberid ///
  bmi underwe normalwe overwe obese ///
  yall yinpat yemer yphis ychiro yambu yopt yassist ///
  ytherapy ydrug ydental yhome yvision yequip ///
  age male white black hispanic asian otrace ///
  marry widow divorce single ///
  eduyears edulhs eduhs edusomec educolle ///
  pincome fincome pov1 pov2 pov3 pov4 pov5 ///
  smoke ///
  northea midwest south west ///
  uninsure privatei medicaid medicare otherins ///
  dentali drugi seatbelt pweight

saveold ../data_final/2012meps_old, replace
save ../data_final/2012meps, replace
outfile using ../data_final/2012meps.csv, wide comma replace
