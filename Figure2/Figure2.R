################################################################################
#Title:  Simulations of PK models
#Author : Itziar Irurzun-Arana
#Created  : 15-Mar-2020
################################################################################
#load libraries
library(mrgsolve)
library(ggplot2)
library(dplyr)

############ PK model for Lopinavir and Ritonavir  ###############
##################################################################

## https://aac.asm.org/content/aac/55/6/2775.full.pdf

modelcode='
$PARAM @annotated
CL_L   :  4.73 : Clearance for Lopinavir (L/h)
V_L    : 55.7 : Central volume for Lopinavir (L)
KA_L   : 0.325 : Absorption rate constant (1/h)
F_L: 1 : bioavailability for Lopinavir
ALAG_L  : 0.875 : Lag time for Lopinavir (h)
CL_R   :  22 : Clearance for Ritonavir (L/h)
V_R    : 177 : Central volume for Ritonavir (L)
KA_R   : 0.74 : Absorption rate constant (1/h)
F_R: 1 : bioavailability for Ritonavir
ALAG_R  : 1.10 : Lag time for Ritonavir (h)
c1: 0.000624: RTV effect coefficient
fuo: 0.01: fraction of unbound concentration
CP_R_mean: 299: Mean concentration of Ritonavir for the 100mg dose

$CMT  @annotated
A1   : Depot compartment Lopinavir
CENT_L : Central compartment Lopinavir
A2   : Depot compartment Ritonavir
CENT_R : Central compartment Ritonavir

$MAIN
double CL_Ri = exp(log(CL_R) + ETA_CL_R);
double V_Ri = exp(log(V_R) + ETA_V_R);
double KA_Ri = exp(log(KA_R)  + ETA_KA_R);
double ALAG_Ri = exp(log(ALAG_R) + ETA_ALAG_R);
double c1i = exp(log(c1) + ETA_c1);
double CL_Li = exp(log(CL_L) + ETA_CL_L);
double V_Li = exp(log(V_L) + ETA_V_L);
double ALAG_Li = exp(log(ALAG_L) + ETA_ALAG_L);
double fu = exp(log(fuo) + ETA_fu);

ALAG_A1 = ALAG_Li;
F_A1 = F_L;
ALAG_A2 = ALAG_Ri;
F_A2 = F_R;

$OMEGA @labels ETA_CL_R ETA_V_R ETA_KA_R ETA_ALAG_R ETA_c1 ETA_CL_L ETA_V_L ETA_ALAG_L ETA_fu
0 0 0 0 0 0 0 0 0

$GLOBAL
#define CP_L ((CENT_L/V_Li)*1000)
#define CP_R ((CENT_R/V_Ri)*1000)
#define CP_L_free (CP_L*fu)

$ODE
dxdt_A1  = -KA_L*A1;
dxdt_A2  = -KA_Ri*A2;
dxdt_CENT_R = KA_Ri*A2  - (CL_Ri/V_Ri)*CENT_R;
dxdt_CENT_L = KA_L*A1  - (CL_Li*exp(-c1*(CP_R-CP_R_mean))/V_Li)*CENT_L;



$CAPTURE @annotated
CP_L : Plasma concentration for Lopinavir (mass/volume)
CP_R : Plasma concentration for Ritonavir (mass/volume)
CP_L_free : Free plasma concentration for Lopinavir (mass/volume)'

model = mcode("Lop_Rit", modelcode)


# Dose and schedule:
e_l = ev(amt = 400, ii = 12, addl = 8, time=0, cmt=1)
e_r = ev(amt = 100, ii = 12, addl = 8, time=0, cmt=3)

# First simulation
PK=as.data.frame(mrgsim(model,events=e_l+e_r,end=7*24))
PK$Regimen="Lopinavir+Ritonavir"

## Simulate Lopinavir alone witout ritonavir
PK_Lop=as.data.frame(mrgsim(model,events=e_l,end=7*24)) 
PK_Lop$Regimen="Lopinavir"


PK_Lop=rbind(PK,PK_Lop)


#pdf("Figure2A.pdf",width=4,height=3)
ggplot(data=PK_Lop,aes(x=time,y=CP_L,col=Regimen,linetype=Regimen))+geom_line(size=1.3)+
  ylab(" Lopinavir plasma concentration (ng/ml)")+xlab("Time (hours)")+
  theme_classic()+  theme(text = element_text(size = 10))+theme(legend.position="top")+
  scale_color_manual(values=c("#0066CC","orange")) 
#dev.off()

## Add variability to the model
omat(model) #variance-covariance matrix

model=omat(model,dmat(0.255, 0.38,1.14,0.487,0.283^2,0.24^2,0.552^2,0.605^2,0.0547))

set.seed(124)
PK_IIV=as.data.frame(mrgsim(model,events=e_l+e_r,end=7*24,nid=1000)) 

#Free Ritonavir:
PK_IIV$CP_R_free=PK_IIV$CP_R*0.01
PK$CP_R_free=PK$CP_R*0.01

# 5% of patients with higher Cmax values:
Cmax1=PK_IIV %>% group_by(ID) %>% summarize(Cmax_L=max(CP_L_free),Cmax_R=max(CP_R_free))

Cmax1=Cmax1[order(Cmax1$Cmax_L,decreasing = T),]

Ribbon=PK_IIV %>% group_by(time) %>% summarize(Cmax_L=max(CP_L_free),Cmax_R=max(CP_R_free),Cmin_L=min(CP_L_free),Cmin_R=min(CP_R_free))

Ribbon$Cmin_L[Ribbon$Cmin_L<0.5]=0.5
Ribbon$Cmin_R[Ribbon$Cmin_R<0.001]=0.001

#pdf("Figure2B_Free_Lopinavir_IIV_log.pdf",width=4.5,height=3.5)
ggplot()+#geom_line(size=1.5,col="gray87")+
  ylab(" Lopinavir free plasma concentration (ng/ml)")+xlab("Time (hours)")+
  geom_ribbon(data=Ribbon,aes(x=time,ymin = Cmin_L,ymax=Cmax_L),fill="gray87",alpha=0.8)+
  geom_hline(yintercept = c(3250.95),linetype="dashed",col="#3399FF",size=1.2)+
  geom_hline(yintercept = c(9810.06),linetype="dashed",col="#FF9933",size=1.2)+
  geom_hline(yintercept = (Cmax1$Cmax_L[50]+Cmax1$Cmax_L[1])/2,linetype="dashed",col="red",size=1.2)+
  geom_line(data=PK,aes(x=time,y=CP_L_free),col="#0066CC",size=1.5)+
  geom_text(aes(x = 150, y =(Cmax1$Cmax_L[50]+Cmax1$Cmax_L[1])/2+80, label = "227.5"),color = "red",size=4)+
  #geom_text(aes(x = 150, y = 4100, label = "IC50 IKr"),color = "#3399FF",size=3)+
  #geom_text(aes(x = 150, y = 12500, label = "IC50 ICaL"),color = "#FF9933",size=3)+
  theme_classic()+  theme(text = element_text(size = 11))+scale_y_continuous(trans="log10",limits=c(0.5,13500))
#dev.off()

Cmax1=Cmax1[order(Cmax1$Cmax_R,decreasing = T),]

#pdf("Figure2C_Free_Ritonavir_IIV_log.pdf",width=4.5,height=3.8)
ggplot()+
  ylab(" Free ritonavir plasma concentration (ng/ml)")+xlab("Time (hours)")+
  geom_ribbon(data=Ribbon,aes(x=time,ymin = Cmin_R,ymax=Cmax_R),fill="gray87",alpha=0.8)+
  geom_hline(yintercept = (Cmax1$Cmax_R[50]+Cmax1$Cmax_R[1])/2,linetype="dashed",col="red",size=1.2)+
  geom_line(data=PK,aes(x=time,y=CP_R_free),col="orangered",size=1.5)+
  geom_hline(yintercept = c(3717.94),linetype="dashed",col="#3399FF",size=1.1)+
  geom_hline(yintercept = c(5931.98),linetype="dashed",col="#FF9933",size=1.1)+
  geom_hline(yintercept = c(5172.82),linetype="dashed",col="#336633",size=1.1)+
  theme_classic()+
  geom_text(aes(x = 150, y =(Cmax1$Cmax_R[50]+Cmax1$Cmax_R[1])/2+10, label = "19"),color = "red",size=3)+
  # geom_text(aes(x = 140, y = 3000, label = "IC50 IKr"),color = "#3399FF",size=3)+
  # geom_text(aes(x = 140, y = 7500, label = "IC50 ICaL"),color = "#FF9933",size=3)+
  # geom_text(aes(x = 170, y = 6500, label = "IC50 INaL"),color = "#336633",size=3)+
  theme(text = element_text(size = 11))+scale_y_continuous(trans="log10",limits=c(0.001,8000))
#dev.off()

###################################################################################################
###################################################################################################

############ PK model for Chloroquine  ###########################
##################################################################

modelcode='
$PARAM @annotated
TVCL   :  54.6 : Clearance (volume/time)
TVV    : 2920 : Central volume (volume)
TVKA   : 0.943 : Absorption rate constant (1/time)
Q    : 47.2 : Inter-compartmental clearance (volume/time)
Vp   : 4700 : Peripheral volume of distribution (volume)
TVF: 1 : bioavailability
ALAG  : 0 : Lag time (time)

$CMT  @annotated
EV   : Extravascular compartment
CENT : Central compartment
PERIPH : Peripheral compartment

$MAIN
double CL = exp(log(TVCL) + ETA_CL);
double V2 = exp(log(TVV)  + ETA_V);
double KA = exp(log(TVKA)  + ETA_KA);
double V3 = exp(log(Vp)  + ETA_Vp);
double F = exp(log(TVF)  + ETA_F);


ALAG_EV = ALAG;
F_EV = F;

$OMEGA @labels ETA_CL ETA_V ETA_KA ETA_Vp ETA_F
0 0 0 0 0

$GLOBAL
#define CP (CENT/V2)

$PKMODEL ncmt = 2, depot = TRUE

$CAPTURE @annotated
CP : Plasma concentration (mass/volume)
'
#Read the model:
mod = mcode("2cmt_ev", modelcode)

#molecular weight: 319.8 g/mol
#unbound fraction in plasma: 0.39

# Administer the drug 500 mg BID for 7 days (prosed regimen for COVID19 patients according to 
# the National Health Commission of the People’s Republic of China)

e1 = ev(amt = 500, ii = 12, addl = 13, time=0, cmt=1)
#mrgsim is the function needed to simulate the model
mrgsim(mod,events=e1,end=10*24)  %>% plot(CP*1000~time) #To have ng/ml 

PK=as.data.frame(mrgsim(mod,events=e1,end=10*24))
#Modify units to ng/ml
PK$CP=PK$CP*1000
#Calculate free plasma concentrations
PK$CP_free=PK$CP*0.39

#Add IIV

mod=omat(mod,dmat(0.084^2, 0,0.662^2,0,0.192^2))

set.seed(124)
PK_IIV=as.data.frame(mrgsim(mod,events=e1,end=10*24,nid=1000)) 

PK_IIV$CP=(PK_IIV$CP*1000)
PK_IIV$CP_free=(PK_IIV$CP*0.39)

Cmax1=PK_IIV %>% group_by(ID) %>% summarize(Cmax=max(CP_free))

Cmax1=Cmax1[order(Cmax1$Cmax,decreasing = T),]

Ribbon=PK_IIV %>% dplyr::group_by(time) %>% dplyr::summarize(Cmax=max(CP_free),Cmin=min(CP_free))

#pdf("Figure2D_Free_CQ_IIV_log.pdf",width=4.6,height=4)
pcq=ggplot()+
  #geom_ribbon(aes(ymin = Cmax1$Cmax[50],ymax=Cmax1$Cmax[1]),fill="red",alpha=0.2)+
  geom_ribbon(data=Ribbon,aes(x=time,ymin = Cmin,ymax=Cmax),fill="gray87",alpha=0.8)+
  #geom_point(size=1,col="gray")+
  geom_line(data=PK,aes(x=time,y=CP_free),col="purple4",size=1.5)+
  geom_hline(yintercept = (Cmax1$Cmax[50]+Cmax1$Cmax[1])/2,linetype="dashed",col="red",size=1.2)+
  geom_hline(yintercept = c(3554.04),linetype="dashed",col="#3399FF",size=1.05)+
  geom_hline(yintercept = c(5465.96),linetype="dashed",col="orchid4",size=1.05)+
  ylab("CQ free plasma concentration (ng/ml)")+xlab("Time (hours)")+
  theme_classic()+theme(text = element_text(size = 14))+
  geom_text(aes(x =230 , y =(Cmax1$Cmax[50]+Cmax1$Cmax[1])/2+95, label = "771 ng/ml"),color = "red",size=5)+
  scale_x_continuous(limits = c(0,245))+
  scale_y_continuous(trans="log10",limits=c(10,5600))
pcq
#dev.off()

###################################################################################################
###################################################################################################

############ PK model for Azythromycin ###########################
##################################################################

#Code for a 3cmt model with extravascular administration:

modelcode='
$PARAM @annotated
CL   : 258 : Clearance (volume/time)
V    : 160 : Central volume (volume)
KA   : 0.88 : Absorption rate constant (1/time)
Q1   : 207 : First inter-compartmental clearance (volume/time)
VP1  : 1190 : First peripheral volume of distribution (volume)
Q2   :  101 : Second inter-compartmental clearance (volume/time)
VP2  : 9721 : Second peripheral volume (volume) 
ALAG : 1.45 : Lag time (time)

$CMT  @annotated
EV   : Extravascular compartment
CENT : Central compartment
PERIPH : First Peripheral compartment
PERIPH2 : Second peripheral compartment

$MAIN
double CL1 = exp(log(CL) + ETA_CL);
double V1 = exp(log(V)  + ETA_V);
double ALAG1 = exp(log(ALAG)  + ETA_ALAG);

ALAG_EV = ALAG1;

$OMEGA @labels ETA_CL ETA_V ETA_ALAG
0 0 0


$GLOBAL
#define CP (CENT/V1)
#define CT (PERIPH/VP1)
#define CT2 (PERIPH2/VP2)

$ODE
dxdt_EV= -KA*EV;
dxdt_CENT = KA*EV - ((CL1+Q1+Q2)/V1)*CENT  + (Q1/VP1)*PERIPH + (Q2/VP2)*PERIPH2;
dxdt_PERIPH = (Q1/V1)*CENT - (Q1/VP1)*PERIPH;
dxdt_PERIPH2 = Q2*CP - Q2*CT2;

$CAPTURE @annotated
CP : Plasma concentration (mass/volume)'

#Read the model:
mod = mcode("3cmt_ev", modelcode)

# Administer the drug 500 mg once daily
e1 = ev(amt = 500, ii = 24, addl = 4, time=0, cmt=1)

#Simulate population PK
PK=as.data.frame(mrgsim(mod,events=e1,end=7*24)) 
#Change units
PK$CP=PK$CP*1000
#CP is already free plasma concentrations according to the model

#Add IIV

mod=omat(mod,dmat(0.293^2, 1.68^2,0.176^2))

set.seed(124)
PK_IIV=as.data.frame(mrgsim(mod,events=e1,end=5*24,nid=1000)) 

PK_IIV$CP=PK_IIV$CP*1000

Cmax1=PK_IIV %>% group_by(ID) %>% summarize(Cmax=max(CP))

Cmax1=Cmax1[order(Cmax1$Cmax,decreasing = T),]

Ribbon=PK_IIV %>% group_by(time) %>% summarize(Cmax=max(CP),Cmin=min(CP))


#pdf("Figure2E_Free_AZ_IIV_log.pdf",width=4.6,height=4)
ggplot()+
  geom_ribbon(data=Ribbon,aes(x=time,ymin = Cmin,ymax=Cmax),fill="gray87",alpha=0.8)+
  geom_hline(yintercept = (Cmax1$Cmax[50]+Cmax1$Cmax[1])/2,linetype="dashed",col="red",size=1.05)+
  geom_line(data=PK,aes(x=time,y=CP),col="darkgoldenrod1",size=1.5)+
  geom_hline(yintercept = c(53026.20),linetype="dashed",col="#3399FF",size=1.05)+
  geom_hline(yintercept = c(66484.24),linetype="dashed",col="deeppink3",size=1.05)+
  geom_hline(yintercept = c(141656.87),linetype="dashed",col="#336633",size=1.05)+
  geom_hline(yintercept = c(352128.12),linetype="dashed",col="#CC9900",size=1.05)+
  ylab(" Unbound Azythromycin \n plasma concentration (ng/ml)")+xlab("Time (hours)")+
  theme_classic()+
  geom_text(aes(x = 110, y =(Cmax1$Cmax[50]+Cmax1$Cmax[1])/2+400, label = "768 ng/ml"),color = "red",size=3)+
  theme(text = element_text(size = 12))+scale_x_continuous(limits = c(0,130))+
  scale_y_continuous(trans="log10",limits=c(0.001,355000),breaks=c(1e-01,1e+01,1e+03,1e+05))
#dev.off()

