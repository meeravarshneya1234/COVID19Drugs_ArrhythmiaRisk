function currents = calculate_currents(p,c,state_vars)
%--------------------------------------------------------------------------
                       %% -- calculate_currents.m -- %%
% Description: calculate currents for a single cell 

% Inputs:
% --> p - [struct] model parameters  
% --> c - [struct] model parameters varied in the population  
% --> state_vars  - [double array] state variables  

% Outputs: 
% --> currents - [struct] 
% -------------------------------------------------------------------------
%%
[V, Nai, Nass, Ki, Kss, Cai, Cass, Cansr, Cajsr, m, hf, hs, j, hsp, jp, mL, hL,...
    hLp, a, iF, iS, ap, iFp, iSp, d, ff, fs, fcaf, fcas, jca, nca, ffp, fcafp,...
    xrf, xrs, xs1, xs2, xk1, Jrelnp, Jrelp, CaMKt] = deal(state_vars{:});

%CaMK constants
KmCaMK=0.15;
aCaMK=0.05;
bCaMK=0.00068;
CaMKo=0.05;
KmCaM=0.0015;
%update CaMK
CaMKb=CaMKo.*(1.0-CaMKt)./(1.0+KmCaM./Cass);
CaMKa=(CaMKb+CaMKt)*p.K_CaMKa;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reversal potentials
ENa=(p.R*p.T/p.F).*log(p.Nao./Nai);
EK=(p.R*p.T/p.F).*log(p.Ko./Ki);
PKNa=0.01833;
EKs=(p.R*p.T/p.F).*log((p.Ko+PKNa*p.Nao)./(Ki+PKNa.*Nai));

%convenient shorthand calculations
vffrt=V*p.F*p.F/(p.R*p.T);
vfrt=V*p.F/(p.R*p.T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate INa
Ahf=0.99;
Ahs=1.0-Ahf;
h=Ahf*hf+Ahs*hs;
hp=Ahf*hf+Ahs*hsp;
fINap=(1.0./(1.0+KmCaMK./CaMKa));
currents.INa=c.G.GNa.*(V-ENa).*m.^3.0.*((1.0-fINap).*h.*j+fINap.*hp.*jp);

%% calculate INaL
fINaLp=(1.0./(1.0+KmCaMK./CaMKa));
currents.INaL=c.G.GNaL.*(V-ENa).*mL.*((1.0-fINaLp).*hL+fINaLp.*hLp);

%% calculate Ito
AiF=1.0./(1.0+exp((V-213.6)./151.2));
AiS=1.0-AiF;
i=AiF.*iF+AiS.*iS;
ip=AiF.*iFp+AiS.*iSp;
fItop=(1.0./(1.0+KmCaMK./CaMKa));
currents.Ito = c.G.Gto.*(V-EK).*((1.0-fItop).*a.*i+fItop.*ap.*ip);

%% calculate ICaL, ICaNa, ICaK
Aff=0.6;
Afs=1.0-Aff;
f=Aff*ff+Afs*fs;
Afcaf=0.3+0.6./(1.0+exp((V-10.0)/10.0));
Afcas=1.0-Afcaf;
fca=Afcaf.*fcaf+Afcas.*fcas;
fp=Aff.*ffp+Afs.*fs;
fcap=Afcaf.*fcafp+Afcas.*fcas;
PhiCaL=4.0*vffrt.*(Cass.*exp(2.0*vfrt)-0.341*p.Cao)./(exp(2.0*vfrt)-1.0);
PhiCaNa=1.0*vffrt.*(0.75*Nass.*exp(1.0*vfrt)-0.75*p.Nao)./(exp(1.0*vfrt)-1.0);
PhiCaK=1.0*vffrt.*(0.75*Kss.*exp(1.0*vfrt)-0.75*p.Ko)./(exp(1.0*vfrt)-1.0);
PCap=1.1*c.G.PCa_;
PCaNa=0.00125*c.G.PCa_;
PCaK=3.574e-4*c.G.PCa_;
PCaNap=0.00125*PCap;
PCaKp=3.574e-4*PCap;
fICaLp=(1.0./(1.0+KmCaMK./CaMKa));
ICa_L=(1.0-fICaLp).*c.G.PCa_.*PhiCaL.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp*PCap.*PhiCaL.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca);
ICaNa=(1.0-fICaLp).*PCaNa.*PhiCaNa.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp*PCaNap.*PhiCaNa.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca);
ICaK=(1.0-fICaLp).*PCaK.*PhiCaK.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp*PCaKp.*PhiCaK.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca);
currents.ICaL = ICa_L + ICaNa + ICaK;

%% calculate IKr
Axrf=1.0./(1.0+exp((V+54.81)/38.21));
Axrs=1.0-Axrf;
xr=Axrf.*xrf+Axrs.*xrs;
rkr=1.0./(1.0+exp((V+55.0)/75.0))*1.0./(1.0+exp((V-10.0)/30.0));
currents.IKr=c.G.GKr_*sqrt(p.Ko/5.4).*xr.*rkr.*(V-EK);

%% calculate IKs
KsCa=1.0+0.6./(1.0+(3.8e-5./Cai).^1.4);
currents.IKs=c.G.GKs_.*KsCa.*xs1.*xs2.*(V-EKs);

%% calculate IK1 
rk1=1.0./(1.0+exp((V+105.8-2.6*p.Ko)./9.493));
currents.IK1 = c.G.GK1.*sqrt(p.Ko).*rk1.*xk1.*(V-EK);

%% calculate INaCa - Main 80% 
kna1=15.0;
kna2=5.0;
kna3=88.12;
kasymm=12.5;
wna=6.0e4;
wca=6.0e4;
wnaca=5.0e3;
kcaon=1.5e6;
kcaoff=5.0e3;
qna=0.5224;
qca=0.1670;
hca=exp((qca*(V)*p.F)./(p.R*p.T));
hna=exp((qna*(V)*p.F)./(p.R*p.T));
h1=1+Nai./kna3.*(1+hna);
h2=(Nai.*hna)./(kna3.*h1);
h3=1.0./h1;
h4=1.0+Nai./kna1.*(1+Nai./kna2);
h5=Nai.*Nai./(h4.*kna1.*kna2);
h6=1.0./h4;
h7=1.0+p.Nao./kna3.*(1.0+1.0./hna);
h8=p.Nao./(kna3.*hna.*h7);
h9=1.0./h7;
h10=kasymm+1.0+p.Nao./kna1.*(1.0+p.Nao./kna2);
h11=p.Nao*p.Nao/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12.*p.Cao.*kcaon;
k2=kcaoff;
k3p=h9.*wca;
k3pp=h8.*wnaca;
k3=k3p+k3pp;
k4p=h3.*wca./hca;
k4pp=h2.*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6.*Cai.*kcaon;
k7=h5.*h2.*wna;
k8=h8.*h11.*wna;
x1=k2.*k4.*(k7+k6)+k5.*k7.*(k2+k3);
x2=k1.*k7.*(k4+k5)+k4.*k6.*(k1+k8);
x3=k1.*k3.*(k7+k6)+k8.*k6.*(k2+k3);
x4=k2.*k8.*(k4+k5)+k3.*k5.*(k1+k8);
E1=x1./(x1+x2+x3+x4);
E2=x2./(x1+x2+x3+x4);
E3=x3./(x1+x2+x3+x4);
E4=x4./(x1+x2+x3+x4);
KmCaAct=150.0e-6;
allo=1.0./(1.0+(KmCaAct./Cai).^2.0);
zna=1.0;
JncxNa=3.0.*(E4.*k7-E1.*k8)+E3.*k4pp-E2.*k3pp;
JncxCa=E2.*k2-E1.*k1;
zca=2.0;

currents.INaCa_i = 0.8.*c.G.Gncx.*allo.*(zna.*JncxNa+zca.*JncxCa);

%% calculate INaCa_ss
h1=1+Nass./kna3.*(1+hna);
h2=(Nass.*hna)./(kna3.*h1);
h3=1.0./h1;
h4=1.0+Nass./kna1.*(1+Nass./kna2);
h5=Nass.*Nass./(h4.*kna1.*kna2);
h6=1.0./h4;
h7=1.0+p.Nao./kna3.*(1.0+1.0./hna);
h8=p.Nao./(kna3.*hna.*h7);
h9=1.0./h7;
h10=kasymm+1.0+p.Nao./kna1.*(1+p.Nao./kna2);
h11=p.Nao*p.Nao/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12*p.Cao*kcaon;
k2=kcaoff;
k3p=h9.*wca;
k3pp=h8.*wnaca;
k3=k3p+k3pp;
k4p=h3.*wca./hca;
k4pp=h2.*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6.*Cass.*kcaon;
k7=h5.*h2.*wna;
k8=h8.*h11.*wna;
x1=k2.*k4.*(k7+k6)+k5.*k7.*(k2+k3);
x2=k1.*k7.*(k4+k5)+k4.*k6.*(k1+k8);
x3=k1.*k3.*(k7+k6)+k8.*k6.*(k2+k3);
x4=k2.*k8.*(k4+k5)+k3.*k5.*(k1+k8);
E1=x1./(x1+x2+x3+x4);
E2=x2./(x1+x2+x3+x4);
E3=x3./(x1+x2+x3+x4);
E4=x4./(x1+x2+x3+x4);
KmCaAct=150.0e-6;
allo=1.0./(1.0+(KmCaAct./Cass).^2.0);
JncxNa=3.0.*(E4.*k7-E1.*k8)+E3.*k4pp-E2.*k3pp;
JncxCa=E2.*k2-E1.*k1;
currents.INaCa_ss =  0.2.*c.G.Gncx.*allo.*(zna.*JncxNa+zca.*JncxCa);

%% calculate INaK
k1p=949.5;
k1m=182.4;
k2p=687.2;
k2m=39.4;
k3p=1899.0;
k3m=79300.0;
k4p=639.0;
k4m=40.0;
Knai0=9.073;
Knao0=27.78;
delta=-0.1550;
Knai=Knai0.*exp((delta.*(V).*p.F)./(3.0*p.R*p.T));
Knao=Knao0.*exp(((1.0-delta).*(V).*p.F)./(3.0*p.R*p.T));
Kki=0.5;
Kko=0.3582;
MgADP=0.05;
MgATP=9.8;
Kmgatp=1.698e-7;
H=1.0e-7;
eP=4.2;
Khp=1.698e-7;
Knap=224.0;
Kxkur=292.0;

P=eP./(1.0+H./Khp+Nai./Knap+Ki./Kxkur);
a1=(k1p.*(Nai./Knai).^3.0)./((1.0+Nai./Knai).^3.0+(1.0+Ki./Kki).^2.0-1.0);
b1=k1m*MgADP;
a2=k2p;
b2=(k2m.*(p.Nao./Knao).^3.0)./((1.0+p.Nao./Knao).^3.0+(1.0+p.Ko./Kko).^2.0-1.0);
a3=(k3p.*(p.Ko./Kko).^2.0)./((1.0+p.Nao./Knao).^3.0+(1.0+p.Ko./Kko).^2.0-1.0);
b3=(k3m.*P.*H)./(1.0+MgATP./Kmgatp);
a4=(k4p.*MgATP./Kmgatp)./(1.0+MgATP./Kmgatp);
b4=(k4m.*(Ki./Kki).^2.0)./((1.0+Nai./Knai).^3.0+(1.0+Ki./Kki).^2.0-1.0);
x1=a4.*a1.*a2+b2.*b4.*b3+a2.*b4.*b3+b3.*a1.*a2;
x2=b2.*b1.*b4+a1.*a2.*a3+a3.*b1.*b4+a2.*a3.*b4;
x3=a2.*a3.*a4+b3.*b2.*b1+b2.*b1.*a4+a3.*a4.*b1;
x4=b4.*b3.*b2+a3.*a4.*a1+b2.*a4.*a1+b3.*b2.*a1;
E1=x1./(x1+x2+x3+x4);
E2=x2./(x1+x2+x3+x4);
E3=x3./(x1+x2+x3+x4);
E4=x4./(x1+x2+x3+x4);
zk=1.0;
JnakNa=3.0.*(E1.*a3-E2.*b3);
JnakK=2.0.*(E4.*b1-E3.*a1);
currents.INaK=c.G.Pnak.*(zna.*JnakNa+zk.*JnakK);

%% calculate IKb
xkb=1.0./(1.0+exp(-(V-14.48)./18.34));
currents.IKb=c.G.GKb.*xkb.*(V-EK);

%% calculate INab
currents.INab=c.G.PNab.*vffrt.*(Nai.*exp(vfrt)-p.Nao)./(exp(vfrt)-1.0);

%% calculate ICab
currents.ICab=c.G.PCab.*4.0.*vffrt.*(Cai.*exp(2.0.*vfrt)-0.341.*p.Cao)./(exp(2.0.*vfrt)-1.0);

%% calculate IpCa
currents.IpCa=c.G.GpCa.*Cai./(0.0005+Cai);

%% calculate ryanodione receptor calcium induced calcium release from the jsr
fJrelp=(1.0./(1.0+KmCaMK./CaMKa));
currents.Jrel= c.G.RyR_total.*((1.0-fJrelp).*Jrelnp+fJrelp.*Jrelp );

%% calculate SERCA and Leak 
%calculate serca pump, ca uptake flux
Jupnp=0.004375.*Cai./(Cai+0.00092);
Jupp=2.75*0.004375.*Cai./(Cai+0.00092-0.00017);
% if strcmp(p.celltype,'epi')==1
%     Jupnp=Jupnp*1.3;
%     Jupp=Jupp*1.3;
% end
fJupp=(1.0./(1.0+KmCaMK./CaMKa));
currents.Jleak = c.G.Leak_total*0.0039375*Cansr/15.0;
currents.Jup= c.G.SERCA_total.*((1.0-fJupp).*Jupnp+fJupp.*Jupp) - currents.Jleak;

%% calculate tranlocation flux
ttr = 100.0 ;
currents.Jtr=c.G.Trans_total.*(Cansr-Cajsr)./ttr;
