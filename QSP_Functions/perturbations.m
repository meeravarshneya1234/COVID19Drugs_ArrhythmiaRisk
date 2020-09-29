function [p,c]= perturbations(c,p,pert)

% Gs 
c.G.GNa  = c.G.GNa  * pert.GNa;
c.G.GNaL = c.G.GNaL * pert.GNaL;
c.G.Gto  = c.G.Gto  * pert.Gto;
c.G.GKr_ = c.G.GKr_ * pert.GKr;
c.G.GKs_ = c.G.GKs_ * pert.GKs;
c.G.GK1  = c.G.GK1  * pert.GK1;
c.G.Gncx = c.G.Gncx * pert.GNCX;
c.G.GKb  = c.G.GKb  * pert.GKb;
c.G.GpCa = c.G.GpCa * pert.GpCa;
c.G.PCa_ = c.G.PCa_ * pert.PCa;
c.G.Pnak = c.G.Pnak * pert.NaK;
c.G.PNab = c.G.PNab * pert.PNab;
c.G.PCab = c.G.PCab * pert.PCab;
c.G.SERCA_total = c.G.SERCA_total * pert.SERCA ;
c.G.RyR_total   = c.G.RyR_total   * pert.RyR ;
c.G.Leak_total  = c.G.Leak_total  * pert.Leak;
c.G.Trans_total = c.G.Trans_total * pert.Trans;

% protocols 
p.Ko = p.Ko * pert.Ko;
p.Inject = pert.Inject;
p.Cao = p.Cao * pert.Cao;
p.Nao = p.Nao * pert.Nao;

% Vs 
c.V.Vd = c.V.Vd + pert.Vd;
c.V.Vm = c.V.Vm + pert.Vm;
c.V.Vh = c.V.Vh + pert.Vh;
c.V.Vj = c.V.Vj + pert.Vj;
c.V.Vhp = c.V.Vhp + pert.Vhp;
c.V.Vjp = c.V.Vjp + pert.Vjp;
c.V.VmL = c.V.VmL + pert.VmL;
c.V.VhL = c.V.VhL + pert.VhL;
c.V.VhLp = c.V.VhLp + pert.VhLp;
c.V.Va = c.V.Va + pert.Va;
c.V.Vi = c.V.Vi + pert.Vi;
c.V.Vap = c.V.Vap + pert.Vap;
c.V.Vip= c.V.Vip + pert.Vip;
c.V.Vf= c.V.Vf + pert.Vf;
c.V.Vfca= c.V.Vfca + pert.Vfca;
c.V.Vjca= c.V.Vjca + pert.Vjca;
c.V.Vfp= c.V.Vfp + pert.Vfp;
c.V.Vfcap= c.V.Vfcap + pert.Vfcap;
c.V.Vxr= c.V.Vxr + pert.Vxr;
c.V.Vxs1= c.V.Vxs1 + pert.Vxs1;
c.V.Vxs2= c.V.Vxs2 + pert.Vxs2;
%c.V.Vxk1= 0;
c.V.Vncx= c.V.Vncx + pert.Vncx;
c.V.Vnak= c.V.Vnak + pert.Vnak;

c.p.ph = c.p.ph * pert.ph;
c.p.pfcapf = c.p.pfcapf * pert.pfcapf;
c.p.pjca = c.p.pjca * pert.pjca;