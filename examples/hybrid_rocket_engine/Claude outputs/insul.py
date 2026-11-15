exec(open('nozzle1d.py').read().split('eq=(cp_eq')[0])
import numpy as np
from scipy.integrate import solve_ivp
Tg_=np.array([300,400,500,600,800,1000,1200,1500,2000,2500,3000.])
cpg_=4184*(0.44391+0.30795e-4*Tg_-0.61257e5/Tg_**2+0.10795e8/Tg_**3)
kg_=130*np.array([1.0,0.90,0.82,0.76,0.67,0.61,0.57,0.51,0.42,0.39,0.37])
cpg=lambda T:np.interp(T,Tg_,cpg_); kg=lambda T:np.interp(T,Tg_,kg_)
rho_p,cp_p,k_p=1250.,1400.,0.25
def run(t_ins,mode,cpset,hc=5000.,Ng=30,Np=12):
    # mode 'replace': insulation eats graphite OD; 'add': aluminium moved outward, same Al thickness
    rgo = ro-t_ins if mode=='replace' else ro
    rio = rgo+t_ins
    rao = (ra_o) if mode=='replace' else ra_o+t_ins
    eg=np.linspace(0,1,Ng+1)**1.6; rfg=ri+(rgo-ri)*eg
    rfp=np.linspace(rgo,rio,Np+1) if t_ins>0 else np.array([rgo])
    rf=np.concatenate([rfg,rfp[1:]]); n=len(rf)-1; rc=0.5*(rf[1:]+rf[:-1]); V=np.pi*(rf[1:]**2-rf[:-1]**2)
    isg=np.arange(n)<Ng; Va=np.pi*(rao**2-rio**2)
    Qin=[0.0]
    def props(T):
        k=np.where(isg,kg(T),k_p); c=np.where(isg,cpg(T),cp_p); rho=np.where(isg,rho_g,rho_p); return k,c,rho
    def rhs(t,Y):
        T=Y[:n];Ta=Y[n]; k,c,rho=props(T); q=np.zeros(n)
        Rf=np.log(rf[1:-1]/rc[:-1])/(2*np.pi*k[:-1])+np.log(rc[1:]/rf[1:-1])/(2*np.pi*k[1:])
        if t_ins>0: Rf[Ng-1]+=1/(hc*2*np.pi*rgo)   # contact at graphite/insulator
        Q=(T[:-1]-T[1:])/Rf; q[:-1]-=Q; q[1:]+=Q
        Qg=0.0
        if t<t_burn:
            Ts=T[0]
            for _ in range(3):
                h,Taw,_=bartz(1.0,False,Ts,*cpset); Rg=1/(h*2*np.pi*ri); Rh=np.log(rc[0]/ri)/(2*np.pi*k[0])
                Qg=(Taw-T[0])/(Rg+Rh); Ts=T[0]+Qg*Rh
            q[0]+=Qg
        Rh=np.log(rio/rc[-1])/(2*np.pi*k[-1]); Rc_=1/(hc*2*np.pi*rio); Ra=np.log((rio+rao)/2/rio)/(2*np.pi*k_a)
        Qa=(T[-1]-Ta)/(Rh+Rc_+Ra); q[-1]-=Qa
        return np.concatenate([q/(rho*c*V),[Qa/(rho_a*cp_a*Va), Qg]])
    s=solve_ivp(rhs,(0,t_end),np.concatenate([np.full(n+1,294.),[0.]]),method='BDF',max_step=0.05,dense_output=True)
    ts=np.linspace(0,t_end,4001); Y=s.sol(ts); Ta=Y[n]-273.15
    i6=np.searchsorted(ts,t_burn)
    t150 = ts[np.argmax(Ta>150)] if (Ta>150).any() else None
    return Ta[i6], Ta.max(), t150, Y[-1][-1]/1e3, (rgo-ri)*1e3
eq=(cp_eq,Pr_eq); fr=(cp_f,Pr_f)
print("mode     t_ins  props   graphite_wall  Q_abs(kJ/m)  Al@6s  Al_soak(40s)  t(Al>150C)")
for cps,nm in [(fr,'frozen'),(eq,'equil')]:
  for mode in ['replace','add']:
    for t in [0,1e-3,2e-3,3e-3,5e-3]:
        a6,am,t150,Q,wg=run(t,mode,cps)
        print(f"{mode:8s} {t*1e3:4.1f}mm {nm:7s} {wg:6.2f}mm   {Q:8.0f}     {a6:6.0f}  {am:6.0f}   {('%.1f s'%t150) if t150 is not None else 'never'}")
