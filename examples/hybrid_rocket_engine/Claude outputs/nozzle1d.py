import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ---------- gas / CEA (P~19 bar, O/F~6.9) ----------
T0=3253.5; pc=2.0e6; cstar=1609.1; mu=9.79e-5; gam=1.153
cp_eq=3929.0; Pr_eq=0.463
gam_f=1.20; R=8314.46/25.54; cp_f=gam_f/(gam_f-1)*R; Pr_f=4*gam_f/(9*gam_f-5)
Dt=0.014; rt=Dt/2; Rc=rt

def mach(Ar,sup):
    f=lambda M:(1/M)*((2/(gam+1))*(1+(gam-1)/2*M*M))**((gam+1)/(2*(gam-1)))-Ar
    return 1.0 if abs(Ar-1)<1e-9 else (brentq(f,1.0001,10) if sup else brentq(f,1e-4,0.9999))

def bartz(Ar,sup,Tw,cp,Pr):
    M=mach(Ar,sup); g=1+(gam-1)/2*M*M
    sigma=1/((0.5*Tw/T0*g+0.5)**0.68*g**0.12)
    h=0.026/Dt**0.2*(mu**0.2*cp/Pr**0.6)*(pc/cstar)**0.8*(Dt/Rc)**0.1*(1/Ar)**0.9*sigma
    Taw=T0*(1+Pr**(1/3)*(gam-1)/2*M*M)/g
    return h,Taw,M

print(f"frozen cp={cp_f:.0f} Pr={Pr_f:.3f};  cp/Pr^0.6: eq={cp_eq/Pr_eq**0.6:.0f}  frozen={cp_f/Pr_f**0.6:.0f}")
print(" station        A/At   M     h_eq(Tw=300) h_eq(Tw=2000) h_fr(300) h_fr(2000)  Taw")
for name,r,sup in [("conv inlet",16.66e-3,False),("mid conv",11.8e-3,False),("throat",7e-3,False),("mid div",10.5e-3,True),("exit",14e-3,True)]:
    Ar=(r/rt)**2
    a=bartz(Ar,sup,300,cp_eq,Pr_eq);b=bartz(Ar,sup,2000,cp_eq,Pr_eq);c=bartz(Ar,sup,300,cp_f,Pr_f);d=bartz(Ar,sup,2000,cp_f,Pr_f)
    print(f" {name:12s} {Ar:5.2f} {a[2]:5.2f}   {a[0]:8.0f}    {b[0]:8.0f}     {c[0]:8.0f}  {d[0]:8.0f}  {a[1]:6.0f}")

# ---------- graphite / aluminium ----------
Tg_=np.array([300,600,1000,1500,2000,3000.]);cpg_=np.array([710,1200,1600,1900,2050,2150.]);kg_=np.array([130,80,60,45,40,38.])
cpg=lambda T:np.interp(T,Tg_,cpg_); kg=lambda T:np.interp(T,Tg_,kg_)
rho_g=1760; rho_a=2700; cp_a=896; k_a=167
ri,ro=7e-3,17.465e-3; ra_o=20.965e-3
t_burn=6.0; t_end=40.0

def run_resolved(hc,props,cpset,N=50,hmult=1.0):
    e=np.linspace(0,1,N+1)**1.6; rf=ri+(ro-ri)*e; rc=0.5*(rf[1:]+rf[:-1])
    V=np.pi*(rf[1:]**2-rf[:-1]**2)          # per metre length
    Va=np.pi*(ra_o**2-ro**2)
    def rhs(t,T):
        Tgr=T[:N]; Ta=T[N]
        k=kg(Tgr) if props else np.full(N,130.); c=cpg(Tgr) if props else np.full(N,710.)
        q=np.zeros(N); 
        # internal faces
        for j in range(N-1):
            Rf=np.log(rf[j+1]/rc[j])/(2*np.pi*k[j])+np.log(rc[j+1]/rf[j+1])/(2*np.pi*k[j+1])
            Q=(Tgr[j]-Tgr[j+1])/Rf; q[j]-=Q; q[j+1]+=Q
        # gas side: surface temperature from half-cell resistance
        if t<t_burn:
            Ts=Tgr[0]
            for _ in range(3):
                h,Taw,_=bartz(1.0,False,Ts,*cpset); h*=hmult
                Rg=1/(h*2*np.pi*ri); Rh=np.log(rc[0]/ri)/(2*np.pi*k[0])
                Q=(Taw-Tgr[0])/(Rg+Rh); Ts=Tgr[0]+Q*Rh
            q[0]+=Q
        # graphite -> contact -> aluminium
        Rh=np.log(ro/rc[-1])/(2*np.pi*k[-1]); Rc_=1/(hc*2*np.pi*ro) if hc else 0.0
        Ra=np.log((ro+ra_o)/2/ro)/(2*np.pi*k_a)
        Q=(Tgr[-1]-Ta)/(Rh+Rc_+Ra); q[-1]-=Q
        return np.append(q/(rho_g*c*V), Q/(rho_a*cp_a*Va))
    s=solve_ivp(rhs,(0,t_end),np.full(N+1,294.),method='BDF',max_step=0.05,dense_output=True)
    return s,rc
def run_lumped(hc,props,cpset,hmult=1.0):
    Vg=np.pi*(ro**2-ri**2); Va=np.pi*(ra_o**2-ro**2); rn=np.sqrt(ri*ro)
    def rhs(t,T):
        Tg,Ta=T; k=kg(Tg) if props else 130.; c=cpg(Tg) if props else 710.
        Qg=0.0
        if t<t_burn:
            h,Taw,_=bartz(1.0,False,Tg,*cpset); Qg=h*hmult*2*np.pi*ri*(Taw-Tg)   # like the current ROM: no half-cell on gas side
        Rh=np.log(ro/rn)/(2*np.pi*k); Rc_=1/(hc*2*np.pi*ro) if hc else 0.0; Ra=np.log((ro+ra_o)/2/ro)/(2*np.pi*k_a)
        Q=(Tg-Ta)/(Rh+Rc_+Ra)
        return [(Qg-Q)/(rho_g*c*Vg),Q/(rho_a*cp_a*Va)]
    return solve_ivp(rhs,(0,t_end),[294.,294.],method='BDF',max_step=0.05,dense_output=True)

def report(label,s,resolved):
    ts=np.linspace(0,t_end,2001); Y=s.sol(ts)
    Ta=Y[-1]; i6=np.searchsorted(ts,t_burn)
    surf=Y[0]
    print(f" {label:44s} T_gr,inner@6s={surf[i6]-273:6.0f}C  T_gr,outer@6s={Y[-2][i6]-273:6.0f}C  Al@6s={Ta[i6]-273:5.0f}C  Al_peak={Ta.max()-273:5.0f}C at t={ts[Ta.argmax()]:4.1f}s")

eq=(cp_eq,Pr_eq); fr=(cp_f,Pr_f)
print("\nThroat station, 1D radial, per unit length (T_aw at M=1, Bartz sigma with live T_w)")
for props in [False,True]:
    for cpset,nm in [(eq,"eq props"),(fr,"frozen props")]:
        for hc in [None,5000.,1000.]:
            tag=f"{'T-dep' if props else 'room-T'} graphite, {nm}, hc={'perfect' if hc is None else int(hc)}"
            report("LUMPED   "+tag,run_lumped(hc,props,cpset),False)
            s,_=run_resolved(hc,props,cpset); report("RESOLVED "+tag,s,True)
        print()
