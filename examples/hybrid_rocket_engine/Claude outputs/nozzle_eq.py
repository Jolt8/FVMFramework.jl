"""Isentropic nozzle expansion: (A) constant-gamma ideal gas, (B1) equilibrium SP sweep, (B2) frozen SP sweep.
Thermo: GRI-Mech 3.0 NASA-7 polynomials. Gibbs equilibrium by element-potential method."""
import numpy as np
from scipy.optimize import fsolve, brentq

Ru=8.314462618; Pref=1e5
# name: (Tmid, low[7], high[7], {element:count}, MW g/mol)
S={
'O'  :(1000,[3.1682671,-3.27931884e-03,6.64306396e-06,-6.12806624e-09,2.11265971e-12,2.91222592e+04,2.05193346],[2.56942078,-8.59741137e-05,4.19484589e-08,-1.00177799e-11,1.22833691e-15,2.92175791e+04,4.78433864],{'O':1},15.999),
'O2' :(1000,[3.78245636,-2.99673416e-03,9.84730201e-06,-9.68129509e-09,3.24372837e-12,-1063.94356,3.65767573],[3.28253784,1.48308754e-03,-7.57966669e-07,2.09470555e-10,-2.16717794e-14,-1088.45772,5.45323129],{'O':2},31.998),
'H'  :(1000,[2.5,7.05332819e-13,-1.99591964e-15,2.30081632e-18,-9.27732332e-22,2.54736599e+04,-0.446682853],[2.50000001,-2.30842973e-11,1.61561948e-14,-4.73515235e-18,4.98197357e-22,2.54736599e+04,-0.446682914],{'H':1},1.008),
'OH' :(1000,[3.99201543,-2.40131752e-03,4.61793841e-06,-3.88113333e-09,1.3641147e-12,3615.08056,-0.103925458],[3.09288767,5.48429716e-04,1.26505228e-07,-8.79461556e-11,1.17412376e-14,3858.657,4.4766961],{'O':1,'H':1},17.007),
'H2' :(1000,[2.34433112,7.98052075e-03,-1.9478151e-05,2.01572094e-08,-7.37611761e-12,-917.935173,0.683010238],[3.3372792,-4.94024731e-05,4.99456778e-07,-1.79566394e-10,2.00255376e-14,-950.158922,-3.20502331],{'H':2},2.016),
'H2O':(1000,[4.19864056,-2.0364341e-03,6.52040211e-06,-5.48797062e-09,1.77197817e-12,-3.02937267e+04,-0.849032208],[3.03399249,2.17691804e-03,-1.64072518e-07,-9.7041987e-11,1.68200992e-14,-3.00042971e+04,4.9667701],{'H':2,'O':1},18.015),
'CO' :(1000,[3.57953347,-6.1035368e-04,1.01681433e-06,9.07005884e-10,-9.04424499e-13,-1.4344086e+04,3.50840928],[2.71518561,2.06252743e-03,-9.98825771e-07,2.30053008e-10,-2.03647716e-14,-1.41518724e+04,7.81868772],{'C':1,'O':1},28.010),
'CO2':(1000,[2.35677352,8.98459677e-03,-7.12356269e-06,2.45919022e-09,-1.43699548e-13,-4.83719697e+04,9.90105222],[3.85746029,4.41437026e-03,-2.21481404e-06,5.23490188e-10,-4.72084164e-14,-4.8759166e+04,2.27163806],{'C':1,'O':2},44.009),
'N2' :(1000,[3.298677,1.4082404e-03,-3.963222e-06,5.641515e-09,-2.444854e-12,-1020.8999,3.950372],[2.92664,1.4879768e-03,-5.68476e-07,1.0097038e-10,-6.753351e-15,-922.7977,5.980528],{'N':2},28.014),
'NO' :(1000,[4.2184763,-4.638976e-03,1.1041022e-05,-9.3361354e-09,2.803577e-12,9844.623,2.2808464],[3.2606056,1.1911043e-03,-4.2917048e-07,6.9457669e-11,-4.0336099e-15,9920.9746,6.3693027],{'N':1,'O':1},30.006),
}
names=list(S); E=['C','H','O','N']
A=np.array([[S[s][3].get(e,0) for s in names] for e in E],float)
MWs=np.array([S[s][4] for s in names])/1000
def coef(s,T): Tm,lo,hi,_,_=S[s]; return np.array(hi if T>=Tm else lo)
def cp_R(T): return np.array([np.dot(coef(s,T)[:5],[1,T,T**2,T**3,T**4]) for s in names])
def h_RT(T): return np.array([(lambda a:a[0]+a[1]*T/2+a[2]*T**2/3+a[3]*T**3/4+a[4]*T**4/5+a[5]/T)(coef(s,T)) for s in names])
def s_R(T): return np.array([(lambda a:a[0]*np.log(T)+a[1]*T+a[2]*T**2/2+a[3]*T**3/3+a[4]*T**4/4+a[6])(coef(s,T)) for s in names])

# propellant elements per kg: N2O oxidizer, (CH2)n fuel, O/F from table row
OF=6.9393939393939394; P0=1918367.3469387756; T0=3253.5120568821453
mox=OF/(1+OF); mf=1/(1+OF)
b=np.array([mf/0.014027, 2*mf/0.014027, mox/0.044013, 2*mox/0.044013])  # mol element / kg

def equil(T,P,guess=[None]):
    g=h_RT(T)-s_R(T)  # g°/RT
    def F(x):
        lam=x[:4]; lnN=x[4]
        n=np.exp(lnN+A.T@lam-g-np.log(P/Pref))
        return np.append(A@n-b, n.sum()-np.exp(lnN))
    x0=guess[0] if guess[0] is not None else np.array([-20,-10,-15,-12,np.log(38.)])
    x,info,ier,msg=fsolve(F,x0,full_output=True,xtol=1e-13)
    if ier!=1 or np.max(np.abs(F(x)))>1e-9:
        for trial in range(30):
            x,info,ier,msg=fsolve(F,x0+np.random.randn(5),full_output=True,xtol=1e-13)
            if ier==1 and np.max(np.abs(F(x)))<1e-9: break
    guess[0]=x
    return np.exp(x[4]+A.T@x[:4]-g-np.log(P/Pref))   # mol/kg

def props(T,P,n):
    N=n.sum(); X=n/N
    h=Ru*T*np.dot(n,h_RT(T))                                    # J/kg
    sm=s_R(T)-np.log(np.maximum(X,1e-300))-np.log(P/Pref)
    s=Ru*np.dot(n,sm)                                            # J/kg/K
    MW=1/N; rho=P*MW/(Ru*T)
    return h,s,rho,MW
def eq_state(T,P):
    n=equil(T,P); return n,*props(T,P,n)

n0,h0,s0,rho0,MW0=eq_state(T0,P0)
Y=n0*MWs; 
print("=== Chamber validation vs your CEA table row (19.18 bar, O/F 6.94, T fixed at 3253.5 K) ===")
print(f"MW: mine {MW0*1000:.2f}  CEA 25.54")
untr=sum(Y[i] for i,s in enumerate(names) if s in ('O','H','OH','NO'))
print(f"untracked (O,H,OH,NO) mass frac: mine {untr:.4f}  CEA 0.0288")
cea={'N2':0.56808,'O2':0.007383,'CO2':0.12187,'CO':0.18141,'H2O':0.11735,'H2':0.003896}
for s,v in cea.items():
    i=names.index(s); print(f"  {s:4s} renormalized Y: mine {Y[i]/(1-untr):.4f}  CEA {v:.4f}")
dT=1.0
h_p=eq_state(T0+dT,P0)[1]; h_m=eq_state(T0-dT,P0)[1]
cp_eq=(h_p-h_m)/(2*dT); cp_fr=Ru*np.dot(n0,cp_R(T0))
print(f"cp_eq: mine {cp_eq:.0f}  CEA 3929   | cp_frozen: {cp_fr:.0f}")

# ---- isentropic sweeps ----
_Tprev=[T0]
def T_at_s_eq(P):
    lo=_Tprev[0]-400; hi=_Tprev[0]+20
    T=brentq(lambda T: eq_state(T,P)[2]-s0, lo, hi, xtol=1e-6)
    _Tprev[0]=T; return T
def T_at_s_fr(P):
    def sfr(T): return props(T,P,n0)[1]
    return brentq(lambda T: sfr(T)-s0, 500, T0+50, xtol=1e-6)
Ps=P0*np.concatenate([np.linspace(0.999,0.62,60),np.linspace(0.62,0.52,80),np.logspace(np.log10(0.52),np.log10(0.025),80)])
res={}
for mode in ["eq","fr"]:
    _Tprev[0]=T0
    out=[]
    for P in Ps:
        if mode=='eq':
            T=T_at_s_eq(P); n,h,s,rho,MW=eq_state(T,P)
        else:
            T=T_at_s_fr(P); h,s,rho,MW=props(T,P,n0)
        u=np.sqrt(2*(h0-h)); out.append((P,T,u,rho*u))
    out=np.array(out); res[mode]=out
    Gt=out[:,3].max(); it=out[:,3].argmax()
    _Tprev[0]=T0
    print(f"\n[{mode}] throat: P/P0={out[it,0]/P0:.4f}  T={out[it,1]:.0f} K  c*={P0/Gt:.1f} m/s" + ("  (CEA eq c* 1609.1)" if mode=='eq' else ""))

# equilibrium gamma_s at chamber
dP=P0*1e-4
_Tprev[0]=T0+5; Tp=T_at_s_eq(P0+dP); rp=eq_state(Tp,P0+dP)[3]; _Tprev[0]=T0+5; Tm=T_at_s_eq(P0-dP); rm=eq_state(Tm,P0-dP)[3]
print(f"gamma_s chamber: mine {np.log((P0+dP)/(P0-dP))/np.log(rp/rm):.4f}  CEA 1.1534")

# ---- approach A: constant gamma ----
g=1.153404265555152
def mach_A(Ar,sup):
    f=lambda M:(1/M)*((2/(g+1))*(1+(g-1)/2*M*M))**((g+1)/(2*(g-1)))-Ar
    return 1.0 if abs(Ar-1)<1e-12 else (brentq(f,1.000001,20) if sup else brentq(f,1e-6,0.999999))
def at_area(mode,Ar,sup):
    out=res[mode]; Gt=out[:,3].max(); it=out[:,3].argmax()
    seg=out[:it+1] if not sup else out[it:]
    ar=Gt/seg[:,3]
    order=np.argsort(ar)
    T=np.interp(np.log(Ar),np.log(ar[order]),seg[order,1]); P=np.interp(np.log(Ar),np.log(ar[order]),seg[order,0]); u=np.interp(np.log(Ar),np.log(ar[order]),seg[order,2])
    return T,P,u
print("\n=== Static T (K) vs area ratio ===   A=const gamma 1.153 | B1=equilibrium | B2=frozen")
print(" branch  A/At     T_A    T_B1   T_B2   |  T_aw_A  T_aw_B1  T_aw_B2   (r=0.774)  | P_A   P_B1 (bar)")
r=0.463**(1/3)
for sup,lst in [(False,[5.66,3,2,1.5,1.25,1.1]),(True,[1.0,1.1,1.25,1.5,2,2.5,3,4])]:
    for Ar in lst:
        M=mach_A(Ar,sup); TA=T0/(1+(g-1)/2*M*M); PA=P0*(TA/T0)**(g/(g-1))
        TB1,PB1,_=at_area('eq',Ar,sup); TB2,PB2,_=at_area('fr',Ar,sup)
        aw=lambda T:T+r*(T0-T)
        print(f" {'sup' if sup else 'sub'}    {Ar:5.2f}  {TA:6.0f} {TB1:6.0f} {TB2:6.0f}   |  {aw(TA):6.0f}  {aw(TB1):6.0f}   {aw(TB2):6.0f}             | {PA/1e5:5.2f} {PB1/1e5:5.2f}")
# exit Isp check
TB1,PB1,uB1=at_area('eq',4.0,True); Gt=res['eq'][:,3].max()
print(f"\nIsp_vac (eq, eps=4): mine {(uB1+PB1*4/Gt)/9.80665:.1f} s   CEA 264.1 s")
