from nozzle1d import report
from nozzle1d import run_resolved
exec(open('nozzle1d.py').read().split('eq=(cp_eq')[0])
eq=(cp_eq,Pr_eq); fr=(cp_f,Pr_f)
import numpy as np
# grid check
for N in [25,50,100]:
    s,_=run_resolved(5000.,True,fr,N=N); report(f"grid N={N}",s,True)
# lower T0 (c* efficiency 0.9 -> T0 ~ 0.81 T_CEA), Bartz x0.5
T0_save=T0
for T0v,hm,lab in [(T0_save,1.0,"T0=CEA, Bartz x1"),(0.81*T0_save,1.0,"T0=0.81 CEA (eta_c*=0.9)"),(0.81*T0_save,0.5,"T0=0.81 CEA, Bartz x0.5")]:
    T0=T0v
    s,_=run_resolved(5000.,True,fr,hmult=hm); report("RESOLVED frozen hc=5000 "+lab,s,True)
