exec(open('nozzle1d.py').read().split('eq=(cp_eq')[0])
eq=(cp_eq,Pr_eq); fr=(cp_f,Pr_f)
for cps,nm in [(fr,'frozen'),(eq,'eq')]:
    s,_=run_resolved(5000.,True,cps,N=1); report(f"1 node WITH gas-side half-cell R {nm}",s,True)
    s,_=run_resolved(5000.,True,cps,N=50); report(f"N=50 reference {nm}",s,True)
# energy check: heat absorbed vs final stored, resolved N=50 frozen
