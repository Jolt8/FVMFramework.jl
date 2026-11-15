exec(open('nozzle1d.py').read().split('eq=(cp_eq')[0])
eq=(cp_eq,Pr_eq); fr=(cp_f,Pr_f)
for cps,nm in [(fr,'frozen'),(eq,'eq')]:
    report(f"LUMPED (1 node) {nm}",run_lumped(5000.,True,cps),False)
    for N in [2,3,5,10,50]:
        s,_=run_resolved(5000.,True,cps,N=N); report(f"radial nodes N={N} {nm}",s,True)
