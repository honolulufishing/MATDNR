function [gen,bus,Lnbr,trsfm,shtc,shtr,vctr,loss,sysdt]=pf(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt)

% Step 1: Converts external to internal bus numbering(non-conseutive to consecutive)
[i2e,gen,bus,Lnbr,trsfm,shtc,shtr,vctr] = ext2int(gen,bus,Lnbr,trsfm,shtc,shtr,vctr);

% Step 2: Intialize and pre-process before using power flow method
[Ybus,Yf,Yt,V0,ref,pv,pq,Sbus,tol,pfmxt] = pfinitial(gen,bus,Lnbr,trsfm,shtc,shtr,sysdt);

% Step 3: Run Newton Raphson Method
[V,converged,iterations] = matnewtonpf(Ybus,V0,ref,pv,pq,Sbus,tol,pfmxt);

% Step 4: Record the Solver Performance Information
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
    ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;

if converged == 0
    sysdt(SUCCESS) = 0;
else sysdt(SUCCESS) = 1;
end

sysdt(PFITER) = iterations;

% Step 5: Bus Solution and Line Flows 
[loss,gen,bus,Lnbr,trsfm,vctr,sysdt] = pfsoln(gen,bus,Lnbr,trsfm,vctr,sysdt,Ybus,Yf,Yt,V,ref,pv,pq,converged,iterations);

% Step 6: Bus Solution and Line Flows 
% [gen,bus,Lnbr,trsfm,shtc,shtr,vctr] = int2ext(i2e, gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt,loss);
