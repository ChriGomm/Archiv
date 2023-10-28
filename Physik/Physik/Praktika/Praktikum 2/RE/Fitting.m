G1 = @(a,f1)a(1)*a(2)*f1./(a(1)^2*f1.^2+a(2)^2*(a(1)*a(3)*f1.^2-1).^2).^(1/2);
g = G1;
a03 = [4 22000 22/10^9];
options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',2000000,'MaxIterations',2000000,'StepTolerance',10^(-14),'FunctionTolerance',10^(-14));
lb = [0 0 0];
ub = [100000 60000 100000];
a = lsqcurvefit(g,a03,f1,Z1,lb,ub,options);
gZ = a(1)*a(2)*f1./(a(1)^2*f1.^2+a(2)^2*(a(1)*a(3)*f1.^2-1).^2).^(1/2);
scatter(f1,Z1)
hold on
plot(f1,gZ)
F = a(1)*a(2)*y/(a(1)^2*y^2+a(2)^2*(a(1)*a(3)*y^2-1)^2)^(1/2);
fr2 = vpa(solve(diff(F)==0,y));