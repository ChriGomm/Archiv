H1 = @(a)sum((Z2-sqrt((a(1)^2+a(2)^2*f2.^2)./(a(1)^2*a(3)^2*f2.^2+(f2.^2*a(2)*a(3)-1).^2))).^2);
a02 = [10000 2 22/10^9];
options = optimset('MaxFunEvals',20000,'MaxIter',20000);
a = fminsearch(H1,a02,options);
Hz = sqrt((a(1)^2+a(2)^2*f2.^2)/(a(1)^2*a(3)^2*f2.^2+(f2.^2*a(2)*a(3)-1).^2));
scatter(f2,Z2)
hold on
scatter(f2,Hz)