lIA3 = log(IA3./10^6);
kUG3 = UG3-IA3./10^6.*MB3;
mdl3 = fitlm(kUG32,lIA32);
IA32 = [IA3(14:15); IA3(17:19)]./10^6;
S1 = 5./MB1./10^6;
S32 = [S3(14:15); S3(17:19)];
A3 = (S32+0.0005.*IA32)./IA32./kUG32;
SUG3 = (0.0005+0.0005.*UG3)+(S3+0.0005.*IA3./10^6).*MB3;
SUG32 = [SUG3(14:15); SUG3(17:19)];
B3 = abs(SUG32.*lIA32./kUG32.^2);
syF3 = (sum(A3)+sum(B3))/5;
Re3 = -10.778478882036925.*kUG3+0.693330509537266;
scatter(kUG3,lIA3)
hold on
plot(kUG3,Re3)
ylabel('${\bar{I}_1}$','interpreter','latex', 'FontWeight','bold')
legend('Measurement Data','Regression')
xlabel('${U_G\;in\;V}$','interpreter','latex', 'FontWeight','bold')
U1=mean(U21);
g1=mean(gam(79:80));
RFe=mean(RFe1(79:80));
LM=mean(LM1(79:80));
phi=mean(phi1(79:80));
I10=mean(I11(79:80));
