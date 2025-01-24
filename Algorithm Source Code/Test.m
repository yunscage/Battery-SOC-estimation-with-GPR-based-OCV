% subplot(2,3,1);
% xlim([0 6500]);
% subplot(2,3,4);
% xlim([0 6500]);
% 
% subplot(2,3,2);
% xlim([0 6500]);
% subplot(2,3,5);
% xlim([0 6500]);

% 
% gprMdl = fitrsvm(train_in',train_out');
% ypred = resubPredict(gprMdl);
% plot(train_out,'k.'); hold on;
% plot(ypred);

L=length(x_soc2);
Q=1e-5;
R=1e-2;
Pk=1e-3;
x_fused=zeros(1,L);
x_fused(1)=0.6;
B=-0.1/(3600*C_cell);
x_minus=0;
for k=1:L-1
    x_minus=x_fused(k)+B*I(k);
    Pk=Pk+Q;
    Kk=Pk/(Pk+R);
    x_fused(k+1)=x_minus+Kk*(x_soc2(k+1)-x_minus);
    Pk=Pk-Kk*Pk;
end

x_ackf3(3,:)=x_fused;

e1=(soc_ref-x_ackf(3,:));
e2=(soc_ref-x_soc2);
e3=(soc_ref-x_ackf3(3,:));
e4=(soc_ref-x_ackf4(3,:));

SE1=(e1.^2);
SE2=(e2.^2);
SE3=(e3.^2);
SE4=(e4.^2);


Plot_view;