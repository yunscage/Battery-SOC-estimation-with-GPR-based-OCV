clear;
% The clear for function  is used to clear the persistent value. 
% When delete the clear of the function of "ACKF,...", the persistent value will be affected by previous run.
% clear ACKF ACKF3 ACKF4 ACKF4_hysis
PreTime=datetime('now');
 %%  选择循环号 cnt_j 温度   cnt_degc
cnt_degc=2;  % 1 2 4 5
cnt_j=1;  %  1 UDDS; 3 LA92; 4 US06 
if cnt_j==2
    disp('Stop!');
end
StandTem=[-20 -10 0 10 25];
disp(['Degc: ' num2str(StandTem(cnt_degc))]);
load('PathName.mat');
load('Network.mat')

% load('GPR.mat');
% % 
% gprmdl=GPR{cnt_degc};
filepath=[Path_1{cnt_degc} Path_2{cnt_degc,cnt_j}];
load(filepath);
Cycle=Path_2{cnt_degc,cnt_j}(20:24);
if Cycle(end)=='_'
    Cycle(end)=' ';
end

%% Pre-Process
I=-meas.Current';
xindexc=ZeroIndex(I,0);
I=I(1:xindexc);

z=meas.Voltage(1:xindexc)';

TemBatt=meas.Battery_Temp_degC(1:xindexc)';
Tem=mean(TemBatt(end-500:end));
C_cell=ECM_Para(0,Tem,'Cap');
soc_ref=(meas.Ah(1:xindexc)')/C_cell+1;
BattCap=C_cell;

% XTest=[z-3.6;I./40;TemBatt./10];%(:,numTimeStepsTrain:end);
XTest=[z-3.7;I./40;TemBatt./20];
GRUnet=resetState(GRUnet);
YPred = predict(GRUnet,XTest);
x_soc2=(YPred+1)./2;
step=length(soc_ref);
t=1:step;
t=t./10;


%% *********************  CKF  ******************* %
n=3;
Q1=diag([1e-3,1e-3,1e-2]);
Q3=diag([1e-3,1e-3,1e-2]);
Q4=Q1;
Q_hys=Q1;
Rec_Q1=zeros(1,step);
Nor_Q=zeros(1,step);
Nor_Q(1)=0.01;
R1=0.1*eye(1);
R3=diag(0.1);
R4=R1;
R_hys=R1;

Xhat=[0;0;0.5];         %初始值
x_ackf=zeros(n,step);    %存储ackf的估计值；
x_ackf(:,1)=Xhat;
x_ackf3=x_ackf;
x_ackf4=x_ackf;
x_ackf_hysis=x_ackf;

y_sim=zeros(1,step);
y_sim(1)=z(1);
y_sim3=y_sim;
y_sim4=y_sim;

Rec_Q1(1)=Q1(3,3);
rec=zeros(1,step-1);
rec2=rec;
ParaRec=[];
% isted=mean(I);
% disp([Cycle 'mean current = ' num2str(isted) 'A']);
MHV=zeros(1,step);
for cnt=1:step
    MHV(cnt)=Tab_MHV(soc_ref(cnt),TemBatt(cnt));
end
[HysisV]= hystVolt(I,MHV,BattCap);
for i=2:step
    [x_ackf(:,i),Q1,R1(i),y_sim(i-1),~] = ACKF(x_ackf(:,i-1),I(:,i-1),z(:,i-1),Q1,R1(i-1),TemBatt(i),BattCap);

    [x_ackf3(:,i),Q3,R3(i),y_sim3(i-1)] = ACKF3(x_ackf3(:,i-1),I(:,i-1),z(:,i-1),Q3,R3(i-1),TemBatt(i),BattCap);

    [x_ackf4(:,i),Q4,R4(i),y_sim3(i-1)] = ACKF4(x_ackf4(:,i-1),I(:,i-1),z(:,i-1),Q4,R4(i-1),TemBatt(i),BattCap);

    [x_ackf_hysis(:,i),Q_hys,R_hys,~] = ACKF4_hysis(x_ackf_hysis(:,i-1),I(:,i-1),z(:,i-1),Q_hys,R_hys,TemBatt(i),HysisV(i-1),BattCap);

 end
y_sim(end)=y_sim(end-1);
y_sim4(:,end)=y_sim4(:,end-1);
%% 

% DNN_H = resetState(DNN_H);
% Ypred2=predict(DNN_H,XTest);
% Ypred2=(Ypred2+1)./2;

L=length(x_soc2);
Q=1e-3;
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

e1=(soc_ref-x_ackf(3,:));
e2=(soc_ref-x_soc2);
e3=(soc_ref-x_ackf3(3,:));
e4=(soc_ref-x_ackf4(3,:));
e5=(soc_ref-x_ackf_hysis(3,:));

SE1=(e1.^2);
SE2=(e2.^2);
SE3=(e3.^2);
SE4=(e4.^2);
SE5=e5.^2;

%% 绘制估计结果
NColor1='#E36255';
Ncolor2='#9DC877';
NColor3='#FFA500';
figure;
set(gcf, 'Position', [200 200 320 280]);
subplot(1,3,1);
hold on; % 保持所有后续绘图在同一图表中
Lengths=601;
plot(t(1:Lengths), soc_ref(1:Lengths) * 100, 'k', 'linewidth', 1);
plot(t(1:Lengths), x_soc2(1:Lengths) * 100, 'linewidth', 0.87, 'color', Ncolor2);
plot(t(1:Lengths), x_ackf(3,1:Lengths) * 100, 'linewidth', 1, 'color', NColor1);
plot(t(1:Lengths), x_ackf3(3,1:Lengths)*100,'linewidth',1,'color',NColor3);
plot(t(1:Lengths), x_ackf4(3,1:Lengths) * 100, 'b', 'linewidth', 1);
plot(t(1:Lengths), x_ackf_hysis(3,1:Lengths) * 100, 'y', 'linewidth', 1);
hold off;
xlim([0,60]);
ylim([40 105]);
subplot(1,3,2);
hold on;
plot(t,soc_ref*100,'k','linewidth',1);
plot(t,x_soc2*100,'linewidth',0.87,'color',Ncolor2);
plot(t,x_ackf(3,:)*100,'linewidth',1,'color',NColor1);

plot(t,x_ackf3(3,:)*100,'linewidth',1,'color',NColor3);
plot(t,x_ackf4(3,:)*100,'b','linewidth',1);
plot(t, x_ackf_hysis(3,:) * 100, 'y', 'linewidth', 1);
% legend('Reference','ACKF','RNN','Prop-2D','Prop-3D');
ylabel('SOC (%)');
xlabel('Time (s)');
title([Cycle,'Degc: ' num2str(StandTem(cnt_degc))]);
leng=length(x_soc2)/10;
xlim([0,leng]);
box on;

subplot(1,3,3);
hold on;
plot(t,e2*100,'linewidth',0.87,'color',Ncolor2);
plot(t,e1*100,'linewidth',1,'color',NColor1);

plot(t,e3*100,'linewidth',1,'color',NColor3);
plot(t,e4*100,'b','linewidth',1);
plot(t,e5*100,'y','linewidth',1);
% plot(t,RSE4,'b');
ylabel('Error (%)');
% legend('ECM','RNN','Prop-2D','Prop-3D');
xlabel('Time (s)');
box on;
xlim([0,leng]);
ylim([-7.5,7.5])
set(gcf, 'position',[200,200,900,200]);

K=200;
RMSE1=100*sqrt(mean(SE1(K:end)));
RMSE2=100*sqrt(mean(SE2(K:end)));
RMSE3=100*sqrt(mean(SE3(K:end)));
RMSE4=100*sqrt(mean(SE4(K:end)));
RMSE5=100*sqrt(mean(SE5(K:end)));
disp(['*********Cycle : ' Cycle ' ********']);
disp([' ACKF      RMSE1 : ' num2str(RMSE1) '%']);
disp([' GRU       RMSE2 : ' num2str(RMSE2) '%']);
disp([' Hybrid   RMSE3 : ' num2str(RMSE3) '%']);
disp([' ACKF_Hysis   RMSE4 : ' num2str(RMSE5) '%']);
disp([' Prop       RMSE4 : ' num2str(RMSE4) '%']);


CostTime=datetime('now')-PreTime;
disp(CostTime);

%% 函数部分
function xindex=ZeroIndex(y0,offset)
% 去除长时间电流为0的数据段
LEN=length(y0);
K_cnt=0;
xindex=LEN;
for index=LEN:-1:ceil(0.5*LEN)
    if (y0(index)==0)
        K_cnt=K_cnt+1;
    else
        if (K_cnt>500)
             xindex=index+offset;
            break;
        end
    end
end
end
