function [X,Q,R,Zhat_k] = ACKF4(X0,u,Z,Q,R,Tem,BattCap)
%ACKF 此处显示有关此函数的摘要
%   Zk        本时刻测量输出
%   X0        上一时刻估计状态
%   u k-1  uk 上一时刻输入
%   Pk  本时刻协方差
n=length(X0);
m=2*n;
w=1/m;
%isACKF=1;
persistent Pk_k_1 Hk u_k_ Q_ index;
if isempty(Pk_k_1)
    Hk=0;
    Pk_k_1=diag([1e-1 1e-3 2e-4]);
    u_k_=0;
    Q_=0;
    index=0;
end
soc_local=X0(3);

Para(1)=ECM_Para(soc_local,Tem,'R0');
Para(2)=ECM_Para(soc_local,Tem,'A1');
Para(3)=ECM_Para(soc_local,Tem,'A2');
Para(4)=ECM_Para(soc_local,Tem,'B1');
Para(5)=ECM_Para(soc_local,Tem,'B2');
Para(6)=Tem;
Para(7)=BattCap;

Xi_k_1=zeros(n,m); 
Xi_k=zeros(n,m); 
Zi_k_=zeros(1,m);       %Zi，k|k-1传播容积点
Xi_k_minus=zeros(n,m);   %X*i，k|k-1  传播容积点

ksi=sqrt(m/2)*[eye(n),-eye(n)];%用与计算ksi的矩阵

%% 预测更新 #1
% [Uk_1,Sk_1,~]=svd(Pk_k_1);     %对CKF算法的协方差矩阵Pk-1进行奇异值分解
% Sk=Uk_1*sqrt(Sk_1);
% for j=1:m
%     Xi_k_1(:,j)=Sk*ksi(:,j)+X0;   %计算求容积点 
% end
% for k=1:m
%     Xi_k_minus(:,k) = State(Xi_k_1(:,k),u_k_,Para);
% %     通过状态方程传播后的容积点
% end
% Xk_minus=0;       
% for j=1:m  
%     Xk_minus=Xk_minus+w*Xi_k_minus(:,j);       %状态预测
% end
% Pk_k_1= Q;
% for j=1:m
%      Pk_k_1= w*(Xi_k_minus(:,j)-Xk_minus)*(Xi_k_minus(:,j)-Xk_minus)'+Pk_k_1;%协方差预测
% end

%% 预测更新 #2
a11=Para(2);
a22=Para(3);

A=diag([a11 a22 1]);%状态转移矩阵
Xk_minus =State(X0,u_k_,Para);    
Pk_k_1=A*Pk_k_1*A'+ Q;

%% 测量更新
[Uk_1,Sk_1,~]=svd(Pk_k_1);%对Pxx进行奇异值分解
Sk=Uk_1*sqrt(Sk_1);
for j=1:m
    Xi_k(:,j)=Sk*ksi(:,j)+Xk_minus;   %计算求容积点 
end

% Xnn=zeros(m,3);
for j=1:m      
    Zi_k_(:,j) = Measure(Xi_k(:,j),u,Para);
    % Xnn(j,1)=Xi_k(3,j);
end
% Xnn(:,2)=Tem/10*ones(m,1);
% Xnn(:,3)=u*ones(m,1)/40;
% Xnn=[Xk_minus(3) Tem/10 u/40];
% Errarr=NN_error_com(Xnn');
% Zhat_k=mean(Errarr);
% if Xk_minus(3)<0.33
%     Zhat_k=0;
% end
Zhat_k=0;
for j=1:m  
    Zhat_k=w*Zi_k_(:,j)+Zhat_k;      %量测预测
end
% u_m=(u+u_k_)/2;
% DR=0.05*tanh(u_m^2/900);
% Pzz=0.5*R+0.5*DR;
Pzz=R;
for j=1:m
    Pzz=Pzz + w*(Zi_k_(:,j)-Zhat_k)*(Zi_k_(:,j)-Zhat_k)';%协方差预测
end
Pxz=zeros(n,1);
for j=1:m
    Pxz=Pxz + w*(Xi_k(:,j)-Xk_minus)*(Zi_k_(:,j)-Zhat_k)';%互协方差预测
end
u_k_=u;
%******更新增益矩阵******%
Kk=Pxz/Pzz;
%******更新状态向量******%
X=Xk_minus+Kk*(Z-Zhat_k);
%******更新状态协方差矩阵******%
Pk_k_1=Pk_k_1-Kk*Pzz*Kk';

%% Adaptive Part
Lw=5;
% R_max=0.5;
R_min=1e-3;

if index>Lw
    tau=-0.01;
    Q_max=9e-4*exp(tau*index)+2e-7;
    Q_min=9e-7*exp(tau*index)+1e-9;

    Weight=0.92;
    Hk=Weight*Hk+(1-Weight)*(Z-Zhat_k)*(Z-Zhat_k)';
    Q_=Weight*Q_+(1-Weight)*Kk*Hk*Kk';
    % 限制方差更新
    NOR_Q=norm(Q_);
    if (NOR_Q>Q_min) && (NOR_Q<=Q_max)
      Q=Q_;
    end
    R_=Hk; 
    for j=1:m
    %R=R+w*(Zi_k_(:,j)-Zhatk_k_1)*(Zi_k_(:,j)-Zhatk_k_1)';
     R_=R_+w*(Zi_k_(:,j)-Z)*(Zi_k_(:,j)-Z)';
    end
    %  R=0.99*R+0.01*R_;
    NOR_R=norm(R_);
    if (NOR_R>=R_min)% && NOR_R<=R_max
      R=R_;
    end
    % end
    
else
    Hk=Hk+1/Lw*(Z-Zhat_k)*(Z-Zhat_k)';
end
index=index+1;

end

function X1 = State(X0,I,Para)
% 此处显示有关此函数的摘要
%   此处显示详细说明
T=0.1;
a11=Para(2);
a22=Para(3);
b1=Para(4);
b2=Para(5);
C_cell=Para(7);

eta=1;
b3=-eta*T/3600/C_cell;
A=diag([a11 a22 1]);%状态转移矩阵
B=[b1;b2;b3];%控制矩阵
X1=A*X0+B*I;
X1=min(X1,[0.8;0.7;1]);
end

function ut = Measure(x,I,Para)
% 此处显示有关此函数的摘要
%   此处显示详细说明
R0=Para(1);
Tem=Para(6);
soc=x(3);
% Tem=round(Tem/5)*5;
% Eocv=NN_ocv([soc,Tem/20,I/40]')+3.6;
Eocv=NN_ocv_3d([soc,Tem/20,I/40]')+3.6;
if soc<0
    Eocv=2.8+3*soc;
end
if soc>=1
    Eocv=4.19+3*(soc-1);
end
if soc<0.11&& soc>0
    Eocv=Tab_OCV(soc,Tem);
end

ut=Eocv-[1 1 0]*x-R0*I;
end