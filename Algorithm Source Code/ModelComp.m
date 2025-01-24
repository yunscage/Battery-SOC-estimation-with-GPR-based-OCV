load('..\DataProcess\PathName.mat');
cnt_degc=1;
cnt_j=4;
filepath=[Path_1{cnt_degc} Path_2{cnt_degc,cnt_j}];
load(filepath);
Cycle=Path_2{cnt_degc,cnt_j}(20:24);
if Cycle(end)=='_'
    Cycle(end)=' ';
end
disp(Cycle);
I=-meas.Current';
xindexc=ZeroIndex(I,0)-500;
I=I(1:xindexc);

z=meas.Voltage(1:xindexc)';

TemBatt=meas.Battery_Temp_degC(1:xindexc)';
Tem=mean(TemBatt(end-500:end));
C_cell=ECM_Para(0,Tem,'Cap');
soc_ref=(meas.Ah(1:xindexc)')/C_cell+1;
BattCap=C_cell;

Leng=length(I);
X0=[0;0;1];
x1=zeros(3,Leng);
x1(:,1)=X0;
x2=x1;
y1=zeros(1,Leng);
y2=y1;
for cnt=2:Leng
    Tem=TemBatt(:,cnt-1);
    soc_local=x1(3,cnt-1);
    Para=GetPara(soc_local,Tem,BattCap);
    % ECM
    x1(:,cnt)=State1(x1(:,cnt-1),I(cnt-1),Para);
    y1(cnt-1)=Measure1(x1(:,cnt-1),I(cnt-1),Para);
    % Comp
    x2(:,cnt)=State1(x2(:,cnt-1),I(cnt-1),Para);
    y2(cnt-1)=Measure2(x1(:,cnt-1),I(cnt-1),Para);
end
y1(Leng)=Measure1(x1(:,Leng),I(Leng),Para);
y2(Leng)=Measure2(x1(:,Leng),I(Leng),Para);
figure;
subplot(2,1,1);
plot(z,'k'); hold on;
plot(y1,'r');
plot(y2,'b');
legend('Measure','ECM','Advanced');
% xlim([1,32000]);
subplot(2,1,2);
plot(z-y1,'r'); hold on;
plot(z-y2,'b');
err=z-y2;
for cnt=1:length(err)
    if abs(err(cnt))>0.05
        err(cnt)=0.5*err(cnt);
    end
end
plot(err,'b');
% xlim([1,32000]);
figure;
% subplot(2,1,2);
binWidth = 0.002;
xlim([-0.1 0.1]); 
hold on;
histogram(z-y2,'FaceAlpha', 0.5, 'BinWidth', binWidth, ...
                      'EdgeColor', 'none','Normalization', 'probability');
histogram(z-y1,'FaceAlpha', 0.5, 'BinWidth', binWidth, ...
                      'EdgeColor', 'none','Normalization', 'probability');
xlabel('Model error (V)');
legend('Comp','ECM')
mean1=mean(z-y1)
mean2=mean(z-y2)
var1=std(z-y1)
var2=std(z-y2)




%% 函数部分
function X1 = State1(X0,I,Para)
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
end

function ut = Measure1(x,I,Para)
% 此处显示有关此函数的摘要
%   此处显示详细说明
R0=Para(1);
Tem=Para(6);
soc=x(3);
Eocv=Tab_OCV(soc,Tem);
ut=Eocv-[1 1 0]*x-R0*I;
end

function ut = Measure2(x,I,Para)
% 此处显示有关此函数的摘要
R0=Para(1);
Tem=Para(6);
soc=x(3);
Eocv=NN_ocv([soc,Tem/20,I/40]')+3.6;
if soc<0
    Eocv=2.8+3*soc;
end
if soc>=1
    Eocv=4.19+3*(soc-1);
end
ut=Eocv-[1 1 0]*x-R0*I;
end

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

function Para=GetPara(soc_local,Tem,BattCap)
    Para(1)=ECM_Para(soc_local,Tem,'R0');
    Para(2)=ECM_Para(soc_local,Tem,'A1');
    Para(3)=ECM_Para(soc_local,Tem,'A2');
    Para(4)=ECM_Para(soc_local,Tem,'B1');
    Para(5)=ECM_Para(soc_local,Tem,'B2');
    Para(6)=Tem;
    Para(7)=BattCap;
end
