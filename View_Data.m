figure
subplot(3,1,1);
plot(meas.Time,meas.Voltage,'r')
hold on;
ylabel('Voltage');
subplot(3,1,2);
plot(meas.Time,meas.Ah,'r'); 
% Ahcont=max(meas.Ah)-min(meas.Ah);
ylabel('Ah');
% plot(meas.Time,meas.Current,'r'); 
% ylabel('Current');
subplot(3,1,3);
plot(meas.Time,meas.Battery_Temp_degC,'r');
hold on;
ylabel('Temperature');
Tem=mean(meas.Battery_Temp_degC);

%% OCV curve get
[AhMIN,K1]=min(meas.Ah);
Ahdsz=meas.Ah(K1);
[AhMAX,K2]=max(meas.Ah);
% 
Discell=meas.Ah(1)-Ahdsz;
soc_ref=1-meas.Ah(1:K1)/Ahdsz;
ocv_dis=meas.Voltage(1:K1);
figure;
plot(soc_ref,ocv_dis,'b');
hold on;
x=0:0.001:1;
y=interp1(soc_ref,ocv_dis,x);
plot(x,y);

% SOC_1=y(1:501);
% SOC_2=y(502:1001);