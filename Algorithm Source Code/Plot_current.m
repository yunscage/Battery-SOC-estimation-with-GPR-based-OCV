cnt_degc=4;
cnt_j=3;
if cnt_j==2
    disp('Stop!');
end
Tem=[-20 -10 0 10 25];
disp(['Degc: ' num2str(Tem(cnt_degc))]);
load('..\DataProcess\PathName.mat');
filepath=[Path_1{cnt_degc} Path_2{cnt_degc,cnt_j}];
load(filepath);
Cycle=Path_2{cnt_degc,cnt_j}(20:24);
if Cycle(end)=='_'
    Cycle(end)=' ';
end

%% Pre-Process
I=-meas.Current';
plot(I);