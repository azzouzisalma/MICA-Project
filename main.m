clc;
close all;
clear all;
%% le spectrogramme de
load('/net/t/sazzouzi/Documents/Traitement_numerique_signal/MICA_project/data/PatientData.mat')
ecg=ecg{7,1};
Fs=200;
% ecg = (ecg(1:42200));
N_fft= 1024 ;
d=2000;
w= hamming(4*Fs);
% group delay introduce by the pom and tampkin algorithm and the causality
% of differentiated filter = 2 samples
delay1= 0+ 2 + 4 ;
% group delay introduce by g1 and g2
delay2= 6;
%% spectrogramms visualisation
[X, f, t ] = stft(ecg,w,d,N_fft,Fs);
[Sx, f, t] = spectro(ecg,w(1:512),d,N_fft,Fs);

imagesc(t,f,10*log(Sx)) % a logarithmique
xlabel('time in s');
ylabel('frequency in hz');
title('spectrogram of an ecg_VF  with windows of duration 4s');

figure();
plot(ecg);
xlim([0 10e2]);
xlabel('time in s');
ylabel('ecg');
title('ecg signal for a data number 2');
% figue();
% spectrogram(ecg , w,1000,N_fft);
%% combinaition of low and high filter
A_low = [1 -2 1]; % denominator 
B_low= [1 0 0 0 0 0 -2 0 0 0 0 0 1];

B_high =[-1 zeros(1,15) 32 -32 zeros(1,14) 1];
A_high= [1 -1];
ecg_filter_low = filter( B_low , A_low,ecg);
ecg_filter_high = filter(B_high , A_high ,ecg);
B_der =Fs/8*[1 2 0 -2 -1];
A_der=[1]; % filtre causal , on a eleminer le z-2 pour le rendre causal , atytention faut prendre le retard en considération
ecg_filter_derivation =filter( B_der , A_der , ecg);
S_seq= abs(ecg_filter_derivation()).^2;

S_MWI= conv(S_seq , 1/(10)*ones(1 ,10));
% subplot(6, 1, 1);
% plot(ecg_filter_low);
% title('signal low-filtred');
% xlim([0 4e3]);

% subplot(6,1,2);
% plot(ecg);
% hold on;
% plot(T_position , ecg(T_position),  '*');
% hold on;
% plot(P_pos , ecg(P_pos), '+');
% title('signal ecg');
% xlabel("temps en s");
% xlim([0 10e2]);

subplot(2,2,1);
plot(ecg_filter_high);
title('signal high and low-filtred');
xlabel("temps en s");
xlim([0 4e2]);

subplot(2,2,2)
plot(ecg_filter_derivation);
title('derivated signal ');
xlim([0 4e2]);
xlabel("temps en s");

subplot(2,2,3);
plot(S_seq);
title('squared signal');
xlim([0 2e2]);
xlabel("temps en s");
%% threasholding 
[peak_1, index]= findpeaks(S_MWI);

%tresholding =(mean(peak_1(1:40)));

tresholding = mean(S_MWI) ;
subplot(2,2,4);
plot(S_MWI);
xlim([0 2*Fs]);
hold on
plot(index , tresholding*ones(length(index)));
title("the output signal of moving window and the tresholding");
xlabel("time in s");
ylabel("amplitude ");


mask = S_MWI > (tresholding);
R_time= zeros(1 ,length(mask));
for i=1:length(mask)-1
    if mask(i)==1
        R_time(1,i)= i/Fs -delay1/Fs;
    end
end

peak= zeros(1, length(mask));
for i=1:length(mask)-1
    if mask(i)==1
        peak(1,i)=S_MWI(i);
    end
end

% hold on;
% plot(T_position , ecg(T_position),  '*');
% hold on;
% % plot(P_pos , ecg(P_pos), '+');
% title('signal ecg');
% xlabel("temps en s");
% xlim([0 1e2]);

[R_values , position1] = findpeaks(peak);
position= position1 - delay1;
R_position= zeros(1,length(position));
for i=1:length(position)
    R_position(1,i)=R_time(position1(i));  % position1 is the position without delay and it is suitable with the index of r_times déja vu
end
% figure()
% plot(ecg);
% hold on;
% plot(position , ecg(position),  '*');
% xlim([0 10e2]);
BPM=60./diff(R_position(1:16)); % we take hust the 16 first values of r position
BPM_moy=mean(BPM);
tachycardia=0;
bradycardia=0;
tach_pos=[];
brad_pos=[];

for i=1:length(BPM)
    if(BPM(i)<60)
        bradycardia=1;
        brad_pos=[brad_pos i];
    end
    if(BPM(i)>100)
        tachycardia=1;
        tach_pos=[tach_pos i];
    end
end



%% Determenation of S and Q
temps= (0:length(peak)-1)/Fs;
D_peak = diff(peak);

BPM= 60./diff(R_position);
 %Q and S wave location
 index=[];                          %index of beginning and last R interval
 for k=1:length(mask)-1
     if mask(k)==0 && mask(k+1)==1
         index=[index k+1];
     end
     if mask(k)~=0 && mask(k+1)==0
         index=[index k];
     end
 end
 
%Q and S wave location
 
Q_pos=[];
S_pos=[];
c=1;
for i=1:2:length(index)-1
    diff_left= diff(peak(index(i):position1(c)));
    diff_right= diff(peak(position1(c):index(i+1)));
    for j=1:length(diff_left)                           %On calcule l'indice du minimum dans la partie gauche et droite de chaque intervalle R
        l_max=1;
        if diff_left(j) > diff_left(l_max)
            l_max=j;
        end
      end
 
    for k=1:length(diff_right)
        r_min=1;
        if diff_right(k) < diff_right(r_min)
            r_min=k;
        end
    end
     l_max=position(c)+ l_max - length(diff_left);
     r_min=r_min + position(c);
     c=c+1;
     Q_pos=[Q_pos l_max];
     S_pos=[S_pos r_min] ;
end
% Q_pos_delay = Q_pos - delay1;
% S_pos_delay = S_pos -delay1;
% figure()
% plot(ecg);
% hold on;
% plot(position , ecg(position),  '*');
% hold on;
% plot(Q_pos(1:50) , ecg(Q_pos(1:50)), '*');
% hold on;
% plot(S_pos , ecg(S_pos), '+');
% hold on;
% plot(T_position, ecg(T_position));
% legend('ecg','R','Q','S');
% xlim([0 10e2]);

%% Detection of T and P position
B_g1= [1 0 0 0 0 0 -1];
A_g1=[1];
B_g2= [1 0 0 0 0 0 0 0 -1];
A_g2=[1 -1];
ecg_g1 = filter(B_g1, A_g1 , abs(ecg));
ecg_g2 = filter(B_g2, A_g2 , ecg_g1);
temps=(0:1000)/Fs ;
figure()
plot(ecg_g2);
title('ecg to detect T');
xlim([0 1000]);
xlabel('temps en s ');
ylabel("ecg");

T_position=[];
extremum=[];
taille=[1];
extremum_value=[];
for i=1:length(position)-1
    d=0;
    RR= ecg_g2(position(i):position(i+1));
    for j=position(i):round(0.7*(position(i+1)-position(i)))+position(i)
        if ecg_g2(j) >0 && ecg_g2(j+1)<0 % tu prend j 
           extremum=[extremum j]; % Tous les maximum locaux
           extremum_value= [extremum_value ecg(j)];
           d=d+1;
        end
    end
    
    taille=[taille length(extremum)];
%     vmax= extremum(1);
%      for k=1:length(extremum)
%          %vmax= extremum(1);
%          if abs(ecg(vmax)) < abs(ecg(extremum(k)))
%              vmax= extremum(k);     
%         end        
%    end
%     
%     T_position=[T_position (vmax)];
%     extremum=[];
        
end
extremum_ecg=ecg(extremum);
for i=1:length(taille)-1
    extremum_local= extremum(taille(i):taille(i+1));
    max_ecg=extremum_local(1);
    for j=1:length(extremum_local)
        if abs(ecg(max_ecg))<abs(ecg(extremum_local(j)))
            max_ecg=extremum_local(j) ;
        end
    end
    T_position=[T_position max_ecg];
end

% P detection

P_pos= [];
for k=1:length(position)-5
    for i=round(0.7*(position(k+1)-position(k)))+position(k):position(k+1)-7
        imax= round(0.7*(position(k+1)-position(k)))+position(k);
        if (ecg(imax)) < (ecg(i))
            imax=i;
        end
    end
    
    P_pos= [P_pos imax];
end

 figure()      
plot(ecg);
% hold on;
% plot(S_pos_delay , ecg(S_pos_delay), '+');
% hold on;
% plot(position, ecg(position), '*' );
% hold on ;
% plot(Q_pos_delay , ecg(Q_pos_delay), '');
% title('signal ecg');
% xlabel("temps en s");
% xlim([0 10e2]);
% legend('ecg','S','R');
figure()
plot(ecg);
hold on;
plot(position , ecg(position),  '*');
hold on;
plot(Q_pos(1:10) , ecg(Q_pos(1:10)), '*');
hold on;
plot(S_pos , ecg(S_pos), '+');
% hold on;
% plot(T_position-delay2, ecg(T_position-delay2),'+');
% hold on;
% plot(P_pos, ecg(P_pos),'+');
legend('ecg','R','Q','S');
title('the QRS positions for a normal person (Data normal1)')
xlim([0 10e2]);
xlabel('time in s');
ylabel('ECG');


figure()
plot(ecg);
hold on;
plot(position , ecg(position),  '*');
hold on;
plot(Q_pos(1:10) , ecg(Q_pos(1:10)), '*');
hold on;
plot(S_pos , ecg(S_pos), '+');
hold on;
plot(T_position-delay2, ecg(T_position-delay2),'+');
hold on;
plot(P_pos, ecg(P_pos),'+');
legend('ecg','R','Q','S','T','P');
title('the PQRST positions for a vf person ')
xlim([0 10e2]);
xlabel('time in s');
ylabel('ECG');
    
%% pathologies
%Tachycardia/Bradycardia
BPM=60./diff(R_position(1:16)); % we take hust the 16 first values of r position
BPM_moy=mean(BPM);
tachycardia=0;
bradycardia=0;
tach_pos=[];
brad_pos=[];

for i=1:length(BPM)
    if(BPM(i)<60)
        bradycardia=1;
        brad_pos=[brad_pos position(i)];
    end
    if(BPM(i)>100)
        tachycardia=1;
        tach_pos=[tach_pos position(i)];
    end
end
figure();
hold on;
plot(ecg);
hold on;
plot(brad_pos,ecg(brad_pos),'*');
hold on;
plot(tach_pos,ecg(tach_pos),'+');
legend('ecg','bradycardia','tachycardia');
title('position of bradycardia and tachycardia on an ecg')
xlim([0 10e3]);
xlabel('time in s');
ylabel('ECG');

%%Ectopic beat
thresholding=0.2425;
ectopic=0;
delta=diff(position);
ectopic_pos=[];                 %give the number of RR interval where ectopic pathologie is present
for i=1:length(delta)-1
    if abs(delta(i+1)-delta(i))>= thresholding
        ectopic=1;
        ectopic_pos=[ectopic_pos position(i)];
    end
end

figure();
hold on;
plot(ecg);
hold on;
plot(ectopic_pos,ecg(ectopic_pos),'*');
legend('ecg','ectopic position');
title('position of ectopic beats on the patient data')
xlim([0 10e3]);
xlabel('time in s');
ylabel('ECG');

figure();
plot(ecg);
xlim([0 10e4]);
xlabel('time in s');
ylabel('ecg');
title('rcg of the patient data 3')
    