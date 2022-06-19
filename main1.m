clear all
clc
close all


%%
% figure
% spectrogram(ecg,hamming(1000),0,2048,Fs,'yaxis')
% figure
% plot(ecg)

%%
 load('/net/t/sazzouzi/Documents/Traitement_numerique_signal/MICA_project/data/ecg_normal_1.mat')
Ts=1/Fs;                                      
Nfft=1024;
[Sx,f,t]=spectro(ecg,hamming(8*Fs),5000,Nfft,Fs);    %window of 8 s 
figure
imagesc(10*log10(Sx))
colorbar
title('spectrogram of normal person')
figure
spectrogram(ecg,hamming(Nfft),511,Nfft,Fs,'yaxis');

%%R wave detection
    
    %Band-pass filter
ecg_filtered_low=filter([1 0 0 0 0 0 -2 0 0 0 0 0 1],[1 -2 1],ecg);
delay1=5; %samples


num_coeff=zeros(1,33); %numerator coefficient of the high-pass filter
num_coeff(1)=-1;
num_coeff(17)=32;
num_coeff(18)=-32;
num_coeff(33)=1;

ecg_filtered_high= filter(num_coeff,[1 -1],ecg_filtered_low);
delay2=16.5;

    %Derivative to provide the QRS complex slope information
 ecg_derivated=filter([1 2 0 -2 -1],[8*Ts],ecg_filtered_high);      %Filter became causal by mutliplying by z^-2(retarde le signal de -2Ts)
 delay3=2;
 
    %square
 s_sq=abs(ecg_derivated.^2);
    
    %Moving Window Integration
w_length=0.1*Fs; 
s_MWI=filter(ones(1,w_length),[1/w_length],s_sq);
delay4=9.5;

delay_total=delay1+delay2+delay3+delay4+2;      %le 2 est pour le filtre d�riv� que l'on rend causal

    %Thresholding
thresholding=mean(s_MWI);
mask=s_MWI>thresholding;                        % Tab that have 1 when the signal is higher than thresholding


% R position
R_interval= zeros(1 ,length(mask));             %Tableau contenant les intervalles R
for i=1:length(mask)
    if mask(i)==1
        R_interval(i)= i*Ts;
    end
end
 
peak= zeros(1, length(mask));               %Tableau contenant les valeurs dans les intervalles R
for i=1:length(mask)
    if mask(i)==1
        peak(i)=s_MWI(i);
    end
end

 index=[];                          %index of beginning and last R interval
 for k=1:length(peak)-1
     if peak(k)==0 && peak(k+1)~=0
         index=[index k+1];
     end
     if peak(k)~=0 && peak(k+1)==0
         index=[index k];
     end
 end

R_pos=[];                                   %R_pos renvoie la position du R 
for k=1:2:length(index)-2
    interval=peak(index(k):index(k+1));
    i_max=1;
    for j=1:length(interval)
        if interval(j)>interval(i_max)
            i_max=j;
        end
    end
    i_max=i_max+index(k);
    
    R_pos=[R_pos i_max];
end
 

%%Plot
 t=(0:length(ecg)-1)*Ts;
 new_t=t+2*Ts;                  %Adjust time after delay due to derivative filter

 figure
 plot(t,ecg_filtered_low)
 title('ecg after low filter')
 xlim([0 5])
 xlabel('time (s)')
 
 figure
 subplot(5,1,1)
 plot(t,ecg)
 title('ecg')
 xlim([0 20])
 

 subplot(5,1,2)
 plot(t,ecg_filtered_high)
 title('ecg after band filter')
 xlim([0 20])
 xlabel('time (s)')
 
 subplot(5,1,3)
 plot(new_t,ecg_derivated)
 title('ecg derivated')
 xlim([0 20])
 xlabel('time (s)')

 
 subplot(5,1,4)
 plot(new_t,s_sq)
 title('ecg squared')
 xlim([0 20])
 xlabel('time (s)')
 
 subplot(5,1,5)
 hold on
 plot(s_MWI)
 plot(thresholding*ones(size(s_MWI)))
 title('ecg after moving window integer')
 xlim([0 1000]) 
 xlabel('sample')
 
 %Q and S wave location

Q_pos=[];
S_pos=[];
c=1;
for i=1:2:length(index)-2
    diff_left= diff(peak(index(i):R_pos(c)));
    diff_right= diff(peak(R_pos(c):index(i+1)));
    for j=1:length(diff_left)                           %On calcule l'indice du maximum dans la partie gauche et l'indice du minimum droite de chaque intervalle R
        l_max=1;
        if diff_left(j)>diff_left(l_max) 
            l_max=j;
        end
    end
    
    for k=1:length(diff_right)
        r_min=1;
        if diff_right(k)<diff_right(r_min)
            r_min=k;
        end
    end
     l_max= l_max + index(i);
     r_min=r_min + R_pos(c);
     c=c+1;
     Q_pos=[Q_pos l_max];
     S_pos=[S_pos r_min] ;
end

R_pos=R_pos -delay_total;%Prend en compte le delay total introduit par chaque filtre
R_pos(1)=[];               %Quand le retard est plus grand que les premieres positions de R on ne les prend pas en compte

index=index-delay_total;
R_interval= R_interval-delay_total;
Q_pos= Q_pos -delay_total;
S_pos=S_pos -delay_total;

%%P and T wave detection
ecg_diff=filter([1 0 0 0 0 0 -1],[1],ecg);
ecg_diff_low=filter([1 0 0 0 0 0 0 0 -1],[1 -1],ecg_diff);
delay_tot= 6;

    %T detection
 T_pos=[];
 i_maximum=[];                                              %T is the first passage in zero in RR interval in derivated ecg
 for k=1:length(R_pos)-1
     for i=R_pos(k):round(0.7*(R_pos(k+1)-R_pos(k))+R_pos(k))       %search all maximum in a R-70%R interval
         if ecg_diff_low(i)>0 && ecg_diff_low(i+1)<0
             i_maximum=[i_maximum i];
         end
     end
     for j=1:length(i_maximum)
         j_max=i_maximum(1);
         if ecg(i_maximum(j))>ecg(j_max)
             j_max=i_maximum(j);
         end
     end
     i_maximum=[];
     T_pos=[T_pos j_max];
 end
 
    %P detection
 P_pos=[];                                                                      %search the highest peak of ecg in the interval 70%R-R
 for k=1:length(R_pos)-1
     for i=round(0.7*(R_pos(k+1)-R_pos(k))+R_pos(k)):R_pos(k+1)
         i_max=round(0.7*(R_pos(k+1)-R_pos(k))+R_pos(k));
         if ecg(i)>ecg(i_max)
             i_max=i;
         end
     end
     P_pos=[P_pos i_max];
 end
 
 T_pos=T_pos -delay_tot;
 P_pos=P_pos -delay_tot;

figure
hold on
plot(ecg)
plot(P_pos,ecg(P_pos),'x')
plot(Q_pos,ecg(Q_pos),'*')
plot(R_pos,ecg(R_pos),'*')
plot(S_pos,ecg(S_pos),'*')
plot(T_pos,ecg(T_pos),'*')
legend('ecg','P','Q','R','S','T')
xlabel('sample'),title('ecg')
xlim([0 500])
hold off

figure
subplot(3,1,1)
plot(ecg)
xlabel('sample'),title('ecg')
xlim([0 3000])

subplot(3,1,2)
plot(ecg_diff)
xlabel('sample'),title('ecg differentiated')
xlim([0 1000])

subplot(3,1,3)
hold on
plot(ecg_diff_low)
plot(zeros(size(ecg_diff_low)))
xlabel('sample'),title('ecg low-pass filtered')
xlim([0 1000])
 


%Tachycardia/Bradycardia
BPM=60./diff(R_pos*Ts);
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


%Ectopic beat
ectopic=0;
delta=diff(R_pos);
ectopic_pos=[];                 %give the number of RR interval where ectopic pathologie is present
for i=1:length(delta)
    if abs(delta(i))>= thresholding
        ectopic=1;
        ectopic_pos=[ectopic_pos i];
    end
end

%Fibrillation
    %Atrial fibrillation
% Nfft=512;
% TF_delta=fft(delta,Nfft);
% k=0:Nfft-1;
% f=(Fs/Nfft)*k;
% figure
% plot(f,fftshift(abs(TF_delta)))

% auto_cov=[];
% for k=1:length(R_pos)-1
%     auto_cov=[auto_cov gamma_(k,R_pos)];
% end




        