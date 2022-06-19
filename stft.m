function [X, f, t ] = stft(x,w,d,N_fft,Fs)% This function computes the stft for m = [0, d, 2d, 3d...]
% This function outputs are:% -> X, which is  a matrix of n_fft lines and M columns
L=length(x);
M=round(L/d) ;
X= zeros(length(w), M);

%    M is the number of elements of m
%    X(i,j) is the value of the spectrogram for time t(i) and frequency f(j)
% -> f, is a column vector of the frequencies (in Hz)
% -> t, is a row vector containing the times of the beginning of the windows
% 
i=1;
for j=1:M
    X(:,j)=x(i:i+length(w)-1);
    i=i+d;
end
% multiplication par la fenetre avec w une colonne
for j=1:M
    X(:, j)= X(:,j).*w;
end

X= fft(X , N_fft);

f=0:Fs/N_fft:Fs*(N_fft-1)/N_fft;
 t= N_fft/(Fs*2):1/Fs:(M+N_fft/2-1)/Fs;
% nbr_column=100;
% t=0:w(0):nbr_column*w(0)-1;

end


