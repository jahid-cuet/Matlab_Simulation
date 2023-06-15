function [Fd,fPPG_IR,fPPG_Red,sPPG_IR,sPPG_Red] = PPG_filter(Fs,N,PPG_IR,PPG_Red)
% linear detrending
nPoly = 1;
dPPG_IR = detrend(PPG_IR,nPoly);
dPPG_Red = detrend(PPG_Red,nPoly);

%4th order bandpass Butterworth filter (0.5 - 5 Hz) 
Fc_low = 0.5;                           %Low cutoff freq, Hz
Fc_high = 5;                            %High cutoff freq, Hz
[b,a] = butter(4,[Fc_low Fc_high]/(Fs/2));
fPPG_IR = filtfilt(b,a,dPPG_IR);        %Filtered IR PPG, no phase shift
fPPG_IR = fPPG_IR-mean(fPPG_IR);        %Shift mean point to 0
fPPG_Red = filtfilt(b,a,dPPG_Red);      %Filtered Red PPG, no phase shift
fPPG_Red = fPPG_Red-mean(fPPG_Red);     %Shift mean point to 0

%Make copies for use in SpO2 calculations, and reference graph
sPPG_IR = fPPG_IR;
sPPG_Red = fPPG_Red;

%Normalize the PPG amplitudes for graphing nicely with ECG
fPPG_IR = (fPPG_IR - min(fPPG_IR))/(max(fPPG_IR)-min(fPPG_IR))*2 - 1;
fPPG_IR = fPPG_IR - mean(fPPG_IR);

%% Apply FFT to the Filtered IR PPG
y = fft(fPPG_IR);                          
P2 = abs(y/N);
P1 = P2(1:round(N/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs/N*(0:(N/2));
[~,I] = max(P1);
Fd = f(I);                              %Dominant frequency
end

