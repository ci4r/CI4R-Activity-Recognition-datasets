function [ ] = datToImage_Anchortech( fNameIn, fNameOut, fnameBin )

fileID = fopen(fNameIn, 'r');
dataArray = textscan(fileID, '%f');
fclose(fileID);
radarData = dataArray{1};
clearvars fileID dataArray ans;
fc = radarData(1); % Center frequency
Tsweep = radarData(2); % Sweep time in ms
Tsweep=Tsweep/1000; %then in sec
NTS = radarData(3); % Number of time samples per sweep
Bw = radarData(4); % FMCW Bandwidth. For FSK, it is frequency step;
% For CW, it is 0.
Data = radarData(5:end); % raw data in I+j*Q format

fs=NTS/Tsweep; % sampling frequency ADC
record_length=length(Data)/NTS*Tsweep; % length of recording in s
nc=record_length/Tsweep; % number of chirps

%% Reshape data into chirps and do range FFT (1st FFT)
Data_time=reshape(Data, [NTS nc]);

%Part taken from Ancortek code for FFT and IIR filtering
tmp = fftshift(fft(Data_time),1);
Data_range(1:NTS/2,:) = tmp(NTS/2+1:NTS,:);
% IIR Notch filter
ns = oddnumber(size(Data_range,2))-1;
Data_range_MTI = zeros(size(Data_range,1),ns);
[b,a] = butter(4, 0.01, 'high');
% [h, f1] = freqz(b, a, ns);
for k=1:size(Data_range,1)
  Data_range_MTI(k,:) = filter(b,a,Data_range(k,:));
end

%% Spectrogram processing for 2nd FFT to get Doppler
% This selects the range bins where we want to calculate the spectrogram
bin_indl = 3;
bin_indu = 60;
%Parameters for spectrograms
MD.PRF=1/Tsweep;
MD.TimeWindowLength = 200;
MD.OverlapFactor = 0.95;
MD.OverlapLength = round(MD.TimeWindowLength*MD.OverlapFactor);
MD.Pad_Factor = 4;
MD.FFTPoints = MD.Pad_Factor*MD.TimeWindowLength;
MD.DopplerBin=MD.PRF/(MD.FFTPoints);
MD.DopplerAxis=-MD.PRF/2:MD.DopplerBin:MD.PRF/2-MD.DopplerBin;
MD.WholeDuration=size(Data_range_MTI,2)/MD.PRF;
MD.NumSegments=floor((size(Data_range_MTI,2)-MD.TimeWindowLength)/floor(MD.TimeWindowLength*(1-MD.OverlapFactor)));

nfft = 2^12;window = 128;noverlap = 100;shift = window - noverlap;
sx = myspecgramnew(Data_range_MTI(19,:),window,nfft,shift); % best bin 19-20 for asl
sx1 = flipud(fftshift(sx,1));
sx_scaled = sx1(1408:3688,:);
sx2 = abs(sx1);
MD.TimeAxis=linspace(0,MD.WholeDuration,size(Data_range_MTI,2));
%% Denoising

% Part 1: Isodata thresholding
sx2 = rescale(sx2,0,255);
ctr = 0;
t = 4;%mean(sx2(:)); % initial threshold
epst = 100;
while (1)
    low_idx = find(sx2(:)<t);
    high_idx = find(sx2(:)>t);
    mH = mean(sx2(high_idx));
    mL = mean(sx2(low_idx));
    prev_t = t;
%     t = (mH+mL)/2;
    t = (mH+prev_t)/2; 
    if abs(t-prev_t) < epst % mH-t
       break
    end
    ctr = ctr+1;
end
sx2(low_idx) = 0;
med = median(median(nonzeros(sx2)));
sx2(sx2 <= med/2) = 0;
% clear t prev_t


% num_parts = 20; % 20
% mean_energy = Energy(sx2)/num_parts;
% parts_energy = zeros(1,num_parts);
% stepsize = floor(size(sx2,1)/num_parts);
% for i=1:num_parts
%     part = sx2((i-1)*stepsize+1:i*stepsize,:);
%     parts_energy(i) = Energy(part);
%     if parts_energy(i) < mean_energy/10
%         sx2((i-1)*stepsize+1:i*stepsize,:) = 0;
%     end
% end 

%%

fig= figure('units','normalized','outerposition',[0 0 .5 .5]);

% imagesc(MD.TimeAxis,MD.DopplerAxis.*3e8/2/5.8e9,20*log10(abs(Data_spec_MTI2))); colormap('jet'); axis xy
colormap(jet(256)); %xlim([1 9])
imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(sx2./max(max(sx2))))); % normalization = /max(max(abs(Data_spec_MTI2)))
axis xy
ylim([-6/(3e8/2/5.8e9) 6/(3e8/2/5.8e9)]);
%ylim([-6 6]);
ylabel('Frequency(Hz)', 'FontSize',9);
clim = get(gca,'CLim');
set(gca, 'CLim', clim(2)+[-45,0]);
% set(gca,'xtick',[],'FontSize',9)
xlabel('Time(s)', 'FontSize',9);
% ylabel('Velocity [m/s]','FontSize',9)
save(strcat(fNameOut(1:end-4),'.mat'),'sx_scaled');
saveas(fig,strcat(fNameOut(1:end-4),'.fig'));
fclose('all');
close all
end