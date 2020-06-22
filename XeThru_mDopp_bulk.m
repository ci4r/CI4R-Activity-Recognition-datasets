function [] = XeThru_mDopp_bulk(filename, fOut)

% Open file.
fid = fopen(filename, 'rb');

ctr = 0;

while(1)
       

    ContentId=fread(fid,1,'uint32');
    Info=fread(fid,1,'uint32');
 
    ctr=ctr+1; 
    Length=fread(fid,1,'uint32');
    Data=fread(fid,182,'float');
    
    if feof(fid)
           break
    end
        
    Datastream(:,ctr)=Data;
    data_length(ctr)=Length;

  
end 
    
frame_start=0.3; %0.3858
frame_stop=20;

%Generate range vector
bin_length = 8 * 1.5e8/23.328e9; % range_decimation_factor * (c/2) / fs.
range_vector = (frame_start-1e-5):bin_length:(frame_stop+1e-5); % +-1e-5 to account for float precision.

for n=1:size(Datastream,2)
        Data=Datastream(:,n);
        iq_mat(:,n) = Data(1:end/2) + 1i*Data(end/2 + 1:end);
end

%% MTI Filter
[b,a]=butter(4, 0.01, 'high'); %  4th order is 24dB/octave slope, 6dB/octave per order of n
                                 % [B,A] = butter(N,Wn, 'high') where N filter order, b (numerator), a (denominator), ...
                                 % highpass, Wn is cutoff freq (half the sample rate)
[h,fl]=freqz(b,a,size(iq_mat,2));
for k=1:size(iq_mat,1)
    Data_RTI_complex_MTIFilt(k,:)=filter(b,a,iq_mat(k,:));
end

PRF=512;

TimeAxis=(1:size(iq_mat,2))/PRF;
DopplerAxis=linspace(-PRF/2,PRF/2,size(iq_mat,2));

nfft = 2^12;window = 150;noverlap = 128;shift = window - noverlap; 
%sx = myspecgramnew(sum(iq_mat(5:50,:)),window,nfft,shift);                  % why row 8?
% sx = myspecgramnew(sum(Data_RTI_complex_MTIFilt(5:90,:)),window,nfft,shift);                  % why row 8?
sx = myspecgramnew(Data_RTI_complex_MTIFilt(27,:),window,nfft,shift); %18 front, 27 corner
sx1 = fftshift(sx,1);
sx_scaled = sx1(1408:3688,:);
sx2 = abs(sx1);
%% Denoising

% Part 1: Isodata thresholding
sx2 = rescale(sx2,0,255);
ctr = 0;
prev_t = 0;
t = mean(sx2(:))+1; % initial threshold
epst = 10;
while (1)
    low_idx = find(sx2(:)<t);
    high_idx = find(sx2(:)>t);
    mH = mean(sx2(high_idx));
    mL = mean(sx2(low_idx));
    prev_t = t;
    t = (mH+mL)/2;
%     t = (mH+prev_t)/2; 
    if abs(t-prev_t) < epst % mH-t
       break
    end
    ctr = ctr+1;
end
sx2(low_idx) = 0;

num_parts = 100; % 20
mean_energy = Energy(sx2)/num_parts;
parts_energy = zeros(1,num_parts);
stepsize = floor(size(sx2,1)/num_parts);
for i=1:num_parts
    part = sx2((i-1)*stepsize+1:i*stepsize,:);
    parts_energy(i) = Energy(part);
    if parts_energy(i) < mean_energy/200
        sx2((i-1)*stepsize+1:i*stepsize,:) = 0;
    end
end 

%% micro-Doppler Spectrogram
fig = figure('units','normalized','outerposition',[0 0 .8 .8]);
colormap(jet(256));  
imagesc(TimeAxis,DopplerAxis,20*log10(sx2./max(max(sx2))));
% doppSignMTI = imagesc(TimeAxis,DopplerAxis,20*log10(abs(sx1)));
axis xy
% set(gca,'FontSize',10)
% title({'Micro Doppler Signature After MTI Filter'; fnameBin});
% xlabel('Time (sec)');
% ylabel('Frequency (Hz)','FontSize',20);
% ax = gca;
% ax.FontSize = 16; 
caxis([-45 0])
axis([0 max(TimeAxis) -PRF/2 PRF/2])
%% Save
save(strcat(fOut(1:end-4),'.mat'),'sx_scaled');
saveas(fig,strcat(fOut(1:end-4),'.fig'));
% F = getframe(gca);
% [img, ~] = frame2im(F);
% imwrite(img,fOut)
% %% Gray
% colormap(gray)
% F = getframe(gca);
% [img, ~] = frame2im(F);
% imwrite(img,fOut_gray)
close all
fclose('all');
end
