clear;

%% Load data
load( 'broadband_hologram_2b' ); % hologram data from paper
[Ny, Nx] = size(mat_holo);

% define normalization area
vec_norm_x = 4:11;    vec_norm_y = 87:94;


%% Parameters for reconstruction
Nd = 17; % number of interferogram points
Nterm = 8; % number of terms
filter_width_x = 1/1; filter_width_y = 1/(Nterm*4);

% resample final image in y to make image ~square, to compensate for 
% 6-times oversampling in y
resample = 1/6; 

% total interferogram length = step width * number of interferogram points
len_inter = 5.2779e-06 * Nd; 
% FFT bin (IR) frequency of the subsampled interferogram (assuming one 
% wrapping, adjust as necessry)
dec_freq = 0.01/(2*len_inter) *( (1:Nd)-1 + Nd )  ; % 
% IR frequency of the individual terms    
holo_freq = dec_freq(2:Nterm+1);



%% Preprocess data
% crop hologram to integer number of interferograms
Ny = floor(Ny/Nd)*Nd;     
mat_holo = mat_holo(1:Ny,:);
ky = Ny / (Nterm*2+1); kx = 0;
Nyresize = floor(Ny*resample);

%% Filtering in Fourier Space
mat_holo_fft = fftshift(fft2(mat_holo));
figure(1); imagesc(abs(mat_holo));
figure(1); imagesc(log10(abs(mat_holo_fft)));

% reconstruct and resize
[Ny, Nx] = size(mat_holo);
[X,Y] = meshgrid(1:Nx, 1:Ny);
sum_window = zeros(Ny,Nx);
mat_reco = zeros(Nyresize,Nx,Nterm);
for ct=1:Nterm
    windowX = Nx/2+ct*kx+1;
    windowY = Ny/2+(ct*ky+1);
    windowW = Nx*filter_width_x/1.0;
    windowH = Ny*filter_width_y/1.0;
    filter1 = 0.4 + 0.6 * cos( (Y-windowY) / windowH * pi*0.6 );
    filter2 = ( 0.4 + 0.6 * cos( (X-windowX) / windowW * pi*0.6 ) );
    mask1 = abs(Y-windowY) < windowH;
    mask2 = abs(X-windowX) < windowW;
    mat_window = filter1.*filter2.*mask1.*mask2;
    sum_window = sum_window + mat_window;
    mat_fft_filt = mat_holo_fft .* mat_window;
    mat_phasecorr= exp( -2*pi*1i*Y/Ny*(ct*ky) ) .* exp( -2*pi*1i*X/Nx*(ct*kx) );
    tmp = ifft2(fftshift(mat_fft_filt)) .* mat_phasecorr;;
    mat_reco(:,:,ct) = imresize( tmp,[Nyresize Nx]);
end

%normalize (e.g. on substrate as specified in vec_norm_y/x)
vec_ref = mean(mean(mat_reco(vec_norm_y,vec_norm_x,:),2),1);
holo_norm = bsxfun( @rdivide, mat_reco, vec_ref) ;


%% Plot data
%figure(2); imagesc(sum_window);

figure(3); clf;
for q=1:Nterm
subplot(2,Nterm,q); imagesc(abs(holo_norm(:,:,q))); caxis([0 1.2]); 
title(holo_freq(q));
subplot(2,Nterm,q+Nterm); imagesc(angle(holo_norm(:,:,q))); caxis([-0.5 1]); %caxis([-1 1]*1);
colormap 'gray';
end
