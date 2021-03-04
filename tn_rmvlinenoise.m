function [data] = tn_rmvlinenoise(data,nf,nharm,srate)

% removed power line artifacts with discrete fourier filter
% data: should be row vector
% nf:  frequency of artifact
% nharm: number of harmonics you want to remove
% srate: sampling frequency

% remove mean
datamean  = nanmean(data);
data      = data-datamean;
nsamples  = length(data); % get num samples

for nfh = nf:nf:(nf*nharm)
    n         = round(floor(nsamples .* nf./srate) * srate./nf);
    sel       = 1:n(1);
    time      = (0:nsamples-1)/srate;
    tmp       = exp(1i*2*pi*nfh*time);            
    ampl      = 2*data(:,sel)/tmp(:,sel);                
    est       = ampl*tmp;                              
    filt      = data - est;                              
    data      = real(filt);
end

data = data + datamean; % add mean back
end