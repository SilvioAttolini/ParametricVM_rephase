function [rephased_diffuse] = getRephased(params,meanDiffuse_MIC, mic_pos, vm_pos)

%fLen=params.Fs;
fLen = length(meanDiffuse_MIC);
freq =linspace(0,params.Fs/2,fLen/2);
freqs =  [freq, flip(freq)];
dists=norm(mic_pos-vm_pos);

omega=2*pi*freqs;
%display("here0");

%FFT segnale input
diffuseFFT(:)=fft(meanDiffuse_MIC,fLen);
%display(size(diffuseFFT));

rephase(:) = exp(1i*dists*omega/params.c);
%display(size(rephase));

rephased_diffuse=ifft(diffuseFFT.*rephase,fLen);
%display(size(rephased_diffuse));
%display("here");

end