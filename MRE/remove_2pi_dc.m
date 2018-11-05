function out = remove_2pi_dc(input,dim)
%out = remove_2pi_dc(input,dim)
%remove 2*pi dc component from the dim dimention
[data,perm,nshifts] = shiftdata(input,dim);
my_ph_uw_fft = ifft(data,[],1);

sc = sign(my_ph_uw_fft(1,:)) .* round(abs(my_ph_uw_fft(1,:))/2/pi);
my_ph_uw_fft(1,:) =  my_ph_uw_fft(1,:) - sc*2*pi;
%   my_ph_uw_fft(:,:,1) = 0;
out = fft(my_ph_uw_fft,[],1);
out = unshiftdata(out,perm,nshifts);

