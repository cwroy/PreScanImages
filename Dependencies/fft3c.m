function res = fft3c(x)
fctr = size(x,1)*size(x,2)*size(x,3);
res=x;%zeros(size(x));
for i=1:size(x,4)
res(:,:,:,i)= 1/sqrt(fctr)*fftshift(fftn(ifftshift(x(:,:,:,i))));
end
end




