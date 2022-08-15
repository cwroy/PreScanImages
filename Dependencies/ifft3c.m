function res = ifft3c(x)
fctr = size(x,1)*size(x,2)*size(x,3);
res=zeros(size(x));
for i=1:size(x,4)
    for j=1:size(x,5)
        for k=1:size(x,6)
            res(:,:,:,i,j,k) = sqrt(fctr)*fftshift(ifftn(ifftshift(x(:,:,:,i,j,k))));
        end
    end
end