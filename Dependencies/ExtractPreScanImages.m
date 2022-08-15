function [IBODY,ISURFACE] = ExtractPreScanImages(Filepath,ReconResolution)
addpath(genpath('Dependencies'))

% Function to calculate individual (body and surface) coil images from prescan data. To use
% this function your protocol should have prescan normalize activated and
% you should save the multi-raid file. Only tested for 3D radial so far.

% Inputs:
%
% Filepath: full path to rawdata file
%i.e. Filepath='C:/meas_XXXX.dat';
%
% ReconResolution: The base resolution you want to interpolate the prescan
% images to

% Output:
% Centered and interpolated body (IBODY) and surface (ISURFACE) coil images

IBODY=[];
ISURFACE=[];
% Modified mapVBVD to include prescan body and surface coil images
twix_obj = mapVBVD_PreScan(Filepath);
while(1)
    % If multi-raid your twix_obj should be a cell
    if iscell(twix_obj)
        
        if ~isfield(twix_obj{1,1},'prescanbody')
            break
        end
        
        % Unsorted k-space from body coil elements
        tempBody=twix_obj{1,1}.prescanbody.unsorted();
        % Unsorted k-space from surface coil elements
        tempSurface=twix_obj{1,1}.prescansurface.unsorted();
        
        % sort k-space for body coil images
        KBODY=zeros(max(twix_obj{1,1}.prescanbody.Par),max(twix_obj{1,1}.prescanbody.Lin),twix_obj{1,1}.prescanbody.NCol,2,'single');
        for iLin=1:twix_obj{1,1}.prescanbody.NLin
            KBODY(twix_obj{1,1}.prescanbody.Par(iLin),twix_obj{1,1}.prescanbody.Lin(iLin),:,:)=tempBody(:,:,iLin);
        end
        
        % sort k-space for surface coil images
        KSURFACE=zeros(max(twix_obj{1,1}.prescansurface.Par),max(twix_obj{1,1}.prescansurface.Lin),twix_obj{1,1}.prescansurface.NCol,twix_obj{1,1}.prescansurface.NCha,'single');
        for iLin=1:twix_obj{1,1}.prescansurface.NLin
            KSURFACE(twix_obj{1,1}.prescansurface.Par(iLin),twix_obj{1,1}.prescansurface.Lin(iLin),:,:)=tempSurface(:,:,iLin);
        end
        
    else
        if ~isfield(twix_obj,'prescanbody')
            break
        end
                % Unsorted k-space from body coil elements
        tempBody=twix_obj.prescanbody.unsorted();
                % Unsorted k-space from surface coil elements
        tempSurface=twix_obj.prescansurface.unsorted();
        % sort k-space for body coil images
        KBODY=zeros(max(twix_obj.prescanbody.Par),max(twix_obj.prescanbody.Lin),twix_obj.prescanbody.NCol,2,'single');
        for iLin=1:twix_obj.prescanbody.NLin
            KBODY(twix_obj.prescanbody.Par(iLin),twix_obj.prescanbody.Lin(iLin),:,:)=tempBody(:,:,iLin);
        end
        
        % sort k-space for surface coil images
        KSURFACE=zeros(max(twix_obj.prescansurface.Par),max(twix_obj.prescansurface.Lin),twix_obj.prescansurface.NCol,twix_obj.prescansurface.NCha,'single');
        for iLin=1:twix_obj.prescansurface.NLin
            KSURFACE(twix_obj.prescansurface.Par(iLin),twix_obj.prescansurface.Lin(iLin),:,:)=tempSurface(:,:,iLin);
        end
        
        
    end
    
    % Very basic realignment of the prescan acquisition to the scan
    % acquisition. This will work well for tranverse 3D data with offsets
    % in x,y,z. Arbitrary orientations will require more complex
    % transformations
   
    
    % Calculate the offset from the center for the image acquisition
    M=(twix_obj{1,length(twix_obj)}.image.slicePos(1:3,1)-twix_obj{1,1}.image.slicePos(1:3,1));%/PixelSize;
    Mx=-M(2);My=M(1);
    Mz=-M(3);

    % Calculate the k-space coordinates of the prescan images
    [kx,ky,kz]=meshgrid(linspace(-1,1,size(KBODY,1))/2,linspace(-1,1,size(KBODY,2))/2,linspace(-1,1,size(KBODY,3))/2);
    

% calculate shift between prescan and scan taking into acount the prescan
% pixel size
    PreScanPz=twix_obj{1,1}.hdr.Config.ReadFoV/twix_obj{1,1}.hdr.Config.BaseResolution;
    PreScanPxy=PreScanPz/twix_obj{1,1}.hdr.Config.PhaseResolution;

    CENTER_OFFSET=single(exp(-2*pi*1i*(kx.*Mx/PreScanPxy+ky.*My/PreScanPxy+kz.*Mz/PreScanPz)));
    
    % reconstruct 3D images of the body and surface coils that have been
    % shifted to the center of the scan FOV
    tempBODY=fft3c(CENTER_OFFSET.*KBODY);    
    tempSURFACE=fft3c(CENTER_OFFSET.*KSURFACE);
    


    % Interpolate the prescan images to the full scan FOV and a desired
    % ReconResolution. Assuming you are using this for coil sensitivity
    % estimates I have found ReconResolution = 48-64 to be good values.
    
    PZ=0:PreScanPz:PreScanPz*2*twix_obj{1,1}.hdr.Config.BaseResolution;PZ=PZ(1:end-1)-twix_obj{1,1}.hdr.Config.ReadFoV;
    PXY=0:PreScanPxy:PreScanPxy*twix_obj{1,1}.hdr.Config.BaseResolution/2;PXY=PXY(1:end-1)-twix_obj{1,1}.hdr.Config.ReadFoV/2;
        
    
    [px,py,pz]=ndgrid(PXY,PXY,PZ);
    S=linspace(-twix_obj{1,length(twix_obj)}.hdr.Config.ReadFoV,twix_obj{1,length(twix_obj)}.hdr.Config.ReadFoV,ReconResolution+1);S=S(1:end-1);
    [sx,sy,sz]=ndgrid(S);
    
    IBODY=zeros(ReconResolution,ReconResolution,ReconResolution,size(tempBODY,4),'single');
    for iCoil=1:size(tempBODY,4)
        F = griddedInterpolant(px,py,pz,(double(tempBODY(:,:,:,iCoil))),'spline');
        IBODY(:,:,:,iCoil)=F(sx,sy,sz);
    end
    
    
    ISURFACE=zeros(ReconResolution,ReconResolution,ReconResolution,size(tempSURFACE,4),'single');
    for iCoil=1:size(tempSURFACE,4)
        F = griddedInterpolant(px,py,pz,(double(tempSURFACE(:,:,:,iCoil))),'spline');
        ISURFACE(:,:,:,iCoil)=F(sx,sy,sz);
    end
    
    
    break
end
end