function [fullsol,enuf,kx] = ReadDumpOutput(filein,nVs)
% Reads data from h5 file and puts into a matrix in the same format as my
% other code mat(ny,nx,nz)
global file num_Lin num_MF;
% file = '/Users/jsquire/Documents/QL_DNS/QL_DNS/Data/CompTest/FullSolution.h5';
file = filein;
num_Lin = nVs(1);
num_MF = nVs(2);

% Get info, work out dimensions and number of time slices
hinfo = hdf5info(file);

slicenames = {};
num_slices = length(hinfo.GroupHierarchy.Groups);
for kk=1:num_slices
    slicenames = [slicenames hinfo.GroupHierarchy.Groups(kk).Name];
end

%kz and ky
ky = hdf5read(file,'/ky');
kz = hdf5read(file, '/kz');
% nz
nz_full_choices = [16 32 48 64 96 128 192 256 384 512]; % Choices for non-dealiased size
nz = length(kz); nxy = length(ky); 
[~,ind] = min(abs(nz*3/2-nz_full_choices));
nz_full = nz_full_choices(ind);
disp(['Using NZ full = ' num2str(nz_full)]);
% Form ky and kx index vectors
ny = diff(find(ky==0,2)); nx= nxy/ny;
kyind = ky/ky(2);
kxind = reshape(repmat([0:(nx-1)/2 -(nx-1)/2:-1],[ny,1]),[nxy,1]);
indkxky = [(0:nxy-1)' kxind kyind]; %Indices as labeled in file
mind = round([reshape(repmat(1:nx,[ny,1]),[nxy,1]) kyind+1]);% Index for matrix 

% Loop through individual slices
disp('Found time-slices:')
disp(slicenames)
fullsol = cell(num_slices,1);
to_get = 1:length(slicenames);
% This is not correct right now! Need to reorder kx values before putting
% them into matrix

dealiasinds = [1:floor(nz_full/3)+1 ceil(2*nz_full/3)+1:nz_full];

for kk = to_get
    % Find all the kx data
    kxdata = hdf5read(file,[slicenames{kk} '/kx']); 
    
    fullsol{kk}.time = slicenames{kk};
    fullsol{kk}.data = zeros(ny,nx,nz_full,num_Lin);
    
    kxtmp = zeros(ny,nx);
    kytmp = reshape(ky,[ny,nx]);
    datatmp = zeros(ny,nx,nz,num_Lin);
    % Get data, unordered in kx
    for jj = 1:size(indkxky,1)
        datatmp(mind(jj,2),mind(jj,1),:,:) = kxkyfetch(slicenames{kk}, indkxky(jj,:));
        kxtmp(mind(jj,2),mind(jj,1)) = kxdata(indkxky(jj,1)+1);% Storing the kx value associated with each element
    end
    % Go through each y slice and reorder
    for yy = 1:ny
        [~,order] = sort(kxtmp(yy,:));
        datatmp(yy,:,:,:) = datatmp(yy,ifftshift(order),:,:);
    end
    % Add in zeros since dealiased version is saved
    fullsol{kk}.data(:,:,1:(nz/2+1/2),:) = datatmp(:,:,1:(nz/2+1/2),:);
    fullsol{kk}.data(:,:,(nz_full-nz/2+3/2):nz_full,:) = datatmp(:,:,(nz/2+3/2):nz,:);

    MFtmp = hdf5read(file,[slicenames{kk} '/Mean']); 
    MFtmp = (MFtmp(1,:)+1i*MFtmp(2,:)).';
    fullsol{kk}.MF = zeros(nz_full,num_MF);
    for mf = 1:num_MF
        os = (mf-1)*nz;
        fullsol{kk}.MF(1:(nz/2+1/2),mf) = MFtmp(1+os:(nz/2+1/2)+os);
        fullsol{kk}.MF((nz_full-nz/2+3/2):nz_full,mf) = MFtmp(os+(nz/2+3/2):nz+os);
    end
    

    
    % Energy
    mfac = 2*ones([ny,nx,nz]); 
    kzf = repmat(reshape(kz.^2,[1,1,nz]),[ny,nx,1]);
    kyf = repmat(kytmp.^2,[1,1,nz]);
    kxf = repmat(kxtmp.^2,[1,1,nz]);
    lap2=kyf.^2+kzf.^2;
    lapF=kxf.^2+kyf.^2+kzf.^2;
    lap2(kyf==0 & kzf==0)=1;
    mfac(kyf==0)=1;
    enu{kk}= lapF./lap2.*abs(fullsol{kk}.data(:,:,dealiasinds,1)).^2+1./lap2.*abs(fullsol{kk}.data(:,:,dealiasinds,2)).^2;
    enu{kk} = mfac.*enu{kk}/(2*(nx+1)^2*(2*ny)^2*nz_full^2);
end

enuf=0;
for kk=to_get
    enuf = enuf+enu{kk};
end
enuf=enuf/length(to_get);
    kx = kxtmp(1,:);
end



function data = kxkyfetch(t, inds)

global file num_Lin;

   % t is string with time, inds (index, kx, ky)
   slicename = [t '/Ind ' num2str(inds(1)) ': kx ' num2str(inds(2)) ' ky ' num2str(inds(3))]; 
   
   data = hdf5read(file,slicename); 
   data = (data(1,:)+1i*data(2,:)).';
   
   data = reshape(data,[length(data)/num_Lin num_Lin]);
   
end


function out=Rifftn(in)
% Takes ifftn on a fourier function that is missing half of its first
% dimension

    dims=size(in);
    parft=ifft(ifft(in,[],2),[],3);
    out=ifft([parft;zeros(1,dims(2),dims(3));conj(parft(end:-1:2,:,:))], [],1);

end


function energy=solEnergy(P,K,Ckl)

energy=zeros(P.nx-1,P.ny/2);
% Calcuates the energy of Ckl solution

mfac=[1;2*ones(P.ny/2-1,1)];

for yyy=1:P.ny/2
    for xxx=1:P.nx-1
        % Various k values and matrices
        ky=K.KY(yyy)
        kxt = K.kx_array(xxx,yyy)
        % Form laplacians
        lap2=ky^2 + K.KZ.^2;
        lapF=ky^2+kxt^2+K.KZ.^2;
        if ky==0
            lap2(1)=1;% Shouldn't be energy in this mode anyway
        end

        Mkl=abs([lapF./lap2;1./lap2;lapF./lap2;1./lap2]);
%         if ts==0
%         disp(['kx = ', num2str(K.kx(xxx)), 'ky = ', num2str(K.ky(yyy))])
%         disp(num2str(mfac(yyy)*sum(Mkl.*diag(Ckl(:,:,xxx,yyy)))));
%         end

        energy(xxx,yyy)=mfac(yyy)*sum(Mkl.*diag(Ckl(:,:,xxx,yyy)));
        
        if kxt==0.0 && ky==0.0
            if sum(sum(abs(Ckl(:,:,xxx,yyy))))>1e-10
                warning('Energy in kx=ky=0 fluctuation tensor')
            end
        end

    end
end
energy=energy/(2*P.nx^2*P.ny^2*P.nz^2);

end
