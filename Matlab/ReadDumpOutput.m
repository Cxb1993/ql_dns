function fullsol = ReadDumpOutput
% Reads data from h5 file and puts into a matrix in the same format as my
% other code mat(ny,nx,nz)
global file num_Lin num_MF;
file = '/Users/jsquire/Documents/QL_DNS/QL_DNS/Data/FullQL_Pm2Noise2_L184/FullSolution.h5';
num_Lin = 4;
num_MF = 4;

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
for kk = to_get
    % Find all the kx data
    kxdata = hdf5read(file,[slicenames{kk} '/kx']); 
    
    fullsol{kk}.time = slicenames{kk};
    fullsol{kk}.data = zeros(ny,nx,nz_full,num_Lin);
    
    kxtmp = zeros(ny,nx);
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
    

end


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
