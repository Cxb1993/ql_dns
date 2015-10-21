function WriteInitialConditions(ICs)
% Writes initial conditions for QL_DNS
% Format is an HDF5 file, same as is saved by QL_DNS itself, save data in
% /t=0 group
global file
QLdns_input = 'SW_test';

fullfolder = ['/Users/jsquire/Documents/QL_DNS/QL_DNS/Data/' QLdns_input];
if ~exist(fullfolder,'dir')
    mkdir(fullfolder);
end
file = [fullfolder '/FullSolution.h5'];
if exist(file, 'file')
    delete(file)
end

[NY,NX,NZ] = size(ICs.u);
NZde = 2*floor(NZ/3)+1;
dealias = [0:NZ/2-1 -NZ/2:-1];nzV = 1:NZ;
dealias = nzV(abs(dealias)<=floor(NZ/3));
if length(dealias)~=NZde
    warning('Lengths dont match up!')
end
NXY = (NX-1)*NY;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Fluctuating variables
ind = 0;
for nx = [1:NX/2 NX/2+2:NX]
    for ny = 1:NY
        
        % Get data from ICs
        kxyslice = reshape(ICs.u(ny,nx,dealias),[NZde,1]);
        kxyslice = [kxyslice; reshape(ICs.zeta(ny,nx,dealias),[NZde,1])];
%         kxyslice = [kxyslice; reshape(ICs.b(ny,nx,dealias),[NZde,1])];
%         kxyslice = [kxyslice; reshape(ICs.eta(ny,nx,dealias),[NZde,1])];
        % Work out kx and index
        if nx<=NX/2
            kxi = nx;
        else
            kxi = nx-NX;
        end
        inds = [ind kxi-1 ny-1];
        % Write to hdf5
        kxkywrite(kxyslice, inds)
        
        ind = ind+1;
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean fields
MFdata = reshape(ICs.Bx(:,:,dealias),[NZde,1]);
MFdata = 0*[MFdata; reshape(ICs.By(:,:,dealias),[NZde,1])];
% if any(strcmp(properties(ICs),'Ux'))
%     MFdata = [MFdata; reshape(ICs.Ux(:,:,dealias),[NZde,1])];
%     MFdata = [MFdata; reshape(ICs.Uy(:,:,dealias),[NZde,1])];
% end
MFdata = [real(MFdata) imag(MFdata)].';
% Write to file
h5create(file,'/t=0/Mean',size(MFdata),'Datatype','single'); 
h5write(file,'/t=0/Mean', MFdata); 



end


function kxkywrite(data, inds)
% Writes data into kx ky slice of hdf5 file

global file;

   % t is string with time, inds (index, kx, ky)
   slicename = ['/t=0/Ind ' num2str(inds(1)) ': kx ' num2str(inds(2)) ' ky ' num2str(inds(3))]; 
   
   data = [real(data) imag(data)].';
   %Create file
   h5create(file,slicename,size(data),'Datatype','single'); 
  h5write(file,slicename, data);   
end