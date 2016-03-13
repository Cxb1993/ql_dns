% MatlabRead.m

n2s = @(str) strrep(num2str(str),'.','');
n2s = @(s) num2str(s);

% Reads in data produced by MRIDSS and plots basic things
base_dir='/Users/jsquire/Documents/QL_DNS/';
data_dir = 'QL_DNS/Data/';

run_dir = ['CompTest' '/'];
% run_dir = 'YousefN500_S2L16/';

MFdim=64; % Dimension in the mean fields
LZ=1;
numMF = 4; % number of mean fields
num_Energy_AM = 4;
num_reynolds = 4;
By=0.0001;

fid_time = fopen([base_dir data_dir run_dir 'time_vec.dat']);
fid_en =  fopen([base_dir data_dir run_dir 'energy.dat']);
fid_AM =  fopen([base_dir data_dir run_dir 'angular_momentum.dat']);
fid_diss =  fopen([base_dir data_dir run_dir 'dissipation.dat']);
fid_MF =  fopen([base_dir data_dir run_dir 'mean_fields.dat']);
fid_rey =  fopen([base_dir data_dir run_dir 'reynolds_stress.dat']);

% Read time vector and simulation length
time = fread(fid_time, inf ,'double');
time = time(1:end);
sim_len= length(time);

figure
% Read in energy and plot
subplot(312)
energy = fread(fid_en,[num_Energy_AM,sim_len],'double');
semilogy(time , sum(energy), 'k',time, energy  );
semilogy(time , energy(1,:));
% title('Energy');
if size(energy,1)==4
    legend('Total','Mean U','Mean B','Fluct u','Fluct b')
else
    legend('Total','Mean U','Fluct u')
end% 
% Angular momentum
AM = fread(fid_AM,[num_Energy_AM,sim_len],'double');
% tAM=AM(1,:)+AM(3,:)-AM(2,:)-AM(4,:);
tAM=AM(1,:)-AM(2,:);
subplot(313);
plot(time,tAM,'k',time, (AM) );
title('AM stress');
if size(energy,1)==4
    legend('Total','Mean U','Mean B','Fluct u','Fluct b')
else
    legend('Total','Mean U','Fluct u')
end
% % Dissipation
% diss = fread(fid_diss,[num_Energy_AM,sim_len],'double');
% subplot(312);
% TdissU = (((1/8000)^3./(diss(1,:)+diss(2,:))).^0.25)*(2*pi*128*(1/3));
% TdissB = (((1/16000)^3./(diss(2,:)+diss(4,:))).^0.25)*(2*pi*128*(1/3));
% plot(time, TdissU, time ,TdissB);plot(time,sum(diss),'k',time, diss,'--' );
% title('Dissipation');
% legend('Total','Mean U','Mean B','Fluct u','Fluct b')
% 
% dt=time(2)-time(1);
% disscheck = [-sum(energy);1.5*cumtrapz(time,tAM);-cumtrapz(time,sum(diss))];
% disscheck = [gradient(sum(energy),time);-(1.5*tAM-sum(diss))];
% plot(time,disscheck,time,sum(disscheck),'k')

%Reynolds stress
subplot(311)
reynolds = fread(fid_rey,[num_reynolds,sim_len],'double');
% reynolds=filter(ones(1,10)/10,1,reynolds')';
plot(time, reynolds);
title('Terms in the dynamo equation');
legend('Bx','By','Ux','Uy')

close gcf
plot(time, energy(3:4,:))
ginput


% Mean fields;
% tmp = fread(fid_MF , inf ,'double');
% sim_len = floor(length(tmp)/(MFdim*numMF));
% tmp=tmp(1:(sim_len*MFdim*numMF));
tmp = fread(fid_MF , MFdim*numMF*sim_len ,'double');
if length(tmp)~=MFdim*numMF*sim_len
    error('Number of MF elements does not match up');
end
tmp = reshape(tmp,[MFdim*numMF,sim_len]);
for ii=0:numMF-1
    MF{ii+1} = tmp(ii*MFdim+1:(ii+1)*MFdim,:);
    MFmax{ii+1} = max(abs(MF{ii+1}));
    MFmean{ii+1} = sqrt(sum(abs(MF{ii+1}).^2)/MFdim);
end
subplot(211);
htime = 10;tf=floor(length(time));
% plot(time(1:htime), log10(MFmean{2}(1:htime)),time(1:htime),
% log10(MFmean{1}(1:htime)),...
%     time(1:htime), log10(MFmean{3}(1:htime)),time(1:htime), log10(MFmean{4}(1:htime)))
% plot( linspace(0,2*pi,MFdim), MF{2}(:,end),linspace(0,2*pi,MFdim), MF{1}(:,end) )
% Remove repeated times - annoying but oh well
% tf = find(time>700);
normalts = [1; find(diff(time)>0)+1];
% normalts=normalts(1:tf);
time = time(normalts);time = time(htime:tf);
zp=0:LZ/MFdim:LZ-LZ/MFdim;
MFByplot= MF{2}(:,normalts);MFByplot=MFByplot(:,htime:tf);
MFBxplot= MF{1}(:,normalts);MFBxplot=MFBxplot(:,htime:tf);
% MFVyplot= MF{4}(:,normalts);MFVyplot=MFVyplot(:,htime:tf);
% MFVxplot= MF{3}(:,normalts);MFVxplot=MFVxplot(:,htime:tf);

la=@(f) f;% log10(abs(f));
imagesc(time, zp,la(MFBxplot))%,20,'Linestyle','none')
xlabel('t');ylabel('z');title('Ux')
colorbar;
% xlim([0 100])
subplot(212)
imagesc(time, zp,la(MFByplot))%,20,'Linestyle','none')
xlabel('t');ylabel('z');title('Uy')
colorbar;
xlim([0 100])
% 
% 
% 
% % Reynolds stresses when full solution is saved
% subplot(212)
% tmp = fread(fid_rey , MFdim*num_reynolds*sim_len ,'double');
% if length(tmp)~=MFdim*num_reynolds*sim_len
%     error('Number of Reynolds elements does not match up');
% end
% tmp = reshape(tmp,[MFdim*num_reynolds,sim_len]);
% kz = (2i*pi/LZ)*[0:MFdim/2-1  -MFdim/2:-1].';
% for ii=0:num_reynolds-1
%     Rey{ii+1} = tmp(ii*MFdim+1:(ii+1)*MFdim,:);
% end
% EMF{1} = fft(Rey{2})./repmat(kz,[1 sim_len]);EMF{1}([1,MFdim/2+1],:)=0;
% EMF{1} = ifft(EMF{1});
% EMF{2} = -fft(Rey{1})./repmat(kz,[1 sim_len]);EMF{2}([1,MFdim/2+1],:)=0;
% EMF{2} = ifft(EMF{2});
% 
% % meant = find(time>1,1);
% % etayx = mean(EMF{2}(2,meant:end))/(kz(2)*By*MFdim/2);
% % etat = mean(EMF{1}(2,meant:end))/(kz(2)*By*MFdim/2);
% % en = mean(energy(3,meant:end))
% % etayx
% % etat
% 
% 
% imagesc(time,zp,EMF{2})
% colorbar
% 
% % % For AnalyzeCorrelations
% % MFs.z= zp;
% % MFs.t=time;
% % MFs.bx = MF{1}(:,normalts);
% % MFs.by=MF{2}(:,normalts);
% % MFs.emfx = EMF{1}(:,normalts);
% % MFs.emfy = EMF{2}(:,normalts);
% % AnalyzeCorrelations(MFs)
% 
fclose(fid_en);
fclose(fid_time);
fclose(fid_AM);
fclose(fid_diss);
fclose(fid_MF);
fclose(fid_rey);

% end
% end


