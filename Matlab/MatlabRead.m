% MatlabRead.m

n2s = @(str) strrep(num2str(str),'.','');

% Reads in data produced by MRIDSS and plots basic things
base_dir='/Users/jsquire/Documents/QL_DNS/';
data_dir = 'QL_DNS/Data/';

run_dir = ['FullQL_Pm2Noise2Rm16000_L0511/'];
% run_dir = ['SW_test/'];

MFdim=128; % Dimension in the mean fields
numMF = 4; % number of mean fields
num_Energy_AM = 4;

fid_time = fopen([base_dir data_dir run_dir 'time_vec.dat']);
fid_en =  fopen([base_dir data_dir run_dir 'energy.dat']);
fid_AM =  fopen([base_dir data_dir run_dir 'angular_momentum.dat']);
fid_diss =  fopen([base_dir data_dir run_dir 'dissipation.dat']);
fid_MF =  fopen([base_dir data_dir run_dir 'mean_fields.dat']);

% Read time vector and simulation length
time = fread(fid_time, inf ,'double');
time = time(1:end-1);
sim_len= length(time);

figure
% Read in energy and plot
subplot(312)
energy = fread(fid_en,[num_Energy_AM,sim_len],'double');
plot(time , sum(energy), 'k',time, energy  );
title('Energy');
legend('Total','Mean U','Mean B','Fluct u','Fluct b')

% Angular momentum
AM = fread(fid_AM,[num_Energy_AM,sim_len],'double');
tAM=AM(1,:)+AM(3,:)-AM(2,:)-AM(4,:);
subplot(313);
plot(time,tAM,'k',time, (AM) );
title('AM stress');
legend('Total','Mean U','Mean B','Fluct u','Fluct b')

% % Dissipation
% diss = fread(fid_diss,[num_Energy_AM,sim_len],'double');
% subplot(312);
% TdissU = (((1/8000)^3./(diss(1,:)+diss(3,:))).^0.25)*(2*pi*128*(1/3));
% TdissB = (((1/16000)^3./(diss(2,:)+diss(4,:))).^0.25)*(2*pi*128*(1/3));
% plot(time, TdissU, time ,TdissB);%plot(time,sum(diss),'k',time, diss,'--' );
% title('Dissipation');
% legend('Total','Mean U','Mean B','Fluct u','Fluct b')

% Mean fields;
% tmp = fread(fid_MF , inf ,'double');
% sim_len = length(tmp)/(MFdim*numMF);
tmp = fread(fid_MF , MFdim*numMF*sim_len ,'double');
if length(tmp)~=MFdim*numMF*sim_len
    error('Number of NF elements does not match up');
end
tmp = reshape(tmp,[MFdim*numMF,sim_len]);
for ii=0:numMF-1
    MF{ii+1} = tmp(ii*MFdim+1:(ii+1)*MFdim,:);
    MFmax{ii+1} = max(abs(MF{ii+1}));
    MFmean{ii+1} = sqrt(sum(abs(MF{ii+1}).^2)/MFdim);
end
subplot(311);
htime = 1;floor(length(time)/2);
% plot(time(1:htime), log10(MFmean{2}(1:htime)),time(1:htime),
% log10(MFmean{1}(1:htime)),...
%     time(1:htime), log10(MFmean{3}(1:htime)),time(1:htime), log10(MFmean{4}(1:htime)))
% plot( linspace(0,2*pi,MFdim), MF{2}(:,end),linspace(0,2*pi,MFdim), MF{1}(:,end) )
% Remove repeated times - annoying but oh well
normalts = find(diff(time)>0);
time = time(normalts);
zp=linspace(0,1,MFdim);
MFplot= MF{2}(:,normalts);
% MFplot = log10(abs(fft(MFplot)).^2);MFplot=MFplot(2:floor(MFdim/3),:);
imagesc(time(htime:end), zp,MFplot(:,htime:end))%,20,'Linestyle','none')
colorbar;
% MFplot= MFmean{2}(normalts);
% plot(time,MFplot.^2)

fclose(fid_en);
fclose(fid_time);
fclose(fid_AM);
fclose(fid_diss);
fclose(fid_MF);


