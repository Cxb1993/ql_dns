% MatlabRead.m

n2s = @(str) strrep(num2str(str),'.','');
n2s = @(s) num2str(s);

% Reads in data produced by MRIDSS and plots basic things
base_dir='/Users/jsquire/Documents/QL_DNS/';
data_dir = 'QL_DNS/Data/';

noiseV = 6:2:50;
RmV = 4000:500:8000;
Bymat = zeros(length(noiseV),length(RmV));
Bxmat = zeros(length(noiseV),length(RmV));

% for nkk = 1:length(noiseV)
% for Rmkk = 1:length(RmV)
% run_dir = ['FullKNoise_Rm6000Pm8Noise4/'];
% run_dir = ['K0to20Noise_Rm' n2s(Rm) 'Pm' n2s(Pm) 'Noise' n2s(noise) '/']; % file name
Lz=8; q=2;
run_dir = ['YousefDynamo_Lz' n2s(Lz) 'q' n2s(q) '/']; % file name
run_dir = 'FullQL_Pm4Noise1_L148/';

MFdim=512; % Dimension in the mean fields
numMF = 4; % number of mean fields
num_Energy_AM = 4;
num_reynolds = 5;

fid_time = fopen([base_dir data_dir run_dir 'time_vec.dat']);
fid_en =  fopen([base_dir data_dir run_dir 'energy.dat']);
fid_AM =  fopen([base_dir data_dir run_dir 'angular_momentum.dat']);
fid_diss =  fopen([base_dir data_dir run_dir 'dissipation.dat']);
fid_MF =  fopen([base_dir data_dir run_dir 'mean_fields.dat']);
fid_rey =  fopen([base_dir data_dir run_dir 'reynolds_stress.dat']);

% Read time vector and simulation length
time = fread(fid_time, inf ,'double');
time = time(1:end-1);
sim_len= length(time);

figure
% Read in energy and plot
subplot(312)
energy = fread(fid_en,[num_Energy_AM,sim_len],'double');
semilogy(time , sum(energy), 'k',time, energy  );
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

% %Reynolds stress
% reynolds = fread(fid_rey,[num_reynolds,sim_len],'double');
% endp=length(time);
% mean(reynolds,2)
% % reynolds=filter(ones(1,100)/10,1,reynolds')';
% plot(time, reynolds(1:3,:),time,reynolds(4:5,:),'LineWidth',1);
% title('Terms in the dynamo equation');
% legend('Shear term','B_y emf','B_y dissipation','B_x emf','B_x dissipation')



% Mean fields;
% tmp = fread(fid_MF , inf ,'double');
% sim_len = length(tmp)/(MFdim*numMF);
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
subplot(311);
htime = 3;tf=floor(length(time));
% plot(time(1:htime), log10(MFmean{2}(1:htime)),time(1:htime),
% log10(MFmean{1}(1:htime)),...
%     time(1:htime), log10(MFmean{3}(1:htime)),time(1:htime), log10(MFmean{4}(1:htime)))
% plot( linspace(0,2*pi,MFdim), MF{2}(:,end),linspace(0,2*pi,MFdim), MF{1}(:,end) )
% Remove repeated times - annoying but oh well
normalts = [1; find(diff(time)>0)+1];
time = time(normalts);time = time(htime:tf);
zp=linspace(0,1,MFdim);
MFyplot= MF{2}(:,normalts);MFyplot=MFyplot(:,htime:tf);
MFxplot= MF{1}(:,normalts);MFxplot=MFxplot(:,htime:tf);
% MFplot = log10(abs(fft(MFplot)).^2);MFplot=MFplot(2:floor(MFdim/3),:);
imagesc(time, zp,MFyplot)%,20,'Linestyle','none')
colorbar;

% subplot(312)
% semilogy(time, sqrt(mean(MFyplot.^2,1)),time, sqrt(mean(MFxplot.^2,1)))
% timebnds = ginput(2);
% tb =[find(time>timebnds(1,1),1) find(time>timebnds(2,1),1)]; sort(tb);
% growthdata = MFmean{2}(normalts); growthdata = growthdata(htime:tf);
% growth_rate = (log(growthdata(tb(2)))-log(growthdata(tb(1))))/(time(tb(2))-time(tb(1)))


fclose(fid_en);
fclose(fid_time);
fclose(fid_AM);
fclose(fid_diss);
fclose(fid_MF);

% end
% end


