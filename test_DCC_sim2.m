%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example script for Dynamic Correlation Toolbox
%
% File created by Martin Lindquist 07/22/14
%
% Makes use of functions from the UCSD_Garch toolbox by Kevin Shepard (Please see license agreement)
%
% Before running this script, begin by adding the DC_toolbox and all its subdirectories to the Matlab path.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create simulated data set 
% tic
% Define data dimensions
K = 10; % dimension of time series
T = 200; % time series length
burnIn = 1000;

% parameters for modelling the volatility of each dimension as a Garch(1,1) process
theta = [0.01 0.05 0.9]; 
% parameters for conditional correlations
para = [0.05 0.93];
% generate data with true DCC model
%       r = T-by-K data matrix
%       H0: K*K*T matrix containing H_1,...,H_T, where H_t is the true conditional covariance matrix at time t. 
%       R0: K*K*T matrix containing R_1,...,R_T, where R_t is the true conditional correlation matrix at time t. 
for rep = 1: 10
    rep
[ r, H0, R0 ] =   generateData( K, T, theta, para, burnIn );
dat = r - repmat(mean(r),T,1); % demean    

% Note the input data has dimensions T-by-p (time by #nodes)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit DCC

% Simple quick version that assumes common dynamic across all nodes. Not
% recommended for problems with larger values of p.

[Ct1 ] = DCCsimple(dat);

% Slower, more accurate version
[Ct2 ] = DCC(dat);
% Ct1 is the dynamic correlation matrix 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit sliding-window correlations
windowsize = 20;
[ Ct3 ] = sliding_window(dat,windowsize);
% toc
% Ct2 is the sliding window correlation matrix 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot some of the results


% figure
% subplot 221
% imagesc(Ct2(:,:,100), [-1 1])          % Plot the conditional correlation matrix at time 100
% colorbar
% title('DCC - conditional correlation at time 100')
% 
% subplot 222
% plot(squeeze(Ct2(1,3,:)))    % Plot the dyanamic correlation between nodes 1 and 3.
% ylim([-0.7 0.7])
% hold
% plot(1:T,squeeze(R0(1,3,:)),'-r')    % Plot the true dyanamic correlation between nodes 1 and 3.
% title('DCC - dynamic correlation between nodes 1 and 3')
% 
% 
% subplot 223
% imagesc(Ct3(:,:,100), [-1 1])          % Plot the conditional correlation matrix at time 100
% colorbar
% title('SWC - conditional correlation at time 10')
% 
% subplot 224
% plot(squeeze(Ct3(1,3,:)))    % Plot the dyanamic correlation between nodes 1 and 3.
% ylim([-0.7 0.7])
% hold
% plot(1:T,squeeze(R0(1,3,:)),'-r')    % Plot the true dyanamic correlation between nodes 1 and 3.
% title('SWC - dynamic correlation between nodes 1 and 3')

for i = 1 : T
    SWC_v(i,:) = mat2vec(squeeze(Ct3(:,:,i)));
    DCC_v(i,:) = mat2vec(squeeze(Ct2(:,:,i)));
    TrueC(i,:) = mat2vec(squeeze(R0(:,:,i)));
end
for i = windowsize+1 : size(SWC_v,1)
    sim_SWC_w(rep,i) = corr(SWC_v(i,:)',TrueC(i,:)');
     sim_DCC_w(rep,i) = corr(DCC_v(i,:)',TrueC(i,:)');
end
for i = 1 : size(SWC_v,2)
    sim_SWC_n(rep,i) = corr(SWC_v(windowsize+1:end,i),TrueC(windowsize+1:end,i));
     sim_DCC_n(rep,i) = corr(DCC_v(windowsize+1:end,i),TrueC(windowsize+1:end,i));
end
end

sim_SWC_w(:,1:windowsize) = [];
sim_DCC_w(:,1:windowsize) = [];