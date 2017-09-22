% Elena Allen, 2/3/15, edited by Oliver Xie 
% Code adapted from Allen et al., "Tracking whole-brain connectivity dynamics in the resting state." Cereb Cortex. 2014 24(3):663-76. 
% Dependency: SimTB toolbox (http://mialab.mrn.org/software/simtb/)
%
% Objective: create a toy simulation with time varying connectivity
% In this example we model 10 components whose connectivity structure varies through 4 discrete states.
% Variables of interest are the state vector ('STATE') and the final TCs ('TC').
% The order and duration (dwell time) of states are controlled by variables 'Sorder' and 'Sdwell'
% The amplitude of unique events (aU) relative to the shared events will  make the modular structure more or less apparent.

%%
clear all; close all; clc
rng(100) % to enable repeated generation of the same simulation


% number of components
nC = 10;

% number of time points
nT = 400;

% TR
TR = 1.5;

% number of different connectivity states
nStates = 4;

% probability of unique events
pU = 0.5;

% amplitude of unique events (relative to module-specific events)
aU = .5;

% probability of state specific events 
pState = .5;

%Module membership for each state
ModMem = zeros(nC,nStates);

% Number of event types (i.e., number of different modules)
nE = 3;

% Modules are formed by sharing of events.  In each state there are up to
% nE different modules. The matrix ModMem stores the membership of each
% component to a module.Note tjat module 2 in state 2 has nothing to do with 
% module 2 in the other states, it's just an index.
% Negative numbers indicate activity is negatively related to the events in
% the given module.
ModMem(1,:) = [2   -2   3    2  ];
ModMem(2,:) = [2   -2   1    2  ];
ModMem(3,:) = [2   -2   3    -2  ];
ModMem(4,:) = [-2  3    3    2  ];
ModMem(5,:) = [-2  -1    2    3  ];
ModMem(6,:) = [-2  3    2    2  ];
ModMem(7,:) = [-2  2    2    1 ];
ModMem(8,:) = [1   2    -1   1   ];
ModMem(9,:) = [1   2    1   -2  ];
ModMem(10,:)= [1   2    1   -2  ];

% Check out the figure below -- should help make the ModMembership clear

%% Create figure of the connectivity matrix for each state
F = figure('color','w','Name', 'sim_neural_connectivity');

for ii = 1:nStates
    subplot(1,nStates,ii)
    CM = zeros(nC,nC);
    for jj = 1:nC
        for kk = 1:nC
            if ModMem(jj,ii) == ModMem(kk,ii)
                CM(jj,kk) = 1;
            elseif abs(ModMem(jj,ii)) == abs(ModMem(kk,ii))
                CM(jj,kk) = -1;
            else
                CM(jj,kk) = 0;
            end
        end
    end
    H = simtb_pcolor(1:nC, 1:nC, .8*CM);
    axis square; 
    axis ij
    set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])%, 'XColor', [1 1 1], 'YColor', [1 1 1])
    c = get(gca, 'Children');
    set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
    title(sprintf('State %d', ii))
end


%% Create the event time courses

% random aspects (different for each component)
eT = rand(nT, nC) < pU;
eT = eT.*sign(rand(nT, nC)-0.5);
eT = eT*aU;

% define the order and time in each state
Sorder = [1 2 3 2]; % state order
Sdwell = [100 100 100 100];
%NOTE: the Sdwell should sum to nT, check here and amend the last partition:
if sum(Sdwell) ~= nT
    Sdwell(end) = nT - sum(Sdwell(1:end-1));
end
Cdwell = cumsum(Sdwell);
Cdwell = [0 Cdwell];
STATE = zeros(1,nT); % state vector
for ii = 1:length(Sorder)
    sIND = Cdwell(ii)+1:Cdwell(ii+1);
    % events related to each module
    e = rand(length(sIND),nE) < pState;
    e = e.*sign(rand(length(sIND), nE)-0.5);
    for cc = 1:nC
        eT(sIND,cc) = eT(sIND,cc) + sign(ModMem(cc,Sorder(ii)))*e(:,abs(ModMem(cc,Sorder(ii))));
    end
    STATE(sIND) = Sorder(ii);
end

% event time series are stored in eT
%% Convolve event TCs
[tc, MDESC, P, PDESC] = simtb_TCsource(eT(:,1), TR, 1);
P =  zeros(7,nC);

P(1,:) = 6;     % delay of response (relative to onset)
P(2,:) = 15;    % delay of undershoot (relative to onset)
P(3,:) = 1;     % dispersion of response
P(4,:) = 1;     % dispersion of undershoot
P(5,:) = 3;     % ratio of response to undershoot
P(6,:) = 0;     % onset (seconds)
P(7,:) = 32;    % length of kernel (seconds)
P = P.*(1+0.01*(randn(7,nC))); 
TC  = zeros(nT,nC);
for cc = 1:nC
    TC(:,cc) = simtb_TCsource(eT(:,cc), TR, 1, P(:,cc)); % all use same HRF
    %TC(:,cc) = simtb_TCsource(eT(:,cc), TR, 1); % different HRFs
end

% Add a little gaussian noise
TC = TC + 0.1*randn(nT,nC);

dat = TC - repmat(mean(TC),nT,1); % demean    

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
CM = zeros(nC,nC,nStates);
% for ii = 1:nStates
%     for jj = 1:nC
%         for kk = 1:nC
%             if ModMem(jj,ii) == ModMem(kk,ii)
%                 CM(jj,kk,ii) = 1;
%             elseif abs(ModMem(jj,ii)) == abs(ModMem(kk,ii))
%                 CM(jj,kk,ii) = -1;
%             else
%                 CM(jj,kk,ii) = 0;
%             end
%         end
%     end
% end


for ii = 1:nStates
    CM(:,:,ii) = corr(TC(Sdwell2(ii):Sdwell2(ii+1)-1,:));
end
% Sorder = [1 5 2 5 3 5 4 5 2 5 1 5 4 5 3]; % state order
% Sdwell = [120 8 120 8 120 8 120 8 120 8 120 8 120 8 120];
Sdwell2 = cumsum([1 Sdwell]);
CM_t = zeros(nC,nC,nT);
for i = 1 :  length(Sorder)-1
    CM_t(:,:,Sdwell2(i):Sdwell2(i+1)-1) = repmat(CM(:,:,Sorder(i)),1,1,Sdwell(i));
end

for i = 1 : nT
    SWC_v(i,:) = mat2vec(squeeze(Ct3(:,:,i)));
    DCC_v(i,:) = mat2vec(squeeze(Ct2(:,:,i)));
    TrueC(i,:) = mat2vec(squeeze(CM_t(:,:,i)));
end
for i = windowsize+1 : size(SWC_v,1)
    sim_SWC_w(i) = corr(SWC_v(i,:)',TrueC(i,:)');
     sim_DCC_w(i) = corr(DCC_v(i,:)',TrueC(i,:)');
end
for i = 1 : size(SWC_v,2)
    sim_SWC_n(i) = corr(SWC_v(windowsize+1:end,i),TrueC(windowsize+1:end,i));
     sim_DCC_n(i) = corr(DCC_v(windowsize+1:end,i),TrueC(windowsize+1:end,i));
end
n = 2;
figure
plot(SWC_v(:,n))
hold on 
plot(DCC_v(:,n),'k')
plot(TrueC(:,n),'r')
hold off
legend('SWC','DCC','GT')


% %% Figure to display the states, TCs, and correlation matrices for each partition
% F=figure('color','w','Name', 'sim_TCs_CorrMatrices'); 
% subplot(4, length(Sorder), 1:length(Sorder))
% plot((0:nT-1)*TR, STATE , 'k', 'Linewidth', 1); axis tight; box off
% ylabel('State')
% set(gca, 'YTick', 1:nStates, 'XTick', Cdwell*TR, 'TickDir', 'out', 'Layer', 'Bottom'); grid on
% 
% subplot(4, length(Sorder), length(Sorder)+1:length(Sorder)*2)
% plot((0:nT-1)*TR, TC, 'LineWidth',0.75);
% xlabel('Time (s)')
% ylabel('Amplitude')
% set(gca, 'TickDir', 'out', 'XTick', Cdwell*TR, 'Xgrid', 'on'); 
% axis tight; box off
% 
% for ii = 1:length(Sorder)
%     sIND = Cdwell(ii)+1:Cdwell(ii+1);
%     if length(sIND)~=8
%         subplot(4,length(Sorder),length(Sorder)*3+ii)
%         temp = corr(TC(sIND,:));
%         H = simtb_pcolor(1:nC, 1:nC, temp);
%         axis square; axis ij 
%         set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])%, 'XColor', [1 1 1], 'YColor', [1 1 1])
%         c = get(gca, 'Children');
%         set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');        
%         text(1.5,-2,sprintf('Partition %d\nState %d', ii, Sorder(ii)), 'Fontsize', 8);
%     end
% end

% wsize = 20;
% % adaptive bandpass filtering
% TR = 1.5; % seconds
% Lcutoff = 1/wsize;
% Hcutoff = 0.18;
% NyqF = (1/TR)/2;
% N = 3;
% [B,A] = butter(N,[Lcutoff Hcutoff]/NyqF);
% for i = 1 : nC
%     TC(:,i) = filtfilt(B,A,TC(:,i));
% end
% 
% %% make the sliding window
% sigma = 3; % for smoothing Gaussian
% %wsize = 22; % for box (44 seconds)
% gw = gaussianwindow(nT,nT/2,sigma);
% b = zeros(nT,1);  b((nT/2 - wsize/2 + 1):(nT/2+wsize/2)) = 1;
% c = conv(gw,b); c = c/max(c); c=c(nT/2+1:end-nT/2 +1);
% 
% A = repmat(c,1,nC);
% 
% %% initialize
% Nwin = nT-wsize;
% FNCdyn = zeros(Nwin,nC,nC);
% 
% 
% Pdyn = zeros(Nwin, nC,nC);
% LAMBDA = zeros(1);
% lambdas = [.1:.03:.40];
% NLL = cell(1,1);
% nRep = 10;
% %%
% 
% speclength = (2^nextpow2(nT))/2+1;
% TR=1.5;
% params.tapers = [3 5]; params.Fs = (1/TR); params.fpass = [0.0, 1/(2*TR)];
% 
% %initialize
% sesInfo.numOfSub=1; % delete later
% sesInfo.numComp=nC; % delete later
% Pall = zeros(1, speclength, nC);
%  
% 
% 
% nTwin = length(find(c > 1e-4));
% for sub = 1:1
% %     tstart = tic;
% %     fprintf('Working on subject %d of %d\n', sub, M)
%     tc = TC;
%     tcwin = zeros(Nwin, nT, nC);
%     tcwin2 = zeros(Nwin, nTwin, nC);
%     %%
%     for ii=1:Nwin
%         Ashift =circshift(A,-nT/2 + wsize/2 + ii);
%         cshift = circshift(c,-nT/2 + wsize/2 + ii);
%         csind = find(cshift > 1e-4)';
% %         dfc = diff(csind);
% %         wrap_pnt = find(dfc>1);
% %         if ~isempty(wrap_pnt) 
% %             cwinind = [csind(wrap_pnt +1:end) csind(1:wrap_pnt)];
% %         else
% %             cwinind = csind;
% %         end       
%         tcwin(ii,:,:) = tc.*Ashift;        
% %         tcwin2(ii,:,:) = tcwin(ii,cwinind,:);
% %     [P,f] = compute_spectrum(tc, TR, 1);
%     end
%          for ii = 1:Nwin
%             a = corr(squeeze(tcwin2(ii,:,:)));
%             FNCdyn(ii,:,:) = a;
%          end
% %     telapsed = toc(tstart);
%    % fprintf('\tTime Elapsed: %0.2f minutes\n', telapsed/60)
%     
% end
% %% Fisher z-score the correlations
% FNCdyn = atanh(FNCdyn);
% FNCdynflat = zeros(Nwin, nC*(nC-1)/2);
% for ii = 1:Nwin
%     FNCdynflat(ii,:) = mat2vec(squeeze(FNCdyn(ii,:,:)));
% end


