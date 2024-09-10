function [alpha_k0, j_pi, Nkey] = sw_keyhole_recon(RBC_k0, ktraj, ind, nps, saveStatus)

%% Sliding Window Keyhole Reconstruction

% calculates alpha_k0
% selects k0 projections for the keyhole reconstruction using a sliding window
% calculates iterative dcf for each keyhole
% finds j_pi

% ------------------------------------------------------------------------%

% Inputs:

% RBC_k0 --- [nspokes, 1], RBC k0 signal
% ktraj --- [1, nspokes x nsamples, 3], k-space trajectory
% ind --- logical array for reconstruction
% nps --- number of samples per spoke
% saveStatus --- 0 = do not save, 1 = save

% ------------------------------------------------------------------------%

% Outputs:

% alpha_k0 --- k0 oscillation amplitude
% S_all --- [mtx_reco, mtx_reco, mtx_reco, Nkey], all keyhole images
% S_min --- [mtx_reco, mtx_reco, mtx_reco], min RBC signal across all
% keyhole images
% j_pi --- [mtx_reco, mtx_reco, mtx_reco], index of keyhole image where
% chosen projections are pi out of phase with the projections for the
% first keyhole image
% NKey --- number of keyhole images

% ------------------------------------------------------------------------%

% Note that this script is missing the actual gridding reconstruction so
% will not calculate S_all or S_min. These are provided in the example data

% ------------------------------------------------------------------------%

% Contributors:

% Guilhem J. Collier
% Jemima H. Pilgrim-Morris
% Mika Takigawa

% ------------------------------------------------------------------------%

%%

if nargin < 4
    nps = 13;
end

if nargin < 5
    saveStatus = 0;
end

%% Calculate alpha_k0

t = 0:0.015:0.015*933;

% peak detection
[~, locs_Pos] = findpeaks(RBC_k0,'MinPeakDistance',40);
[~, locs_Neg] = findpeaks(-RBC_k0,'MinPeakDistance',40);

figure
set(gca,'LineWidth',2,'FontName','Arial','FontSize',14);
plot(t,RBC_k0+1,LineWidth=2)
hold on
plot(t(locs_Pos),RBC_k0(locs_Pos)+1,'rv','MarkerFaceColor','r')
plot(t(locs_Neg),RBC_k0(locs_Neg)+1,'rs','MarkerFaceColor','b')
title('Peak Detect')
xlabel('Time (s)')
ylabel('RBC Signal (a.u.)')
set(gcf,'Color','w');

% choose high SNR peaks in the first 7s to calculate alpha_k0
locs_Pos_good = locs_Pos(2:8);
locs_Neg_good = locs_Neg(2:8);

figure
set(gca,'LineWidth',2,'FontName','Arial','FontSize',14);
plot(t,RBC_k0+1,LineWidth=2)
hold on
plot(t(locs_Pos_good),RBC_k0(locs_Pos_good)+1,'rv','MarkerFaceColor','r')
plot(t(locs_Neg_good),RBC_k0(locs_Neg_good)+1,'rs','MarkerFaceColor','b')
alpha_k0 = 100*(mean(RBC_k0(locs_Pos_good)+1)-mean(RBC_k0(locs_Neg_good)+1));
xlabel('Time (s)')
ylabel('RBC Signal (a.u.)')
titleInfo = sprintf('alpha_{k0} = %.2f %%', alpha_k0);
title(titleInfo, FontSize=14);
set(gcf,'Color','w');
if saveStatus == 1
    print('RBC_Osc_k0','-r300','-dpng');
end

%% Sliding Window

% choose peaks to include in sliding window
locs_Pos2 = locs_Pos(2:end-1);
locs_Neg2 = locs_Neg(2:end);
nPos = size(locs_Pos2,1);

% number of points per peak (so that each key contains 20% of the spokes)
ppeakPos = round(size(RBC_k0,1)/(5*nPos));

% pdiff = points between adjacent maxima
% calculate number of keyhole images, Nkey
pdiff = []; 
for i = 1:nPos-1
    pdiff = [pdiff (locs_Pos2(i+1) - locs_Pos2(i))];
end

pcycle = mean(pdiff);
Nkey = round(pcycle);

%% Sliding Window Recon
% commented out because we have removed the actual recon code
% but here for reference :)

% nk = size(kdata,1)*nps;
% mtx_reco = size(Mask,1);

% j = 1;
% 
% for c = 0:Nkey-1
% 
%     ykey = zeros(size(kdata));
%     key = zeros(size(kdata));
%     if rem(ppeakPos,2) == 0
%         for i = 1:nPos-1
%             a = locs_Pos2(i) - (ppeakPos-2*(rem(i,2)))/2;
%             b = locs_Pos2(i) + (ppeakPos-2*(rem(i+1,2)))/2;
%             ykey((a+c):(b+c)) = kdata((a+c):(b+c));
%             key((a+c):(b+c)) = 1;
%         end
%     else
%         for i = 1:nPos-1
%             a = locs_Pos2(i) - (ppeakPos-1)/2;
%             b = locs_Pos2(i) + (ppeakPos-1)/2;
%             ykey((a+c):(b+c)) = kdata((a+c):(b+c));
%             key((a+c):(b+c) )= 1;
%         end
%     end
% 
%     key = logical(key);
%     
%     % generate ind for the keys + keyholes
%     key_pt = 6;
% 
%     ind_key = ind;
%     for i = 3:(key_pt+2)
%         ind_key(21:end,i) = key;
%     end
%     ind_key = logical(ind_key);
%     
%     % generate indices for keys
%     indices_key = [];
%     for i = 0:size(key)-1
%         x = (i*13)+1;
%         if key(i+1) == 1 
%             indices_key = [indices_key x:(x+key_pt-1)];
%         end
%     end
%     
%     % generate indices for keyholes
%     hole_ind = [];
%     for x = 0:((nk/13)-1)
%         hole_ind = [hole_ind (key_pt+1:13)+x*13];
%     end
%     
%     % combine and reorder
%     index_key = sort([hole_ind indices_key]);
%     
%     
%     % k-space trajectories for the keys
%     k_key = k(:,index_key,:);
%     
%     % iterative dcf
%     dcf_key = calc_sdc3(k_key,mtx_reco,15);
% 
%     % IMAGE RECON HERE %
%     % each iteration in the loop reconstructs S_key
%     % S_all = S_key for every key
%     % S_all = (mtx_reco, mtx_reco, mtx_reco, Nkey)
% 
% 
% %     if j == 1
% %         S_min = S_key;
% %     else
% %         S_min = min(S_min, S_key);
% %     end
% 
%     j = j + 1;
%  end

%% Find j_pi
      
% keys through which to look for j_pi
startKey = round(Nkey/4);
endKey = round(Nkey*(3/4));

% log y values
yPointLog = zeros(size(locs_Pos2,1), endKey-startKey);

% calculate difference between each centre point of the iteration and
% its corresponding actual negative peak
dNegPeak = zeros(size(locs_Pos2,1), endKey-startKey);

%total discrepancy with real value per key
discKeyTotal = []; 
    
z=1;
for w=startKey:endKey
    loc2y = locs_Pos2+w;
    if loc2y(end) > size(RBC_k0,1)
        loc2y(end) = [];
        yPointLog(end) = [];
    end
        
    yPointLog(:,z) = RBC_k0(loc2y)+1;
    dNegPeak(:,z) = abs(RBC_k0(locs_Neg2)-RBC_k0(loc2y));
    discKeyTotal = [discKeyTotal sum(dNegPeak(:,z))];
    z=z+1;
            
    %image produced     
    figure
    set(gca,'LineWidth',2,'FontName','Arial','FontSize',14);
    plot(t,RBC_k0+1)
    hold on
    plot(t(loc2y),RBC_k0(loc2y)+1, '.', MarkerSize=20);
    hold off
    xlabel('Time (s)', FontSize=14)
    ylabel('RBC Signal (a.u.)')
    title(['j = ' num2str(w)]);
    set(gcf,'Color','w');
end

leastDisc = find(discKeyTotal == min(discKeyTotal));
j_pi = (startKey+leastDisc)-1;

% end of function    
end 