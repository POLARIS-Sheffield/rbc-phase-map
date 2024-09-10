function outTable = sw_phase_map(S_all, S_min, j_pi, Nkey, Mask, saveStatus)

%% Sliding Window Keyhole Phase Mapping

% Calculates and plots RBC oscillation amplitude and phase

% ------------------------------------------------------------------------%

% Inputs:

% Calculated from sw_keyhole_recon.m:
% S_all --- [mtx_reco, mtx_reco, mtx_reco, Nkey], all keyhole images
% S_min --- [mtx_reco, mtx_reco, mtx_reco], min RBC signal across all
% keyhole images
% j_pi --- [mtx_reco, mtx_reco, mtx_reco], index of keyhole image where
% chosen projections are pi out of phase with the projections for the
% first keyhole image
% NKey --- number of keyhole images

% Other:
% Mask --- lung mask from dissolved 129Xe recon
% saveStatus --- 0 = do not save, 1 = save

% ------------------------------------------------------------------------%

% Outputs:

% outTable --- summary of results

% ------------------------------------------------------------------------%

% Contributors:

% Jemima H. Pilgrim-Morris

% ------------------------------------------------------------------------%

%%

if nargin < 6
    saveStatus = 0;
end

%% Find S_max and j_max

    mtx_reco = size(Mask,1);
    S_all_smooth = smoothdata(S_all,4,"gaussian",10);
    
    for xx = 1:mtx_reco
        for yy = 1:mtx_reco
            for zz = 1:mtx_reco
                % for every pixel
                [val_max, key_max] = max(S_all_smooth (xx,yy,zz,:));  % find the max value across all keys
                list_peaks = S_all_smooth (xx,yy,zz,:)>=val_max*0.9;  % find values above 0.9*max value
                q = [];
                for i = 1:size(S_all_smooth ,4)
                    if key_max-i > 20
                        if list_peaks(i) == 1
                            q = [q i];
                        end
                    end
                end
                if size(q,1) >= 1
                    [val_max, q_max] = max(S_all_smooth (xx,yy,zz,q));
                    key_max = q(q_max);
                end
                S_max(xx,yy,zz) = val_max;
                j_max(xx,yy,zz) = key_max;
            end
        end
    end
    
    %% For plotting

    s = 8;
    e = 23;
    
    custom_map = [[0,0,0];[0, 0.2, 0.4];colMapGen([0 0.33 0.33],[0 0.27 0.49],Nkey-3);[0,0.2,0.4]];
    
    %% Calculate phase map
    
    % Convert j_max to phase
    j_max = j_max.*Mask;   
    phase_map = zeros(mtx_reco,mtx_reco,mtx_reco);
    
    for xx = 1:mtx_reco
        for yy = 1:mtx_reco
            for zz = 1:mtx_reco
                if Mask(xx,yy,zz) == 1
                    if j_max(xx,yy,zz) > j_pi
                        phase_map(xx,yy,zz) = pi*(j_max(xx,yy,zz)-Nkey-1)/(Nkey-j_pi);
                    else
                        phase_map(xx,yy,zz) = pi*(j_max(xx,yy,zz)-1)/(j_pi-1);
                    end
                else
                    phase_map(xx,yy,zz) = -pi -0.01;
                end
            end
        end
    end

    rmsPhase = rms(phase_map(Mask));
    sdPhase = std(phase_map(Mask));
    cvPhase = sdPhase/rmsPhase;
    medPhase = median(phase_map(Mask));
    minPhase = min(phase_map(Mask));
    maxPhase = max(phase_map(Mask));
    iqrPhase = iqr(phase_map(Mask));
    
    % plot phase
    figure
    ax = gca;
    ax.FontSize = 15;
    montage(phase_map, Indices=s:e, DisplayRange=[-pi-0.01, pi], ThumbnailSize=[]);
    c = colorbar;
    c.Label.String = ('\phi (rad)');
    c.Ticks = [-pi, 0, pi];
    c.TickLabels = {'-\pi', '0', '\pi'};
    c.FontSize = 15;
    titleInfo = sprintf('Median = %.2f rad, IQR = %.2f rad', medPhase, iqrPhase);
    title(titleInfo);
    colormap(gca,custom_map)
    set(gcf, 'Position', get(0, 'Screensize'));
    if saveStatus == 1
        print('RBC_Osc_Map_SW_phase','-r300','-dpng');
    end
    
    
    %% Calculate alpha_SW
    
    % normalise pixelwise by S_mean
    
    S_diff = S_max - S_min;
    S_mean = mean(S_all, 4);
    alpha_SW = 100*Mask.*S_diff./S_mean;
    
    mRBCOsc=mean(alpha_SW(Mask));
    sdRBCOsc=std(alpha_SW(Mask));
    medRBCOsc=median(alpha_SW(Mask));
    maxRBCOsc=max(alpha_SW(Mask));
    minRBCOsc=min(alpha_SW(Mask));
    iqrRBCOsc=iqr(alpha_SW(Mask));
    cvRBCOsc = sdRBCOsc/mRBCOsc;
    perc_highRBCOsc = 100*sum(alpha_SW(Mask) > mRBCOsc+sdRBCOsc, 'all')/size(alpha_SW(Mask),1);
    perc_lowRBCOsc = 100*sum(alpha_SW(Mask) < mRBCOsc-sdRBCOsc, 'all')/size(alpha_SW(Mask),1);

    % plot alpha_SW
    figure
    ax = gca;
    ax.FontSize = 15; 
    montage(alpha_SW, Indices=s:e, DisplayRange=[0, 50], ThumbnailSize=[]);
    c = colorbar;
    c.Label.String = ('\alpha_{SW} (%)');
    c.FontSize = 15;
    titleInfo = sprintf('Median = %.2f %%, IQR = %.2f %%', medRBCOsc, iqrRBCOsc);
    title(titleInfo);
    colormap("hot")
    set(gcf, 'Position', get(0, 'Screensize'));
    if saveStatus == 1
        print('RBC_Osc_Map_SW','-r300','-dpng');
    end

    outTable = table(mRBCOsc, sdRBCOsc, cvRBCOsc, medRBCOsc, minRBCOsc, maxRBCOsc, iqrRBCOsc, perc_highRBCOsc, perc_lowRBCOsc,  rmsPhase, sdPhase, cvPhase, medPhase, minPhase, maxPhase, iqrPhase);

    % save the data
    if saveStatus == 1
        save('sw_phase_recon')
        % save as a table 
        writetable(outTable,'RBCOscMap_SW_RESULTS.xlsx'); 
    end

% end of function
end