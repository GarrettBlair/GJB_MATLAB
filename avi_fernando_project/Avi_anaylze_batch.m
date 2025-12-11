% toppath = 'G:/.shortcut-targets-by-id/1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY/';
% msDir =[toppath 'miniscope_v4/Mov_Grid/ms_data/'];
% exp = 'Mov_Grid';
% dirs = {[toppath 'miniscope_v4/Mov_Grid/A_baseline_6mov/camkII_mpfc/2023_09_12/15_49_27/' ],...
% [toppath 'miniscope_v4/Mov_Grid/A_shocks_12mov/camkII_mpfc/2023_09_13/15_16_35/' ],...
% [toppath 'miniscope_v4/Mov_Grid/A_retrieval_24mov/camkII_mpfc/2023_09_15/15_14_12/' ],...
% [toppath 'miniscope_v4/Mov_Grid/B_baseline_6mov/camkII_mpfc/2023_09_12/16_16_30/' ],...
% [toppath 'miniscope_v4/Mov_Grid/B_shocks_12mov/camkII_mpfc/2023_09_13/15_43_47/' ],...
% [toppath 'miniscope_v4/Mov_Grid/B_retrieval_24mov/camkII_mpfc/2023_09_15/15_43_03/' ],...
% [toppath 'miniscope_v4/Mov_Grid/C_baseline_6mov/camkII_mpfc/2023_09_12/16_42_31/' ],...
% [toppath 'miniscope_v4/Mov_Grid/C_shocks_12mov/camkII_mpfc/2023_09_13/16_11_58/' ],...
% [toppath 'miniscope_v4/Mov_Grid/C_retrieval_24mov/camkII_mpfc/2023_09_15/16_15_02/' ],...
% [toppath 'miniscope_v4/Mov_Grid/D_baseline_6mov/camkII_mpfc/2023_09_25/16_11_42/' ],...
% [toppath 'miniscope_v4/Mov_Grid/D_shocks_12mov/camkII_mpfc/2023_09_26/15_51_39/' ],...
% [toppath 'miniscope_v4/Mov_Grid/D_retrieval_24mov/camkII_mpfc/2023_09_28/15_34_33/' ],...
% [toppath 'miniscope_v4/Mov_Grid/E_baseline_6mov/camkII_mpfc/2023_09_25/16_42_58/' ],...
% [toppath 'miniscope_v4/Mov_Grid/E_shocks_12mov/camkII_mpfc/2023_09_26/16_15_43/' ],...
% [toppath 'miniscope_v4/Mov_Grid/E_retrieval_24mov/camkII_mpfc/2023_09_28/16_05_37/' ],...
% [toppath 'miniscope_v4/Mov_Grid/F_baseline_6mov/camkII_mpfc/2023_09_12/17_07_17/' ],...
% [toppath 'miniscope_v4/Mov_Grid/F_shocks_12mov/camkII_mpfc/2023_09_13/16_38_17/' ],...
% [toppath 'miniscope_v4/Mov_Grid/F_retrieval_24mov/camkII_mpfc/2023_09_15/16_42_22/' ],...
% [toppath 'miniscope_v4/Mov_Grid/G_baseline_6mov/camkII_mpfc/2023_09_25/17_08_13/' ],...
% [toppath 'miniscope_v4/Mov_Grid/G_shocks_12mov/camkII_mpfc/2023_09_26/17_00_14/' ],...
% [toppath 'miniscope_v4/Mov_Grid/G_retrieval_24mov/camkII_mpfc/2023_09_28/16_33_40/' ]};

if false
toppath = 'G:/.shortcut-targets-by-id/1BLCj3dIx9Q2KP6SQL4j-8Vo5CNepB8fz/moving_grid_no_shock/';
msDir =[toppath 'ms_data/'];
exp = 'Mov_Grid';
dirs = {[toppath 'NA_hab/camkII_mpfc/2025_03_10/13_10_14/'],...
[toppath 'NA_shock/camkII_mpfc/2025_03_11/14_43_12/'],...
[toppath 'NA_shock/camkII_mpfc/2025_03_11/14_58_07/'],...
[toppath 'NA_retrieval/camkII_mpfc/2025_03_13/11_25_13/'],...
[toppath 'NA_retrieval_day2/camkII_mpfc/2025_03_14/15_07_29/'],...
[toppath 'NA_retrieval_day2b/camkII_mpfc/2025_03_14/15_30_39/'],...
[toppath 'NB_hab/camkII_mpfc/2025_03_10/13_28_04/'],...
[toppath 'NB_shock/camkII_mpfc/2025_03_11/15_20_00/'],...
[toppath 'NB_retrieval/camkII_mpfc/2025_03_13/11_54_56/'],...
[toppath 'NB_retrieval_day2/camkII_mpfc/2025_03_14/15_59_17/']};
else
    toppath = 'H:/My Drive/';
msDir =[toppath 'ms_data/'];
exp = 'noshock';
dirs = {   'H:/My Drive/1_noshock_mpfc/camkII_mpfc/2025_07_28/17_38_07/';...
            'H:/My Drive/1_noshock_mpfc_acquisition/camkII_mpfc/2025_07_29/14_31_41/';...
            'H:/My Drive/1_noshock_mpfc_retrieval/camkII_mpfc/2025_07_31/11_18_12/';...
            'H:/My Drive/2_noshock_mpfc/camkII_mpfc/2025_07_28/18_01_19/';...
            'H:/My Drive/2_noshock_mpfc_acquisition/camkII_mpfc/2025_07_29/14_49_02/';...
            'H:/My Drive/2_noshock_mpfc_retrieval/camkII_mpfc/2025_07_31/11_42_52/';...
            'H:/My Drive/3_noshock_mpfc/camkII_mpfc/2025_07_28/18_19_24/';...
            'H:/My Drive/3_noshock_mpfc_acquisition/camkII_mpfc/2025_07_29/15_06_52/';...
            'H:/My Drive/3_noshock_mpfc_retrieval/camkII_mpfc/2025_07_31/12_04_31/';...
            'H:/My Drive/4_noshock_mpfc/camkII_mpfc/2025_07_28/18_36_22/';...
            'H:/My Drive/4_noshock_mpfc_acquisition/camkII_mpfc/2025_07_29/15_22_01/';...
            'H:/My Drive/4_noshock_mpfc_retrieval/camkII_mpfc/2025_07_31/12_26_12/';...
            'H:/My Drive/5_noshock_mpfc/camkII_mpfc/2025_07_28/18_53_25/';...
            'H:/My Drive/5_noshock_mpfc_acquisition/camkII_mpfc/2025_07_29/15_38_40/';...
            'H:/My Drive/5_noshock_mpfc_retrieval/camkII_mpfc/2025_07_31/12_51_08/';...
            'H:/My Drive/6_noshock_mpfc/camkII_mpfc/2025_07_28/19_16_39/';...
            'H:/My Drive/6_noshock_mpfc_acquisition/camkII_mpfc/2025_07_29/15_59_33/';...
            'H:/My Drive/6_noshock_mpfc_retrieval/camkII_mpfc/2025_07_31/13_15_02/';...
            'H:/My Drive/7_noshock_mpfc/camkII_mpfc/2025_07_28/19_35_56/';...
            'H:/My Drive/7_noshock_mpfc_acquisition/camkII_mpfc/2025_07_29/16_15_23/';...
            'H:/My Drive/7_noshock_mpfc_retrieval/camkII_mpfc/2025_07_31/13_39_40/';...
            'H:/My Drive/NA_noshock_mpfc/camkII_mpfc/2025_03_10/13_10_14/';...
            'H:/My Drive/NA_noshock_mpfc_acquisition/camkII_mpfc/2025_03_11/14_58_07/';...
            'H:/My Drive/NA_noshock_mpfc_retrieval/camkII_mpfc/2025_03_13/11_25_13/';...
            'H:/My Drive/NB_noshock_mpfc/camkII_mpfc/2025_03_10/13_28_04/';...
            'H:/My Drive/NB_noshock_mpfc_acquisition/camkII_mpfc/2025_03_11/15_20_00/';...
            'H:/My Drive/NB_noshock_mpfc_retrieval/camkII_mpfc/2025_03_13/11_54_56/'};
    
savedir = 'F:/Fernando/noshock_mpfc/';

end

behavcam = 'My_WebCam/';
msCam = 'My_V4_Miniscope/';

% extractbehav = true;
% if  extractbehav==true % extract behavior
% for i = 1%:length(dirs)
%     saveName = [dirs{i} behavcam 'extracted_pos.mat'];
%     if isfile(saveName) == false
%         fprintf('Run %s\n', saveName)
%         extract_mouse_pos_MovGrid(dirs{i}, saveName);
%     else
%         fprintf('~ Skip %s\n', saveName)
%     end
% end
% end
%%

extractbehav = true;
params = cell(length(dirs), 1);
parpool(8);
if  extractbehav==true % extract behavior
for i = 1:length(dirs)
    saveName = [dirs{i} behavcam 'extracted_pos.mat'];
    if isfile(saveName) == false
        vid = VideoReader( [[dirs{i} behavcam] '0.avi']);
        im = uint8(vid.read(1));
        [cropX, cropY, ~] = fernando_helper_click_behavROI(im);
        [lastX, lastY, lastfnum] = fernando_helper_click_mini(im, cropX, cropY, 1);
        params{i} = [];
        params{i}.cropX = cropX;
        params{i}.cropY = cropY;
        params{i}.lastX = lastX;
        params{i}.lastY = lastY;
        params{i}.lastfnum = lastfnum;
    else
%         fprintf('~ Skip %s\n', saveName)
    end
end
for i = 13%1:length(dirs)
    saveName = [dirs{i} behavcam 'extracted_pos.mat'];
    if true % isfile(saveName) == false
        fprintf('Run %s\n', saveName)
        extract_mouse_pos_MovGrid([dirs{i} behavcam], [], false, params{i});
    else
        fprintf('~ Skip %s\n', saveName)
    end
end
end

%%
altDir = [toppath 'miniscope_v4\Mov_Grid\ms_data_decode\'];
seps = '/';
analysis_ver = 1.0;
analyze_too = false;

timebins = 5; % for seconds pre/post PETH
% xbins = [0:.05:1];
% ybins = [0:.06:.18]; % should be 3 y bins

% im = imread("G:\.shortcut-targets-by-id\1mQqRDGRSlffLsw4MnHKbIatTEeGRdoaY\miniscope_v4\Mov_Grid\D_baseline_6mov\camkII_mpfc\2023_09_25\16_11_42\My_WebCam\behavCam.tiff", 1);
% cmperpix = 19.5/(165-27);
cmperpix = 17/(135-29);

for sessLoop = 1:length(dirs)
    behavName = [dirs{sessLoop} behavcam 'extracted_pos.mat'];
    s = strfind(behavName, seps);
    sessName = behavName(s(5)+1:s(6)-1);
    animal = sessName(1:2);
    s = strfind(sessName, '_');
    if contains(sessName, 'hab')
        nmoves = 6;
    elseif contains(sessName, 'shocks')
        nmoves = 12;
    elseif contains(sessName, 'retriev')
        nmoves = 24;
    end
%     nmoves = str2double( sessName(s(2)+1:end-3) );
    fprintf('%s\n', sessName)

    msdataName = [msDir sessName '.mat'];
    exptFileName = [dirs{sessLoop} msCam sessName '_analysis.mat'];
    alt_FileName = [altDir sessName '_analysis.mat'];
    checker = isfile(behavName) & isfile(msdataName) & (~isfile(exptFileName) | ~isfile(alt_FileName));
    msName = [dirs{sessLoop} msCam 'extracted_pos.mat'];
    if checker == true % run analysis
        fprintf('\t>Run %s\n', sessName)
        tempbehav   = load(behavName);
        tempms        = load(msdataName);
        cmn_YrA = tempms.ms.neuron.YrA;
        cmn_C   = tempms.ms.neuron.C;
        ms_time = tempms.ms.timestamps;
        ms_dt    = tempms.ms.dt;
        frameNum = tempms.ms.frameNum;


        barx = tempbehav.barpos(tempbehav.beh.frameNum+1); % was zero indexed
        behav_dt = tempbehav.beh.dt(tempbehav.beh.frameNum+1);
        behav_time = tempbehav.beh.timestamps;
        mousex = tempbehav.mouse.x;
        mousey = tempbehav.mouse.y;
        ledx = tempbehav.led.x;
        ledy = tempbehav.led.y;
        dx = ledx-mousex;
        % dy = ledy - mousey;
        dy = ledy-mousey;
        % dists = sqrt( (dx.^2) + (dy.^2) );
        ang =  rad2deg( atan2(dy, dx) );
        if max(mousex)<max(mousey)
            warning('hmmm')
            return
        end
        behav_scalar = 1/cmperpix; % max([mousex;mousey]);
        mousex = mousex./behav_scalar;
        mousey = mousey./behav_scalar;
        ledx = ledx./behav_scalar;
        ledy = ledy./behav_scalar;
        barx = barx./behav_scalar;

        sbarx = movmean(barx, 5);

        barspd = [0; abs(diff(sbarx))]./median(behav_dt);
        
        ms_mx  = interp1(behav_time, mousex, ms_time, 'linear', 'extrap');
        ms_my  = interp1(behav_time, mousey, ms_time, 'linear', 'extrap');
        ms_ledx  = interp1(behav_time, ledx, ms_time, 'linear', 'extrap');
        ms_ledy  = interp1(behav_time, ledy, ms_time, 'linear', 'extrap');
        ms_ang = interp1(behav_time, ang,    ms_time, 'linear', 'extrap');
        ms_barx= interp1(behav_time, barx,   ms_time, 'linear', 'extrap');
        ms_barspd= interp1(behav_time, barspd,   ms_time, 'linear', 'extrap');
        % ms_bardist = abs(ms_mx - ms_barx);
        %%
        

        % [xmap, ~, xbin] = histcounts(ms_mx, xbins);
        % [ymap, ~, ybin] = histcounts(ms_mx, ybins);
        % sbarx = movmean(ms_barx, 5);
        % bar_spd  = [0; abs(diff(sbarx))];
        bar_spd = ms_barspd;
        bar_spd = (bar_spd./max(bar_spd)) ;
        bar_still = bar_spd<.2;
        bar_slow  = bar_spd>=.2;
        bar_fast = bar_spd>=.7;
        barstart = find(bar_still(1:end-1)==1 & bar_still(2:end)==0);
        barstop  = find(bar_still(1:end-1)==0 & bar_still(2:end)==1);
        if length(barstart) - length(barstop) == 1;
            warning('incomplete start-stop')
            barstart = barstart(1:end-1)

            nmoves_a = length(barstart);
        end



        if (length(barstart)~=nmoves && length(barstop)~=nmoves) || any(barstop < barstart)
            warning('start stop mismatch')
        end

        t = ms_time./1000;
        start_id  = zeros(1, length(ms_barx));
        stop_id   = zeros(1, length(ms_barx));
        nb = floor(timebins/median(ms_dt)); % 5 sec win should be 98
        for i = 1:nmoves_a
            % can use time but has small chance to be slightly different length
            % t1 = t(barstart(i))  - timebins; 
            % ind1 = find(min(abs(t1-t)) == (abs(t1 - t)));
            % t1 = t(barstart(i))  + timebins; 
            % ind2 = find(min(abs(t1-t)) == (abs(t1 - t)));

            % just use set bins based on sampling rate
            ind1 = barstart(i) - nb;
            ind2 = barstart(i) + nb;
            start_id(ind1:ind2) = i;

            ind1 = barstop(i) - nb;
            ind2 = barstop(i) + nb;
            stop_id(ind1:ind2)  = i;
        end

        %%
        craw = cmn_C + cmn_YrA;
        % craw = zscore(craw, 1, 2);
        spks = tempms.ms.spks;
        [nsegs, nframes] = size(spks);
        
        sessionDir = dirs{sessLoop};
        saveVars = {'spks','craw','ms_time', 'ms_dt', 'frameNum', 'ms_mx', 'ms_my', 'ms_ang', 'ms_barx',...
            'ms_ledx', 'ms_ledy', 'ms_barspd', 'barstart', 'barstop', 'start_id', 'stop_id', 'nmoves',...
            'sessName', 'animal','sessionDir'};
        fprintf('Saving: %s\n\t-%s\n\t-%s\n\n', sessName, exptFileName, alt_FileName)
        save(exptFileName, saveVars{:})
        save(alt_FileName, saveVars{:})
        %%
        figure;
        set(gcf, 'Name', sessName)
        subplot_tight(3,2,1, [.05 .05])
        plot(ms_mx, ms_my, 'k')
        axis([-10 100 0 18])
        subplot_tight(3,2,2, [.05 .05])
        degsz = 15;
        polarhistogram(deg2rad(ms_ang), deg2rad([-degsz/2:degsz:360-degsz/2]));
        % rose(tempms.ms.ori.yaw, [-pi:pi/12:pi])
        subplot_tight(3,2,5:6, [.05 .05])
        imagesc(spks>0)
        subplot_tight(3,2,3:4, [.05 .05]); hold on
        plot(ms_time./1000, ms_barx, 'm')
        plot(ms_time./1000, ms_mx, 'k')
        plot(ms_time./1000, ms_ledx, 'r')
        ylim([-10 100])

        drawnow()
        %%
        % spks = zscore(spks, 1, 2);
        % cdff = diff(craw, 1, 2);
        % cdff= cat(2, zeros(nsegs,1), cdff);
        %%%%%%%%%%%%% Placefields to check
        % pfields = NaN(nsegs, length(xbins)-1);
        % xocc = NaN(1, length(xbins)-1);
        % for j = 1:length(xbins)-1
        %     thisbin = xbin==j;
        %     xocc(j) = sum(abs(ms_dt(thisbin)));
        % 
        %     pfields(:, j) = sum(craw(:, thisbin), 2) ./ xocc(j);
        % 
        % end
        % [peakrate, ind] = max(pfields, [], 2);
        % pfields = normalize_rows(pfields);
        % 
        % [~, peakord] = sort(ind, 'ascend');
        % pfields = pfields(peakord,:);
        %%%%%%%%%%%%% PETH of grid movement
        if analyze_too == true
        %%
        spks = tempms.ms.spks;% cmn_YrA+cmn_C;
        % spks = cmn_YrA+cmn_C;
        % spks = tempms.ms.neuron.S_matw;% cmn_YrA+cmn_C;
        spks = normalize_rows(spks);
        grid_start_peth = NaN(nsegs, nb*2+1, nmoves);
        grid_stop_peth  = NaN(nsegs, nb*2+1, nmoves);
        temp = [];
        for i = 1:nmoves
            % i
            thisbin = start_id==i;
            xocc = ms_dt(thisbin)';
            occmat = 1;%ones(nsegs,1)*xocc;
            fr = spks(:, thisbin) ./ occmat;
            grid_start_peth(:, :, i) = fr;
            temp = cat(2, temp, fr);

            thisbin = stop_id==i;
            xocc = 1;%ms_dt(thisbin)';
            occmat = ones(nsegs,1)*xocc;
            fr = spks(:, thisbin) ./ occmat;
            grid_stop_peth(:, :, i)  = fr;
        end
        % grid_move_peth = mean(cat(3, grid_start_peth, grid_stop_peth), 3, "omitmissing");
        grid_move_peth = mean(grid_start_peth, 3, "omitmissing");
        % grid_move_peth = mean(grid_stop_peth, 3, "omitmissing");
        [peakrate_peth, ind] = max(grid_move_peth, [], 2);
        grid_move_peth = normalize_rows(grid_move_peth);
        [~, peakord] = sort(ind, 'ascend');
        grid_move_peth = grid_move_peth(peakord,:);
        temp = temp(peakord,:);
        clf
        imagesc(grid_move_peth)
        hold on; plot([nb+1, nb+1], [.5 nsegs+.5], 'w--', 'LineWidth', 2)
        end
    else
        fprintf('~ Skip %s\n', sessName)
    end
end



%%
































