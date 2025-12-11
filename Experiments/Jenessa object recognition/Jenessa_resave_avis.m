clear
close all


% in_dir = 'H:\Jenessa\ObjRec_MiniPilot\Videos\';
% out_dir = 'H:\Jenessa\ObjRec_MiniPilot\Videos\videos_blank_removed\';

in_dir = 'I:\Jenessa\L2L 2025\ObjRec\Retention\';
out_dir = 'I:\Jenessa\L2L 2025\ObjRec\Retention\videos_blank_removed\';


% first load all of the tracker info files
txts = dir([in_dir '*_TrackerVideo_info.txt']);
frames2cut = 2;
% parameter for interactively drawing a mask (4 point polygon via gui)
% default - false
use_mask = false;

% whether to interpolate the bad frames
% default - true
interpolate_avis = 'nearest';

for tnum = 1:length(txts)
    fname = [in_dir txts(tnum).name];
    d = readtable(fname, 'Delimiter', '\t');
    aviname = d.Var1;
    f1 = d.Var2;
    f2 = d.Var3;
    
    save_name = [out_dir aviname{1}];
    if ~isfile(save_name)
        %%
        fprintf('RUNNING:\n%s \n\t%s...\n', out_dir, aviname{1})
            vidname = [in_dir aviname{1}];
            [Y, vt, badf] = load_tracker_avi(vidname, [], [], [], interpolate_avis);
            %%%%%%%%%%%%%%%%%%%%%%% excluding last frame because of errors
            %%%%%%%%%%%%%%%%%%%%%%% occasionally in sleap with bad last
            %%%%%%%%%%%%%%%%%%%%%%% frame
            Y    = Y(:,:,1:end-frames2cut);
            vt   = vt(1:end-frames2cut);
            badf = badf(1:end-frames2cut);
            
            fnum = 1:size(Y,3);
            T = table(fnum', vt', ~badf, 'VariableNames', {'FrameNumber', 'Time', 'isGood'});
            save_output(save_name, Y, T)            
            %%%%%%%%%%%%%%%%%%%%%%%
%         for vidnum = 1:length(aviname)
%             %%
%             vidname = [in_dir aviname{vidnum}];
%             [Y, vt, badf] = load_tracker_avi(vidname, [], [], [], interpolate_avis);
%             %%%%%%%%%%%%%%%%%%%%%%%
%             Y    = Y(:,:,1:end-1);
%             vt   = vt(1:end-1);
%             badf = badf(1:end-1);
%             %%%%%%%%%%%%%%%%%%%%%%%
%             transistions = vt>=(max(vt)*.95);
%             dt = median(abs(diff(vt)));
%             if vidnum==1
%                 im = squeeze(Y(:,:,1));
%                 if use_mask == true
%                     [xGrid, yGrid] = meshgrid(1:size(im,2), 1:size(im,1));  % x: cols, y: rows
%                     figure();
%                     imagesc(im);
%                     [px, py] = ginput(4);
%                     mask = ~inpolygon(xGrid, yGrid, px, py);
%                 else
%                     mask = im<-1;
%                 end
%                 Yall = Y;
%                 vid_timestamps = vt;
%                 badframes = badf';
%                 alltrans = transistions;
% 
%             else
%                 Yall = cat(3, Yall, Y);
%                 vid_timestamps = cat(2, vid_timestamps, vid_timestamps(end)+dt+vt);
%                 alltrans = cat(2, alltrans, transistions);
%                 badframes = cat(2, badframes, badf');
%             end
%         end
%         fnum = 1:size(Yall,3);
%         T = table(fnum', vid_timestamps', ~badframes', 'VariableNames', {'FrameNumber', 'Time', 'isGood'});
%         for i = 1:size(Yall,3)
%             im = squeeze(Yall(:,:,i));
%             im2 = im;
%             im2(mask) = 255;
%             Yall(:,:,i) = im2;
%         end
        %%
%         save_output(save_name, Yall, T)
    else
        fprintf('%s \n\t%s video found, skipping...\n', out_dir, aviname{1})
    end
end

% %%
% figure; colormap bone
% alltrans2 = conv(alltrans, ones(100,1), 'same')>0;
% clrs = cumsum(cumsum(diff(alltrans2)));
% clrs = [clrs(1), clrs];
% clrs = round(100*(clrs./max(clrs)));
% clrmap = jet(100);
% 
% fs = find(alltrans2);
% for i = 1:length(fs)
%     %%
%     clf
%     ind = fs(i);
%     im = squeeze(Yall(:,:,ind));
%     imagesc(im, [0 255]);
%     hold on
%     plot([0 0], [0 size(im,1)*(clrs(ind)/100)], 'LineWidth', 5, 'Color', clrmap(clrs(ind)+1,:))
%     axis image off
%     title(sprintf('t = %3.1f sec', vid_timestamps(ind)/1000))
%     drawnow();
% end


%%
function save_output(save_name, Y, vt)
fprintf('Done in %d sec\n', round(toc()))
tic
fprintf('2. Writing output to:\n\t\t%s...\n\t', save_name)
v = VideoWriter(save_name);
open(v)
for ii = 1:size(Y, 3)
    v.writeVideo(squeeze(Y(:,:,ii)));
end
close(v)
fprintf('Done in %d sec\n', round(toc()))
aviend = strfind(save_name, '.avi');
tsName = [save_name(1:aviend-1) '_timeStamps.csv'];
fprintf('3. Timestamps saved to:\n\t\t%s\n', tsName)
%     csvwrite(tsName, vt); % doesn't have enough precision
if istable(vt)
    writetable(vt, tsName)
else
    dlmwrite(tsName, vt, 'precision', '%.3f')
end

end