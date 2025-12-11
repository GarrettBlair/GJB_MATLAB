function [Y, vt, bads] = load_tracker_avi(fname, save_name, stopframe, room_tracking_fname, interpolate_bads)
fprintf('1. Loading video:\n\t\t%s...\n\t', fname)
v = VideoReader(fname);

if isempty(stopframe)
    stopframe=Inf;
end

if ~isempty(room_tracking_fname)
    [room , ~]  = read_APA_csv(room_tracking_fname, []);
    % [arena, arena_params] = read_APA_csv(arena_tracking_fname, []);
    vt = room.timestamps; %   NaN(nf,1);
    nf = length(vt);
else
%     warning('timestamps assumed')
    nf = ceil(v.Duration*v.FrameRate);
    dt = 1000/v.FrameRate;
    vt = dt*linspace(0,nf-1, nf);
end

Y = uint8(false(v.Height, v.Width, nf));
bads = true(nf,1);

idx = 0;
tic

% v.open
while v.hasFrame && idx<=stopframe% && v.CurrentTime<timetoget
    im = rgb2gray(v.readFrame);
    idx = idx+1;
    if any(im(:)>1) % 'empty frames will be sometimes banded 0's and 1's
        Y(:,:,idx) = im;
        bads(idx) = false;
    end
end
if nf~=idx && stopframe==Inf
    if abs(idx-nf)>1
        disp(idx-nf)
        warning('dimension mismatch!')
    end
    if nf<idx
        Y = Y(:,:,1:nf);
        bads = bads(1:nf);
        vt = vt(1:nf);
    else
        Y = Y(:,:,1:idx);
        bads = bads(1:idx);
        vt = vt(1:idx);
    end
end


nearest_interp = true;
if any(bads)
    if strcmp(interpolate_bads, 'none'); % false
        Y    = Y(:,:,bads==false);
        vt = vt(bads==false);
        bads = bads(bads==false);
    else
        badidx = find(bads);
        goodidx = find(~bads);
        for ii = 1:length(badidx)
            %%
            thisidx = badidx(ii);
            prevgood = goodidx(find(goodidx<badidx(ii), 1, 'last'));
            nextgood = goodidx(find(goodidx>badidx(ii), 1, 'first'));
            if isempty(prevgood)
                im1 = squeeze(Y(:, :, nextgood))*0;
                im2 = squeeze(Y(:, :, nextgood));
                im3 = squeeze(Y(:, :, nextgood));
            elseif isempty(prevgood)
                im1 = squeeze(Y(:, :, prevgood));
                im2 = squeeze(Y(:, :, prevgood))*0;
                im3 = squeeze(Y(:, :, prevgood));
            else
                im1 = Y(:, :, prevgood);
                im2 = Y(:, :, nextgood);
                w1 = thisidx-prevgood;
                w2 = nextgood-thisidx;
                if strcmp(interpolate_bads, 'nearest'); % 
                    if w1<=w2
                        im3 = im1;
                    else
                        im3 = im2;
                    end
                elseif strcmp(interpolate_bads, 'linear'); % 
                    totaldist = w1+w2;
                    if totaldist > 2
                        warning('check distances')
                        disp(thisidx)
                        disp(totaldist)
                        disp(save_name)
                    end
                    im1_w = double(im1).*(1-(w1/totaldist));
                    im2_w = double(im2).*(1-(w2/totaldist));
                    im3 = im1_w + im2_w;
                    im3(im3<0) = 0;
                    im3(im3>255) = 255;
                    im3 = uint8(im3);
                else
                    error('unkown interpolation')
                end
            end
            if false % ii<5
                figure(thisidx)
                title(thisidx)
                imagesc(cat(2, im1, im3, im2), [0 255])
                axis image
            end
            % burn in an id tag to the last 6 pixels to signify it was
            % replaced
            im3(end-5:end) = [0 255 0 255 0 255];
            Y(:, :, badidx(ii)) = im3;
%             bads(badidx(ii)) = false;
        end
    end
end



fprintf('Done in %d sec\n', round(toc()))
if ~isempty(save_name)
    tic
    fprintf('2. Writing output to:\n\t\t%s...\n\t', save_name)
    v = VideoWriter(save_name);
    open(v)
    for ii = 1:length(vt)
        v.writeVideo(squeeze(Y(:,:,ii)));
    end
    close(v)
    fprintf('Done in %d sec\n', round(toc()))
    aviend = strfind(save_name, '.avi');
    tsName = [save_name(1:aviend-1) '_timeStamps.csv'];
    fprintf('3. Timestamps saved to:\n\t\t%s\n', tsName)
    %     csvwrite(tsName, vt); % doesn't have enough precision
    dlmwrite(tsName, vt, 'precision', '%.3f')
end
end