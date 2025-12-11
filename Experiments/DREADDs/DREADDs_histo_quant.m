%%
clear 
ddir = "\\sshfs.r\garrettb@monk.cns.nyu.edu\f\fentonlab\RAWDATA\CaImage\GarrettBlair\GJB Histo\DREADDs rats slide scan\";
names = {'24508' '24509' '24520' '24521'};
slides = {  'slide1_20250104' 'slide2_20250104' 'slide3_20250104' 'slide4_20250104' 'slide5_20250104'...
            'slide1_20250105' 'slide2_20250105' 'slide3_20250105' 'slide4_20250105' 'slide5_20250105'};

max_images = 12;

laststr = {'.tif' '.btf'};
% laststr = {'.thumbnail.jpg'};
nr=4; nc=6;
fnames  = cell(length(names), nr*nc);
badfns  = cell(length(names), nr*nc);
imnames = cell(length(names), nr*nc);
images  = cell(length(names), nr*nc);
images2  = cell(length(names), nr*nc);
num_idx  = zeros(length(names),1);
%
for nameloop = 1:length(names)
    idx = 0;%num_idx(nameloop);
    for slideloop = 1:length(slides)
        for imloop = 1:max_images
            for lastloop = 1:length(laststr)
                fn = sprintf('%s%s\\%s\\_%02d%s', ddir, names{nameloop}, slides{slideloop}, ...
                    imloop, laststr{lastloop});
                if isfile(fn)
                    idx = idx+1;
                    imname = sprintf('%s_%s_%02d%s', names{nameloop}, slides{slideloop}(1:6), ...
                        imloop, laststr{lastloop});
                    
                    fnames{ nameloop, idx} = fn;
                    imnames{nameloop, idx} = imname;
                    fprintf('%s\n', imname)
                    try
                        im = imread(fn, 1);
                        im2 = imread(fn, 2);
                        images{nameloop, idx} = imresize(im, .2);
                        images2{nameloop, idx} = imresize(im2, .2);
                    catch
                        warning(sprintf('\nUnable to read %s\n', fn))
                        badfns{ nameloop, idx} = fn;
                    end % TRY
                end % IF is file
            end % FOR file endings .tif .btf
        end % FOR image number loop
    end % FOR slide loop
    num_idx(nameloop) = idx;
end % FOR animal loop
%%
mv = .1;
for nameloop = 1:length(names)
    idx=0;
    figure(nameloop); clf; colormap magma
for imloop = 1:size(images,2)
    if ~isempty(images{nameloop, imloop})
        idx = idx+1;
        subplot_tight(nr,nc, idx, [.01 .01])
        name = imnames{nameloop, idx};
        im   = images{nameloop, idx}*10;
        im2   = images2{nameloop, idx}*0;
        im_color = cat(3, im, im*0, im2);
%         im_color(im_color>mv) = mv;
%         im_color = im_color./mv; 
        image(im_color)
        axis image off
%         imagesc(im2)
%         fprintf('\n%s', imname)
    end
end    
end
%%
for nameloop = 1:length(names)
    idx=0;
    figure(nameloop+100); clf; colormap magma
for imloop = 1:size(images,2)
    if ~isempty(images{nameloop, imloop})
        idx = idx+1;
        subplot(nr,nc, idx)
        name = imnames{nameloop, idx};
        im   = images{nameloop, idx}*0;
        im2   = images2{nameloop, idx}*20;
        im_color = cat(3, im, im*0, im2);
%         im_color(im_color>mv) = mv;
%         im_color = im_color./mv; 
        image(im_color)
%         imagesc(im2)
%         fprintf('\n%s', imname)
    end
end    
end