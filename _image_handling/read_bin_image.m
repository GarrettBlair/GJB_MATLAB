function [ms] = read_bin_image(pathname, sessname, spatial_ds, temporal_ds)
d1 = 480;% height
d2 = 752; % width
write_size = (d1*d2+8); % 8 bits reserved for beginning of frame [255 255 255 255] flag
search_size  = write_size*1;
image_data_start = 9; % imaging data starts on the 9th bit

filename = [pathname '\downsampled_msCam.tiff'];
cd(pathname)
ImageFile = fopen([pathname '\msImage.bin'], 'r');

fseek(ImageFile, write_size*-1, 'eof');
temppp = fread(ImageFile, 8, 'uint8');
totalFrames = sum(temppp(5:8)'.*[2^0 2^8 2^16 2^24]);
frameNum = NaN(totalFrames,1);
ms.spatialDownSample = spatial_ds;
ms.temporalDownSample = temporal_ds;

fTIF = Fast_Tiff_Write(filename);

ind = 0;
fseek(ImageFile, 0, 'bof');
for frame = 2:ms.temporalDownSample:totalFrames % first frame is blank it seems?
    ind = ind+1;
    start_ind = write_size*(frame-1);
%     bbb = aaa + search_size;
    fseek(ImageFile, start_ind, 'bof');
    imData = fread(ImageFile, search_size, 'uint8');

    frameVec = imData(5:8); % frame number written in bytes [5:8] (8-bit)
    frameNum(ind) = sum(frameVec'.*[2^0 2^8 2^16 2^24]);% + 1;

    temp = reshape(imData(image_data_start:end), d2, d1)';
    temp = uint8(imresize(temp, 1/ms.spatialDownSample));
    if mod(frame, floor(totalFrames*.1))==0
        fprintf('\nFrame %5.0d / %5.0d done...', frame, totalFrames);
    end
    fTIF.WriteIMG(temp'); % much faster method
%    imwrite(temp./254,  filename, 'tiff', 'writemode', 'append', 'compression', 'none'); % 255 is reserved for beginning of frame, so 254 is max value
end
fTIF.close;

frameNum = frameNum(~isnan(frameNum));

fclose(ImageFile);
ms.dirName = pathname;
ms.frameNum = frameNum;
ms.numFrames = length(ms.frameNum);
ms.vidNum = ones(ms.numFrames,1);
ms.height = 480/ms.spatialDownSample;
ms.width = 752/ms.spatialDownSample;
% ms.Experiment = "A16"
% ms.camNumber = 0
ms.time = (ms.frameNum-1)*(ms.temporalDownSample/30);
% ms.maxBufferUsed = 1;
ms.dateNum = NaN;
ms.columnCorrection = 0;
ms.columnCorrectionOffset = 0;


% ms.alignmentROI = [50; 50; ms.width-50; ms.height-50];
if ~isfield(ms, 'alignmentROI')
    numROIs=0;
    userInput = 'Y';
    refFrame = imread(filename, 1);
    imshow(uint8(refFrame))
    hold on
    while (strcmp(userInput,'Y'))
        numROIs = numROIs+1;
        display(['Select ROI #' num2str(numROIs)])
        rect = getrect();
        rect(3) = rect(3) - mod(rect(3),2);
        rect(4) = rect(4) - mod(rect(4),2);
        
        ms.alignmentROI(:,numROIs) = rect; %uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
        rectangle('Position',rect,'LineWidth',2);
        userInput = upper(input('Select another ROI? (Y/N)','s'));
    end
end

fprintf('\nRaw frames written!\nBeginning motion correction...\n');

plotting = true;
ms = msAlignmentFFT_parallel_tiff(ms, filename, plotting);
ms.selectedAlignment = 1;

ms = msMeanFrame_parallel_tiff(ms, filename, 5); % 5 - downsampling for average frame, min frame, and max frame

ms.tiff_filename = ['C:\Users\BlairLab-01\Documents\MATLAB\AnalysisData\A36_shock_data\' sessname '_alignedCam'];

msWriteMeanTiff_andNE_tiff(ms, filename, [1 ms.numFrames], ms.tiff_filename, 1, 1, 5, 30, false);
ms.tiff_filename = ['C:\Users\BlairLab-01\Documents\MATLAB\AnalysisData\A36_shock_data\' sessname '_alignedCam_raw.tiff'];
end
