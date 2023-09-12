function [Y] = load_tiffstack(fileName, inds)
% fileName = "E:\RecordingData\GarrettBlair\PKCZ_imaging\mHPC23454\2023_08_12\16_55_19_RET8\HPC_miniscope1\dff_small_close4_gaussblur4.tiff";
tiffInfo = imfinfo(fileName);  % Get the TIFF file information
num_frames = numel(tiffInfo);    % Get the number of images in the file
if nargin < 2
    inds = 1:num_frames;
else
    num_frames = length(inds);
end
Y = uint8(zeros(tiffInfo(1).Height, tiffInfo(1).Width, num_frames)); 
for iFrame = 1:length(inds)
  Y(:,:,iFrame) = double(imread(fileName,'Index',iFrame,'Info',tiffInfo));
end
