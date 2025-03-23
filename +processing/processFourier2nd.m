function wrappedPhase = processFourier2nd(hologram, centerX, centerY, widthRec, heightRec, selectRegion_HCN, maxIntensity)
% processFourier: thực hiện các phép biến đổi Fourier và chọn vùng bậc +1

if nargin < 6 || isempty(maxIntensity) || isempty(selectRegion_HCN)
    selectRegion_HCN = true;
    maxIntensity = 100000;
end

%% Chuyển ảnh hologram sang grayscale
hologramGray = myConvGrayScale(hologram);  % Chuyển đổi ảnh sang grayscale

[numRows, numCols] = size(hologramGray);  % Lấy kích thước ảnh (chiều cao và chiều rộng)

%% Biến đổi Fourier và chọn ROI
fourierTransform = fftshift(fft2(hologramGray));  % Thực hiện biến đổi Fourier và dịch chuyển zero-frequency về trung tâm

% % Hiển thị ảnh đã biến đổi Fourier
% figure;
% imshow(abs(fourierTransform), [1, maxIntensity]);
% title('Biến đổi Fourier của Hologram');

%% Xác định +1 của phổ
if selectRegion_HCN
    % Trích xuất nội dung ROI (bậc 1)
    xRec = centerX - round(widthRec / 2);
    yRec = centerY - round(heightRec / 2);
    roiContent = fourierTransform(yRec:yRec + heightRec - 1, xRec:xRec + widthRec - 1);
else
    error('Hàm này chỉ hỗ trợ định dạng hình chữ nhật (HCN).');
end

newFourierTransform = zeros(size(fourierTransform));

if selectRegion_HCN
    % Tính toán vị trí trung tâm mới
    newCenterX = round((numRows - heightRec) / 2);  % Tọa độ x của tâm
    newCenterY = round((numCols - widthRec) / 2); % Tọa độ y của tâm

    newFourierTransform(newCenterX:newCenterX + heightRec - 1, newCenterY:newCenterY + widthRec - 1) = ...
        roiContent;
end

finalPhase = ifft2(ifftshift(newFourierTransform));
wrappedPhase = angle(finalPhase);

% Hiển thị newFourierTransform
figure;
imshow(abs(newFourierTransform), [1, maxIntensity]);
title('Ảnh sau khi dịch vào tâm');

wrappedPhase = wrappedPhase';
% figure();
% mesh(wrappedPhase);
% title('Pha chưa Unwrapping');

end

%% Hàm chuyển sang grayscale
function output = myConvGrayScale(inputImage)
    if size(inputImage, 3) > 1
        inputImage = rgb2gray(inputImage);
    end
    output = double(inputImage); 
end