function [wrappedPhase, centerX, centerY, widthRec, heightRec] = processFourier(hologram, selectRegion_HCN, maxIntensity)
% processFourier: thực hiện các phép biến đổi Fourier và chọn vùng bậc +1

if nargin < 2 || isempty(maxIntensity) || isempty(selectRegion_HCN)
    selectRegion_HCN = true;
    maxIntensity = 100000;
end

%% Chuyển ảnh hologram sang grayscale
hologramGray = myConvGrayScale(hologram);  % Chuyển đổi ảnh sang grayscale

[numRows, numCols] = size(hologramGray);  % Lấy kích thước ảnh (chiều cao và chiều rộng)

%% Biến đổi Fourier và chọn ROI
fourierTransform = fftshift(fft2(hologramGray));  % Thực hiện biến đổi Fourier và dịch chuyển zero-frequency về trung tâm

% Hiển thị ảnh đã biến đổi Fourier
figure;
imshow(abs(fourierTransform), [1, maxIntensity]);
title('Biến đổi Fourier của Hologram');

%% Xác định +1 của phổ
if selectRegion_HCN
    [centerX, centerY, widthRec, heightRec] = myDrawRec();
    % Tính tọa độ góc trên bên trái của hình chữ nhật
    xRec = centerX - round(widthRec / 2);
    yRec = centerY - round(heightRec / 2);
    % Trích xuất nội dung ROI (bậc 1)
    roiContent = fourierTransform(yRec:yRec + heightRec - 1, xRec:xRec + widthRec - 1);
else
    error('Hàm này chỉ hỗ trợ định dạng hình chữ nhật (HCN).');
end

newFourierTransform = zeros(size(fourierTransform));

if selectRegion_HCN
    % Tính toán vị trí trung tâm mới
    newCenterX = round((numRows - heightRec) / 2);  % Tọa độ x của tâm mới
    newCenterY = round((numCols - widthRec) / 2); % Tọa độ y của tâm mới

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
figure();
mesh(wrappedPhase);
title('Pha chưa Unwrapping');

end

%% hàm vẽ hình chữ nhật
function [centerX, centerY, widthRec, heightRec] = myDrawRec()
% Vẽ hình chữ nhật để chọn ROI
roi = drawrectangle();

% Lấy tọa độ tâm ban đầu
centerRec = round([roi.Position(1) + roi.Position(3)/2, roi.Position(2) + roi.Position(4)/2]);

% Vẽ dấu cộng tại tâm
hold on;
centerMarker = plot(centerRec(1), centerRec(2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
hold off;

% Thiết lập callback để cập nhật tâm trong quá trình di chuyển hình chữ nhật
addlistener(roi, 'MovingROI', @(src, evt) updateCenterRectangle(src, centerMarker));

wait(roi);  % Double-click to confirm the ROI

% Get the position of the rectangle [x, y, width, height]
pos = round(roi.Position);
xRec = pos(1);
yRec = pos(2);
widthRec = pos(3);
heightRec = pos(4);

% Tính tọa độ tâm
centerX = xRec + widthRec / 2;
centerY = yRec + heightRec / 2;
end

% Hàm cập nhật tâm khi di chuyển ROI
function updateCenterRectangle(roi, centerMarker)
% Cập nhật vị trí của dấu cộng theo tọa độ tâm mới
centerX = roi.Position(1) + roi.Position(3)/2;
centerY = roi.Position(2) + roi.Position(4)/2;
centerMarker.XData = centerX;
centerMarker.YData = centerY;
drawnow;
end

%% hàm vẽ đường tròn
function [centerCir_X, centerCir_Y, radiusCir] = myDrawCir()
% Vẽ vòng tròn để chọn ROI
roi = drawcircle();

% Lấy tọa độ ban đầu của tâm và bán kính
centerCir = round(roi.Center);

% Vẽ dấu cộng tại tâm
hold on;
centerMarker = plot(centerCir(1), centerCir(2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
hold off;

% Thiết lập callback để cập nhật tâm trong lúc di chuyển vòng tròn
addlistener(roi, 'MovingROI', @(src, evt) updateCenter(src, centerMarker));

wait(roi);
centerCir_X = round(roi.Center(1));  % Tọa độ x của trung tâm hình tròn
centerCir_Y = round(roi.Center(2));  % Tọa độ y của trung tâm hình tròn
radiusCir = round(roi.Radius);
end

function updateCenter(roi, centerMarker)
% Cập nhật vị trí của dấu cộng theo tọa độ tâm mới
centerMarker.XData = roi.Center(1);
centerMarker.YData = roi.Center(2);
drawnow;
end

%% Hàm chuyển sang grayscale
function output = myConvGrayScale(inputImage)
    if size(inputImage, 3) > 1
        inputImage = rgb2gray(inputImage);
    end
    output = double(inputImage); 
end