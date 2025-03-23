%% start

%% Xóa dữ liệu trong workspace và đóng tất cả các figure
clear;
close all;
clc;

%% Additional file (thêm file ảnh cần thiết - nếu nằm ngoài thư mục source code)
filePath = 'C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Ảnh mẫu';
addpath(filePath);


%%
% Đường dẫn tới thư mục chứa các ảnh
image_folder = 'C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Ảnh mẫu';
[images, num_images] = processing.load_folder_images(image_folder);


%% Biến toàn cục (tham số đầu vào các hàm)
DPD = 25;
he_so = 1;
poly_order = 3;
maxIntensity =100000 ;  % Cường độ tối đa để hiển thị ảnh
%%
wrappedPhase = cell(1, num_images);
unwrapped_Phase = cell(1, num_images);

%% processFourier: thực hiện các phép biến đổi Fourier và chọn vùng bậc +1

%%

for i = 1:num_images
%% Chọn thuật toán unwrapping
if i ==1 
[wrappedPhase{1}, centerX, centerY, widthRec, heightRec] = processing.processFourier(images{1});
else
wrappedPhase{i} = processing.processFourier2nd(images{i}, centerX, centerY, widthRec, heightRec);
end
% Ví dụ:
%   unwrapped_Phase = unwrapPhase(wrappedPhase, 'ls', 'dct'); % LS với DCT
%   unwrapped_Phase = unwrapPhase(wrappedPhase, 'tie', 'fft'); % TIE với FFT
%   unwrapped_Phase = unwrapPhase(wrappedPhase, 'linh'); % Phương pháp của a Linh
%   unwrapped_Phase = unwrapPhase(wrappedPhase, '2dweight'); % 2D weighted phase unwrapping
%
% Nếu không nhập methodGroup và methodType, mặc định sử dụng 'ls' với 'dct'.

methodGroup = 'poisson';
methodType ='';
unwrapped_Phase{i} = unwrapping.unwrapPhase(wrappedPhase{i}, methodGroup);


end
figure;
surf(unwrapped_Phase{1});
title("anh 1");
figure;
surf(unwrapped_Phase{2});
title("anh 2");
%%
sigma = 0.5;
 % Khởi tạo danh sách để lưu trữ các ảnh giữ lại
    kept_images = cell(1, num_images);
    j = 1;
     % Vòng lặp để xử lý các ảnh
    for i = 1:num_images-1 
        % Tính toán sự khác biệt pha
        phase_diff = abs(unwrapped_Phase{i+1} - unwrapped_Phase{i});   

        % Kiểm tra giá trị sigma
        if sum(phase_diff(:) > sigma)
            kept_images{j} = unwrapped_Phase{i};
            kept_images{j+1} = unwrapped_Phase{i+1};
            j = j+1;
        end
    end

%%
wavelength = 633e-9;
reconSurface = (unwrapped_Phase .* wavelength .* he_so) / (4*pi);
reconSurface = reconSurface * 10^6;     

offSet = 10;      
reconSurface = reconSurface(offSet:end-offSet,offSet:end-offSet);  % cắt/chọn vùng để vẽ đồ thị

temp_r = reconSurface;
%% Chuyển đổi đơn vị (sang nanomet nếu cần)
[reconSurface, dimensional] = processing.postProcess.myConvertUnit(reconSurface);

%% Xác định chiều vân (ngang/dọc)
detectFringe = processing.postProcess.detectFringeSobel(reconSurface);
disp(detectFringe);
if strcmpi(detectFringe, 'vân ngang')
    reconSurface = rot90(reconSurface); % Xoay 90 độ theo chiều dương
end

%% Post Processing
imagesc(reconSurface);
hold on;    
title('Mặt phẳng pha'); % Đặt tiêu đề cho hình ảnh

% Vẽ đường thằng cắt ngang
positionLine = processing.postProcess.myDrawLine();  

crossLine = processing.postProcess.myCrossSection(reconSurface, positionLine);

% tinh đường trung bình mean_line
meanLine = processing.postProcess.myMeanLine(crossLine, poly_order);

%% Tính toán độ nhám 2D
[Ra, Ra_line] = roughness.myCalcRa(crossLine, meanLine, DPD);

Rz = roughness.myCalcRz(crossLine, meanLine);

%% Tính toán độ nhám 3D
Sz=  roughness.myCalcSz(reconSurface,poly_order);
% Độ nhám trung bình Sa
Sa = roughness.myCalcSa(reconSurface, poly_order);
% Độ nhám Sq
Sq = roughness.myCalcSq(reconSurface, poly_order);

%% Tổng hợp thông số độ nhám
surfaceParams = {reconSurface, Ra, Rz, Sa, Sq, Sz, Ra_line, positionLine,...
                                    crossLine, meanLine, dimensional, DPD}; 

%% Vẽ đồ thị 2D
visualization.display2D(surfaceParams);

%% Tạo figure 3D với đơn vị thích hợp
enableDisplayCrossSection3D = true;     %bật tắt mặt cắt ngang ảnh 3D

visualization.display3D(surfaceParams, true);


%% end