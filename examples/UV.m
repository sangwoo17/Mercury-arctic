% UV intensity

% === 1. 파일 경로 설정 ===
filename = 'UV.he5'; % 여기에 파일명 입력

% === 2. HDF5 구조 확인 ===
% 한번 구조 확인 (변수 이름 확인용)
h5disp(filename);

% === 3. 필요한 변수 읽기 ===
% 예시: 'UVB_Irradiance_305nm' 라고 가정 (구체적 이름은 h5disp로 확인 필요)
% 위도, 경도 변수도 확인해서 읽기

uv_data = h5read(filename, '/HDFEOS/GRIDS/OMI_UV_Grid/Data Fields/UVB_Irradiance_305nm');
latitude = h5read(filename, '/HDFEOS/GRIDS/OMI_UV_Grid/Data Fields/Latitude');
longitude = h5read(filename, '/HDFEOS/GRIDS/OMI_UV_Grid/Data Fields/Longitude');

% === 4. 데이터 전처리 (필요시 NaN 처리) ===
uv_data(uv_data < 0) = NaN; % 음수 값 제거 (센서 오류 값 보정용)

% === 5. 10도 bin 설정 ===
lat_bins = -90:10:90; % bin 엣지 설정 (-90~90도까지)

% === 6. 각 bin 별 평균 계산 ===
uv_mean_per_bin = zeros(length(lat_bins)-1, 1);

for i = 1:length(lat_bins)-1
    % 현재 bin 범위
    lat_min = lat_bins(i);
    lat_max = lat_bins(i+1);
    
    % 해당 위도 범위에 해당하는 인덱스 찾기
    lat_idx = (latitude >= lat_min) & (latitude < lat_max);
    
    % 해당 bin에서 UV 데이터 추출
    uv_in_bin = uv_data(lat_idx);
    
    % 평균 계산 (NaN 무시)
    uv_mean_per_bin(i) = mean(uv_in_bin, 'omitnan');
end

% === 7. 결과 출력 ===
fprintf('Latitude Bin\tMean UV Intensity (W/m²)\n');
for i = 1:length(uv_mean_per_bin)
    fprintf('%d to %d\t\t%.3f\n', lat_bins(i), lat_bins(i+1), uv_mean_per_bin(i));
end

% === 8. Optional: 그래프로 보기 ===
figure;
bar(lat_bins(1:end-1)+5, uv_mean_per_bin); % 중앙값에 plot
xlabel('Latitude (deg)');
ylabel('Mean UV Intensity at 305nm (W/m²)');
title('Latitude-binned Mean UV Intensity at 305nm');
grid on;