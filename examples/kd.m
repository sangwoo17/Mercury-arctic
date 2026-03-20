% 파일명 설정
filename = 'Kd490.nc'; % 여기 파일명 수정하세요

% NetCDF 변수명 설정
lat = ncread(filename, 'lat'); % 또는 'latitude' → NetCDF 파일에 따라 이름 확인 필요
lon = ncread(filename, 'lon'); % 또는 'longitude'

% 원하는 변수명 설정 (예: 'Kd_490', 'Kd_325', 'Kd_380', 'chlor_a', 등)
var = ncread(filename, 'Kd_490');

% 확인: var의 dimension (lon x lat) 인지 (lat x lon) 인지 확인 필요
% 보통 MODIS/SeaDAS 제품은 (lon, lat) 구조 → 위에서는 lon 1st, lat 2nd로 가정

% Latitude bin 정의
lat_bins = 0:10:90; % → 0,10,20,...,90 → 총 9개 구간 (0-10, ..., 80-90)

% 결과 저장용 배열
mean_values = nan(1, length(lat_bins)-1);

% 루프 돌면서 각 구간별 평균 계산
for i = 1:length(lat_bins)-1
    % 현재 위도 구간
    lat_min = lat_bins(i);
    lat_max = lat_bins(i+1);
    
    % 해당 위도 index 찾기
    lat_idx = find(lat >= lat_min & lat < lat_max);
    
    % 예외 처리: 만약 해당 구간에 index가 없다면 넘어가기
    if isempty(lat_idx)
        mean_values(i) = NaN;
        continue;
    end
    
    % 변수에서 해당 lat_idx slice 추출
    % → 전체 longitude 에 대해 해당 latitude slice 사용
    % var 구조가 (lon x lat) 인 경우:
    var_slice = var(:, lat_idx);
    
    % NaN 제외하고 평균 계산
   mean_values(i) = mean(var_slice(:), 'omitnan');
    
    % 결과 출력 (optional)
    fprintf('Latitude %d-%d deg: mean = %.4f\n', lat_min, lat_max, mean_values(i));
end

% 최종 결과 표시
disp('Summary of mean values by latitude bin:')
for i = 1:length(mean_values)
    fprintf('%2d-%2d deg: %.4f\n', lat_bins(i), lat_bins(i+1), mean_values(i));
end
