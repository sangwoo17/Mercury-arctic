% === 옵션 늘린 Nelder-Mead 피팅 예제 ===
tbl = readtable('your_data.xlsx');    % 1열 t, 2열 c(t)
t = tbl{:,1};  c = tbl{:,2};

p0 = [10 1 5 280];                    % 초기 추정치

% ---------- fminsearch 옵션 ----------
opts = optimset('fminsearch');        % 기본 템플릿
opts.Display      = 'iter';           % 진행 상황 보고 (필요 없으면 'off')
opts.MaxFunEvals  = 5e4;              % 함수 평가 한도 ↑
opts.MaxIter      = 5e4;              % 반복 한도 ↑
opts.TolX         = 1e-8;             % 파라미터 수렴 허용오차
opts.TolFun       = 1e-8;             % 함수값 수렴 허용오차

% ---------- 피팅 실행 ----------
[p_est,fval,exitflag,output] = ...
    fminsearch(@(p) ssq(p,t,c), p0, opts);

% ---------- 결과 ----------
fprintf('Exitflag %d | SSE %.6g | Evals %d\n', ...
        exitflag, fval, output.funcCount);
disp(p_est)                              % [k1 k2 k3 A0]

c_fit = predict_c(p_est,t);
plot(t,c,'o', t,c_fit,'-','LineWidth',1.5);
legend('Data','Fit'); xlabel('Time'); ylabel('c(t)');

% -------- 내부 함수들 --------
function sse = ssq(p,t,c)
    if any(p<=0) || any(~isfinite(p)), sse = Inf; return; end
    y = predict_c(p,t); 
    if any(~isfinite(y)), sse = Inf; return; end
    sse = sum((c-y).^2);
end

function y = predict_c(p,t)          % 모델 (양수 파라미터 가정)
    k1=p(1); k2=p(2); k3=p(3); A0=p(4);
    sum_r  = k1+k2+k3;  prod_r = k1*k3;
    disc   = sum_r^2-4*prod_r;  if disc<=0, y=NaN(size(t)); return; end
    r1 = (sum_r+sqrt(disc))/2;  r2 = (sum_r-sqrt(disc))/2;
    y  = A0.*(1 + (r2/(r1-r2))*exp(-r1*t) + (r1/(r2-r1))*exp(-r2*t));
end
