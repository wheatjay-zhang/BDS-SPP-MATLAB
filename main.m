% BeiDou Single Point Positioning
% Author: WheatJay
% Version: 1.1
% Description: BDS SPP based on MATLAB

%% ===================== 0. 填入文件 =====================
obsFile = '';  % 在此填入观测文件 eg: ALIC00AUS_R_20240320000_01D_30S_MO.rnx
navFile = '';  % 在此填入导航文件 eg: BRDC00IGS_R_20240320000_01D_MN.rnx

%% ===================== 1. 观测数据 =====================
approx_coords = approx_position(obsFile);
obsData = rinexread(obsFile);
obs = timetable2table(obsData.BeiDou);
obs = obs(~isnan(obs.C2I) & ~isnan(obs.C6I), :);
obs = sortrows(obs, 'Time');
obs.SatelliteID = compose('C%02d', obs.SatelliteID);

%% ===================== 2. 导航数据 =====================
navData = rinexread(navFile);
nav = navData.BeiDou;

satList = unique(nav.SatelliteID);
satStr  = compose('C%02d', satList);
numSats = numel(satList);
t0 = dateshift(obs.Time(1), 'start', 'day');   % 当天00:00:00 第一个观测时刻 
t_in = (weekday(t0) - 1) * 86400;
step = 30;
%% ===================== 3. IF组合伪距 =====================
f1 = 1561.098; f3 = 1268.520;
a = f1^2/(f1^2 - f3^2);  b = f3^2/(f1^2 - f3^2);   c = 299792458;

P_IF = NaN(height(obs),1);

for i = 1:height(obs)
    satNum = str2double(obs.SatelliteID{i}(2:end));
    
    nav_sat = nav(nav.SatelliteID==satNum & nav.SatH1==0,:);
    if isempty(nav_sat), continue; end
    
    t_rx = seconds(obs.Time(i) - t0) + t_in;
    curDT = obs.Time(i);
    
    valid = nav_sat.Time <= curDT;    % 用 Time 匹配
    if ~any(valid), continue; end
    
    nav_valid = nav_sat(valid,:);
    dt = seconds(curDT - nav_valid.Time);
    
    [~,k ] = min(dt);
    tgd = nav_valid.TGD1(k);
    
    C2 = obs.C2I(i) - tgd*c;
    P_IF(i) = a*C2 - b*obs.C6I(i);  %tgd改正
end

obs.P_IF = P_IF;
obs = rmmissing(obs, 'DataVariables','P_IF');

%% ===================== 4. 构建历元矩阵 =====================
t_vec = t_in:step:t_in+86400-step;
numEpochs = numel(t_vec);
obsSec = seconds(obs.Time - t0) + t_in;
epochIdx = round((obsSec - t_in)/step) + 1;
Pmat = NaN(numSats, numEpochs);

for i = 1:height(obs)
    col = find(strcmp(obs.SatelliteID{i}, satStr));
    ep = epochIdx(i);
    
    if ~isempty(col) && ep>=1 && ep<=numEpochs
        Pmat(col, ep) = obs.P_IF(i);
    end
end

%% ===================== 5. 卫星坐标计算 =====================
rows = cell(0,7);
for ep = 1:numEpochs
    t_rx = t_vec(ep);
    curDT = t0 + seconds(t_rx - t_in);
    
    for s = 1:numSats
        satNum = satList(s);
        satID  = satStr{s};
        
        nav_sat = nav(nav.SatelliteID==satNum & nav.SatH1==0,:);
        if isempty(nav_sat), continue; end
        
        valid = nav_sat.Time <= curDT;        % Time匹配
        if ~any(valid), continue; end
        
        nav_valid = nav_sat(valid,:);
        dt = seconds(curDT - nav_valid.Time);
        
        [dt_min,k] = min(dt);
        if dt_min > 7200, continue; end        % 星历有效性限制
        
        e_raw = nav_valid(k,:);
        e = struct( ...
            'toe', e_raw.Toe, ...
            'af0', e_raw.SVClockBias, ...
            'af1', e_raw.SVClockDrift, ...
            'af2', e_raw.SVClockDriftRate, ...
            'sqrtA', e_raw.sqrtA, ...
            'Delta_n', e_raw.Delta_n, ...
            'M0', e_raw.M0, ...
            'e', e_raw.Eccentricity, ...
            'i0', e_raw.i0, ...
            'IDOT', e_raw.IDOT, ...
            'omega', e_raw.omega, ...
            'Cuc', e_raw.Cuc, ...
            'Cus', e_raw.Cus, ...
            'Crc', e_raw.Crc, ...
            'Crs', e_raw.Crs, ...
            'Cic', e_raw.Cic, ...
            'Cis', e_raw.Cis, ...
            'omega0', e_raw.OMEGA0, ...
            'omegaDot', e_raw.OMEGA_DOT ...
        );
        
        dt_toe = t_rx - e.toe;
        V_ts = e.af0 + e.af1*dt_toe + e.af2*dt_toe^2;
        P = Pmat(s,ep);
        if isnan(P), continue; end
        
        t_emit = t_rx - P/c - V_ts;
        if satNum <= 5 || satNum >= 59
            [X,Y,Z,dt_rel] = GEO_STATE_tsv(t_emit, e);
        else
            [X,Y,Z,dt_rel] = MEO_STATE_tsv(t_emit, e);
        end
        V_ts = V_ts + dt_rel;

        rows(end+1,:) = {satID, t_rx, curDT, X, Y, Z, V_ts}; %#ok
    end
end
XYZ_Table = cell2table(rows, 'VariableNames', {'SatelliteID','t','Time','X','Y','Z','V_ts'});
writetable(XYZ_Table,'AAAbeidou_coords.csv');

%% ===================== 6. 最小二乘 =====================
[results_table, removed_table] = least_squares(XYZ_Table, obs, approx_coords);
writetable(results_table,'AAAresults.csv');
writetable(removed_table,'AAAremoved_satellites.csv');

%% ===================== 7. 绘图 ==========================
t = results_table.times;
XYZ = results_table.final_coords;
dev = XYZ - approx_coords;

figure('Position', [100 100 1000 500]); 
plot(t, dev(:,1), 'r.-', 'DisplayName', 'X方向', 'LineWidth', 1.2);
hold on;
plot(t, dev(:,2), 'g.-', 'DisplayName', 'Y方向', 'LineWidth', 1.2);
plot(t, dev(:,3), 'b.-', 'DisplayName', 'Z方向', 'LineWidth', 1.2);
hold off;
grid on;xlabel('时间');ylabel('偏差 (m)');title('SPP定位结果与接收机近似坐标的差值');
legend('Location', 'best');

%% ===================== 8. 本地函数 =======================
function [Xk,Yk,Zk,dt_rel] = GEO_STATE_tsv(tsv,GEO_eph)
mu = 3.986005e14;
c = 299792458;
radv = 7.2921150e-05;  %地球旋转速率（CGCS2000下）
GM = 3.986004418e+14;  %地球标准引力常数（CGCS2000下）

toe = GEO_eph.toe;
sqrtA = GEO_eph.sqrtA;
Delta_n = GEO_eph.Delta_n;
M0 = GEO_eph.M0;
e = GEO_eph.e;
omega = GEO_eph.omega;
Cuc = GEO_eph.Cuc;
Cus = GEO_eph.Cus;
Crc = GEO_eph.Crc;
Crs = GEO_eph.Crs;
Cic = GEO_eph.Cic;
Cis = GEO_eph.Cis;
IDOT = GEO_eph.IDOT;
i0 = GEO_eph.i0;
omega0 = GEO_eph.omega0;
omegaDot = GEO_eph.omegaDot;
tk = tsv - toe - 14;%计算规划时间(观测历元到参考历元的时间差) 

A=(sqrtA)^2;      %计算长半轴
n = sqrt(GM/A^3) + Delta_n;%改正平均角速度
Mk=M0+n*tk;       %计算平近点角Mk
Mk=mod(Mk,2*pi);  %将计算结果的角度限制在[0, 2*pi)中

%迭代计算偏近点角E，由 Mk = E - esinE，其中Mk已知，迭代求解:Ei=M+esin(Ei-1)
E0 = Mk;
while 1
    Ei = Mk + e*sin(E0);
    delt_E = Ei - E0;
    delt_E = rem(delt_E,2*pi);%rem函数将结果限制在（-2pi,2pi）中
    if abs(delt_E) < 1e-12%abs取绝对值
        E = Ei;
        break;  %绝对值小于阈值时，将当时的Ei赋给E，跳出
    else
        E0 = Ei;
    end
end

dt_rel = -2 * sqrt(mu) * e * sqrtA * sin(E) / c^2;  %相对论效应
vk = atan2(sqrt(1-e^2)*sin(E),cos(E)-e);   %计算真近点角vk
phi_k = vk+omega;  %计算纬度幅角参数

delt_uk = Cus*sin(2*phi_k) + Cuc*cos(2*phi_k);   %纬度幅角改正项
delt_rk = Crs*sin(2*phi_k) + Crc*cos(2*phi_k);   %径向改正项
delt_ik = Cis*sin(2*phi_k) + Cic*cos(2*phi_k);   %轨道倾角改正项

uk = phi_k + delt_uk;           %计算改正后的纬度幅角
rk = A*(1-e*cos(E)) + delt_rk;  %计算改正后的径向
ik = i0 + IDOT*tk + delt_ik;    %计算改正后的轨道倾角

omega_k = omega0 + omegaDot*tk - radv*toe;  %计算历元升交点赤经（惯性系中 因为是GEO卫星）
omega_k = mod(omega_k,2*pi);

xk = rk*cos(uk);     %卫星在轨道平面内的坐标
yk = rk*sin(uk);

X_GK = xk*cos(omega_k) - yk*cos(ik)*sin(omega_k);     %计算GEO卫星在自定义坐标系下的坐标
Y_GK = xk*sin(omega_k) + yk*cos(ik)*cos(omega_k);
Z_GK = yk*sin(ik);

R_x = [1,0,0;0,cosd(-5),sind(-5);0,-sind(-5),cosd(-5)];   %旋转矩阵
R_z = [cos(radv*tk),sin(radv*tk),0;-sin(radv*tk),cos(radv*tk),0;0,0,1];

matrixB = [X_GK; Y_GK; Z_GK];
matrixA = R_z*R_x*matrixB;      %卫星在CGCS2000中的坐标
Xk=matrixA(1);  Yk=matrixA(2);  Zk=matrixA(3);
end

function [Xk,Yk,Zk,dt_rel] = MEO_STATE_tsv(tsv,MEO_eph)
mu = 3.986005e14;
c = 299792458;
GM = 3.986004418e+14;%地球标准引力常数（CGCS2000下）
radv = 7.2921150e-05;%地球旋转速率（CGCS2000下）

sqrtA = MEO_eph.sqrtA;
toe = MEO_eph.toe;
Delta_n = MEO_eph.Delta_n;
M0 = MEO_eph.M0;
e = MEO_eph.e;
i0 = MEO_eph.i0;
IDOT = MEO_eph.IDOT;
omega = MEO_eph.omega;
Cuc = MEO_eph.Cuc;
Cus = MEO_eph.Cus;
Crc = MEO_eph.Crc;
Crs = MEO_eph.Crs;
Cic = MEO_eph.Cic;
Cis = MEO_eph.Cis;
omega0 = MEO_eph.omega0;
omegaDot = MEO_eph.omegaDot;
tk = tsv - toe-14 ; 

A = (sqrtA)^2;
n = sqrt(GM/A^3) + Delta_n;
Mk = M0 + n*tk;
Mk=mod(Mk,2*pi);

E0 = Mk;
while 1
    Ei = Mk + e*sin(E0);
    delt_E = Ei - E0;
    delt_E = rem(delt_E,2*pi);%rem函数将结果限制在（-2pi,2pi）中
    if abs(delt_E) < 1e-12%abs取绝对值
        E = Ei;
        break;%绝对值小于这个阈值时，将当时的Ei赋给E，跳出循环
    else
        E0 = Ei;
    end
end

dt_rel = -2 * sqrt(mu) / c^2 * e * sqrtA * sin(E);
vk = atan2(sqrt(1-e^2)*sin(E),cos(E)-e);
phi_k = vk + omega;

delt_uk = Cus*sin(2*phi_k) + Cuc*cos(2*phi_k);
delt_rk = Crs*sin(2*phi_k) + Crc*cos(2*phi_k);
delt_ik = Cis*sin(2*phi_k) + Cic*cos(2*phi_k);

uk = phi_k + delt_uk;
rk = A*(1-e*cos(E)) + delt_rk;
ik = i0 + IDOT*tk + delt_ik;

xk = rk*cos(uk);%在轨道平面的坐标
yk = rk*sin(uk);

omega_k = omega0 + (omegaDot-radv)*tk - radv*toe;%升交点赤经（此过程与GEO卫星不同）
omega_k = mod(omega_k,2*pi);

Xk = xk*cos(omega_k) - yk*cos(ik)*sin(omega_k);%在CGCS2000坐标中的坐标
Yk = xk*sin(omega_k) + yk*cos(ik)*cos(omega_k);
Zk = yk*sin(ik);
end

function [results_table, removed_table] = least_squares(sat_table, obs_table, approx_coords)
    c = 299792458;
    approx_coords_initial = approx_coords;
    merged = innerjoin(sat_table, obs_table, 'Keys', {'Time', 'SatelliteID'});% 合并卫星坐标与观测值
    times  = unique(merged.Time);

    final_coords = zeros(length(times), 3);
    iter_counts  = zeros(length(times), 1);
    removed_log  =  cell(length(times), 1);
 
    max_iterations = 10;
    threshold = 0.01;
    for k = 1:length(times)
        t = times(k);
        idx = merged.Time == t;
        data_t = merged(idx, :);
    
        s_coords = [data_t.X, data_t.Y, data_t.Z];
        P_IF = data_t.P_IF;
        V_ts = data_t.V_ts;
    
        nSat = size(s_coords, 1);
        dt=0;
        % 迭代求解（无剔除）
        for iter = 1:max_iterations
        [B, L, el] = calculate_B_L(s_coords, P_IF, approx_coords, V_ts,dt);
        P = diag(sin(el).^2);     
        x_hat = (B'*P*B)\(B'*P*L);
            if norm(x_hat(1:3)) < threshold && abs(x_hat(4)) < 1e-6 
                break;
            end
        approx_coords = approx_coords + x_hat(1:3)';
        dt = dt + x_hat(4)/c;
        end
        v = B * x_hat + L;
        v_table = table(data_t.SatelliteID, v, 'VariableNames', {'SatelliteID', 'Residual_m'});
    
        % 剔除残差大于100m的卫星（循环直到满足条件）
        removed_this_epoch = {};
        while max(abs(v_table.Residual_m)) > 100 && nSat >= 4
            [~, idx_max] = max(abs(v_table.Residual_m));
            sat_to_remove = char(v_table.SatelliteID(idx_max));
            removed_this_epoch{end+1} = sat_to_remove;%#ok
        
            keep_idx = ~strcmp(data_t.SatelliteID, sat_to_remove);
            data_t = data_t(keep_idx, :);
        
            s_coords = [data_t.X, data_t.Y, data_t.Z];
            P_IF = data_t.P_IF;
            V_ts = data_t.V_ts;
            nSat = size(s_coords, 1);
            approx_coords = approx_coords_initial;
            dt = 0 ;

            for iter = 1:max_iterations   % 重新迭代
            [B, L, el] = calculate_B_L(s_coords, P_IF, approx_coords, V_ts,dt);
            P = diag(sin(el).^2);     
            x_hat = (B'*P*B)\(B'*P*L);
               if norm(x_hat(1:3)) < threshold && abs(x_hat(4)) < 1e-6
                 break;
               end
            approx_coords = approx_coords + x_hat(1:3)';
            dt = dt + x_hat(4)/c;
            end
            v = B * x_hat + L;
            v_table = table(data_t.SatelliteID, v, 'VariableNames', {'SatelliteID', 'Residual_m'});
        end

        removed_log{k} = removed_this_epoch;
        final_coords(k, :) = approx_coords;
        iter_counts(k) = iter;
    end
    results_table = table(times, final_coords, iter_counts);% 构建结果表（与原有格式一致）

    has_removal = ~cellfun(@isempty, removed_log);    % 构建剔除记录表
    if any(has_removal)
        removed_times = times(has_removal);
        removed_list = removed_log(has_removal);
        removed_str = cellfun(@(x) strjoin(x, ', '), removed_list, 'UniformOutput', false);
        removed_table = table(removed_times, removed_str, 'VariableNames', {'Time', 'RemovedSatellites'});
    else
        removed_table = table();
    end
end

function [B, L, el] = calculate_B_L(s_coords, P_IF, approx_coords, V_ts, dt)
    c = 299792458;             
    radv = 7.2921150e-5;          % 地球自转角速度 (rad/s)
    nSat = size(s_coords, 1);

    rou0_values   = zeros(nSat, 1);
    direction_cos = zeros(nSat, 3);
    trop_values   = zeros(nSat, 1);
    el            = zeros(nSat, 1);   %  预分配高度角输出向量

    tol_tau = 1e-12;
    max_iter_tau = 10;

    % 将近似坐标转换为大地坐标（纬度、经度、椭球高）
    sat_lla = ecef2lla(approx_coords);
    lat = deg2rad(sat_lla(1));   % 纬度 (rad)
    lon = deg2rad(sat_lla(2));   % 经度 (rad)
    h_s = sat_lla(3);            % 高度 (m)

    % 标准气象参数（海平面）
    P0 = 1013.25;   % 气压 (hPa)
    T0 = 288.15;    % 气温 (K)
    e0 = 11.69;     % 水汽压 (hPa)

    % ===== 气象参数高度改正 =====
    if h_s > 11000
        e_s = 0;
    else
        e_s = e0 * (1 - 0.0065 * h_s / T0)^4;
    end
    P_s = P0 * (1 - 0.0065 * h_s / T0)^5.2561;
    T_s = T0 - 0.0068 * h_s;

    % ===== 天顶对流层延迟 (Saastamoinen 模型) =====
    denom = 1 - 0.00266 * cos(2*lat) - 0.00028 * h_s/1000;
    ZHD = (0.002277 * P_s) / denom;
    ZWD = 0.002277 * (1255 / T_s + 0.05) * e_s / denom;
    ZTD = ZHD + ZWD;   % 天顶总延迟 (m)

    % ===== ENU 旋转矩阵（固定）=====
    R_neu = [    -sin(lon),           cos(lon),        0 ;
        -sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat) ;
         cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat)];

    % 主循环：处理每一颗卫星
    for i = 1:nSat
        rcv = approx_coords(:);
        sat_r = s_coords(i, :)';

        % ===== 地球自转改正（迭代计算信号传播时间）=====
        tau = 0;
        for k = 1:max_iter_tau
            theta = radv * tau;
            Rz = [ cos(theta), sin(theta), 0;
                  -sin(theta), cos(theta), 0;
                   0,          0,          1];
            sat_e = Rz * sat_r;            % 改正后的卫星位置

            rho = norm(rcv - sat_e);
            if isnan(rho) || rho < 1e6
                error('rho异常，检查输入数据');
            end

            tau_new = rho / c;
            if abs(tau_new - tau) < tol_tau
                tau = tau_new;%#ok
                break;
            end
            tau = tau_new;
        end
        rou0_values(i) = rho;

        % ===== 方向余弦 =====
        dx = sat_e - rcv;
        direction_cos(i,:) = (dx / rho)';

        % ===== 计算高度角（弧度）=====
        neu = R_neu * dx;
        elev = atan2(neu(3), sqrt(neu(1)^2 + neu(2)^2));   % 当前卫星高度角
        el(i) = elev;                                      % 存储至输出向量

        % ===== 对流层斜延迟（简单映射函数 1/sin(elev)）=====
        if elev <= deg2rad(5)           % 低于 5° 时舍弃对流层改正
            trop = 0;
        else
            trop = ZTD / sin(elev);
        end
        trop_values(i) = trop;
    end

    % ===== 构建设计矩阵 B =====
    B = [-direction_cos, ones(nSat, 1)];

    % ===== 计算观测残差向量 L =====
    L = P_IF - rou0_values + c*(V_ts - dt) - trop_values;
end

function coords = approx_position(filename)   %读取接收机近似坐标
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件: %s', filename);
    end
    coords = [0,0,0];
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, 'END OF HEADER')
            break;
        end
        if contains(line, 'APPROX POSITION XYZ')
            parts = strsplit(strtrim(line));
            coords = [str2double(parts{1}), str2double(parts{2}), str2double(parts{3})];
            break;
        end
    end
    fclose(fid);
end
