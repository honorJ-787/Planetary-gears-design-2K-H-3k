function results = NW_cal(app,ic,m,np,errd,zmin,zmax,gaNW,gamma)
% 预设结果数组
max_results = 60000;
A_store = zeros(max_results,9);
B_store = zeros(max_results,9);
A_count = 0;
B_count = 0;
if isempty(ic) || ic <= 1
    error('ic 必须提供且 > 1')
end
haa = 1.0;  hga = 1.0;
% 读取gamma参数
if exist('app','var') && ~isempty(app)
    try
        %
        if (nargin < 8 || isempty(gaNW)) && isprop(app,'gaNW') && ~isempty(app.gaNW.Value)
            gaNW = app.gaNW.Value;
        end
        if (nargin < 9 || isempty(gamma)) && ~isempty(app) && isprop(app,'gaNW') && ~isempty(app.gaNW.Value)
            gamma = app.gaNW.Value;
        end
        % 齿顶高控件
        if isprop(app,'haStar_a') && ~isempty(app.haStar_a); haa = app.haStar_a; end
        if isprop(app,'haStar_g') && ~isempty(app.haStar_g); hga = app.haStar_g; end
    catch
    end
end
%两种方法独立计算
useMethod1 = isfinite(gamma) && (gamma > 1) && (gamma < 2);
useMethod2 = true;
%方法1：径向轮廓法
if useMethod1
    try
        if ic==1
            error('无效参数：ic 不能等于 1。');
        end
        % gamma=(1+2i1)*(ic-1-i1)/((ic-1)*(i1+1)) 的等价整理
        delta = ((ic-1)*(gamma-2)+1)^2 - 8*(ic-1)*(gamma-1);
        if delta < 0
            error('无解');
        end
        r1 = ( -((ic-1)*(gamma-2)+1) + sqrt(delta) ) / 4;
        r2 = ( -((ic-1)*(gamma-2)+1) - sqrt(delta) ) / 4;
        cand = [r1, r2];
        mask = isreal(cand) & (cand > 0) & (abs(ic.*cand - cand + ic - 1) > eps);
        pos_roots = cand(mask);
        if isempty(pos_roots)
            error('无正实根满足方程与分母非零条件。');
        end

        % 逐个正实根
        for i1 = pos_roots(:).'
            for za = zmin:zmax
                % 邻接上界
                imax = 2*(hga-(haa+hga)*sin(pi/np)-za) / (za*(sin(pi/np)-1));
                i1_eff = i1;
                if isfinite(imax) && i1_eff > imax
                    i1_eff = imax;
                    if i1_eff <= 0
                        continue;
                    end
                end

                % zg≈i1_eff*za 的近邻整数
                zg0 = i1_eff*za;
                zgc = round(zg0);
                if abs(zg0 - zgc) < 1e-12
                    ran_zg = zgc;
                else
                    ran_zg = (floor(zg0)-1):(ceil(zg0)+1);
                    ran_zg = ran_zg(ran_zg>0);
                    if isempty(ran_zg)
                        ran_zg = max(1,floor(zgc-1)):max(1,floor(zgc+2));
                        ran_zg = ran_zg(ran_zg>0);
                    end
                end

                for zg = ran_zg
                    i2 = (ic - 1)/i1_eff;% ic = 1 + i1*i2

                    for zf = zmin:zg
                        zb0 = i2*zf;
                        zbc = round(zb0);
                        if abs(zb0 - zbc) < 1e-12
                            ran_zb = zbc;
                        else
                            ran_zb = (floor(zb0)-1):(ceil(zb0)+1);
                            ran_zb = ran_zb(ran_zb>0);
                            if isempty(ran_zb)
                                ran_zb = max(1,floor(zbc-1)):max(1,floor(zbc+2));
                                ran_zb = ran_zb(ran_zb>0);
                            end
                        end

                        for zb = ran_zb
                            % 同心条件
                            if abs( za*(1+i1_eff) - zf*(i2-1) ) > 1e-6
                                continue;
                            end

                            % 均布
                            s = gcd(zg, zf);
                            q = (za*zf + zb*zg)/(s*np);
                            if mod(za*zf + zb*zg, s*np) == 0
                                beta = 360/np;
                            else
                                qp   = round(q);
                                beta = 360/(za+zb)*qp;
                                % 邻接再校核(同NGW)
                                if ( (za + zg + haa + hga)*sin(pi/np) ) <= (zg + haa)
                                    continue;
                                end
                            end

                            % 安装双缝
                            theta0 = 360/np;
                            gapA   = beta;
                            gapB   = 360 - (np - 2)*theta0 - gapA;
                            if gapB <= 0
                                continue;
                            end
                            gaps = theta0*ones(1, np);
                            gaps(1) = gapA;
                            gaps(1 + floor(np/2)) = gapB;
                            angles = mod(cumsum([0, gaps(1:end-1)]), 360);

                            % 传动比 & 误差
                            i0    = 1 + (zg*zb)/(za*zf);
                            errv  = abs(i0 - ic)/ic;
                            if errv > errd
                                continue;
                            end
                            a = m/2*(za + zg);
                            i_ag = zg/za;
                            i_bf = zb/zf;

                            % 记录 A
                            A_count = A_count + 1;
                            if A_count <= size(A_store,1)
                                i_act = i0;
                                A_store(A_count,:) = [za, zg, zf, zb,a,i_act,errv*100, i_ag, i_bf];
                            end
                        end % zb
                    end % zf
                end % zg
            end % za
        end % i1
    catch ME
        warning('[NW-方法1] 已跳过');
    end
else
    % γ 不在 (1,2)，跳过方法1
end
% 方法2：手册枚举
if useMethod2
    try
        p0 = ic - 1;
        for A = 1.2:0.01:4.2
            B = p0/A;
            if B < 2.4 || B > 4.8 || abs(B-1) < 1e-12
                continue;
            end
            for K = [4,5,6,7]
                za = K*np;
                if za<=zmin
                    continue;
                end
                zf0 = za * ((A+1)/(B-1));    % 同心条件（非角度变位）
                for zf = max(zmin,floor(zf0)-1) : min(zmax,ceil(zf0)+1)
                    if zf<=0, continue; end
                    zg0 = A*za;
                    for zg = max(zmin,floor(zg0)-1) : min(zmax,ceil(zg0)+1)
                        if zg<=0, continue; end
                        zb0 = B*zf;
                        for zb = max(zmin,floor(zb0)-1) : min(zmax,ceil(zb0)+1)
                            if zb<=0, continue; end

                            i0 = 1 + (zg*zb)/(za*zf);
                            if abs(i0 - ic) > 1e-9
                                continue;
                            end
                            % 均布（计算第一排即可）
                            s = gcd(zg, zf);
                            q = (za*zf + zb*zg)/(s*np);
                            if mod(za*zf + zb*zg, s*np) == 0
                                beta = 360/np;
                            else
                                qp   = round(q);
                                beta = 360/(za+zb)*qp;
                                if ( (za + zg + haa + hga)*sin(pi/np) ) <= (zg + haa)
                                    continue;
                                end
                            end
                            i_ag = zg/za;
                            i_bf = zb/zf;
                            % 安装双缝
                            theta0 = 360/np;
                            gapA   = beta;
                            gapB   = 360 - (np - 2)*theta0 - gapA;
                            if gapB <= 0
                                continue;
                            end
                            gaps = theta0*ones(1, np);
                            gaps(1) = gapA;
                            gaps(1 + floor(np/2)) = gapB;
                            angles = mod(cumsum([0, gaps(1:end-1)]), 360);

                            % 记录 B
                            B_count = B_count + 1;
                            if B_count <= size(B_store,1)
                                i_act = i0;
                                errv  = abs(i0 - ic)/ic;
                                if errv > errd
                                    continue;
                                end
                                a=m/2*(za + zg);
                                B_store(B_count,:) = [za, zg, zf, zb,a, i_act, errv*100, i_ag, i_bf];
                            end
                        end % zb
                    end % zg
                end % zf
            end % K
        end % A
    catch ME
    end
end
% 合并返回
if A_count>0 || B_count>0
    results = [A_store(1:A_count,:); B_store(1:B_count,:)];
    results(:,1:4) = round(results(:,1:4));%去除重复结果
    results = sortrows(results, [1 2 3 4 7]);
    [~, ia] = unique(results(:,1:4), 'rows', 'stable');
    results = results(ia, :);
else
    results = [];
end
end
