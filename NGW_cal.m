 % NGW非变位
        function results = NGW_cal(~,i_req, m, n_p, errd)
            z_min = 17;
            z_max = 210;
            % 预设结果数组
            max_results = 1000;
            temp_results = zeros(max_results, 7);
            result_count = 0; 
            % 遍历算法
            for z_a = z_min:z_max
                % 装配系数
                q_ideal = (i_req * z_a) / n_p;
                % q 候选值：floor, round, ceil，附近 +-1
                q_candidates=unique([floor(q_ideal-1), floor(q_ideal), round(q_ideal), ceil(q_ideal), ceil(q_ideal+1)]);
                q_candidates=q_candidates(q_candidates >= 1);
                for q=q_candidates
                    % 根据装配条件计算内齿圈齿数
                    z_b = n_p * q - z_a;
                    if z_b <= 0
                        continue;
                    end
                    
                    % 同心条件：z_g = (z_b - z_a)/2 为正整数
                    if mod((z_b - z_a), 2) ~= 0
                        continue;
                    end
                    z_g = (z_b - z_a) / 2;
                    if z_g <= 0
                        continue;
                    end
                    
                    % 实际传动比
                    i_act = 1 + double(z_b) / double(z_a);
                    errv = abs(i_act - i_req) / i_req;
                    if errv > errd
                        continue;
                    end
                    
                    % 根切
                    min_teeth_allowed = 17;
                    if z_a < min_teeth_allowed || z_b < min_teeth_allowed || z_g < min_teeth_allowed
                        continue;
                    end
                    
                    % 计算中心距
                    a = m / 2 * (z_a + z_g);
                    
                    % 保存结果
                    result_count = result_count + 1;
                    if result_count <= size(temp_results, 1)
                        temp_results(result_count, :) = [z_a, z_g, z_b, a, q, i_act, errv * 100];
                    end
                end
            end
            % 返回有效结果（APP）
            if result_count > 0
                results = temp_results(1:result_count, :);
            else
                results = [];
            end
        end