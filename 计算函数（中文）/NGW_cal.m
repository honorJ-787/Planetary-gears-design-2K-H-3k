% NGW非变位,%%标记算法在历遍时剔除多余数据
function results = NGW_cal(app,ic,m,np,errd,zmin,zmax)%单独使用将app改为~
%zmin=17;
%zmax=210;
% 预设结果数组
max_results  = 60000;
temp_results = zeros(max_results, 6);
result_count = 0;
if isempty(ic) || ic <= 1
    error('ic 必须提供且 > 1')
end
%齿轮参数
%daa=(za+haa)*m;
%dag=(zg+hgahaa)*m;
%dbg=(zb-hba)*m;%斜齿相同
haa=1.0;   %这里可以修改
hga=1.0;   %app中可以直接由用户输入的值返回
if isempty(errd)
    errd = 0;
end
% 以最小太阳轮开始历遍
for za=zmin:zmax
    %邻接条件，%%用以先判断最大行星轮或最大传动比,理想均布状态
    imax=2*(hga-(haa+hga)*sin(pi/np)-za)/(za*(sin(pi/np)-1));
    if ~isempty(ic)
        fprintf('理论最大传动比%.6f\n',imax)
    else
        if abs(ic-imax)>1e-6%最小间隙
            fprintf('需要的传动比过大，请减小或减少行星轮个数。当前允许最大传动比%.6f\n',imax)
        else
            zb0=(imax-1)*za;%由最大齿数找上界
            %zb1=floor(zb0);%向下取
            zb1=floor(zb0)+1;%加一法
        end
    end
    %传动比条件，注意先处理误差0情况
    if isempty(errd)%未输入默认为0
        errd = 0;
    end
    if errd==0
        zb2=(ic-1)*za;
        if abs(zb2-round(zb2))>1e-12
            continue;%下一个za
        end
        ran=round(zb2);%后续只校核该整数
    else
        zb0=(ic-1)*za;
        %保证存在上界
        if ~exist('zb1','var') || isempty(zb1)
            zb1 = inf;
        end
        ran=(floor(zb0)-2):(ceil(zb0)+2);
        ran=ran(ran>=zmin & ran<=zmax);%裁剪边界
        ran=ran(ran<zb1);%当冗余较大用最大传动比控制上界
    end
    for zb=ran%在ran范围内对zb历遍
        i=1+zb/za;%%剔除不符合误差的zb
        errc=(i-ic)/ic;
        errc=abs(errc);
        if errc>errd
            continue;
        end
        %同心条件
        ba=zb-za;
        if mod(ba,2)~=0
            continue;
        end
        zg=ba/2;
        if zg<zmin
            continue;
        end
        %装配条件
        q=(za+zb)/np;%装配系数
        if mod(za+zb,np)==0%判断是否为整数，均布才为整数
            beta=360/(za+zb)*q;
            angles=(0:np-1)*(360/np);
            gaps=ones(1,np)*(360/np);
        else
            qp=round(q);%取整
            beta=360/(za+zb)*qp;
            if qp-(zb-za)/np>=1e-3%需要再次校核邻接条件
                if (za+zg+haa+hga)*sin(pi/np)<=zg+haa
                    continue;
                end
            end
        end
        %计算安装角度
        theta0=360/np;
        gapA=beta;
        gapB=360-(np-2)*theta0-gapA;
        if gapB<=0
            continue;
        end
        gaps=theta0*ones(1,np);
        idx1=1;
        idx2=1+floor(np/2);
        gaps(idx1)=gapA;
        gaps(idx2)=gapB;
        a=m/2*(za+zg);
        % 保存结果
        result_count=result_count + 1;
        if result_count <= size(temp_results, 1)
            temp_results(result_count, :) = [za, zg, zb, a, i, errc * 100];
        end
    end
    angles = cumsum([0, gaps(1:end-1)]);
    angles = mod(angles, 360);
    continue

end%历遍完成

% 返回有效结果（APP）
if result_count > 0
    results = temp_results(1:result_count, :);
else
    results = [];
end
end


