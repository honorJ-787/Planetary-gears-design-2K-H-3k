function results = NGW_cal(app,ic,m,np,errd,zmin,zmax,haax,hagx,jx)

% 预设结果数组
max_results  = 60000;
temp_results = zeros(max_results, 6);
result_count = 0;

%预处理
if isempty(errd)
    errd = 0;
end
if isempty(ic) || ic <= 1
    error('ic 必须提供且 > 1')
end
if isempty(np)
    np = 1;
end

% 以最小太阳轮开始历遍
for za=zmin:zmax
    gaps=[];%新一轮循环，重置

    %传动比条件
    zb0=(ic-1)*za;
    zb1=ceil(zb0);
    if zb1>=zmax
        continue;
    else
    end
    ran=(round(zb0-2)):(round(zb0+2));
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
        a=m*(za+zg);
        %已获得za,zb,zg所有数据，开始校核
        %邻接条件
        L=a*(sin(2*pi/np)/sin(pi/2-pi/np));%正弦定理计算行星轮之间的中心距
        dag=zg*m+2*hagx;%行星轮齿顶圆
        delta=L-dag;%控制间隙
        if np>2%np=1,2，理论上无限
            if delta<jx
                error('传动比过大，请减小。或减少行星轮个数。')
            end
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
                if (za+zg+haax+hagx)*sin(pi/np)<=zg+haax
                    continue;
                end
            end
        end

        %计算安装角度(APP)
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

        % 保存结果
        result_count=result_count + 1;
        if result_count <= size(temp_results, 1)
            temp_results(result_count, :) = [za, zg, zb, a, i, errc * 100];
        end
    end
    angles = cumsum([0, gaps(1:end-1)]);

    %若不存在解，gaps不被定义，以均布的形式兜底
    if ~exist('gaps','var') || isempty(gaps)
        gaps = ones(1,np)*(360/np);
    end
    angles = mod(angles, 360);
    continue

end%历遍完成

% 返回有效结果(app)
if result_count > 0
    results = temp_results(1:result_count, :);
else
    results = [];
end
end
