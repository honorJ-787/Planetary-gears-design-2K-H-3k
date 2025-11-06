function results = WW_cal(app,ic,m,np,errd,zmin,zmax,K,e,haax,hagx,jx)%单独使用将app改为~
%WW需要以50为界分类讨论
%results = WW_cal(0,15,2,3,0.1,15,230,1,1,1.5,1.5)
%WW轮系可视为一般传动链
%预设结果数组
max_results  = 60000;
temp_results = zeros(max_results, 9);
result_count = 0;
results = [];
%前处理
if isempty(ic) || ic <= 1
    error('ic 必须提供且 > 1')
end
if isempty(np)
    np=1;
end
if isempty(errd)
    errd = 0;
end
ic0=abs(ic);
%sig00ma和传动比关系选取
if      ic0 >= 2500
    sig0 = 1;
elseif  ic0 >= 1000
    sig0 = 2;
elseif  ic0 >= 400
    sig0 = 3;
elseif  ic0 >= 100
    sig0 = 4:6;      % 定义为一个整数范围
elseif  ic0 >= 50
    sig0 = 7:10;     % 定义为一个整数范围
else
end

if ic<=50%小传动比
    if K==1&&e==1%凑使za 为np 的倍数加1,即满足装配条件
        for i=1:200
            max0=max(zmin,i*np+1);
            min0=min(zmax,i*np+1);
            for za=max0:min0%步长
                if ic>0
                    zb=za-1;
                    zf=zb;
                    zg=zf-1;
                end
                if ic<0
                    zb=za+1;
                    zf=zb;
                    zg=zf+1;
                end
                if any([zg,zb,zf] < zmin | [zg,zb,zf] > zmax)
                    continue;
                end
                i=za*zf/(za*zf-zb*zg);
                errc=abs((i-ic)/ic);
                if errc > errd
                    continue;
                end
                a = m/2*(za + zg);
                result_count = result_count + 1;
                if result_count <= max_results
                    temp_results(result_count, :) = [za, zg, zf, zb, a, i, errc*100, zg/za, zb/zf];
                end
            end
        end

    else
        for za=zmin:zmax
            zb=za-e;
            zf0=(e*za-e^2)/(za/ic-e);
            ran=(floor(zf0)):(ceil(zf0));
            for zf=ran
                zg=zf-e;
                a=m*(za+zg)/2;
                %邻接条件，%%用以先判断最大行星轮或最大传动比,理想均布状态
                if np>2%np=1,2，理论上无限
                    L=a*(sin(2*pi/np)/sin(pi/2-pi/np));%正弦定理计算行星轮之间的中心距
                    dag=m*(zg+2*hagx);%行星轮齿顶圆
                    delta=L-dag;%控制间隙
                    if delta<jx
                        error('传动比过大，请减小。或减少行星轮个数。')
                    end
                end
                i=za*zf/(za*zf-zb*zg);
                errc=abs((i-ic)/ic);
                if errc > errd
                    continue;
                end
                a = m/2*(za + zg);
                result_count = result_count + 1;
                if result_count <= size(temp_results, 1)
                    temp_results(result_count, :) = [za, zg, zf, zb, a, i, errc * 100, zg/za, zb/zf];
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

                angles = cumsum([0, gaps(1:end-1)]);
                angles = mod(angles, 360);
            end
            continue;
        end
    end
end

if ic>50%大传动比，np至多为2，np=1相当于一般的外啮合传动链系统
    if np>2
        error ('np不大于2,否则行星轮碰撞');
    end
    %依次带sigma计算
    for sig=sig0
        %允许不满足装配条件，以获得最多的组合类型
        za0=sqrt(sig*ic+((sig-1)/2)^2-(sig-1)/2);
        za1=((floor(za0)):(ceil(za0)));
        for za=za1
            zf=za+sig-1;
            zg=za+sig;
            zb=zf-sig;
            a=m/2*(za+zg);
            i=za*zf/(za*zf-zb*zg);
            errc=abs((i-ic)/ic);
            if errc>errd
                continue;
            end
            result_count = result_count + 1;
            if result_count <= size(temp_results, 1)
                temp_results(result_count, :) = [za, zg,zf, zb, a, i, errc * 100, zg/za, zb/zf];
            end
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
