% NGW non-profile-shifted, %% mark: the algorithm prunes redundant data during iteration
function results = NGW_cal(app,ic,m,np,errd,zmin,zmax)% For standalone use, change 'app' to '~'
% zmin=17;
% zmax=210;
% Preallocate result array
max_results  = 60000;
temp_results = zeros(max_results, 6);
result_count = 0;
if isempty(ic) || ic <= 1
    error('ic Must provide and > 1')
end
% Gear parameters
% daa=(za+haa)*m;
% dag=(zg+hgahaa)*m;
% dbg=(zb-hba)*m; % same for helical gears
haa=1.0;   % can be modified here
hga=1.0;   % in the app this can be obtained directly from user input
if isempty(errd)
    errd = 0;
end
% Iterate starting from the minimum sun gear
for za=zmin:zmax
    % Adjacency condition, %% first to determine max planet count or max ratio under ideal uniform distribution
    imax=2*(hga-(haa+hga)*sin(pi/np)-za)/(za*(sin(pi/np)-1));
    if ~isempty(ic)
        fprintf('Theoretical maximum transmission ratio%.6f\n',imax)
    else
        if abs(ic-imax)>1e-6 % minimum clearance
            fprintf('The required transmission ratio is too large, please reduce or reduce the number of planetary gears. The current maximum transmission ratio allowed is%.6f\n',imax)
        else
            zb0=(imax-1)*za; % get upper bound from maximum tooth count
            % zb1=floor(zb0); % round down
            zb1=floor(zb0)+1; % plus-one method
        end
    end
    % Ratio condition; handle zero-error case first
    if isempty(errd) % default to 0 if not provided
        errd = 0;
    end
    if errd==0
        zb2=(ic-1)*za;
        if abs(zb2-round(zb2))>1e-12
            continue; % next za
        end
        ran=round(zb2); % only check this integer later
    else
        zb0=(ic-1)*za;
        % ensure an upper bound exists
        if ~exist('zb1','var') || isempty(zb1)
            zb1 = inf;
        end
        ran=(floor(zb0)-2):(ceil(zb0)+2);
        ran=ran(ran>=zmin & ran<=zmax); % clip to bounds
        ran=ran(ran<zb1); % when redundancy is large, use max ratio to cap upper bound
    end
    for zb=ran % iterate zb within 'ran' range
        i=1+zb/za; %% discard zb not meeting error limit
        errc=(i-ic)/ic;
        errc=abs(errc);
        if errc>errd
            continue;
        end
        % Concentricity condition
        ba=zb-za;
        if mod(ba,2)~=0
            continue;
        end
        zg=ba/2;
        if zg<zmin
            continue;
        end
        % Assembly condition
        q=(za+zb)/np; % assembly factor
        if mod(za+zb,np)==0 % check divisibility; uniform only when divisible
            beta=360/(za+zb)*q;
            angles=(0:np-1)*(360/np);
            gaps=ones(1,np)*(360/np);
        else
            qp=round(q); % round to integer
            beta=360/(za+zb)*qp;
            if qp-(zb-za)/np>=1e-3 % need to recheck adjacency
                if (za+zg+haa+hga)*sin(pi/np)<=zg+haa
                    continue;
                end
            end
        end
        % Compute installation angles
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
        % Save result
        result_count=result_count + 1;
        if result_count <= size(temp_results, 1)
            temp_results(result_count, :) = [za, zg, zb, a, i, errc * 100];
        end
    end
    angles = cumsum([0, gaps(1:end-1)]);
    angles = mod(angles, 360);
    continue

end % iteration finished

% Return valid results (APP)
if result_count > 0
    results = temp_results(1:result_count, :);
else
    results = [];
end
end
