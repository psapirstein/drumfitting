function [TR, TT, ER, best_inlierPerc, best_closepts, t] = ICPinsideFitting(q,p,its,BaseDeltaXY,BaseDeltaZ,BaseTheta,CloseThreshold,CloseWeight)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual implementation

%Show the progress visually if on
showviz = 0;

% Allocate vector for RMS of errors in every iteration.
t = zeros(its+1,1);

% Start timer
tic;

Np = size(p,2);

% Transformed data point cloud
pt = p;

% Allocate vector for RMS of errors in every iteration.
ER = zeros(its+1,1);

% Initialize total transform vector(s) and rotation matric(es).
TT = zeros(3,1, its+1);
TR = repmat(eye(3,3), [1,1, its+1]);

kdOBJ = KDTreeSearcher(transpose(q));

t(1) = toc;

% Go into main iteration loop
for k=1:its
    DeltaXY = BaseDeltaXY*(1 + its-k)/its;
    DeltaZ = BaseDeltaZ*(1 + its-k)/its;
    Theta = BaseTheta*(1 + its-k)/its;
    
    % Do matching
    [match, distances] = match_kDtree(q,pt,kdOBJ);
    
    q_norm = sqrt(sum(q(:,match).^2,1));
    pt_norm = sqrt(sum(pt.^2,1));
    inliers = q_norm > pt_norm;
    closePts = distances < CloseThreshold;
    E = (CloseWeight*Np-sum(inliers))/sum(closePts);
    Einside = sum(inliers)/Np;
    
    if k == 1
        ER(k) = E;
        best_inlierPerc = Einside;
        best_closepts = sum(closePts);
    end
    
    ER(k+1) = E;
    Total_R = TR(:,:,k);
    Total_T = TT(:,:,k);
    for n = randperm(12)
        T = [0;0;0];
        R = eye(3);
        if n == 1
            T = [-DeltaXY;0;0];
        elseif n == 2
            T = [DeltaXY;0;0];
        elseif n == 3
            T = [0;-DeltaXY;0];
        elseif n == 4
            T = [0;DeltaXY;0];
        elseif n == 5
            T = [0;0;-DeltaZ];
        elseif n == 6
            T = [0;0;DeltaZ];
        elseif n == 7
            R = Xrot(-Theta);
        elseif n == 8
            R = Xrot(Theta);
        elseif n == 9
            R = Yrot(-Theta);
        elseif n == 10
            R = Yrot(Theta);
        elseif n == 11
            R = Zrot(-Theta);
        elseif n == 12
            R = Zrot(Theta);
        end
        % Add to the total transformation
        R = R*Total_R;
        T = R*Total_T+T;

        % Apply last transformation
        pt(1:3,:) = R * p(1:3,:) + repmat(T, 1, Np);

        % Do matching
        [match, distances] = match_kDtree(q,pt,kdOBJ);

        q_norm = sqrt(sum(q(:,match).^2,1));
        pt_norm = sqrt(sum(pt.^2,1));
        inliers = q_norm > pt_norm;
        closePts = distances < CloseThreshold;
        E = (CloseWeight*Np-sum(inliers))/sum(closePts);
        Einside = sum(inliers)/Np;

        if E <= ER(k+1) && CloseWeight*Einside > best_inlierPerc
            Total_R = R;
            Total_T = T;
            ER(k+1) = E;
            if Einside > best_inlierPerc
                best_inlierPerc = Einside;
            end
            best_closepts = sum(closePts);
        end
    end
    
    % Add to the total transformation
    TR(:,:,k+1) = Total_R;
    TT(:,:,k+1) = Total_T;
    % Apply last transformation
    pt(1:3,:) = TR(:,:,k+1) * p(1:3,:) + repmat(TT(:,:,k+1), 1, Np);
    t(k+1) = toc;
end

if showviz > 0
    figure(2);
    drawscatter3('Inside fit',q(:,1:8:end),pt(:,~inliers),pt(:,inliers),pt(:,closePts(1:8:end)));
    drawnow;
end

ER = ER(end);
TR = TR(:,:,end);
TT = TT(:,:,end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [match, mindist] = match_kDtree(~, p, kdOBJ)
[match, mindist] = knnsearch(kdOBJ,transpose(p));
match = transpose(match);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = Xrot(degree)
R = [1,0,0;0,cosd(degree),-sind(degree);0,sind(degree),cosd(degree)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = Yrot(degree)
R = [cosd(degree),0,-sind(degree);0,1,0;sind(degree),0,cosd(degree)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = Zrot(degree)
R = [cosd(degree),-sind(degree),0;sind(degree),cosd(degree),0;0,0,1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
