function [Parameters] = ICPFitDrumToModel(model_name, scan_name, scan_path)
startTime = clock;
%Open the scan file
[V,Vc,objLines]=readobj(scan_name,scan_path);

%Each scan should be rotated so the top is up on the Z axis, at least within
%a few degrees of the correct orientation, although the algorithm will check
%for negative tapers in case the scan is upside down. The scan does not
%need to be centered on any particular coordinate.

S = V - repmat(mean(V,1),size(V,1),1);
duration = etime(clock, startTime); nextTime = clock;
fprintf(' in %.1fs. Opening model...',duration);

[VI,~,~]=readobj(model_name,'');
I = VI - repmat(mean(VI,1),size(VI,1),1);
normI = normalizeModel(I);

%Set the sampling ratio for the model and scan
model_sub = ceil((length(normI(:,1))+10000)/40000);
scan_sub = ceil((length(S(:,1))-10000)/20000);
if length(S(1:scan_sub:end,1)) > 20000
    scan_sub = scan_sub+1;
end
if length(S(1:scan_sub:end,1)) < 15000
    scan_sub = scan_sub-1;
end
if scan_sub < 1
    scan_sub = 1;
end

duration = etime(clock, nextTime); nextTime = clock;
fprintf(' in %.1fs.\n Model, Scan sample ratios = %d, %d; sampled / original scan points = %d / %d\n', duration, model_sub, scan_sub, length(V(1:scan_sub:end,1)), length(V(:,1)));

%----------------- Stage 1a: coarse fitting ------------------%
%Set initial parameters for the model from the manually oriented scan.
ScanInitH = max(S(:,3))-min(S(:,3));
ScanInitR = max(max(S(:,1))-min(S(:,1)),max(S(:,2))-min(S(:,2)))/2;
Parameters = [ScanInitH, ScanInitR, ScanInitR, ScanInitR, 0, Inf, ScanInitH];

%Preliminary search through a wide range of radii and orientations
Rlbound = ScanInitR; %It's impossible for the R to be much smaller than the model dimensions
Rubound = Rlbound*2; %At first, search a wide range of radii for possible matches
deltaR = ScanInitR/5;
height = ScanInitH;

icpIts = 30; %Iterations for inside ICP fitting
icpBaseR = 8; %Initial rotation amount for ICP
icpDist = 0.020; %Distance threshold for counting 'close' points
E_best = Inf;
Ein_best = 0;
closeWeight = 1.30; %Weight for proximity, and also instructs to not accept a fit if inside points are this percentage lower than a previous fit
stepDown = 1.1; %This is a multiplier for delta; larger means wider ranges in the next search

sub_fac = 3; %Additional subsampling of the data for speed
Zrotation = 90; %Arbitrary rotation step applied to the scan
rotz = 5; %# times to apply the rotation to the scan
cpucores = 6;
its = 5;

for it = 1:its
    fprintf('\nTrying radii from %.3f to %.3f, step %.3f', Rlbound, Rubound, deltaR);
    icpBaseT = deltaR/2; %Initial translation amount for ICP
    
    R_par = repmat(eye(3,3),[1,1,cpucores]);
    T_par = zeros(3,1,cpucores);
    Rad_par = zeros(1,cpucores);
    E_par = zeros(1,cpucores);
    Ein_par = zeros(1,cpucores);
    close_par = zeros(1,cpucores);
    S_par = S(1:scan_sub*sub_fac:end,:);
    
    parfor i = 1:cpucores
        Rad_par(i) = Rlbound + deltaR*(i-1);
        
        rotTotal_R = repmat(eye(3,3),[1,1,rotz]);
        rotTotal_T = zeros(3,1,rotz);
        rotE = zeros(1,rotz);
        rotEin = zeros(1,rotz);
        rotEclose = zeros(1,rotz);
        
        [Itx, ~, ~] = transformModel(normI,height,Rad_par(i),0);
        It_sub = Itx(1:model_sub*sub_fac:end,:)';
        
        for j = 1:rotz
            Total_R = R_par(:,:,i);
            Total_T = T_par(:,:,i);
                        
            Zrot = (j-1)*Zrotation;
            %Rotate the model by Zrot degrees steps to help ICP
            R = [cosd(Zrot),-sind(Zrot),0;sind(Zrot),cosd(Zrot),0;0,0,1];
            T = [0;0;0];
            if it < 4
                TransL = 2*deltaR*(Rad_par(i)/ScanInitR)^2;
                T = [TransL;TransL;0]; %Shift scan offcenter during early passes
            end
            Stest = R*S_par' + repmat(T,1,size(S_par,1));
            Stest = Stest';
            Total_R = R*Total_R;
            Total_T = R*Total_T+T;
            %At first do regular ICP to get the scan near the model
            if it < 3
                [R, T, ~] = ICPreg(It_sub,Stest',ceil(icpIts/3),'Matching','kDtree');
                Stest = R*Stest' + repmat(T,1,size(Stest,1));
                Stest = Stest';
                Total_R = R*Total_R;
                Total_T = R*Total_T+T;
            end
            %Then move the scan inside the model as best it can
            [R, T, rotE(j), rotEin(j), rotEclose(j)] = ICPinsideFitting(It_sub,Stest',icpIts,icpBaseT,icpBaseT,icpBaseR,icpDist,closeWeight);
            rotTotal_R(:,:,j) = R*Total_R;
            rotTotal_T(:,:,j) = R*Total_T+T;
        end
        
        %For each radius, pick the parameters from the rotation with the best E
        [~,idx] = min(rotE);
        E_par(i) = rotE(idx);
        R_par(:,:,i) = rotTotal_R(:,:,idx);
        T_par(:,:,i) = rotTotal_T(:,:,idx);
        Ein_par(i) = rotEin(idx);
        close_par(i) = rotEclose(idx);
    end
            
    [~,idx] = min(E_par);
    %Since this is already weighted to prefer close points, it rejects
    %solutions that increase outliers significantly
    if min(E_par) < E_best && (Ein_par(idx)*closeWeight) > Ein_best
        E_best = min(E_par);
        if Ein_par(idx) > Ein_best
            Ein_best = Ein_par(idx);
        end
        %Save the best radius & E value; apply the transformation
        Parameters = [-1, -1, Rad_par(idx), -1, 0, E_best, -1];
        R = R_par(:,:,idx); T = T_par(:,:,idx);
        S = R*S' + repmat(T,1,size(S,1));
        S = S';

        %Recreate the model/scans for the graphic
        [It, ~, ~] = transformModel(normI,height,Rad_par(idx),0);
        titleStr = strcat('Coarse fitting of ',scan_name(1:end-4));
        [SvizI, SvizO, ~] = prepfig(It(1:model_sub*sub_fac:end,:)',S(1:scan_sub*sub_fac:end,:)',0.001);
        drawscatter2(titleStr,It(1:model_sub*sub_fac*2:end,:)',SvizO(1:2:end,:)',SvizI(1:3:end,:)',0);
        drawnow;

        fprintf('\n Best E = %.2f, R = %.3f; Close pts = %d / %d%s inside', E_best, Rad_par(idx), close_par(idx), ceil(100*Ein_par(idx)), '%');
    end
    icpBaseR = icpBaseR / 1.5;
    Zrotation = Zrotation / 4.4;
    rotz = rotz - 1;
    if rotz < 1
        rotz = 1;
    end
    
    %Average next height range including close / inside point maxima...
    [~,idx2] = max(close_par);
    [~,idx3] = max(Ein_par);
    Rbest = mean([Parameters(3),Rad_par(idx2),Rad_par(idx3)]);
    
    [Rlbound, Rubound] = setBounds(Rbest,stepDown*deltaR,Parameters(3));
    deltaR = (Rubound-Rlbound)/(cpucores-1);
end

duration = etime(clock, nextTime)/60; nextTime = clock;
fprintf('\nStage 1a done in %.1f mins: E = %.3f, H = %.4f, R = %.4f\n', duration, Parameters(6), ScanInitH, Parameters(3));

%--------------- Stage 1b: height fitting ------------------%
%Now that the model should be oriented, fit the height more carefully
E_best = Inf;
Ein_best = 0;
closeWeight = 1.10; %Allow some points to be outside the model during height test

icpIts = 12;
icpDist = 0.005;

Hubound = max(S(:,3))-min(S(:,3));
Hlbound = Hubound*0.80;
its = 5;
stepDown = 1.1; %Step down by this percentage of previous delta each iteration

fprintf('\nStage 1b, trying heights from %.3f to %.3f', Hlbound, Hubound);
for k = 1:its
    deltaH = (Hubound-Hlbound)/(cpucores-1);
    icpBaseT = 1.5*deltaH;
    icpBaseR = 2 - 1.5*k/its;
    
    H_par = zeros(1,cpucores);
    R_par = repmat(eye(3,3),[1,1,cpucores]);
    T_par = zeros(3,1,cpucores);
    close_par = zeros(1,cpucores);
    E_par = zeros(1,cpucores);
    Ein_par = zeros(1,cpucores);
    Stest = S(1:scan_sub:end,:);
    rad = Parameters(3);
    
    %Do most processing in parallel
    parfor i = 1:cpucores
        H_par(i) = Hlbound+(i-1)*deltaH;
        [It, ~, ~] = transformModel(normI,H_par(i),rad,0);
        [R, T, E_par(i), Ein_par(i), close_par(i)] = ICPinsideFitting(It(1:model_sub:end,:)',Stest',icpIts,icpBaseT/4,icpBaseT,icpBaseR,icpDist,closeWeight);
        R_par(:,:,i) = R*R_par(:,:,i);
        T_par(:,:,i) = R*T_par(:,:,i)+T;
    end
    
    [~,idx] = min(E_par);
    if min(E_par) < E_best && (Ein_par(idx)*closeWeight) > Ein_best
        E_best = min(E_par);
        if Ein_par(idx) > Ein_best
            Ein_best = Ein_par(idx);
        end
        %Save the parameters and apply the transformation to the scan
        Parameters = [H_par(idx), -1, rad, -1, 0, E_best, -1];
        R = R_par(:,:,idx); T = T_par(:,:,idx);
        S = R*S' + repmat(T,1,size(S,1));
        S = S';
        
        %Recreate the inside/outside points for the graphic
        [It, ~, ~] = transformModel(normI,H_par(idx),rad,0);
        [SvizI, SvizO, SvizCl] = prepfig(It(1:model_sub:end,:)',S(1:scan_sub:end,:)',icpDist);
        titleStr = strcat('Fitting ',scan_name(1:end-4),' with E = ',num2str(E_best,'%.3f'));
        drawscatter2(titleStr,It(1:model_sub*4:end,:)',SvizO(1:2:end,:)',SvizI(1:3:end,:)',SvizCl(1:2:end,:)');
        drawnow;
        
        fprintf('\n Best E = %.2f, H = %.4f (dH = %.4f); Close pts = %d / %d%s inside.', E_best, H_par(idx), deltaH, close_par(idx), ceil(100*Ein_par(idx)), '%');
    end
    
    %Average next height range including close / inside point maxima...
    [~,idx2] = max(close_par);
    [~,idx3] = max(Ein_par);
    Hbest = mean([Parameters(1),H_par(idx2),H_par(idx3)]);
    [Hlbound, Hubound] = setBounds(Hbest,stepDown*deltaH,Parameters(1));
end

duration = etime(clock, nextTime)/60;
fprintf('\nStage 1b done in %.1f mins: E = %.3f, H = %.4f, R = %.4f\n\nFinal pass:', duration, Parameters(6), Parameters(1), Parameters(3));

%----------------- Stage 2 ------------------%
%Fitting inside the model
%Sample more of the scan, and especially the model, for proximity tests
model_sub = ceil(model_sub/3);
scan_sub = ceil(scan_sub/2);

ht = Parameters(1); %Height is now fixed, so it won't be altered again
RUfirst = Rubound + 11*deltaR; %Expand the radius search to compensate for tapering
RLfirst = Rlbound - 11*deltaR;
firstTU = 3.5; %Setting the initial range of tapers for which to test
firstTL = -1.5;
TU = firstTU;
TL = firstTL;
deltaT = 1;
bestT = 0;
closestTaper = 0; %These variables will help select tapers for each iteration

E_best = Inf;
closeWeight = 1.05; %During ICP, weight more heavily toward insiders over proximity
icpIts = 4;
icpDist = 0.003; %Only very close points are counted
its = 3;
radits = ceil(0.5*ht/RLfirst); %Search for finer intervals later (start at 1 for normal cases; 2-3 for skinny cols)
stepDown = 1.1;

%Save some processing time by eliminating interior points
[It, ~, ~] = transformModel(normI,ht,Parameters(3),0);
kdOBJ = KDTreeSearcher(It);
[~, distances] = knnsearch(kdOBJ,S);
XYdist = sqrt(sum(S(:,1:2).^2,2));
if radits > 1
    distcutoff = 0.10; %Keep more points if it's tall & skinny
else
    distcutoff = 0.05;
end
exteriorIdx = distances < distcutoff & XYdist > (0.50*max(XYdist));
Sext = S(exteriorIdx,:);

subplot(1,1,1); %Reset to a single pane figure
%Final loops: Tapers and finer radius tests
for k = 1:its
    %Iterations for taper searching
    fprintf(' Trying tapers %.2f-%.2f, step %.2f, with radii %.3f-%.3f\n', TL, TU, deltaT, RLfirst, RUfirst);
    for taper = TL:deltaT:TU
        %Set the initial radius bounds for each taper
        RU = RUfirst;
        RL = RLfirst;
        for kk = 1:radits
            %Step down through increasingly narrow radii
            DR = (RU-RL)/(cpucores-1);
            icpBaseT = DR/2;
            icpBaseTZ = DR/20; %Allow only very small Z-translations, now that the height is fixed
            icpBaseR = 1.5 - kk/its;
            
            Rad_par = zeros(1,cpucores);
            R_par = repmat(eye(3,3),[1,1,cpucores]);
            T_par = zeros(3,1,cpucores);
            close_par = zeros(1,cpucores);
            E_par = zeros(1,cpucores);
            Ein_par = zeros(1,cpucores);
            Stest = Sext(1:scan_sub:end,:);
            
            %Do most processing in parallel
            parfor i = 1:cpucores
                Rad_par(i) = RL + (i-1)*DR;
                [It, ~, ~] = transformModel(normI,ht,Rad_par(i),taper);
                [R, T, E_par(i), Ein_par(i), close_par(i)] = ICPinsideFitting(It(1:model_sub:end,:)',Stest',icpIts,icpBaseT,icpBaseTZ,icpBaseR,icpDist,closeWeight);
                R_par(:,:,i) = R*R_par(:,:,i);
                T_par(:,:,i) = R*T_par(:,:,i)+T;
            end
            
            [~,idx] = min(E_par);
            if min(E_par) < E_best %Drop the minimum inside point requirement [&& (Ein*closeWeight) > Ein_best]
                E_best = min(E_par);
                %Apply transformation to the scan and its subsampled version
                R = R_par(:,:,idx); T = T_par(:,:,idx);
                S = R*S' + repmat(T,1,size(S,1));
                S = S';
                Sext = R*Sext' + repmat(T,1,size(Sext,1));
                Sext = Sext';
                maxht = max(S(:,3))-min(S(:,3)); %Max height of the reoriented scan
                                
                %Recreate the inside/outside points for the graphic
                [It, topR, botR] = transformModel(normI,ht,Rad_par(idx),taper);
                [SvizI, SvizO, SvizCl] = prepfig(It(1:model_sub:end,:)',S(1:scan_sub:end,:)',icpDist);
                
                titleStr = strcat('Fitting ',scan_name(1:end-4),' with E = ',num2str(E_best,'%.3f'));
                drawscatter1(titleStr,It(1:model_sub*6:end,:)',SvizO(1:2:end,:)',SvizI(1:3:end,:)',SvizCl(1:2:end,:)');
                filename = ['ModelFitting',scan_name(1:end-4)];
                savefig(filename);
                drawnow;
                
                Parameters = [Hbest, topR, Rad_par(idx), botR, taper, E_best, maxht];
                fprintf('  Best E = %.2f (dR = %.4f), MR = %.4f, T = %.2f; Close pts = %d / %d%s inside\n', E_best, DR, Rad_par(idx), taper, close_par(idx), ceil(100*Ein_par(idx)), '%');
            end
            
            %Set the next radius averaging the most close & inlier pts with the best E
            [~,idx2] = max(close_par);
            [~,idx3] = max(Ein_par);
            Rbest = mean([Parameters(3),Rad_par(idx2),Rad_par(idx3)]);
            [RL, RU] = setBounds(Rbest,stepDown*DR,Parameters(3));

            %For setting the next taper center, also save the values with
            %the closest points, and the most points inside the model
            if close > closestTaper
                closestTaper = close;
                bestT = taper;
            end
        end
    end
    %Adjust the starting radius bounds for every taper iteration
    RUfirst = Parameters(3) + 2*deltaR;
    RLfirst = Parameters(3) - 2*deltaR;
    radits = radits + 1;
    
    bestT = (bestT + Parameters(5))/2;
    [TL, TU] = setBounds(bestT,stepDown*deltaT,Parameters(5));
    deltaT = (TU-TL)/4;
end

output_filename = [scan_name(1:end-4),'centered',scan_name(end-3:end)];
writeobj(output_filename,S,Vc,objLines);

duration = etime(clock, startTime)/60;
fprintf('\nParameters for %s:\n E = %.3f: H = %.4f (max %.4f), TR = %.4f, BR = %.4f, T = %.2f', scan_name(1:end-4), Parameters(6), Parameters(1), Parameters(7), Parameters(2), Parameters(4), Parameters(5));
fprintf('\n Job completed in %.1f mins.\n\n', duration);
end

function [model,top_radius,bottom_radius] = transformModel(normModel, height, mid_radius, taper)

% Find the top and bottom radius given the height and tapering angle
top_radius = mid_radius - (tand(taper)*height)/2;
bottom_radius = mid_radius + (tand(taper)*height)/2;

% Apply the desired height, top radius, and bottom radius
model = normModel;
model(:,1:2) = model(:,1:2).*repmat((0.5-model(:,3))*bottom_radius + (0.5+model(:,3))*top_radius,1,2);
model(:,3) = model(:,3)*height;

end

function [normModel] = normalizeModel(model)
% Normalize the ideal model to have unit height, and unit top and bottom radius
top_points = model(abs(model(:,3)-max(model(:,3)))<0.01,:);
top_radius = max((max(top_points(:,1)) - min(top_points(:,1)))/2, (max(top_points(:,2)) - min(top_points(:,2)))/2);
bottom_points = model(abs(model(:,3)-min(model(:,3)))<0.01,:);
bottom_radius = max((max(bottom_points(:,1)) - min(bottom_points(:,1)))/2, (max(bottom_points(:,1)) - min(bottom_points(:,1)))/2);
height = max(model(:,3)) - min(model(:,3));
model(:,3) = model(:,3)/height;
model(:,1:2) = model(:,1:2)./repmat((0.5-model(:,3))*bottom_radius + (0.5+model(:,3))*top_radius,1,2);
normModel = model;
end

function [BL, BU] = setBounds(newCtr, span, overlapCtr)
    BL = newCtr - span; BU = newCtr + span; adjustment = 0;
    if BU < overlapCtr
        adjustment = overlapCtr - BU;
    elseif BL > overlapCtr
        adjustment = overlapCtr - BL;
    end
    BU = BU + 1.1*adjustment; BL = BL + 1.1*adjustment;
end

function [pInside,pOutside,pClose] = prepfig(q,p,CloseThreshold)
    kdOBJ = KDTreeSearcher(q');
    [match, distances] = knnsearch(kdOBJ,p');
    
    q_norm = sqrt(sum(q(:,match).^2,1));
    p_norm = sqrt(sum(p.^2,1));
    inliers = q_norm > p_norm;
    outliers = q_norm < p_norm;
    closePts = distances < CloseThreshold;
    
    pClose = p(:,closePts)';
    pInside = p(:,inliers)';
    pOutside = p(:,outliers)';
    %Remove the close points from inside & outside for clarity
    pInside = setdiff(pInside,pClose,'rows');
    pOutside = setdiff(pOutside,pClose,'rows');
end

function drawscatter1(titleStr,mod,out,in,close)
    legendStr = plotScatter(mod,out,in,close);
    view(2);
    axis equal;
    legend(legendStr,'Location','northwest');
    title(titleStr);
end

function drawscatter2(titleStr,mod,out,in,close)
    Xmax = 1.1*max(mod(1,:));
    Xmin = 1.1*min(mod(1,:));
    Ymax = 1.1*max(mod(2,:));
    Ymin = 1.1*min(mod(2,:));
    Zmax = 1.1*max(mod(3,:));
    Zmin = 1.1*min(mod(3,:));
    Xspan = Xmax-Xmin;
    Yspan = Ymax-Ymin;
    Zspan = Zmax-Zmin;

    Padding = 0.03;
    TopBorder = 0.05;
    BotH = (1.00-3*Padding-TopBorder) / (Yspan/Zspan+1);
    BotW = Xspan * BotH / Zspan;
    TopH = (1.00-3*Padding-TopBorder) - BotH;
    TopW = Xspan * TopH / Yspan;
    TopLpos = Padding + (1-Padding*2-TopW)/2;
    BotLpos = Padding + (1-Padding*2-BotW)/2;

    p1 = subplot('Position',[TopLpos, Padding*3+BotH, TopW, TopH]);
    legendStr = plotScatter(mod,out,in,close);
    view(2);
    legend(legendStr,'Location','westoutside');
    title(titleStr);
    axis equal;
    p1.XLim = [Xmin Xmax];
    p1.YLim = [Ymin Ymax];
    p1.ZLim = [Zmin Zmax];
    p1.XGrid = 'off';
    p1.YGrid = 'off';

    p2 = subplot('Position',[BotLpos, Padding, BotW, BotH]);
    plotScatter(mod,out,in,close);
    view([0 0]);
    axis equal;
    p2.XLim = [Xmin Xmax];
    p2.YLim = [Ymin Ymax];
    p2.ZLim = [Zmin Zmax];
    p2.XGrid = 'off';
    p2.YGrid = 'off';
end

function legendStr = plotScatter(model,outliers,inliers,closepts)
    legendidx = 1;
    legendStr = {};
    if model ~= 0
        scatter3(model(1,:)',model(2,:)',model(3,:)',1,[0 0 0.3],'filled');
        legendStr{legendidx} = 'Ideal Model';
        legendidx = legendidx + 1;
    end
    hold on;
    if outliers ~= 0
        scatter3(outliers(1,:)',outliers(2,:)',outliers(3,:)',3,'r','filled');
        legendStr{legendidx} = 'Outside';
        legendidx = legendidx + 1;
    end
    if inliers ~= 0
        scatter3(inliers(1,:)',inliers(2,:)',inliers(3,:)',1,[0.9 0.3 0.7],'filled');
        legendStr{legendidx} = 'Inside';
        legendidx = legendidx + 1;
    end
    if closepts ~= 0
        scatter3(closepts(1,:)',closepts(2,:)',closepts(3,:)',2,'g','filled');
        legendStr{legendidx} = 'Close';
    end
    hold off;
end