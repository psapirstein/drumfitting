function [Parameters] = GEOMFitDrumByRadius(modelreg, modeldeep, modelflat, scan_name, scan_path, nFlutes)
startTime = clock;
[S,Vc,objLines]=readobj(scan_name,scan_path);

duration = etime(clock, startTime); nextTime = clock;
fprintf(' in %.1fs. Scan points = %d\n',duration, size(S,1));

%----------------- Geometric fitting ------------------%
%This assumes the model has already been approximately centered using the 
%inside-ICP fitting, or manually. Usually only 2 iterations are needed to
%refine the centering and rotation of the scan, but more are implemented here
%both to ensure the solution is stable, and also because on the third step,
%the scan is rotated so the thickest flute is on the Y-axis.
%Note that the scan must be rotated so that the arrises align with the
%X axis (so in top view, there are arrises to the right and left).

its = 10;
for it = 1:its
    if it > 1 %Once it has been estimated, apply the transformation
        Rx = [1, 0, 0; 0, cosd(XaxRot), -sind(XaxRot); 0, sind(XaxRot), cosd(XaxRot)];
        Ry = [cosd(YaxRot), 0, sind(YaxRot); 0, 1, 0; -sind(YaxRot), 0, cosd(YaxRot)];
        Rz = [cosd(ZaxRot), -sind(ZaxRot), 0; sind(ZaxRot), cosd(ZaxRot), 0; 0, 0, 1];
        S = Rx*Ry*Rz*S' + repmat([XTrans;YTrans;ZTrans],1,size(S,1));
        S = S';
    end
    
    %More visual feedback is off by default, because this slows processing.
    %By default, scans are divided into 480 radial wedges; this may need to be
    %varied depending on the resolution of the scan.
    withFigure = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nb = 480;
    
    %First, convert the scan into radial coordinates, with points divided into buckets
    [bPs, ZTrans, Height, hiRad] = GEOMnormalizeRadially(S,Nb);
    if it == its
        withFigure = 1; clf; hold on;
    end
    %Examine the sides, trying to filter out points in damaged areas
    distIn = hiRad / 500;
    [coeffs, Tmedian, Rmedian, figUnwrapped] = GEOMgetfilteredCoeffs(Height,bPs,Nb,distIn,withFigure);
    if it == its %Save the figure from the last iteration
        figure(figUnwrapped);
        filename = [scan_name(1:end-4),'-Unwrapped.fig'];
        savefig(filename);
    end
    %Search for local maxima and minima, identifying arrises and flutes
    withFigure = 1; clf; hold on;
    [arrismaxP, fluteminP, ZaxRot, figFlute] = GEOMmeasureFlutes(nFlutes,Nb,coeffs,withFigure);
    
    RmedianArris = median(nonzeros(arrismaxP(:,2)));
    RmedianFlute = median(nonzeros(fluteminP(:,2)));
    RhiArris = prctile(nonzeros(arrismaxP(:,2)),98);
    RloArris = prctile(nonzeros(arrismaxP(:,2)),5); %Assume more low arrises, prone to breaks
    RhiFlute = prctile(nonzeros(fluteminP(:,2)),98);
    RloFlute = prctile(nonzeros(fluteminP(:,2)),2);
    
    if withFigure
        for div = 1:nFlutes
            if arrismaxP(div,2) > RhiArris
                scatter(arrismaxP(div,1),arrismaxP(div,2),80,'mo');
            end
            if arrismaxP(div,2) < RloArris
                scatter(arrismaxP(div,1),arrismaxP(div,2),80,'go');
            end
            if fluteminP(div,2) > RhiFlute
                scatter(fluteminP(div,1),fluteminP(div,2),80,'mo');
            end
            if fluteminP(div,2) < RloFlute
                scatter(fluteminP(div,1),fluteminP(div,2),80,'go');
            end
        end
    end
    arrismaxP(arrismaxP(:,2) > RhiArris,2) = 0; %Delete outliers
    arrismaxP(arrismaxP(:,2) < RloArris,2) = 0;
    fluteminP(fluteminP(:,2) > RhiFlute,2) = 0;
    fluteminP(fluteminP(:,2) < RloFlute,2) = 0;

    if withFigure
        x1 = linspace(1, Nb, 2);
        y1 = RmedianArris+x1*0;
        y2 = RmedianFlute+x1*0;
        plot(x1,y1,'g:','LineWidth',0.05);
        plot(x1,y2,'g:','LineWidth',0.05);
        title(strcat('Analytical section of ',scan_name(1:end-4)));
        
        hold off;
        axis manual;
        set(gca,'xlim',[0 360]);
        set(gca,'ylim',[Rmedian-0.1 Rmedian+0.1]);
        drawnow;
        
        if it == its %Save the figure from the last iteration
            figure(figFlute);
            filename = [scan_name(1:end-4),'-MidRadSection.fig'];
            savefig(filename);
        end
    end
    
    %Estimate the rotation to homogenize the taper, and the translation to
    %homogenize the radii of the arrises and flute centers.
    XaxRot = 0; YaxRot = 0;
    XTrans = 0; YTrans = 0;
    divs = nFlutes*2; %Now sample taper angles to check if scan is tilted
    for d2 = 1:divs %Sample two slopes per arris, one on either side
        lo_idx = 1+floor(Nb*(d2-1)/divs);
        hi_idx = (ceil(Nb*d2/divs));
        slopes = coeffs(lo_idx:1:hi_idx,1);
        if sum(slopes) ~= 0
            theta = median(atand(nonzeros(slopes)))-Tmedian;
            phi = 360*(d2-0.5)/divs; %Angle at the middle of this group of tapers
            XaxRot = XaxRot + 2*sind(phi)*theta/divs; %Separate X/Y directions of tapers around drum
            YaxRot = YaxRot - 2*cosd(phi)*theta/divs;
        end

        Arrisfac = 2/divs; %Weight fitting toward the flutes
        Flutefac = 4/divs;
        phi = 360*(d2-1)/divs;
        if mod(d2,2)
            if arrismaxP(ceil(d2/2),2) > 0
                Roffset = arrismaxP(ceil(d2/2),2)-RmedianArris;
                XTrans = XTrans - Arrisfac*cosd(phi)*Roffset;
                YTrans = YTrans - Arrisfac*sind(phi)*Roffset;
            end
        else
            if fluteminP(ceil(d2/2),2) > 0
                Roffset = fluteminP(ceil(d2/2),2)-RmedianFlute;
                XTrans = XTrans - Flutefac*cosd(phi)*Roffset;
                YTrans = YTrans - Flutefac*sind(phi)*Roffset;
            end
        end
    end
    if abs(XTrans) + abs(YTrans) > 0.03 %Discard any single translation greater than 3 cm
        fprintf('Warning: forced to discard a bad translation due to inaccurately marked points.\n');
        XTrans = 0; YTrans = 0;
    end
    fprintf('Xtilt / Ytilt / Zrot = %.3f / %.3f / %.3f; Xtrans / Ytrans / Ztrans = %.4f / %.4f / %.4f\n',XaxRot,YaxRot,ZaxRot,XTrans,YTrans,ZTrans);
    
    %Estimate the depth of the fluting
    fluteDepths = zeros(nFlutes,1);
    for div = 1:nFlutes
        if fluteminP(div,2) > 0
            nextdiv = div + 1;
            if nextdiv > nFlutes
                nextdiv = 1;
            end
            if arrismaxP(div,2) > 0 || arrismaxP(nextdiv,2) > 0 %Measure depth relative to nearby arris
                if arrismaxP(div,2) > 0
                    Rlocal = arrismaxP(div,2);
                    if arrismaxP(nextdiv,2) > 0 %Average if both arrises are intact
                        Rlocal = (Rlocal + arrismaxP(nextdiv,2))/2;
                    end
                else
                    Rlocal = arrismaxP(nextdiv,2);
                end
            else %No arrises measured nearby
                Rlocal = RmedianArris;
            end
            fluteDepths(div) = Rlocal - fluteminP(div,2);
        end
    end
    estOuterR = (RmedianArris+RhiArris*2)/3; %Put this closer to the maximum than the median
    avChordFluteD = estOuterR*cosd(180/nFlutes)-RmedianFlute;
    
    if it == 3 && sum(fluteminP(:,2) > 0) > 8 %Reorient the model if there are at least 8 flutes identified
        [~,maxFidx] = max(fluteminP(:,2));
        Zmax = (maxFidx-1)*360/nFlutes;
        ZaxRot = ZaxRot + 90 - Zmax; %Rotate so the thickest flute is at the top
    end
end

duration = etime(clock, nextTime)/60;
fprintf('Fitting done in %.1f mins. Est. outer radius = %.2f cm; av. flute depth (radial / chord) = %.2f / %.2f cm\n', ...
    duration, estOuterR*100,median(nonzeros(fluteDepths))*100,avChordFluteD*100);

%Measure the inside/outside points for the final graphic
modelfile = modelreg;
if avChordFluteD < 0.008
    modelfile = modelflat;
elseif avChordFluteD > 0.018
    modelfile = modeldeep;
end
[VI,~,~]=readobj(modelfile,'');
modelI = VI - repmat(mean(VI,1),size(VI,1),1);
icpDist = 0.003; %Lower the threshold for counting points as close to the model
normI = normalizeModel(modelI);

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

%Resize the ideal model to the estimated parameters, and estimate the error
%between the fitted scan to the model. Then display and save a figure.
[It, topR, botR] = transformModel(normI,Height,estOuterR,-Tmedian);
[SvizI, SvizO, SvizCl] = prepfig(It',S',icpDist);
E = (length(S(:,1))-length(SvizI(:,1)))/length(SvizCl(:,1));

clf; hold on;
titleStr = strcat('Fitting ',scan_name(1:end-4),' with E = ',num2str(E,'%.2f'));
drawscatter1(titleStr,It(1:model_sub*3:end,:)',SvizO(1:scan_sub*2:end,:)',SvizI(1:scan_sub*3:end,:)',SvizCl(1:scan_sub:end,:)');
filename = [scan_name(1:end-4),'-Fitting.fig'];
savefig(filename);
hold off;
drawnow;

%Save the estimated parameters of the full drum, and a table of values for
%the estimated values of the individual flutes and arrises. Save a centered
%copy of the scan with only the vertices updated to the new orientation.
fluteFilename = [scan_name(1:end-4),'-Fluting.csv'];
fFluting = fopen(fluteFilename,'w');
fprintf(fFluting, 'Count, Arris Ang. (deg), Arris Rad. (cm), Flute Ang., Flute Rad., Radial Flute depth\n');
for div = 1:nFlutes
    fprintf(fFluting, '%d, %.2f, %.2f, %.2f, %.2f, %.2f\n', div, arrismaxP(div,1), 100*arrismaxP(div,2), fluteminP(div,1), 100*fluteminP(div,2), 100*fluteDepths(div));
end
fclose(fFluting);

Parameters = [Height,topR,estOuterR,botR,-Tmedian,avChordFluteD,E];
output_filename = [scan_name(1:end-4),'centered',scan_name(end-3:end)];
writeobj(output_filename,S,Vc,objLines);

duration = etime(clock, startTime)/60;
fprintf('\nJob completed in %.1f mins.\n\n', duration);
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

function [model,top_radius,bottom_radius] = transformModel(normModel, height, mid_radius, taper)
% Find the top and bottom radius given the height and tapering angle
top_radius = mid_radius - (tand(taper)*height)/2;
bottom_radius = mid_radius + (tand(taper)*height)/2;

% Apply the desired height, top radius, and bottom radius
model = normModel;
model(:,1:2) = model(:,1:2).*repmat((0.5-model(:,3))*bottom_radius + (0.5+model(:,3))*top_radius,1,2);
model(:,3) = model(:,3)*height;
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

function drawscatter1(titleStr,model,outliers,inliers,closepts)
    if model ~= 0
        scatter3(model(1,:)',model(2,:)',model(3,:)',1,[0 0 0.3],'filled','DisplayName','Ideal Model');
        legend('-DynamicLegend');
    end
    hold on;
    if outliers ~= 0
        scatter3(outliers(1,:)',outliers(2,:)',outliers(3,:)',3,'r','filled','DisplayName','Outside');
        legend('-DynamicLegend');
    end
    if inliers ~= 0
        scatter3(inliers(1,:)',inliers(2,:)',inliers(3,:)',1,[0.9 0.3 0.7],'filled','DisplayName','Inside');
        legend('-DynamicLegend');
    end
    if closepts ~= 0
        scatter3(closepts(1,:)',closepts(2,:)',closepts(3,:)',2,'g','filled','DisplayName','Close');
        legend('-DynamicLegend');
    end
    hold off;
    
    view(2);
    axis equal;
    legend('Location','northwest');
    title(titleStr);
end