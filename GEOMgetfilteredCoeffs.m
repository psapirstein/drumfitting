function [coeffs,Tmedian,Rmedian,fig] = GEOMgetfilteredCoeffs(Ht,bucketPts,Nb,distIn,withFigure)
    fig = 0;
    minHcutoff = 0.1*Ht;
    minH = -0.5*Ht;
    maxH = 0.5*Ht;
    if withFigure
        fig = gcf;
        hold off;
    end
    
    TPts_cut = 0;
    TPts = size(vertcat(bucketPts{:}),1);
    max_its = 10;
    for its = 1:max_its
        coeffs = zeros(Nb,2);
        weights = zeros(Nb,1);
        %Perc refines the filtering as more points are removed
        perc = 35*its/max_its;
        for k = 1:Nb
            if ~isempty(bucketPts{k,:,:})
                bps = [bucketPts{k,:,:}];
                [coeffs(k,1), coeffs(k,2), ~, ~, in_idx] = ...
                    GEOMransacLineFitting(bps(:,3),bps(:,2), distIn, 25+perc, round(75+perc*6));
                inHspan = max(bps(in_idx,3)) - min(bps(in_idx,3));
                %fprintf('\ndebug: c1/2: %.3f / %.3f; span %.1f cm H; %d inliers', coeffs(k,1), coeffs(k,2), inHspan*100, length(in_idx));
                if length(in_idx) > 6 && inHspan > minHcutoff
                    weights(k) = length(in_idx) * inHspan;
                end
            end
        end
        %Determine which slices to include, weighted by the total height span 
        %of the fitted points, the number of fitted points, and the stage
        %of the filtering.
        wtcs = coeffs(weights > prctile(weights,80-perc*2),:);
        
        %Get the parameters for tapers and radius for the longer groups
        angles = atand(nonzeros(wtcs(:,1)));
        Tmedian = median(angles);
        Tlow = prctile(angles, 30-perc/1.5);
        Thigh = prctile(angles, 70+perc/1.5);
        Rmedian = median(nonzeros(wtcs(:,2)));
        
        %After angle adjustments, discard points inside the fitted lines
        ItPts_cut = 0;
        all_cutbs = zeros(0,2);
        for k = 1:Nb
            if coeffs(k,1) ~= 0
                adjusted_cs = 0;
                ang = atand(coeffs(k,1));
                newbs = [bucketPts{k,:,:}];
                %Refit the line to the points if the results of the RANSAC
                %estimate for angle are doubtful
                if ang < Tlow || ang > Thigh
                    if ang < Tlow
                        ang = tand(Tlow);
                    else
                        ang = tand(Thigh);
                    end
                    avR = mean(newbs(:,2));
                    avH = mean(newbs(:,3));
                    newRmid = avR - ang*avH;
                    coeffs(k,1) = ang;
                    coeffs(k,2) = newRmid;
                    adjusted_cs = 1;
                end
                
                %Delete points left of the line by increasingly small amounts
                pdists = DistPtLineSigned(newbs(:,2:3),coeffs(k,1),coeffs(k,2));
                Rdiff = distIn * (max_its/its)^1.5 / 2;
                cutbs = newbs(pdists < -Rdiff,:);
                newbs = newbs(pdists >= -Rdiff,:);
                if size(newbs,1) < 8 || max(newbs(:,3))-min(newbs(:,3)) < minHcutoff
                    cutbs = [bucketPts{k,:,:}];
                    bucketPts{k} = []; %Discard if few / narrow point range
                    newbs = 0;
                else
                    bucketPts{k,:} = newbs;
                end
                ItPts_cut = ItPts_cut + size(cutbs,1);
                
                if withFigure
                    Xshift = 7.2*k/Nb;
                    if cutbs ~= 0
                        cuttmp = cutbs(:,2:3) + repmat([Xshift 0],size(cutbs,1),1);
                        all_cutbs = vertcat(all_cutbs, cuttmp);
                    end
                    if its == max_its && size(newbs,1) > 9
                        ptcolor = [(k/Nb) 0 (Nb-k)/Nb];
                        linecolor = 'm';
                        if adjusted_cs > 0
                            linecolor = 'g';
                        end
                        hold on;
                        shiftedPts = newbs(:,2:3) + repmat([Xshift 0],size(newbs,1),1);
                        scatter(shiftedPts(1:8:end,1),shiftedPts(1:8:end,2),3,ptcolor,'filled');
                        
                        Hfit = linspace(minH, maxH, 2);
                        Rfit = polyval([coeffs(k,1), coeffs(k,2)+Xshift], Hfit);
                        plot(Rfit,Hfit,linecolor,'LineWidth',0.1);
                    end
                end
            end
        end
        TPts_cut = TPts_cut + ItPts_cut;
        %fprintf('Pass %d: Removed %d points as breaks.\n',its,ItPts_cut);
        
        if withFigure
            hold on;
            scatter(all_cutbs(1:25:end,1),all_cutbs(1:25:end,2),3,'c','filled');
            xlim([0 Xshift+2*Rmedian]);
            ylim([minH maxH]);
            drawnow;
        end
    end
    
    fprintf('Filtered a total of %d points of %d as breaks.\n', TPts_cut, TPts); 
end

function [ds] = DistPtLineSigned(pts,m,b)
    V1 = [b 0];
    V2 = [0 -b/m];
    ds = zeros(size(pts,1),1);
    for i = 1:size(pts,1)
        d = det([[pts(i,1), pts(i,2)]-V1; V2-V1]) / norm(V2-V1);
        ds(i) = -d*abs(m)/m; %Points left of the line by X will be negative
    end
end