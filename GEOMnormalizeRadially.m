function [bucketPts,ZTrans,Height,hiRad] = GEOMnormalizeRadially(S,Nb)
    step = 360/Nb;
    
    %Convert the scan points into polar cooordinates
    Sangles = atan2d(-S(:,2),-S(:,1)); %Set X=-180 at X=1/Y=0 (not X=-1) 
    Sangles = Sangles+repmat(180,size(Sangles,1),1); %Make them all positive, so 270 = X=0/Y=-1
    Sradii = sqrt(sum(S(:,1:2).^2,2));
    Sheights = S(:,3);
    
    maxRad = max(Sradii);
    minH = min(Sheights);
    maxH = max(Sheights);
    
    %If H is very low, just sample the outermost 15% of points by radius
    %If H is very high (ie, a whole column shaft), keep most of the points
    radCutoff = 0.85 - 0.1*(maxH-minH)/maxRad;
    
    Snorm = horzcat(Sangles,Sradii,Sheights);
    Snorm = sortrows(Snorm,1);
    ExtPtIdx = Snorm(:,2) > radCutoff*maxRad & Snorm(:,3) ...
        > minH+0.01*(maxH-minH) & Snorm(:,3) < maxH-0.01*(maxH-minH);
    ExtPts = Snorm(ExtPtIdx,:);
    NpE = size(ExtPts,1);
    
    %Heights: taken from a sample of the points no more than 5% of H from the top and bottom
    newHmax = prctile(Snorm(Snorm(:,3) > maxH-0.05*(maxH-minH),3),99.5); %Clip the top %0.5 of the remaining points
    newHmin = prctile(Snorm(Snorm(:,3) < minH+0.05*(maxH-minH),3),0.5);
    newHctr = newHmax+newHmin;
    ZTrans = -newHctr/2;
    Height = newHmax-newHmin; %This will be slightly lower than scan point height
	hiRad = prctile(ExtPts(:,2),98); %Find 98% percentile radius of exterior points
    
    %Divvy up by buckets; shifted so first bucket is 1/2 step +/- 0 degrees
    BucketId = zeros(NpE,1);
    for i = 1:NpE
        ptAng = ExtPts(i,1);
        %for the first & last half-steps (at 0/360 deg)
        if ptAng/step < step/2 || ptAng/step >= Nb-step/2
            BucketId(i) = 1;
        else
            BucketId(i) = ceil(ptAng/step+step/2);
        end
    end
    bCounts = zeros(Nb,1);
    for j = 1:Nb
        bCounts(j) = sum(BucketId(:) == j);
    end
    bucketPts = cell(Nb,1);
    minHcutoff = 0.1*(maxH-minH);
    
    %This discards any buckets without enough data for further analysis.
    %All buckets must have at least 15 points, preserving a height at least
    %10% of the full height of the drum.
    for k = 1:Nb
        if bCounts(k) > 15
            bps = ExtPts((BucketId == k),1:3);
            %Shift the points on the Z-axis so they are centered by H
            bps = bps + repmat([0,0,ZTrans], size(bps,1), 1);
            if max(bps(:,3))-min(bps(:,3)) > minHcutoff
                bucketPts{k} = bps;
            end
        end
    end
end