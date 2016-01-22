function [arrismaxP, fluteminP, ZaxRot, fig] = GEOMmeasureFlutes(nFlutes,Nb,coeffs,withFigure)
    fig = 0;
    step = 360/Nb;
    breadth = 4; %Search +/- breadth degrees per arris
    
    [c1Pt, c2Av] = signedCurvature(1:Nb,coeffs(:,2));
    
    cuspIdx = c2Av > prctile(c2Av(c2Av > 0),25) & c1Pt > 0; %Filtered set for cusps
    cminIdx = c1Pt ~= 0 & c1Pt < prctile(c1Pt(c1Pt > 0),10) & c1Pt > prctile(c1Pt(c1Pt < 0),10); %Unfiltered set for low points
    
    offsets = repmat(-10,nFlutes,1); %Stored as degrees, -10 = nil
    cuspPts = zeros(nFlutes,2);
    for div = 1:nFlutes %Search for pattern of arrises, starting w/1 (far right, 0 deg)
        ab = ceil(1+(div-1)*360/nFlutes/step); %Expected bucket for each arris
        aw = ceil(breadth/step);
        rangeAb = ab-aw:ab+aw;
        rangeAb(rangeAb < 1) = rangeAb(rangeAb < 1) + Nb;
        rangeAb(rangeAb > Nb) = rangeAb(rangeAb > Nb) - Nb;
        
        localcuspIdx = cuspIdx(rangeAb);
        if sum(localcuspIdx) > 0 %Evaluate if potential cusps were found
            maxR = max(coeffs(rangeAb,2)); %Radius is just the highest in the area
            [maxRidx,~] = find(coeffs(:,2) == maxR);
            c2tmp = c2Av(rangeAb); %Get the C2 curvatures for this range
            localC2s = c2tmp(localcuspIdx); %Get curvatures at the cusps
            [~,hiCidx] = ismember(localC2s(localC2s >= mean(localC2s)),c2Av); %And the indices of the taller ones
            hiCidx(hiCidx > Nb-aw) = hiCidx(hiCidx > Nb-aw) - Nb; %Deal with the first flute, which spans high and low indices
            degree = step*(sum(hiCidx)/length(hiCidx)-1); %Convert to degrees

            %This will fail at breaks, so check for discontinuities
            rangeLo = floor(sum(hiCidx)/length(hiCidx));
            rangeHi = ceil(sum(hiCidx)/length(hiCidx));
            rangeLo = rangeLo-4:rangeLo-1;
            rangeHi = rangeHi+1:rangeHi+4;
            rangeLo(rangeLo < 1) = rangeLo(rangeLo < 1) + Nb;
            rangeHi(rangeHi > Nb) = rangeHi(rangeHi > Nb) - Nb;
            rangeHi(rangeHi < 1) = rangeHi(rangeHi < 1) + Nb;
            if min(abs(c1Pt(rangeLo))) == 0 || min(abs(c1Pt(rangeHi))) == 0
                if div == 1 && maxRidx > Nb-aw
                    maxRidx = maxRidx - Nb;
                end
                degree = (degree + 2*(maxRidx-1)*step)/3;
            end
            offsets(div) = degree - (div-1)*360/nFlutes;
            if degree < 0
                degree = degree + 360;
            end
            cuspPts(div,:) = [degree,maxR];
        end
    end

    if withFigure
        fig = gcf;
        xlim = linspace(0,359,Nb);
        
        slopePoints = scatter(xlim(coeffs(:,1) ~= 0),coeffs(coeffs(:,1) ~=0,1)+0.5, 3,'g','filled');        
        sectionPoints = scatter(xlim(coeffs(:,2) ~= 0),coeffs(coeffs(:,2) ~= 0,2), 6,'r','filled');
    end
    
    %In the breaks, some cusps should be removed that are unlikely to be
    %arrises. First are those that are very low, and then those that are
    %not aligned well with the others
    offset = sum(offsets(offsets > -10))/sum(offsets > -10);
    cutLowCusps = 0;
    RarrisMin = 3*median(nonzeros(cuspPts(:,2))) - 2*prctile(nonzeros(cuspPts(:,2)),95);
    for div = 1:nFlutes
        cuspOffset = offsets(div) - offset;
        deleteit = 0;
        if cuspOffset > 3 %Delete cusps significantly offset from the others
            fprintf('Cut cusp at arris %d due to offset aberrant by %.2f.\n',div,cuspOffset);
            deleteit = 1;
        elseif cuspPts(div,2) < RarrisMin %Delete very low cusps
            cutLowCusps = cutLowCusps + 1;
            deleteit = 1;
        end
        if deleteit
            if withFigure
                cuspPoint = scatter(cuspPts(div,1),cuspPts(div,2),30,'rd');
            end
            offsets(div) = -10;
            cuspPts(div,:) = [0, 0];
        end
    end
    if cutLowCusps > 0
        plural = ['','s'];
        fprintf('Cut %d cusp%s due to low radius.\n',cutLowCusps,plural(cutLowCusps > 1));
    end
    offset = sum(offsets(offsets > -10))/sum(offsets > -10);
    ZaxRot = -offset;
    
    %Now fit the flutes with parabolas
    flutePF = zeros(nFlutes,3);
    fluteminP = zeros(nFlutes,2);
    fittedFluteLine = 0;
    
    lastflutePts = []; %Need to save this for next loop
    for div = 1:nFlutes
        fo = offset;
        if offsets(div) > -10 %Average the general offset with local offset at this cusp
            fo = (fo + offsets(div) - offset)/2;
        end
        nextdiv = div+1;
        if nextdiv > nFlutes
            nextdiv = 1;
        end
        if offsets(nextdiv) > -10
            fo = (fo + offsets(nextdiv) - offset)/2;
        end
        fb = ceil(1+(fo+(div-0.5)*360/nFlutes)/step);
        fw = ceil((0.5*360/nFlutes-2)/step); %cut out last 2 degrees from sides
        if cuspPts(div,2) > 0 && fb-fw < 1+cuspPts(div,1)/step
            fb = ceil((fb + 1+(cuspPts(div,1) + 180/nFlutes)/step)/2);
            fw = floor((fw + fb - 1-cuspPts(div,1)/step)/2);
        end %But adjust the center if the range overlaps with a cusp
        if cuspPts(nextdiv,2) > 0 && fb+fw > 1+cuspPts(nextdiv,1)/step
            fb = floor((fb + 1+(cuspPts(nextdiv,1) - 180/nFlutes)/step)/2);
            fw = floor((fw + 1+cuspPts(nextdiv,1)/step - fb)/2);
        end

        rangeFb = fb-fw:fb+fw; %Make sure this doens't go out of bounds from the circle
        rangeFb(rangeFb < 1) = rangeFb(rangeFb < 1) + Nb;
        rangeFb(rangeFb > Nb) = rangeFb(rangeFb > Nb) - Nb;
        
        if isnan(rangeFb)
            flutePts = 0;
            rangeFb = 1:10;
        else
            flutePts = horzcat(rangeFb',coeffs(rangeFb,2)); %Use b for X units
            localC1s = cminIdx(rangeFb); %Index of points with good curvature
            flutePts = flutePts(localC1s,:); %Points left for fitting
        end
        if div == nFlutes
            lastflutePts = flutePts;
        end
        if length(flutePts) > 0.4*length(rangeFb) %Only evaluate when enough of the arc is there
            %Do 2nd order polynomial fit
            flutePF(div,:) = polyfit(flutePts(:,1),flutePts(:,2),2);
            if flutePF(div,1) > 0 %Only accept fits with concave flutes
                if withFigure
                    xLineP = linspace(fb-fw-3,fb+fw+3,50);
                    yP = polyval(flutePF(div,:),xLineP);
                    fittedFluteLine = plot((xLineP-1).*step,yP,'b:');
                end

                xMinP = -0.5*flutePF(div,2)/flutePF(div,1);
                yMinP = flutePF(div,1)*xMinP^2 + flutePF(div,2)*xMinP + flutePF(div,3);
                if abs(xMinP-fb) < fw/3 %These should be pretty close
                    fluteminP(div,:) = [(xMinP-1)*step, yMinP]; %Save Xmin as degrees
                    if withFigure
                        cuspPoint = scatter(fluteminP(div,1),yMinP,30,'kx');
                    end
                end
            else
                fprintf('Warning: discarded convex fitting of flute %d.\n',div);
            end
        end
    end
    %Go back, look for intersections between them, and average arrises
    %with previous analysis, including total offset
    arrismaxP = zeros(nFlutes,2);
    nOffsets = 0;
    finalOffset = 0;
    for div = 1:nFlutes
        adeg = (div-1)*360/nFlutes;
        if cuspPts(div,2) > 0 %Copy over the arrises from curvature analysis
            arrismaxP(div,:) = cuspPts(div,:);
        end

        if fluteminP(div,1) > 0
            finalOffset = finalOffset + fluteminP(div,1) - 180/nFlutes - adeg;
            nOffsets = nOffsets + 1;
        end

        lastdiv = div-1;
        if lastdiv < 1
            lastdiv = nFlutes;
        end
        if flutePF(div,3) > 0 && flutePF(lastdiv,3) > 0
            if lastdiv == nFlutes %Redo the polynomial fitting
                lastflutePts(:,1) = lastflutePts(:,1) - Nb; %Shift these to left of origin
                flutePF(lastdiv,:) = polyfit(lastflutePts(:,1),lastflutePts(:,2),2);
            end
            %estimate intersection of the two fitted parabolas
            aP = flutePF(div,1)-flutePF(lastdiv,1);
            bP = flutePF(div,2)-flutePF(lastdiv,2);
            cP = flutePF(div,3)-flutePF(lastdiv,3);
            sqR = sqrt(bP^2 - 4*aP*cP);
            xMaxP = 0.5*(-bP - sqR)/aP;
            xMaxP0 = 0.5*(-bP + sqR)/aP;
            if xMaxP0 < xMaxP && xMaxP0 > -0.5*flutePF(lastdiv,2)/flutePF(lastdiv,1);
                xMaxP = xMaxP0;
            end
            yMaxP = flutePF(div,1)*xMaxP^2 + flutePF(div,2)*xMaxP + flutePF(div,3);

            if isreal(xMaxP) && isreal(yMaxP)
                arrismaxP(div,:) = [(xMaxP-1)*step, yMaxP];
                if cuspPts(div,2) > 0 %Average with curvature analysis data
                    arrismaxP(div,:) = [(arrismaxP(div,1)*2+cuspPts(div,1))/3, ...
                        (arrismaxP(div,2)+cuspPts(div,2)*2)/3];
                end
                finalOffset = finalOffset + arrismaxP(div,1) - adeg;
                nOffsets = nOffsets + 1;
            else
                fprintf('Warning: discarded non-intersecting estimate at arris %d\n',div);
            end

            if lastdiv == nFlutes %Return things to how they were
                lastflutePts(:,1) = lastflutePts(:,1) + Nb;
                flutePF(lastdiv,:) = polyfit(lastflutePts(:,1),lastflutePts(:,2),2);
            end
        end

        if withFigure && arrismaxP(div,2) > 0
            cuspPoint = scatter(arrismaxP(div,1),arrismaxP(div,2),30,'kx');
        end
    end
    
    if withFigure
        if fittedFluteLine > 0
            legend([slopePoints, sectionPoints, cuspPoint,fittedFluteLine],'Slopes','Midsection','Cusp','Fitted flute');
        else
            legend([slopePoints, sectionPoints, cuspPoint],'Slopes','Midsection','Cusp');
        end
    end
    
    if nOffsets > 0
        ZaxRot = -finalOffset/nOffsets;
        if abs(ZaxRot) > 360/nFlutes/5 %Any rotation greater than 1/5 of a flute is discarded as suspect
            fprintf('Warning: forced to discard a bad Z rotation due to inaccurately marked points.\n');
            ZaxRot = 0;
        end
    end
end

function [curv1Pt, curv2Av] = signedCurvature(xPts,yPts)
    Np = length(xPts);
    curv1Pt = zeros(Np,1);
    curv2Av = zeros(Np,1);
    
    for b = 1:Np
        LoR = b-4:b-2; MidR = b-1:b+1; HiR = b+2:b+4;
        LoR(LoR < 1) = LoR(LoR < 1) + Np;
        MidR(MidR < 1) = MidR(MidR < 1) + Np;
        MidR(MidR > Np) = MidR(MidR > Np) - Np;
        HiR(HiR > Np) = HiR(HiR > Np) - Np;

        xAv2s = [xPts(LoR(2)); xPts(MidR(2)); xPts(HiR(2))];
        yAv2s = [mean(yPts(LoR)); mean(yPts(MidR)); mean(yPts(HiR))];

        xPt1s = xPts(MidR);
        yPt1s = yPts(MidR);
        
        if min(yPt1s) > 0
            curv1Pt(b) = -2*((xPt1s(2)-xPt1s(1)).*(yPt1s(3)-yPt1s(1))-(xPt1s(3)-xPt1s(1)).*(yPt1s(2)-yPt1s(1))) ./ ...
            sqrt(((xPt1s(2)-xPt1s(1)).^2+(yPt1s(2)-yPt1s(1)).^2)*((xPt1s(3)-xPt1s(1)).^2+(yPt1s(3)-yPt1s(1)).^2)*((xPt1s(3)-xPt1s(2)).^2+(yPt1s(3)-yPt1s(2)).^2));
        end
        
        if min(yAv2s) > 0
            curv2Av(b) = -2*((xAv2s(2)-xAv2s(1)).*(yAv2s(3)-yAv2s(1))-(xAv2s(3)-xAv2s(1)).*(yAv2s(2)-yAv2s(1))) ./ ...
            sqrt(((xAv2s(2)-xAv2s(1)).^2+(yAv2s(2)-yAv2s(1)).^2)*((xAv2s(3)-xAv2s(1)).^2+(yAv2s(3)-yAv2s(1)).^2)*((xAv2s(3)-xAv2s(2)).^2+(yAv2s(3)-yAv2s(2)).^2));
        end
    end
    curv1Pt = curv1Pt./max(curv1Pt); %Set max to 1; <0 is concave/no data
    curv2Av = curv2Av./max(curv2Av);
end
