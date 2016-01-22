function [TR, TT, ER, outliers, t] = ICPreg(q,p,varargin)
% Perform the Iterative Closest Point algorithm on three dimensional point
% clouds.
%
% [TR, TT] = icp(q,p)   returns the rotation matrix TR and translation
% vector TT that minimizes the distances from (TR * p + TT) to q.
% p is a 3xm matrix and q is a 3xn matrix.
%
% [TR, TT] = icp(q,p,k)   forces the algorithm to make k iterations
% exactly. The default is 10 iterations.
%
% [TR, TT, ER] = icp(q,p,k)   also returns the RMS of errors for k
% iterations in a (k+1)x1 vector. ER(0) is the initial error.
%
% [TR, TT, ER, t] = icp(q,p,k)   also returns the calculation times per
% iteration in a (k+1)x1 vector. t(0) is the time consumed for preprocessing.
%
% EdgeRejection
%       {false} | true
%       If EdgeRejection is true, point matches to edge vertices of q are
%       ignored. Requires that boundary points of q are specified using
%       Boundary or that a triangulation matrix for q is provided.
%
% Matching
%       {bruteForce} | Delaunay | kDtree
%       Specifies how point matching should be done.
%       bruteForce is usually the slowest and kDtree is the fastest.
%       Note that the kDtree option is depends on the Statistics Toolbox
%       v. 7.3 or higher.
%
% Weight
%       {@(match)ones(1,m)} | Function handle
%       For point or plane minimization, a function handle to a weighting
%       function can be provided. The weighting function will be called
%       with one argument, a 1xm vector that specifies point pairs by
%       indexing into q. The weighting function should return a 1xm vector
%       of weights for every point pair.
%
% WorstRejection
%       {0} | scalar in ]0; 1[
%       Reject a given percentage of the worst point pairs, based on their
%       Euclidean distance.
%
% Martin Kjer and Jakob Wilm, Technical University of Denmark, 2012

% Use the inputParser class to validate input arguments.
inp = inputParser;

inp.addRequired('q', @(x)isreal(x) && size(x,1) >= 3);
inp.addRequired('p', @(x)isreal(x) && size(x,1) >= 3);

inp.addOptional('iter', 10, @(x)x > 0 && x < 10^5);

inp.addParameter('EdgeRejection', false, @(x)islogical(x));

validMatching = {'bruteForce','Delaunay','kDtree'};
inp.addParameter('Matching', 'bruteForce', @(x)any(strcmpi(x,validMatching)));

inp.addParameter('Weight', @(x)ones(1,length(x)));

inp.addParameter('WorstRejection', 0, @(x)isscalar(x) && x > 0 && x < 1);

inp.parse(q,p,varargin{:});
arg = inp.Results;
clear('inp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual implementation

% Allocate vector for RMS of errors in every iteration.
t = zeros(arg.iter+1,1);

% Start timer
tic;

Np = size(p,2);

% Transformed data point cloud
pt = p;

% Allocate vector for RMS of errors in every iteration.
ER = zeros(arg.iter+1,1);

% Initialize temporary transform vector and matrix.
T = zeros(3,1);
R = eye(3,3);

% Initialize total transform vector(s) and rotation matric(es).
TT = zeros(3,1, arg.iter+1);
TR = repmat(eye(3,3), [1,1, arg.iter+1]);

% If Matching == 'kDtree', a kD tree should be built (req. Stat. TB >= 7.3)
if strcmp(arg.Matching, 'kDtree')
    kdOBJ = KDTreeSearcher(transpose(q));
end

t(1) = toc;

outliers = [];

% Go into main iteration loop
for k=1:arg.iter
    
    % Do matching
    switch arg.Matching
        case 'bruteForce'
            [match, mindist] = match_bruteForce(q,pt);
        case 'Delaunay'
            [match, mindist] = match_Delaunay(q,pt,DT);
        case 'kDtree'
            [match, mindist] = match_kDtree(q,pt,kdOBJ);
    end
    p_idx = true(1, Np);
    q_idx = match;
    
    
    % If worst matches should be rejected
    if arg.WorstRejection
        edge = round((1-arg.WorstRejection)*sum(p_idx));
        pairs = find(p_idx);
        [~, idx] = sort(mindist);
        outliers = pairs(idx(edge:end));
        p_idx(outliers) = false;
        q_idx = match(p_idx);
        mindist = mindist(p_idx);
    end
    
    if k == 1
        ER(k) = sqrt(sum(mindist.^2)/length(mindist));
    end
    
    [R,T] = eq_point(q(1:3,q_idx),pt(1:3,p_idx), ones(1,length(p_idx)));
%{
    scatter3(pt(1,p_idx)',pt(2,p_idx)',pt(3,p_idx)',2,'b','filled');
    hold on;
    scatter3(q(1,q_idx)',q(2,q_idx)',q(3,q_idx)',2,'r','filled');
    axis equal;
    hold off;
    view(2);
    drawnow;
%}
    % Add to the total transformation
    TR(:,:,k+1) = R*TR(:,:,k);
    TT(:,:,k+1) = R*TT(:,:,k)+T;
    
    % Apply last transformation
    pt(1:3,:) = TR(:,:,k+1) * p(1:3,:) + repmat(TT(:,:,k+1), 1, Np);
    
    % Root mean of objective function
    ER(k+1) = rms_error(q(:,q_idx), pt(:,p_idx));
    
    t(k+1) = toc;
end

TR = TR(:,:,end);
TT = TT(:,:,end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [match, mindist] = match_bruteForce(q, p)
m = size(p,2);
n = size(q,2);
match = zeros(1,m);
mindist = zeros(1,m);
for ki=1:m
    d=zeros(1,n);
    for ti=1:3
        d=d+(q(ti,:)-p(ti,ki)).^2;
    end
    [mindist(ki),match(ki)]=min(d);
end

mindist = sqrt(mindist);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [match, mindist] = match_Delaunay(q, p, DT)
match = transpose(nearestNeighbor(DT, transpose(p)));
mindist = sqrt(sum((p-q(:,match)).^2,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [match, mindist] = match_kDtree(~, p, kdOBJ)
[match, mindist] = knnsearch(kdOBJ,transpose(p));
match = transpose(match);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R,T] = eq_point(q,p,weights)

m = size(p,2);
n = size(q,2);

% normalize weights
weights = weights ./ sum(weights);

% find data centroid and deviations from centroid
q_bar = q * transpose(weights);
q_mark = q - repmat(q_bar, 1, n);
% Apply weights
q_mark = q_mark .* repmat(weights, 3, 1);

% find data centroid and deviations from centroid
p_bar = p * transpose(weights);
p_mark = p - repmat(p_bar, 1, m);
% Apply weights
%p_mark = p_mark .* repmat(weights, 3, 1);

N = p_mark*transpose(q_mark); % taking points of q in matched order

[U,~,V] = svd(N); % singular value decomposition

R = V*diag([1 1 det(U*V')])*transpose(U);

T = q_bar - R*p_bar;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ER = rms_error(p1,p2)
dsq = sum(power(p1 - p2, 2),1);
ER = sqrt(mean(dsq));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
