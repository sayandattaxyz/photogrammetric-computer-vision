clc;
%hard-coding the points that we found out by ginput in program "imagePoints.m"
imagePoints=1.0e+03 *[
    0.2232    0.7066;
    0.3852    0.9345;
    0.7330    1.5612;
    0.8139    1.7531;
    0.8859    0.9075;
    1.0268    0.6797;

];

%hard-coding the points that we measured manually with scale (unit: mm)
objectPoints = [83 155 0;53 112 0; 88 0 7.5;48 0 40; 0 114 83; 0 159 113];

%creating condition matrices
imageTransformation = conditionedMatrix(imagePoints);
objectTransformation = conditionedMatrix(objectPoints);

%converting image and object points to homogeneous points by adding 1 in
%the last row
homogeneousImagePts = [[imagePoints]';ones(1,6)];
homogeneousObjectPts = [[objectPoints]';ones(1,6)];

%conditioning points
conditionedImagePts = (imageTransformation * homogeneousImagePts)';
conditionedObjectPts = (objectTransformation * homogeneousObjectPts)';

%creating design matrix
A = designMatrix(conditionedImagePts, conditionedObjectPts );
disp(A)

%using svd to find Projection matrix
[U,D,V] = svd(A);
SVD = V(:,end);
%reshaping p
projectionMatrix = reshape(SVD,[4,3])';
projectionMatrix = inv(imageTransformation)*projectionMatrix*objectTransformation;

% normalising projection matrix
projectionMatrix = projectionMatrix/projectionMatrix(end,end);
disp("Projection Matrix");
disp(projectionMatrix)

%calculating RQ decomposition of P
lambda = 1/norm(projectionMatrix(3,1:3));
M = projectionMatrix(:,1:3);
if det(M) > 0
    lambda = lambda;
else
    lambda=-lambda
end
M = lambda .* M;
[R,Q] = qr(inv(M));
R = inv(R);
Q = inv(Q)
[U,D,V] = svd(projectionMatrix);
c = V(:,end);

%calculating the exterior values
projection_center = c/c(4,1)
omegaDegree = atand(R(3,2)/R(3,3))
phiDegree = -asind(R(3,1))
kappaDegree = atand(R(2,1)/R(1,1))


%calculating the interior values
principal_distance = Q(1,1);
skewDegree = acotd(-Q(1,2)/Q(1,1))
principal_point = [Q(1,3),Q(2,3)]
principle_distance= Q(1,1)
aspect_ratio = Q(2,2)/Q(1,1)


%defining the funtions

%function for finding the design matrix
function [blankMatrix] = designMatrix(imagePoints, objectPoints)   
    blankMatrix = zeros(2*length(imagePoints(:,1)),12);  %defining a blank matrix for storing the values
    for i = 1:length(imagePoints)
        x = objectPoints(i,:);
        u = imagePoints(i,1);
        v = imagePoints(i,2);
        w = imagePoints(i,3);
        blankMatrix(2*i-1,1:4) = - w*x;
        blankMatrix(2*i,5:8) = -w*x;
        blankMatrix(2*i,9:12) = v*x;
        blankMatrix(2*i-1,9:12) = u*x;
     end
end

% function for conditioning
function [conditioning] = conditionedMatrix(point)
    ctrd = zeros(1,length(point(1,:)));
    %calculating the centroids of the points
    for i = 1:length(point)
        for k = 1:length(point(1,:))
            ctrd(1,k) = ctrd(1,k)+point(i,k);
        end
    end
    ctrd = ctrd/length(point);
    
    meanDist = zeros(1,length(ctrd));
    for i = 1:length(point)
        point(i,:) = point(i,:) - ctrd(1,:);
        for j = 1:length(ctrd)
            meanDist(1,j) = meanDist(1,j)+ abs(point(i,j));        
        end
    end
    meanDist = meanDist/length(point);
    conditioning = zeros(length(ctrd)+1,length(ctrd)+1);
    for i = 1:length(ctrd)
        conditioning(i,i) = 1/meanDist(1,i);
        conditioning(i,length(ctrd)+1) = - ctrd(1,i)/meanDist(1,i);
    end
    conditioning(length(ctrd)+1,length(ctrd)+1) = 1;
end
