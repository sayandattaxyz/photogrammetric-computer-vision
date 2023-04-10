clc;

 %  getting image points from the given bh.dat file                                  
A = readDATFile('bh.dat',  '%f%f%f%f', 4);          
x1 = A(1:2, :); 
x2 = A(3:4, :);
o = ones(1,size(x1,2));
X1 = vertcat(x1, o);
X2 = vertcat(x2, o);

                      
  

disp ("Task 1a:")
%fundamental matrix
F0 = fundamentalMatrix(X1, X2)

%finding the projection matrix 
[mat1,mat2]=projectionMatrix(F0);

%finding k1
k1=mat1(:,1:3)

%finding k2
k2=calibrationMatrixK2(mat2)


disp ("Task 1b:")
%finding essential matrix
essential_matrix= k1*F0*k2

disp ("Task 1c:")
fourFoldAmbiguity(essential_matrix,k1,k2,x1,x2)

singularEssentialMatrix=singularityEnforce(essential_matrix);

disp("Task 1d:")
disp("Epipolar line computed and plotted")
%calculating the epipolar lines
epipolarLineX1 = epipolarLine(singularEssentialMatrix', X2);
epipolarLineX2 = epipolarLine(singularEssentialMatrix, X1);

%ploting the epiploar lines
epipolarLineDrawing(epipolarLineX1, "white.jpg")
epipolarLineDrawing(epipolarLineX2, "white.jpg")

disp("Task 2a:")
%finding the geometric error
errorUnOptimsed = geometricError(singularEssentialMatrix, X1', X2')

disp("Task 2b:")
%optimization of the essential matrix using Levenberg-Marquart technique
optimizedE=optimization(singularEssentialMatrix,X1,X2)

disp("Task 2d:")
%finding the geometric error from the optimized essential matrix 
singularEssentialMatrixOtp=singularityEnforce(optimizedE);
errorAfterOptimization = geometricError(singularEssentialMatrixOtp, X1', X2')

%finding difference between the two error
errorDiff=errorUnOptimsed-errorAfterOptimization


%===================================
%definig the functions

%function for reading the .dat file (referenced from assignment 5)
function points = readDATFile(datFileName, lineFormat, length )
  fh = fopen(datFileName, 'r');
  points = fscanf(fh, lineFormat, [length inf]);
  fclose(fh);
end

% function for conditioning the matrix
function transformationMatrix = conditionedMatrix(points)       %logic referenced from assignment 4
  translationX = mean(points(:,1));                             %translating the x-coordinates to centroid
  translationY = mean(points(:,2));                             %translating the y-coordinates to centroid
  scaleX = mean(abs(points(:,1)-translationX));                 %scaling the x-coordinates 
  scaleY = mean(abs(points(:,2)-translationY));                 %scaling the y-coordinates
  transformationMatrix = [1/scaleX 0    -translationX/scaleX;   %defining the Transformation Matrix
        0   1/scaleY  -translationY/scaleY;
        0   0      1   ];
end


%function for finding the Fundamental Matrix referenced from assigment 5
function F = fundamentalMatrix(points1, points2)   
  t1 = conditionedMatrix(points1); 
  transformation1 = t1 * points1;   
  t2 = conditionedMatrix(points2); 
  transformation2 = t2 * points2;     
  A = designMatrix(transformation1, transformation2);          % design matrix creation    
  f = SVD(A);                                                  
  F = t2' * reshape(f, 3, 3)' * t1;                            % Unconditioning
end

% function for calculating Design matrix referenced assignment 
function A = designMatrix(x1, x2)     
  A = [];
  for i = 1 : 1383
    At = [x1(1,i)*x2(1,i) , x1(2,i)*x2(1,i) , x2(1,i) , x1(1,i)*x2(2,i), x1(2,i)*x2(2,i), x2(2,i) , x1(1,i) , x1(2,i) ,1 ]; 
    A = vertcat(A,At);
  end
end


%function for finding singular value decomposition
function x = SVD(A)    
  [U, D, V] = svd(A);                                          
  x = V(:, end);                                              
end

%function for finding projection matrix
function [mat1,mat2]=projectionMatrix(F)  
mat1=[eye(3,3),zeros(3,1)];
e=SVD(F);
skewMatrix=[  0    -e(3)   e(2)                                 % getting Skew matrix                                                                              
    e(3)    0    -e(1)                                  
   -e(2)   e(1)   0  ];                                  
mat2=[skewMatrix*F+[e e e],e];                                  % second Projection Matrix
end


%function for finding linear triangulation
function value=linearTriangulation(x1,x2,p1,p2) 
S=[];                                                    % Finding Design Matrix S
value=[];
for i=1:size(x1,2)                                       
S=[x1(1,i)*p1(3,:)-p1(1,:)
   x1(2,i)*p1(3,:)-p1(2,:)
   x2(1,i)*p2(3,:)-p2(1,:)
   x2(2,i)*p2(3,:)-p2(2,:)];

[U D V]=svd(S);                                                  % svd of S

value=[value,V(:,size(V,2))./V(4,size(V,2))];                    % getting Object Points
end
end

%function for finding geometric error
function errorFundamentalMatrix = geometricError(fundamentalMatrix, image1Points, image2Points)
 errorFundamentalMatrix = 0;
 image1Points;
 image2Points;
  for i = 1:1383
    imageLine1 = fundamentalMatrix' * image2Points(i, :)';
    imageLine2 = fundamentalMatrix * image1Points(i, :)';
    upErrorFundamentalMatrix = power(image2Points(i, :) * imageLine2, 2);
    downErrorFundamentalMatrix = power(imageLine1(1, 1), 2) + power(image1Points(2, 1), 2) + power(image2Points(1, 1), 2) + power(image2Points(2, 1), 2); 
    errorFundamentalMatrix =errorFundamentalMatrix + upErrorFundamentalMatrix/downErrorFundamentalMatrix;
  end
end

%function for calculating epipolar line
function line = epipolarLine(fundamentalMatrix, imagePoints)
  line = [];
  for i = 1:1383
    line = [
      line, fundamentalMatrix' * imagePoints(:, i);
    ];
  end
end

% function for drawing the epipolar line
function epipolarLineDrawing(line, imageName)
  img = imread(imageName);                                   % reading the image
  figure; 
  imshow(img);                  
  for i = 1:1383
    hline(line(:, i));          %using function hline that has already been provided to draw blue line on the images with the homogeneous co-ordinate 
  end  
 end

%function for calibration matrix k2
function k2=calibrationMatrixK2(mat2)
lambda = 1/norm(mat2(3,1:3));
M = mat2(:,1:3);
if det(M) > 0
    lambda = lambda;
else
    lambda=-lambda;
end
M = lambda .* M;
[R,Q] = qr(inv(M));
R = inv(R);
k2 = inv(Q);
end

%function for singularity enforcement of the essential matrix
function singularEssentialMatrix=singularityEnforce(essential_matrix)
[U, D, V] = svd(essential_matrix);
d = diag(D);
singularEssentialMatrix = U * diag([(d(1)+d(2))/2,(d(1)+d(2))/2, 0]) * V';
end

%non-linear optimization by Levenberg-Marquart technique
%code referenced using matlab documentation for "lsqnonlin"
%reference: https://de.mathworks.com/help/optim/ug/lsqnonlin.html
function optimizedE=optimization(singularEssentialMatrix,X1,X2)
fun = @(singularEssentialMatrix) geometricError(singularEssentialMatrix, X1', X2');
optimizedE = lsqnonlin(fun, singularEssentialMatrix);
end


function fourFoldAmbiguity(essential_matrix,k1,k2,x1,x2)
[U, D, V] = svd(essential_matrix);

% Define the possible rotation and translation matrices
ax = [0, -1, 0; 1, 0, 0; 0, 0, 1];
az = [0, 1, 0; -1, 0, 0; 0, 0, 0];
R1 = U * ax * V';
R2 = U * ax' * V';
R3 = U * az * U';
R4 = U * az' * U';

% Define the translation vector
t1 = U(:, 3);
t2 = -U(:, 3);
t3 = U(:, 3);
t4 = -U(:, 3);

% Define the projection matrices
P1 = k1 * [eye(3), zeros(3, 1)];
P2_1 = k2 * [R1, t1];
P2_2 = k2 * [R2, t2];
P2_3 = k2 * [R3, t3];
P2_4 = k2 * [R4, t4];

% Triangulate the 3D points using each projection matrix
lt1 = linearTriangulation(x1', x2', P1, P2_1);
lt2 = linearTriangulation(x1', x2', P1, P2_2);
lt3 = linearTriangulation(x1', x2', P1, P2_3);
lt4 = linearTriangulation(x1', x2', P1, P2_4);

% Select the solution with positive z-coordinate x(3)
if lt1(3) > 0
    X = lt1;
    P2 = P2_1;
    R = R1;
    t = t1;
elseif lt2(3) > 0
    X = lt2;
    P2 = P2_2;
    R = R2;
    t = t2;
elseif lt3(3) > 0
    X = lt3;
    P2 = P2_3;
    R = R3;
    t = t3;
elseif lt4(3) > 0
    X = lt4;
    P2 = P2_4;
    R = R4;
    t = t4;
else
    error('solution not valid');
end

% Display the selected solution
fprintf('Rotation matrix:\n');
disp(R);
fprintf('Translation vector:\n');
disp(t);
end