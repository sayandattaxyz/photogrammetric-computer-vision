clc;

                                  
A = readDATFile('bh.dat',  '%f%f%f%f', 4);          %  getting image points from the given bh.dat file 
x1 = A(1:2, :); 
x2 = A(3:4, :);
o = ones(1,size(x1,2));
x1 = vertcat(x1, o);
x2 = vertcat(x2, o);                                               
fMatrix = fundamentalMatrix(x1, x2);                      % fundamental matrix
[mat1,mat2]=projectionMatrix(fMatrix);                    % projection matrices
Xp=linearTriangulation(x1,x2,mat1,mat2); 
spatialObjectCordVisualization(Xp);
 
                    
 A= readDATFile('pp.dat', '%f%f%f%f%f%f%f', 7)       % getting image points from the given bh.dat file
 b1 = A(1:2, :);
 b2 = A(3:4, :);
 b3 = A(5:7, :);
 b1(3,:) = 1;
 b2(3,:) = 1;
 b3(4,:) = 1;
objectPonts=linearTriangulation(b1,b2,mat1,mat2);              
H=homography3D(objectPonts,b3);                            % homography for the object point
H_=H*Xp;                                                % applying Transformation H to all projective object points Xp1
transformationH=[];

for i = 1:size(H_,2)
    transformationH = [transformationH,H_(:,i)./H_(4,i)];
end
spatialObjectCordVisualization(transformationH);


%defining the functions

%function for reading the .dat file
function points = readDATFile(datFileName, lineFormat, length )
  fh = fopen(datFileName, 'r');
  points = fscanf(fh, lineFormat, [length inf]);
  fclose(fh);
end


function transformationMatrix = conditionedMatrix(points)       %logic referenced from assignment4
  translationX = mean(points(:,1));                             %translating the x-coordinates to centroid
  translationY = mean(points(:,2));                             %translating the y-coordinates to centroid
  scaleX = mean(abs(points(:,1)-translationX));                 %scaling the x-coordinates 
  scaleY = mean(abs(points(:,2)-translationY));                 %scaling the y-coordinates
  transformationMatrix = [1/scaleX 0    -translationX/scaleX;   %defining the Transformation Matrix
        0   1/scaleY  -translationY/scaleY;
        0   0      1   ];
end


function transformationMatrix3D = conditionedMatrix3D(points)          %logic referenced from assignment4
  translationX = mean(points(:,1));                                    %translating the x-coordinates to centroid
  translationY = mean(points(:,2));                                    %translating the y-coordinates to centroid
  translationZ = mean(points(:,3));                                    %translating the z-coordinates to centroid
  scaleX = mean(abs(points(:,1)-translationX));                        %scaling the x-coordinates
  scaleY = mean(abs(points(:,2)-translationY));                        %scaling the y-coordinates
  scaleZ = mean(abs(points(:,3)-translationZ));                        %scaling the x-coordinates
  transformationMatrix3D = [ 1/scaleX, 0, 0, -translationX/scaleX;     %defining the Transformation Matrix
  0, 1/scaleY, 0, -translationY/scaleY;
  0, 0, 1/scaleZ, -translationZ/scaleZ;
  0, 0, 0, 1 ];
end

% Fundamental Matrix
function F = fundamentalMatrix(points1, points2)   
  t1 = conditionedMatrix(points1); 
  transformation1 = t1 * points1;   
  t2 = conditionedMatrix(points2); 
  transformation2 = t2 * points2;     
  A = designMatrix(transformation1, transformation2);          % design matrix creation    
  f = SVD(A);                                                  
  F = t2' * reshape(f, 3, 3)' * t1;                            % Unconditioning
end

function H =homography3D(x1 , x2)
conditionedMatrix3D1 = conditionedMatrix3D(x1);                         % conditioning the points
conditionedMatrix3D2 = conditionedMatrix3D(x2); 
x1 = conditionedMatrix3D1 * x1;
x2 = conditionedMatrix3D2 * x2;                       
A=[];                                                                  
for i=1:5                                                                   % design matrix creation                                        
A  =[A ; -x2(4,i)*x1(:,i)'     0 0 0 0 0 0 0 0         x2(1,i)*x1(:,i)';
         0 0 0 0         -x2(4,i)*x1(:,i)'   0 0 0 0         x2(2,i)*x1(:,i)';
         0 0 0 0 0 0 0 0       -x2(4,i)*x1(:,i)'  x2(3,i)*x1(:,i)'    ];
end
[U, D, V] = svd(A);                                                         % svd
h = V(:,end);
H = inv(conditionedMatrix3D2) * reshape(h, 4, 4)' * conditionedMatrix3D1     %finding fundamental martrix by reconditioning
end

%function for finding visualization of spatial object coordinates
function spatialObjectCordVisualization(points)
  figure; 
  scatter3(points(1,:), points(2,:), points(3,:), 10, 'filled');
  axis square; 
  view(32, 75);
end

function A = designMatrix(x1, x2)     % Perspective projection of Design matrix
  A = [];
  for i = 1 : size(x1, 2)
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
mat1=[eye(3,3),zeros(3,1)]
e=SVD(F');
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



