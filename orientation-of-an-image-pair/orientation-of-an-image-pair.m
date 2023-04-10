clc;

image1 = imread("image1.png");     % reading image1
imshow(image1);
[x1,y1]=ginput(8);                 % selecting 8 points from the image1

image2 = imread("image2.png");     % reading image2
imshow(image2);
[x2,y2]=ginput(8);                 % selecting 8 points from image2

%converting the image1 and image2 points into homogeneous points by adding
%1 in the last row
image1Points=[x1,y1,ones(8,1)];
size(image1Points);
image2Points=[x2,y2,ones(8,1)];

% conditioning image1 points
transformationImage1= conditionedMatrix(image1Points); 
conditionedImage1 = transformationImage1*image1Points'

% conditioning image2 points
transformationImage2= conditionedMatrix(image2Points); 
conditionedImage2 = transformationImage2*image2Points'

%creating the design matrix
A= designMatrix (conditionedImage1,conditionedImage2)

%finding the fundamental matrix
[U, D, V]=svd(A)
f=V(: , end)
fundamentalMatrix= transformationImage2' * reshape(f,3,3)'* transformationImage1   %finding fundamental martrix by reconditioning 


%finding singular fundamental matrix
[U, D, V] = svd(fundamentalMatrix);
d = diag(D);
singularFundamentalMatrix = U * diag([d(1), d(2), 0]) * V'

% finding epipolar line
epipolarLineImage1 = epipolarLine(singularFundamentalMatrix', image2Points')
epipolarLineImage2 = epipolarLine(singularFundamentalMatrix, image1Points');


% finding geometric error
errorFundamentalMatrix = geometricError(singularFundamentalMatrix, image1Points, image2Points)

% drawing the epipolar line
epipolarLineDrawing(epipolarLineImage1, "image1.png");
epipolarLineDrawing(epipolarLineImage1, "image2.png");


%defining the functions

%defining the function for Transformation Matrix
function transformationMatrix = conditionedMatrix(points)
  translationX = mean(points(:,1));                             %translating the x-coordinates to centroid
  translationY = mean(points(:,2));                             %translating the y-coordinates to centroid
  scaleX = mean(abs(points(:,1)-translationX));                 %scaling the x-coordinates 
  scaleY = mean(abs(points(:,2)-translationY));                 %scaling the y-coordinates
  transformationMatrix = [1/scaleX 0    -translationX/scaleX;   %defining the Transformation Matrix
        0   1/scaleY  -translationY/scaleY;
        0   0      1   ];
end

%defining the function for Design Matrix
function A = designMatrix (pointsImage1,pointsImage2)
  A = [];
  for i = 1:8     % iterating 1 to 8 for all the points             
     A = [A; 
          pointsImage1(1,i).*pointsImage2(1,i),  pointsImage1(2,i).*pointsImage2(1,i) , pointsImage2(1,i),  pointsImage1(1,i).*pointsImage2(2,i),  pointsImage1(2,i).*pointsImage2(2,i) , pointsImage2(2,i),  pointsImage1(1,i), pointsImage1(2,i), 1 ];
       
   end
end

%function for finding out the epipolar line
function line = epipolarLine(fundamentalMatrix, imagePoints)
  line = [];
  for i = 1:8
    line = [
      line, fundamentalMatrix * imagePoints(:, i);
    ];
  end
end

%defining the function for geometric error finder 
function errorFundamentalMatrix = geometricError(fundamentalMatrix, image1Points, image2Points)
 errorFundamentalMatrix = 0;
 image1Points;
 image2Points;
  for i = 1:8
    imageLine1 = fundamentalMatrix' * image2Points(i, :)';
    imageLine2 = fundamentalMatrix * image1Points(i, :)';
    upErrorFundamentalMatrix = power(image2Points(i, :) * imageLine2, 2);
    downErrorFundamentalMatrix = power(imageLine1(1, 1), 2) + power(image1Points(2, 1), 2) + power(image2Points(1, 1), 2) + power(image2Points(2, 1), 2); 
    errorFundamentalMatrix =errorFundamentalMatrix + upErrorFundamentalMatrix/downErrorFundamentalMatrix;
  end
end


%defining the function for drawing epipolar line
function epipolarLineDrawing(line, imageName)
  img = imread(imageName);                                   % reading the image
  figure; 
  imshow(img);                  
  for i = 1:8 
    hline(line(:, i));          %using function hline that has already been provided to draw blue line on the images with the homogeneous co-ordinate 
  end  
  
end

%test
