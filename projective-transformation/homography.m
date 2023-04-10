clear;
image1 = imread('imagein1.jpeg');
image2 = imread('imagein2.jpeg');
image3 = imread('imagein3.jpeg');
ftm12 = ftm(image1,image2);
imageStich12 = geokor(ftm12,image1,image2);
imshow(imageStich12);
ftm23 = ftm(image3, imageStich12);
imageStich123 = geokor(ftm23,image3,imageStich12);
imshow(imageStich123);


function H_Final = ftm(image1, image2) 
% points from image
  imshow(image1);
  [image1X, image1Y] = ginput(4);   %acquiring 4 Points from image1
  imshow(image2);
  [image2X, image2Y] = ginput(4);   %acquiring 4 Points from image2
 
% 2D homography computation ==>

% translation

translationPoint_Source = [ mean(image1X) ; mean(image1Y)] ;       %source image points
translationPoint_Dest = [ mean(image2X) ; mean(image2Y)] ;  %destination image points
  
  translationMatrix_Source = [ 1  0  -translationPoint_Source(1,1);  %defining the translation matrix for the source
         0  1  -translationPoint_Source(2,1);  
         0  0           1    ] ;

  translationMatrix_Dest = [ 1  0  -translationPoint_Dest(1,1);   %defining the translation matrix for the destination
         0  1 -translationPoint_Dest(2,1);
         0  0               1    ] ;

% scaling 

s1 = mean (abs((image1X)-translationMatrix_Source(1)));
s2 = mean (abs((image1Y)-translationMatrix_Source(2)));
s11 = mean (abs((image2X)-translationMatrix_Dest(1)));
s22 = mean (abs((image2Y)-translationMatrix_Dest(2)));

  scalingMatrixS = [1/s1  0  0;        %defining the scaling matrix
         0   1/s2 0;
         0    0    1] ;

  scalingMatrixD = [1/s11  0   0;     %defining the scaling matrix
         0   1/s22 0;
         0    0    1] ;

% coordinate Transformation

transformationS = scalingMatrixS * translationMatrix_Source;
transformationD = scalingMatrixD * translationMatrix_Dest;

% conditioned Coordinates

%source points
 xConditioned1 = transformationS * [image1X(1,1) ; image1Y(1,1) ;1]; 
 xConditioned2 = transformationS * [image1X(2,1) ; image1Y(2,1) ;1]; 
 xConditioned3 = transformationS * [image1X(3,1) ; image1Y(3,1) ;1]; 
 xConditioned4 = transformationS * [image1X(4,1) ; image1Y(4,1) ;1]; 

%destination points
 xConditioned5 = transformationD * [image2X(1,1) ; image2Y(1,1) ;1]; 
 xConditioned6 = transformationD * [image2X(2,1) ; image2Y(2,1) ;1]; %Destination point
 xConditioned7 = transformationD * [image2X(3,1) ; image2Y(3,1) ;1]; %Destination point
 xConditioned8 = transformationD * [image2X(4,1) ; image2Y(4,1) ;1]; %Destination point

% design matrix

A = [-xConditioned1'.*1 0 0 0    (xConditioned5(1,1).*xConditioned1)' ;
     0 0 0    -xConditioned1'.*1 (xConditioned5(2,1).*xConditioned1)' ;
     -xConditioned2'.*1 0 0 0    (xConditioned6(1,1).*xConditioned2)' ;
     0 0 0    -xConditioned2'.*1 (xConditioned6(2,1).*xConditioned2)' ;
     -xConditioned3'.*1 0 0 0    (xConditioned7(1,1).*xConditioned3)' ;
     0 0 0    -xConditioned3'.*1 (xConditioned7(2,1).*xConditioned3)' ;
     -xConditioned4'.*1 0 0 0    (xConditioned8(1,1).*xConditioned4)' ;
     0 0 0    -xConditioned4'.*1 (xConditioned8(2,1).*xConditioned4)' ;
     0 0 0    0 0 0     0 0 0              ]; 
 
% SVD

[U,S,V] = svd (A) ;  
solution_vector = V(:,9) ;         
H_reshape = reshape (V(:,9),3,3);
H_hat = H_reshape';

% Reverse conditioning

H_untransformed = inv(transformationD) * H_hat * transformationS;
H_Final = H_untransformed ./ H_untransformed(3,3)        %the final Homography matrix

end
%Pre-defined geokor function
function i = geokor (H, f, g)
%        ====================
	[h1, w1 d] = size(f); f = double(f);                     % Image dimensions
	[h2 w2 d] = size(g); g = double(g);

         % Transformation of image corners to determine the resulting size
	cp = norm2(H * [1 1 w1 w1; 1 h1 1 h1; 1 1 1 1]);
	Xpr = min([cp(1, :) 0]) : max([cp(1, :) w2]);              % min x : max x
	Ypr = min([cp(2, :) 0]) : max([cp(2, :) h2]);              % min y : max y
	[Xp, Yp] = ndgrid(Xpr, Ypr);
	[wp hp] = size(Xp);                                           % = size(Yp)
	X = norm2(H \ [Xp(:) Yp(:) ones(wp*hp, 1)]');    % Indirect transformation
	
	xI = reshape(X(1, :), wp, hp)';           % Resampling of intensity values 
	yI = reshape(X(2, :), wp, hp)';
	f2 = zeros(hp, wp); i = f2;
	for k = 1 : d
		f2(1:h1, 1:w1, k) = f(:, :, k);
		i(:, :, k) = interp2(f2(:, :, k), xI, yI);    % bilinear interpolation           
	end
                      % Copy the original image g in the rectified image f
	off = -round([min([cp(1, :) 0]) min([cp(2, :) 0])]);
	Index = find(sum(g, 3) > 9); 
	for k = 1 : d
    		iPart = i(1+off(2) : h2+off(2), 1+off(1) : w2+off(1), k);
    		fChannel = g(:, :, k);
    		iPart(Index) = fChannel(Index);
    		i(1+off(2) : h2+off(2), 1+off(1) : w2+off(1), k) = iPart;
	end
	i(~isfinite(i)) = 0;
	i = uint8(i);

end

function n = norm2(x)
	for i = 1 : 3
		n(i,:) = x(i,:) ./ x(3,:);	
	end

end




