img=imread('image.png');
imshow(img); 
[x,y]=ginput(6);
img_point=[x,y];
disp(img_point)