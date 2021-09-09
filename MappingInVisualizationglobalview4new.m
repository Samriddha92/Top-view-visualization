clear all;
close all;

% Load soccer images.
%buildingDir = fullfile('D:', 'SAMRIDDHA' , 'TFOCS-master-panoramicView', 'TFOCS-master' ,  'examples' , 'demos' ,'soccer6');
%D:\SAMRIDDHA\Research Activities\Top View visualization

[filename,pathname,index]=uigetfile('*.mp4','Enter video filename');
messivideo=VideoReader(filename);
Folder = 'D:\SAMRIDDHA\Research_Activities\panoramic view\Top View visualization\soccer7\';
Top_view_coordinate = {};


%% homography of panorama and visualization

temp =xlsread('temp.xlsx');
%temp =xlsread('temp1.xlsx');
 temp3=[];
ntemp =xlsread('ntemp.xlsx');
% ntemp =xlsread('ntemp1.xlsx');

%% homography of top view and global view
   temp1 =xlsread('temp1.xlsx');
   ntemp1 =xlsread('ntemp1.xlsx');
 %% image read
 I2= imread('visualization.jpg');
 [m n]=size(I2);
 
 I4 = imread('global_view.jpg');
 I4 = rgb2gray(I4);
 
 %% input point selection 
 
 I1= imread('panorama3.jpg');
% I1= imread('global_view.jpg');
 I2= imread('visualization.jpg');
%     figure, imagesc(I1)
%     p=input('Enter number of points ');
%     
%     
%     
%     for j=1:p
%         
%      % str=input('enter object ID ');
%         temp1=ginput(1);
%         temp1=temp1';         
%         temp =[temp temp1];   
%         j
%        end
%          hold on
%         % plot(temp(1,:),temp(2,:))
%          scatter(temp(1,:),temp(2,:))
%          hold off
% %        temp = [ temp [ i j str  round(x(:))' round(y(:))']' ];
% 
%     
%     I2=imread('visualization.jpg');
%      imagesc(I2)
%     q=input('Enter number of points ');
%     for j=1:q
%         
%      % str=input('enter object ID ');
%          ntemp1=ginput(1);
%          ntemp1=ntemp1';       
%          ntemp =[ntemp  ntemp1];  
%          
%        end
%          hold on
%          scatter(ntemp(:,1),ntemp(:,2))
%          hold off
% %        temp = [ temp [ i j str  round(x(:))' round(y(:))']' ];
       
    
    
    %% homography from panorama to top view
    
    correp= homography_solve(temp, ntemp);

    %% homography from top view to global view
    correp1= homography_solve(ntemp1,temp1);










for i = 532:3:541 %802 %532%:3:965%3:1000   %70:3:100
   
    b = read(messivideo,i); 
     b1 = read(messivideo,i);
 %   imshow(b)
 % frames = read(vid, iFrame);
 %FILENAME = ['C:\Edgeresults\', filename1];
%imwrite(b,Folder,'jpg')
% imwrite(b, fullfile(Folder, sprintf('%1d.jpg', 2)));

  %% read the images from the folder
  buildingDir = fullfile('D:', 'SAMRIDDHA' , 'Research_Activities', 'panoramic view', 'Top View visualization' ,'soccer7');
  buildingScene = imageDatastore(buildingDir);
  
  %% Display images to be stitched
  %  montage(buildingScene.Files);
   
    % Read the first image from the image set.
    I = readimage(buildingScene, 1);
    
    % Initialize features for I(1)
    grayImage = rgb2gray(I);
    points = detectSURFFeatures(grayImage);
    
%      %% show key points
%     imshow(I); hold on;
%     plot(points.selectStrongest(10));


    [features, points] = extractFeatures(grayImage, points);
    
    % Initialize all the transforms to the identity matrix. Note that the
% projective transfor m is used here because the building images are fairly
% close to the camera. Had the scene been captured from a further distance,
% an affine transform would suffice.
numImages = numel(buildingScene.Files);
tforms(numImages) = projective2d(eye(3));

% Initialize variable to hold image sizes.
imageSize = zeros(numImages,2);

% Iterate over remaining image pairs
for n = 2:numImages
    
    % Store points and features for I(n-1).
    pointsPrevious = points;
    featuresPrevious = features;
        
    % Read I(n).
    I = readimage(buildingScene, n);
    
    % Convert image to grayscale.
    grayImage = rgb2gray(I);    
    
    % Save image size.
    imageSize(n,:) = size(grayImage);
    
    % Detect and extract SURF features for I(n).
    points = detectSURFFeatures(grayImage);    
%     %% show key points
%     imshow(I); hold on;
%     plot(points.selectStrongest(10));
    
    
    
    [features, points] = extractFeatures(grayImage, points);
  
    % Find correspondences between I(n) and I(n-1).
    indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);
       
    matchedPoints = points(indexPairs(:,1), :);
    matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);        
    
    % Estimate the transformation between I(n) and I(n-1).
    tforms(n) = estimateGeometricTransform(matchedPoints, matchedPointsPrev,...
        'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);
    
   % Compute T(n) * T(n-1) * ... * T(1)
   tforms(n).T = tforms(n).T * tforms(n-1).T; 
end



%% given frame and panorama correspondence 
% frame = imwarp(b,tforms(2));
% figure, imshow(frame);





%% show frames
 
%AA=xlsread('mu959.xlsx');
AA=xlsread('MUvsEVfrom532to802.xlsx');
% AA=[AA(:,1) AA(:,6) AA(:,7) AA(:,2) AA(:,3) AA(:,4) AA(:,5)];
% AA=[AA(:,1) AA(:,2) AA(:,3) AA(:,4) AA(:,4)+AA(:,6) AA(:,5) AA(:,5)+AA(:,7)];

%% all player general id
% AA=xlsread('MUvsEV.xlsx');
% AA=[AA(:,1) AA(:,2) AA(:,3) AA(:,4) AA(:,6) AA(:,5) AA(:,7)];
% AAA=AA(find(AA(:,1)==i), :);

%AA=[AA(:,1) AA(:,7) AA(:,2) AA(:,3) AA(:,4) AA(:,5)];
%AA=[AA(:,1) AA(:,6) AA(:,7) AA(:,2) AA(:,4) AA(:,3) AA(:,5)];
AAA=AA(find(AA(:,1)==i), :);
playerId=AAA(:,3);
%playerId= num2str(playerId);
AA1=[AAA(:,3) AAA(:,4) AAA(:,6)];
AA2= [AAA(:,5) AAA(:,7)];
[m n]= size(AA2);
AA3=[AA2;ones(1,n)];
AA4= AA2;
B=AA(find(AA(:,1)==i), :);   

 U1=[];   

 
     for j=1:size(find(B(:,1)==i),1);
         
%         U  = transformPointsForward(tforms(2),AA4(j,1),AA4(j,2));
%        U1= [U1;U];
          
         box=[B(j,4),B(j,6),abs(B(j,5)-B(j,4)),abs(B(j,6)-B(j,7))];
         
         %% New groundtruth
        % box=[B(j,3),B(j,4),abs(B(j,3)-B(j,5)),abs(B(j,4)-B(j,6))];
         
        
          b = insertShape(b, 'rectangle', box, 'LineWidth', 3,'color',[0 0 0]); 
          
           position=[B(j,4) B(j,6)];
           
           %% new groundtruth
         %  position=[B(j,3) B(j,4)];
           
           value=num2str(B(j,3));
            b = insertText(b,position,value,'AnchorPoint','LeftBottom','FontSize',19);
             position1 = [50 50];
           value1=num2str(i);
           b = insertText(b,position1,value1,'AnchorPoint','LeftBottom','FontSize',19);
    %      imshow(b)
         % rectangle('Position', box,'Linewidth',1,'EdgeColor', [10j/220 5j/220 20j/220]);
           
          
           
    %        pause(0.02);
        % imshow(b);
         

     end
% end


% %% given frame and panorama correspondence 
% frame = imwarp(b,tforms(2));
% figure, imshow(frame);

U  = transformPointsForward(tforms(2),AA4);

%% plotting points in panorama
% figure,imagesc(imread('panorama2.jpg'));
% 
% hold on
% scatter(U(:,1),U(:,2),'filled');
% hold off

%% given frame and visualization correspondence
    AA5= U';
    dummy= ones(1,m);
    AA5=[AA5;dummy];
     newpt= homography_transform(AA5,correp);
     newpt1=newpt';
     Top_view_coordinate{i} = newpt1;
     
%%  collection point top view
    



    %% top view and global view correspondence 
     gpoint=newpt;
     gpoint1=[gpoint;dummy];
     gpoint2= homography_transform(gpoint1,correp1);
     gpoint2=gpoint2';
     
I3=I2;
    for ii=1:m
        position=[newpt1(ii,1) newpt1(ii,2)];
       %position=[50.5704 80.9451]
       text= num2str(playerId(ii));  
       %text= num2str(2); 
       %b = insertText(b,position,value,'AnchorPoint','LeftBottom','FontSize',19);
     I3= insertText(I3,position,text,'AnchorPoint','LeftBottom','FontSize',18);
   %  imshow(I3)
         
    end
    I5=I4;
   for ii=1:m
   position_globl=[gpoint2(ii,1) gpoint2(ii,2)]
   textg= num2str(playerId(ii));
    I5= insertText(I5, position_globl,textg,'AnchorPoint','LeftBottom','FontSize',18);
   %I4= imresize(I4,1.2);
   end
   % Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    figure(1);
    subplot(2,3,1);
    imshow(b1);
    title('Fig. (a): Input video ','fontsize', 18)
    subplot(2,3,2);
    imshow(b);
    title('Fig. (b): Player detection with unique id ','fontsize', 18)    
    subplot(2,2,3);     
    imshow(I5);
    title('Fig. (c): Static camera visualization','fontsize', 18)   
       
        
    subplot(2,2,4);
    
    
    
    
    imshow(I3);
    hold on
    plot(newpt(1,:),newpt(2,:),'or')
     title('Fig. (d): Top view visualization','fontsize', 18) 
    hold off 
  
   
%  fname = 'D:\SAMRIDDHA\Research_Activities\panoramic view\Top View visualization\Topviewvideos1';
%         filename= num2str(i);
%         saveas(gca, fullfile(fname, filename), 'jpeg');
%  
end



 






















    
    
    