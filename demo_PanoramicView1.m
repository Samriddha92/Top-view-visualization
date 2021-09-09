clear all;
close all;

%% Matrix completion demo
% This short demo shows how to use TFOCS to perform nuclear norm minimization.
% Nuclear norm minimization is used for recovering all the entries of 
% a partially observed low-rank matrix.

%% Setup a problem
%rng(234923);    % for reproducible results

% for football ground
N1   = 400;       % the matrix is N1 x N2
N2   = 400*3;       % the matrix is N1 x N2
n1  = 400;        % sample image size  
r   = 400;        % the rank of the matrix
%df  = N1*N2; % - r;  % degrees of freedom of a N x N rank r matrix
%df  = 2*N*r - r^2;  % degrees of freedom of a N x N rank r matrix
%nSamples    = 3*df; % number of observed entries
nSamples    = 0.60*N1*N2; % number of observed entries



%nSamples    = df/1.2;%df/1.2; % number of observed entries

% %% vase
% N1   = 100;       % the matrix is N1 x N2
% N2   = 200;       % the matrix is N1 x N2
% n1  = 100;        % sample image size  
% r   = 70;        % the rank of the matrix
% df  = 2*N1*N1 - r;  % degrees of freedom of a N x N rank r matrix
% 
% nSamples    = df; % number of observed entries



%nSamples    = 3*n1*n1-1000; % number of observed entries





% for soccer
 iMax    = 5;
 
 
 i1= imread('f1.jpg');
 i1=rgb2gray(i1);
i1 = imresize(i1, [n1, n1]);
%figure, imshow(i1)

%imshow(i1)
%i1=double(i1);
%i1=i1(:);
%imshow(i1)
%[a1 b1]=size(i1);
%  plot(svd(i1), '.-');
%  axis tight;

i2= imread('f2.jpg');
i2=rgb2gray(i2);
i2 = imresize(i2, [n1, n1]);
%figure, imshow(i2)

%i2=double(i2);
%i2=i2(:);
%imshow(i2)

i3= imread('f3.jpg');
i3=rgb2gray(i3);
i3 = imresize(i3, [n1, n1]);
%figure, imshow(i3)

%imshow(i3)
%i3=double(i3);
%i3=i3(:);



I = [i1 i2 i3];

%I=rgb2gray(imread('embedding.jpg'));

%%

% J=imread('vase.png');
% J=rgb2gray(J);
% sz1 = size(J,1);
% sz2 = size(J,2);
% J1 = imcrop(J,[0 0 sz1 sz2]);
% J2 = imcrop(J,[0 (sz2/2)+1 sz1-1 sz2-1]);
% figure,imshow(J1);
% figure,imshow(J2);


% %% for vase
% 
%  i1= imread('JJ.jpg');
% i1 = imresize(i1, [n1, n1]);
% figure, imshow(i1)
% %i1=rgb2gray(i1);
% %i1=double(i1);
% %i1=i1(:);
% %imshow(i1)
% 
% 
% i2= imread('JJ1.jpg');
% i2 = imresize(i2, [n1, n1]);
% figure, imshow(i2)
%I = [i1 i2];









 

%imshow(I)
% I= I(:);


%%procesing for algorithm

 I=double(I);
X=I;
% X       = randi(iMax,N,r)*randi(iMax,r,N); % Our target matrix



%%
% Now suppose we only see a few entries of X. Let "Omega" be the set
% of observed entries



rPerm   = randperm(N1*N2); % use "randsample" if you have the stats toolbox
 
%omega   = sort( rPerm(1:nSamples) );
omega   = rPerm(1:nSamples);
%%
% Print out the observed matrix in a nice format.
% The "NaN" entries represent unobserved values. The goal
% of this demo is to find out what those values are!

Y = nan(N1,N2);
Y(omega) = X(omega);
%Y(omega) = I(omega);
disp('The "NaN" entries represent unobserved values');
Y=uint8(Y);
%disp(Y)
imshow(Y)
%% Matrix completion via TFOCS


observations = X(omega);    % the observed entries
%observations = I(omega);    % the observed entries
mu           = .001;        % smoothing parameter

% The solver runs in seconds
tic
Xk = solver_sNuclearBP( {N1,N2,omega}, observations, mu );
toc

%%
% To display the recovered matrix, let's round it to the nearest
% .0001 so that it displays nicely:
disp('Recovered matrix (rounding to nearest .0001):')
%Xk=Xk+;
Xk=uint8(Xk);
Xk = Xk(:,:,[1 1 1]);
figure, imshow(Xk);
%colormap

%disp( round(Xk*10000)/10000 )
%figure, imshow( uint8(round(Xk*10000)/10000) )

% %% and for reference, here is the original matrix:
% disp('Original matrix:')
% %disp( X )
% figure, imshow(X)
% % The relative error (without the rounding) is quite low:
% fprintf('Relative error, no rounding: %.8f%%\n', norm(X-Xk,'fro')/norm(X,'fro')*100 );



