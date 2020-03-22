% 
% Very short description of the files in this folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% mixturefit.m 
% Standard Gaussian Mixture Model fit to the data. I got this Matlab code from Todd Leen.
% 
% pnormal.m
% used by mixturefit.m 
%
% plotgmmpdf.m
% A (very stupid) function that plots the GMM pdf. Only works in 2 dimensions. 
% It inputs the GMM parameters along with the boundaries and the sampling density of the pdf to be plotted.
% 
% pc_projection_gmm.m
% This is the main routine that projects any point x into the feature space 
% onto the principal surface of target dimensionality. 
% - Inputs the GMM parameters, and x. 
% - Outputs x* - principal curve projection of x.
%
% pc_projection_multidim.m
% This is the main pc projection routine based on KDE. inputs the data
% samples and the kernel size, outputs the principal curve projection of
% the data
%
% gmmpdf.m
% - Inputs the GMM parameters and x
% - Outputs p(x), g(x), H(x), SI(x) (and also optionally another vector for an
% intermediate step that I use in the projection)
% 
%%%%%%%%%%%%%%%

clear all
close all

% generating data
% dataflag = 1   A single Gaussian in 2D
% dataflag = 2   Gaussian mixture with 3 components in 2D (*-type intersection)
% dataflag = 3   Gaussian mixture with 3 components in 2D (T-type intersection)
dataflag = 3;

if dataflag == 1,
N1=200;cluster1=randn(2,N1); Sc1=[3 0;0 1]; Rot1=[cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)];mean1=[0; 0];
data=Sc1*cluster1+repmat(mean1,1,N1);
end
if dataflag == 2,
x=8;y=4;z=20;t=15;
a=randn(2,100); a(1,:)=4*a(1,:);R=[sqrt(3)/2 1/2; -1/2 sqrt(3)/2];a1=R*a;
a=randn(2,100);a(1,:)=4*a(1,:);a1=[x*ones(1,100);y*ones(1,100)]+a1;
R2=[sqrt(3)/2 -1/2; 1/2 sqrt(3)/2];
a2=R2*a;a2=[-x*ones(1,100);y*ones(1,100)]+a2;
a=randn(2,100);a(1,:)=4*a(1,:);R3=[0 1;1 0 ];a3=R3*a;
a3=[0*ones(1,100);t*ones(1,100)]+a3;
data = [a1,a2,a3];
end
if dataflag == 3,
N1=200;cluster1=randn(2,N1); Sc1=[3 0;0 1]; Rot1=[cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)];mean1=[0; 0];
N2=100;cluster2=randn(2,N2); Sc2=[2 0;0 1]; Rot2=[cos(pi/4) sin(pi/4);-sin(pi/4) cos(pi/4)];mean2=[5; 5];
N3=100;cluster3=randn(2,N3); Sc3=[3 0;0 1]; Rot3=[cos(-pi/4) sin(-pi/4);-sin(-pi/4) cos(-pi/4)];mean3=[-5; 5];
data=[Rot1*Sc1*cluster1+repmat(mean1,1,N1) Rot2*Sc2*cluster2+repmat(mean2,1,N2) Rot3*Sc3*cluster3+repmat(mean3,1,N3)];
end

plot(data(1,:),data(2,:),'.'),axis equal;

% fitting GMM to the data 

% assuming we know the number of components
if dataflag==1,ncomp=1;else,ncomp=3;end

constraint=0; %which means we are looking for full covariance
eps=1e-6; %default value
[alphas,means,mcovs]=mixturefit(data',ncomp,constraint,eps);

% let's see how it looks like
scale=trace(cov(data'));
xmin=min(data(1,:))-scale/10;ymin=min(data(2,:))-scale/10;
xmax=max(data(1,:))+scale/10;ymax=max(data(2,:))+scale/10;
density=10;
plotgmmpdf(xmin,xmax,ymin,ymax,density,alphas,means,mcovs);

% now project the data points onto the principal curve
targetdim = 1; % the target dimensionality is 1 for principal curve. 
               % for any principal surface of different dimensionality just 
               % change this
data_in = data;% initialize on a grid, or onto the data points themselves
pc_projection = pc_projection_GMM(means,alphas,mcovs,targetdim,data_in)

figure,plot(data(1,:),data(2,:),'.');hold on,
title('PC projection with GMM pdf')
plot(pc_projection(1,:),pc_projection(2,:),'.r');axis equal;

% now let's see what happens with a KDE... 
% targetdim = 1 as before, we don't want to find a principal surface etc.
% but the principal curve.
% 
kernel_sigma = 1.5; % play with this parameter... I handpicked the one over
                    % here but there are many ways to select a good kernel 
                    % bandwidth. See my thesis for references...

pc_projection=pc_project_multidim(data,data_in,kernel_sigma,targetdim)
pc_projection = pc_projection';
figure,plot(data(1,:),data(2,:),'.');hold on,
title('PC projection with KDE pdf')
plot(pc_projection(1,:),pc_projection(2,:),'.r');axis equal;

keyboard;