clc;clear all;close all;
I=imread('testsynf.bmp');
% % make image a double and normalize
I=I(:,:,1);
I = double(I);
mx = max(I(:));
mn = min(I(:));
I = (I-mn)/(mx-mn);
% define anisotropic diffusion parameters
nither=2;
alpha=16;
J=7;
 blk0=3;
 bls=(blk0-1)/2;
% get area of uniform speckle

    figure, imshow(I,[],'InitialMagnification','fit');
    rect = getrect;
    close;
counter=0;
f=zeros(1,J*nither*6);

ither=1;
mwait=waitbar(0,char('Manipulating Wavelet coefficients at iteration: ',int2str(ither)));
for ither=1:nither

  [Faf, Fsf] = FSfarras;
  [af, sf] = dualfilt1;
  w = cplxdual2D(I, J, Faf, af);

 adw=w;


%  for k=1:J
%      for rim=1:2
%          for np=1:2
%             for o=1:3
%            S= w{k}{rim}{np}{o};     
%            S2=imcrop(S,rect/(2^J));
%            S3=S3+numel(S2);
%            lambda1=lambda1+sum(S2(:));
%             end
%          end
%      end
%  end
% lambda=alpha*lambda1/S3;
%   err = x - y;
%   max(max(abs(err)))
lambda1=zeros(2,3);
for k=1:J
     blk=(2^(k-1))*(blk0-1)+1;
     for np=1:2
         
            for o=1:3
                
            adw{k}{1}{np}{o}=  w{k}{1}{np}{o}; 
            adw{k}{2}{np}{o}=  w{k}{2}{np}{o};
            [abs ,phase ,Ion]=imReal(w{k}{1}{np}{o},w{k}{2}{np}{o});
            D = abs;
% z=mean(D(:));
            z=conv2(D,ones(blk)/(blk^2),'same');
             
            N=D./(eps+z);
            
        if k==1
             S= D;    
             S1=imcrop(S,rect/(2^k));
             S2=conv2(S1,ones(blk)/(blk^2),'same');
             S3=S1./S2;
             counter=counter+1;
             lambda1(np,o)=alpha*mean(S3(:));
             lambda=lambda1(np,o);
        end
             
            if (k>=2)
          
              counter=counter+1;
            
              lambda=lambda1(np,o)/(sqrt(2)^k);
%               pause(0.5)
            end
            f(counter)=lambda;
          
%             P=zeros(size(D));
%             for i=1:size(D,1)
%                 for j=1:size(D,2)
%                     if (N(i,j)>0&&N(i,j)<=lambda)
%                         P(i,j)=1.5*exp(-(N(i,j).^2-(lambda)^2)./((lambda^2)*(1+lambda^2)));
%                     elseif (N(i,j)>lambda)
%                         P(i,j)=1.8*exp(-3.315./((N(i,j)./lambda).^4));
%                     elseif (N(i,j)<=0)
%                         P(i,j)=0;
%                     end
%                 end
%             end
            
            P11=(N<=lambda);
            P12=(N>0);
            P13=P12.*P11;
            P1=P13.*(0.4*exp(-(N.^2-(lambda)^2)./(eps+(lambda^2)*(1+lambda^2))));
            P22=(N>lambda);
            P2=P22.*(1.8*exp(-3.315./(eps+(N./(eps+lambda)).^4)));  
            P3=-1*((N<=0)-1);
            P=(P1+P2).*P3;
%             D=D.*P;
%             D2=D.*exp(1i.*phase);
%        
      
            adw{k}{1}{np}{o}=w{k}{1}{np}{o}.*P;
            adw{k}{2}{np}{o}=w{k}{2}{np}{o}.*P;
            waitbar(counter/(J*nither*6),mwait,char('Manipulating Wavelet coefficients at iteration: ',int2str(ither)));
            end
      end   
end
 I = icplxdual2D(adw, J, Fsf, sf);
%  A = double(I);
%                  mx = max(A(:));
%                  mn = min(A(:));
%                  I = (A-mn)/(mx-mn);
end
close(mwait);
 figure,

 
  y=I;
%  miny=min(y(:));
%  maxy=max(y(:));
%  y=(y-miny)/(maxy-miny);
% y = uint8(round(y.*255));
% mx = max(I(:));
% mn = min(I(:));
% I = (I-mn)/(mx-mn);
imshow(y);
% grayImage=I;
% lowPass =Io;
% size(Io)
% size(I)
% noiseOnlyImage = im2double(grayImage(:,:,1)) - im2double(lowPass(:,:,1));
% noiseVar = var(noiseOnlyImage);
% size(noiseVar)
% min(noiseVar)
% max(noiseVar)
figure,
title('lambda');
X = linspace(0,J*6*nither,J*6*nither)';
stem(X,f,'filled');