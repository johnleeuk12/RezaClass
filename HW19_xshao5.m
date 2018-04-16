% This HW is to use SVD PCA to analyze the illumination imgages. 
% Using SVD on any data set, eiganvectors corresponding to Principal 
% components can be found through the left sigular vector, and 
% eiganvalues can be found through the sigular values. The function 
% of below at the bottom of the main fucntion. For the yale face data 
% set, 2 PCs are found through the PCA funtion. The mean face of 
% each individual shows an average of all the images, and 
% illumination is linear. For the first PC, it captures the light 
% source changing intensity at the front of the face. The second PC 
% captures the light source changing from the right and the
% left side of the face. 
function HW19_xshao5
for j=1:3 %repeat the process for the 3 individuals
    imgvect=[]; %vector for all the images 
    for i=1:64 %uses all 64 of the images through the load image file
        img=loadimage(j,i); 
        imgvect=[imgvect img(:)];
    end
    [m n]=size(img);
    [mu,Ud,Y]=pcaMS(imgvect,2); %Run SVD PCA on all the images
    meanface=reshape(mu,[m,n]); %reshape the image for mean face
    %Do the same for the 2 eiganfaces
    eigenface1=reshape(Ud(:,1),m,n); a=min(Ud(:,1)); b=max(Ud(:,1));
    eigenface2=reshape(Ud(:,2),m,n);c=min(Ud(:,2));d=max(Ud(:,2));
    %use the variances along PC1 and PC2, and project images along PC1 
    %and PC2
    y1=-std(Y(1,:)):0.2*std(Y(1,:)):std(Y(1,:));
    X1=mu+Ud(:,1)*y1;
    y2=-std(Y(2,:)):0.2*std(Y(2,:)):std(Y(2,:));
    X2=mu+Ud(:,2)*y2;
    
    %plot all the needed information for each individual
    figure();
    subplot(2,2,1);
    imshow(uint8(meanface)); 
    title(['Mean Face of Individual' num2str(j,'%i')]);
    subplot(2,2,2);
    imshow(eigenface1,[a b]);
    title(['EiganFaces U1 of Individual' num2str(j,'%i')]);
    subplot(2,2,3);
    imshow(eigenface2,[c d]);
    title(['EiganFaces U2 of Individual' num2str(j,'%i')]);
    
    figure();
    x0=10;
    y0=5;
    width=1000;
    height=100;
    set(gcf,'units','points','position',[x0,y0,width,height])
    for i=1:length(y1)
        pic=reshape(X1(:,i),m,n);
        subplot(1,length(y1),i);
        imshow(uint8(pic));
    end
    a = axes;
    t1 = title(['faces captured of First PC of Individual' num2str(j,'%i')]);
    a.Visible = 'off'; % set(a,'Visible','off');
    t1.Visible = 'on'; % set(t1,'Visible','on')
    
    figure()
    x0=10;
    y0=5;
    width=1800;
    height=100;
    set(gcf,'units','points','position',[x0,y0,width,height])
    for i=1:length(y2)
        pic=reshape(X2(:,i),m,n);
        subplot(1,length(y2),i);
        imshow(uint8(pic));
    end
    a = axes;
    t1 = title(['faces captured of Second PC of Individual' num2str(j,'%i')]);
    a.Visible = 'off'; % set(a,'Visible','off');
    t1.Visible = 'on'; % set(t1,'Visible','on')
end
end
function [mu,Ud,Y]=pcaMS(X,d)
mu=mean(X')';
X_tilta=X-mu;
[U,S,V]=svd(X_tilta,'econ');
Ud=U(:,1:d);
Y=Ud'*X_tilta;
end