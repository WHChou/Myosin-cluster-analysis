function [M2,MT,Mrej] = mpretrack_wc(img, featsize, barint, barrg, barcc, IdivRg, masscut, Imin, field, draw)
 
% Used to find the parameters (minimum intensity, maximum Rg, maximum eccentricity, feature
% fize, masscut, Imin). These parameters should then be passed to mpretrack.
%
% Note: feature2D requires the Image processing toolbox due to a call to
% imdilate in the localmax subfunction. Use feature2D_nodilate for an
% alternative which works almost as well.
%
% INPUTS :
% im - The image that contains particles
% featsize - The size of the feature you want to find.
% barint - The minimum intensity you want to accept.
% barrg - The maximum Rg squared you want to accept.
% barcc - The maximum eccentricity you want to accept.
% IdivRg - minimum ratio of Intensity/pixel to be accepted (integrated
%           intensity / Rg squared of feature)
% Imin - (optional) the minimum intensity for a pixel to be considered as a potential
%           feature. Defaults to top 30%; if most pixels are 0, discard
%           bottom 20% instead
% masscut - (optional) the masscut parameter for feature2D to remove false positives
%           before rifining the position to speed up the code.
% field - (optional) set to 0 or 1 if image is actually odd or even field of an interlaced 
%           image. All the masks will then be constructed with a 2:1 aspect ratio. Otherwise 
%           set to 2 for progressive scan cameras. Defaults to 2.
% draw - (optional) set to 1 if you want to show plot of images and particles. Defaults to 0.
%
% commented out:
% inv - A logical for inverting the image (1 inverts, 0 doesn't) to look
%           for dark features instead of bright ones.
%
% Also, the program looks for the files "fov#_times.mat" for the "time"
% variable and the images files "fov#\fov#_####.tif" from the basepath.
%
% OUTPUTS 
% 
% M2 - All the features found from calling feature2D
% MT - All the features from feature2D which were accepted given the
% criteria from the inputs
%
% The following appear on the command line as output:
%
% - Number of features found (prior to eliminating anything).
% - A matrix with the positions in columns 1-2 and the parameters int,rg
% and cc in columns 3-5. (accepted features first, then rejected features)
% - Min I, Max Rg, Max Cc these values help you sum up quickly what is in
% the matrix that appears. 
% - A plot showing the image and dots for rejected features and circles for
% accepted features
%
% REVISION HISTORY
% written by Paul Fournier and Vincent Pelletier (Maria Kilfoil's group),
% latest revision 10/18/07
% 10/26/07 Vincent -- commented out the Inv keyword, added a ratio of
% Iint to Rg parameter
% 12/21/07 Maria - added optional field parameter
% 08/04/21 Wen-hung - removed reading image in folder; directly input image
% matrix; commented out draw circle
% 04/08/22 Wen-hung - output all 6 columns from feature2D. Make additional
% function input to control whether to show particles plot. Plotted circles
% now reflect real featsize

if nargin <10, draw = 0; end
if nargin < 9, field = 2; end
if nargin < 8, Imin = 0; end
if nargin < 7, masscut = 0; end

d=0;
MT=[];
Mrej=[];


%     if inv == 1
%         img2=img;
%         img = 255-img;
%     end

    M = feature2D(img,1,featsize,masscut,Imin,field,0);

    a = length(M(:,1));
    
    M2=M;
     for i=1:a
         if ((M(i,5)>barcc))
%              disp('Fail CC')
            Mrej=[Mrej; M(i,:)];
             M(i,1:6)=0;
%          end

         elseif ((M(i,4)>barrg))
%              disp('Fail RG')
            Mrej=[Mrej; M(i,:)];
             M(i,1:6)=0;
%          end

         elseif ((M(i,3)<barint))
%              disp('Fail Int')
            Mrej=[Mrej; M(i,:)];
             M(i,1:6)=0;
         elseif((M(i,3)/M(i,4)<IdivRg))
%              disp('Fail IdivRg')
            Mrej=[Mrej; M(i,:)];
             M(i,1:6)=0;
         end    
     end
    
%    Deleting the zero rows

    M=M(M(:,1)~=0,:);

    a = length(M(:,1));

    MT(d+1:a+d, 1:6)=M(1:a,1:6);   
%     MT(d+1:a+d, 6)=frame;
%     MT(d+1:a+d, 7)=time(frame);

    if draw
        figure;
         %imagesc(img,[0,255]),colormap(gray);
        imshow(uint8(img)),colormap gray;
        %imshow(uint8(bpass(img, 1, 7)), [])
        hold on

        % Making a circle the size of the feature around each feature.

        theta = 0:0.001:2*pi;
        for c = 1:length(M(:,1))
            cx = M(c,1) + featsize*cos(theta);%*2;
            cy = M(c,2) + featsize*sin(theta);%*2;
            plot(cx,cy,'g-','linewidth',1.5)
        end
        if( ~isempty(Mrej)>0 )
             plot( Mrej(:,1), Mrej(:,2), 'r.' );
        end
    end
%     if inv == 1
%         figure;
%         imagesc(img2,[0,255]),colormap(gray);
%         hold on
% 
%         % Making a circle the size of the feature around each feature.
%         for c = 1:length(M(:,1))
%             cx = M(c,1) + featsize*cos(theta);
%             cy = M(c,2) + featsize*sin(theta);
%             plot(cx,cy,'g-','linewidth',1.5)
%         end
%         if( ~isempty(Mrej)>0 )
%             plot( Mrej(:,1), Mrej(:,2), 'r.' );
%         end
%     end
%    axis equal;
    
    format short g
%     disp(M)
    disp(['Kept : ' num2str(size(M,1))])
%     disp(Mrej)
    disp(['Minimum Intensity : ' num2str(min(M(:,3)))])
    disp(['Maximum Rg : ' num2str(max(M(:,4)))])
    disp(['Maximum Eccentricity : ' num2str(max(M(:,5)))])
    clear img;
    clear R;
    clear pic;
    clear X;
    clear t;
    clear i;
    clear j;



