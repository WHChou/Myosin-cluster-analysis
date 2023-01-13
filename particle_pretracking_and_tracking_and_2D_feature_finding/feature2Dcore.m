function r = feature2Dcore(img, lambda, w, xcoor, ycoor)
%%%%%%%%%%%%%%%
% Given a list of speckle positions, this script takes the core codes from 
% feature2D to calculate sub-pixel speckle position, integrated mass, 
% eccentricity and radius of gyration.
%
% Input: 
%       img: (nx,ny) array which presumably contains some features worth finding
%       lambda: length scale of noise to be filtered out, in pixels; typically 1
%       w: a parameter which should be a little greater than the radius of the 
%           largest features in the image.
%       xcoor, ycoor: Speckle position from qFSM (Note that qFSM outputs as (y, x))
% Output: 
%       r(:,1):	the x centroid positions, in pixels.
% 		r(:,2): the y centroid positions, in pixels. 
% 		r(:,3): integrated brightness of the features. ("mass")
% 		r(:,4): the square of the radius of gyration of the features.
% 		    (second moment of the "mass" distribution, where mass=intensity)
% 		r(:,5): eccentricity, which should be zero for circularly symmetric features and 
%                   order one for very elongated images.
%       r(:,6): integrated brightness of the features from the raw image. ("mass")
% Wen-hung Chou 2022.02.21
%%%%%%%%%%%%%%%
    extent=2*w+1;
    sub_img = img(ycoor-10:ycoor+10, xcoor-10:xcoor+10);  % Crop out image centered around particle
    image = bpass(sub_img,lambda,w);
    a=image;
    nmax=length(xcoor);         
    
    x = 11; y = 11;  % Particle position always around center
    m=zeros(length(x),1);
    m2=zeros(length(x),1);
    xl = x - fix(extent/2);
    xh = xl + extent -1;        

    %       Set up some masks
    rsq = rsqd( extent,extent );
    t = thetarr( extent );

    mask = le(rsq,(extent/2)^2);  
    mask2 = ones(1,extent)'*[1:extent];
    mask2 = mask2.*mask;           
    mask3= (rsq.*mask) + (1/6);
    cen = (extent-1)/2 +1;           
    % cmask = vpa(cos(sym('2')*t)).*mask;  
    % ultra high presision, since Matlab and IDL differ here
    cmask = cos(2*t).*mask;
    % smask = vpa(sin(sym('2')*t)).*mask;
    smask = sin(2*t).*mask;
    cmask(cen,cen) = 0.0;
    smask(cen,cen) = 0.0;  

    suba = zeros(extent, extent, nmax);
    subImg = zeros(extent, extent, nmax);
    xmask = mask2;
    ymask = mask2';
    yl = y - fix(extent/2);
    yh = yl + extent - 1;           
    yscale = 1;
    ycen = cen;                   

    %	Setup some result arrays
    xc = zeros(nmax,1);
    yc = zeros(nmax,1);
    rg = zeros(nmax,1);
    e  = zeros(nmax,1);
    %	Calculate feature centers & estimate mass
    for i=1:nmax,
        xc(i) = sum(sum(double(a(fix(yl(i)):fix(yh(i)),fix(xl(i)):fix(xh(i)))).*xmask));  
        yc(i) = sum(sum(double(a(fix(yl(i)):fix(yh(i)),fix(xl(i)):fix(xh(i)))).*ymask));
        m(i) = sum(sum(double(a(fix(yl(i)):fix(yh(i)),fix(xl(i)):fix(xh(i)))).*mask)); 
    end

    %	Correct for the 'offset' of the centroid masks
    xc = xc./m - ((extent+1)/2);             
    yc = (yc./m - (extent+1)/2)/yscale;  
    %	Update the positions and correct for the width of the 'border'
    x = x + xc - 0*fix(extent/2);                                           %%%??????? 0*
    y = ( y + yc - 0*fix(extent/2) ) * yscale;
    x2=x;
    y2=y;
    
    for i = 1:nmax
        suba(:,:,i) = fracshift( double(a(fix(yl(i)):fix(yh(i)),fix(xl(i)):fix(xh(i)))), -xc(i) , -yc(i) );
        m(i) = sum(sum(( suba(:,:,i).*mask )));             % mass
        rg(i) = (sum(sum( suba(:,:,i).*mask3 ))) / m(i);    % squared radius of gyration
        tmp = sqrt(( (sum(sum( suba(:,:,i).*cmask )))^2 ) +( (sum(sum( suba(:,:,i).*smask )))^2 )); 
        tmp2 = (m(i)-suba(cen,ycen,i)+1e-6);
        e(i) = tmp/tmp2;    
        subImg(:,:,i) = fracshift( double(sub_img(fix(yl(i)):fix(yh(i)),fix(xl(i)):fix(xh(i)))), -xc(i) , -yc(i) );
        m2(i) = sum(sum(( subImg(:,:,i).*mask ))); 

        xc(i) = sum(sum(double(suba(:,:,i)).*xmask));  
        yc(i) = sum(sum(double(suba(:,:,i)).*ymask));
    end
    xc = xc./m - ((extent+1)/2);                                           %%% Why correct again??
    yc = (yc./m - (extent+1)/2)/yscale;  %get mass center
    x3 = x2 + xc - 0*fix(extent/2);
    y3 = ( y2 + yc - 0*fix(extent/2) ) * yscale;
    
    x4 = x3+xcoor-11; % Sub-pixel correction for particle position (translate back to original coordinate)
    y4 = y3+ycoor-11; 
    r = [x4 ,y4, m, rg, e, m2];
end