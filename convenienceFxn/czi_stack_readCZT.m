function stack = czi_stack_readCZT(file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Uses bio-formats package to load .czi files   %%%
%%%   and separate out different channels. Store  %%%
%%%   all channels into a struct. Each channel is %%%
%%%   contains a 4D array (x, y, z, t).           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    info = bfopen(file);
    imdata = info{1,1};
    omeMeta = info{1,4};
    chnum = omeMeta.getPixelsSizeC(0).getValue();
    znum = omeMeta.getPixelsSizeZ(0).getValue();
    tnum = omeMeta.getPixelsSizeT(0).getValue();
    px = omeMeta.getPixelsSizeX(0).getValue();
    py = omeMeta.getPixelsSizeY(0).getValue();
    
    stack = struct();
    for i = 1:chnum
        varName = ['c' int2str(i)];
        imzt = zeros(py, px, znum, tnum);
        for j = 1:tnum
            for k = 1:znum
                planeNum = (j-1)*chnum*znum+(k-1)*chnum+i;
                imzt(:,:,k,j) = imdata{planeNum, 1};
            end
        end
        stack.(varName) = imzt;
    end     
end