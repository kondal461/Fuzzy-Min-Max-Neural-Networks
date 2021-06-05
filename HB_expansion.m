function z = HB_expansion( inx,iny )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    n_feat = size(inx,2);
    vj=inx(1,:);
    wj=inx(2,:);
    vk=iny(1,:);
    wk=iny(2,:);
    
    z1=0; 
    for i=1:n_feat
       z1= sum(max(wj(i),vk(i))-min(vj(i),vk(i)));
    end
    z=str2num(sprintf('%.3f',z1));
    %z1=sum(max(x,[],1)-min(y,[],1));
end

