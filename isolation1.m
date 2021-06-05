function z = isolation1(x,y,n_fea)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

n=size(y,2);
y1=zeros(n,n_fea);
for j=1:n
    vj=x{1}(1,:);
    wj=x{1}(2,:);
    vk=y{j}(1,:);
    wk=y{j}(2,:);
    if(vj==wj)
        y1(j,:)=1;
    elseif(vk==wk)
        y1(j,:)=1;
    elseif(vk~=wk)
        for i=1:n_fea
            if((wj(i)<vk(i))||(wk(i)<vj(i)))
                y1(j,i)=1;
            else
                y1(j,i)=0;
            end
        end
    end
end
z=min(max(y1'));
end

