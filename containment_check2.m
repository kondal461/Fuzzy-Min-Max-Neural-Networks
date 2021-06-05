function z1 = containment_check2( g5,g6)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
vj=g5(1,:);
wj=g5(2,:);
vk=g6(1,:);
wk=g6(2,:);
n_fea=size(vj,2);
y=zeros(1,n_fea);
y1=zeros(1,n_fea);
for i=1:n_fea
    if((vj(i)<vk(i))&&(vk(i)<wk(i))&&(wk(i)<wj(i)))
        y(i)=1;
    else
        y(i)=0;
    end
end
for i=1:n_fea
    if((vk(i)<vj(i))&&(vj(i)<wj(i))&&(wj(i)<wk(i)))
        y1(i)=1;
    else
        y1(i)=0;
    end
end
a=max(y');
b=max(y1');
if(a==b)
    
    z1=1;
else
    z1=0;
end

end

