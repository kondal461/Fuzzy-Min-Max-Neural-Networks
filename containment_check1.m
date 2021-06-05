function z1 = containment_check1( g5,g6)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
vj=g5(1,:);
wj=g5(2,:);
vk=g6(1,:);
wk=g6(2,:);
n_fea=size(vj,2);
y=zeros(1,n_fea);
for i=1:n_fea
        if((vj(i)<vk(i))&&(vk(i)<wk(i))&&(wk(i)<wj(i)))
            y(i)=1;
        else
            y(i)=0;
        end
end
a=max(y');
z1=max(a);
end

