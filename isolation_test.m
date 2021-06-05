
function [ b1 ] = isolation_test( x,y )

%n1=size(x,2);
%n2=size(y,2);

n_fea=size(x,2);    %{n1}(1,:),2);
% x=zeros(2,n_fea);
% y=zeros(2,n_fea);
y1=zeros(1,n_fea);  %n2,n_fea);
vj= x(1,:);         %x{n1}(1,:);
wj= x(2,:);         %x{n1}(2,:);
vk= y(1,:);         %y{n2}(1,:);
wk= y(2,:);         %y{n2}(2,:);
if(vj==wj)
    b1=0;
%     b1=str2double(sprintf('%.3f',b1));
else
    for i=1:n_fea
        if((wj(i)<vk(i))||(wk(i)<vj(i)))
            y1(i)=1;
        else
            y1(i)=0;
        end
    end
    b1=max(y1);
%     b1=str2double(sprintf('%.3f',b1));
end
end

