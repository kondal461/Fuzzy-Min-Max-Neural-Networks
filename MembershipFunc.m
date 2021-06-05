function [ bj ] = MembershipFunc( Ah,Vj,Wj,gamma )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n_fea=size(Ah,2);
B=zeros(1,n_fea);
if (Ah==inf(1,n_fea))|(Vj==inf*ones(1,n_fea))|(Wj==inf(1,n_fea))
    bj=0;
else
    for i=1:n_fea
        B(i)=min((1-myramp(Ah(i)-Wj(i),gamma)),(1-myramp(Vj(i)-Ah(i),gamma)));
    end
    B1=min(B);
    bj=str2num(sprintf('%.4f',B1));
end
end

