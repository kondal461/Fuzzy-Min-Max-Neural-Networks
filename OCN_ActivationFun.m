function [ dp1,dp2 ] = OCN_ActivationFun( Ah,x,y,Vj,Wj,gamma,n_fea)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
V1=x(1,:);
W1=x(2,:);
V2=y(1,:);
W2=y(2,:);
bj_ocn=MembershipFunc(Ah,Vj,Wj,gamma);
for p=1
    for i=1:n_fea
        summ(i)=max((Ah(i)/W1(i)),(V1(i)/Ah(i)));
    end
    summ=sum(summ);
    dp1=mystep((bj_ocn-1)*(-1+(1/n_fea)*summ));
end
for p=2
    for i=1:n_fea
        summ(i)=max((Ah(i)/W2(i)),(V2(i)/Ah(i)));
    end
    summ=sum(summ);
    dp2=mystep((bj_ocn-1)*(-1+(1/n_fea)*summ));
end

end

