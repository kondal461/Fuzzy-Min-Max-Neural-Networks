function [ ej_ccn ] = CCN_ActivationFun( Ah,Vj,Wj,gamma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

bj_ccn=MembershipFunc(Ah,Vj,Wj,gamma);
ej_ccn=-1*(mystep(bj_ccn-1));

end

