function output=myramp(x,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if((x*y)>1)
        output=1;
    elseif (0<=(x*y))&&((x*y)<=1)
        output=x*y;
    else
        output=0;
    end

end

