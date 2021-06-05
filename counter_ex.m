 N=1000;
c=[2 3 5 7 11 13];
counter=0;
ctr=1;
for x=1:N
    if sum(rem(x,c))==0
        counter=counter+1;
        P(ctr)=x;
        ctr=ctr+1;
    else
        continue;
    end
end