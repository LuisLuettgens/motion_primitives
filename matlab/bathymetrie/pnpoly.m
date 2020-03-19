function within=pnpoly(vertx,verty,testx,testy)

%
% Copyright (c) 1970-2003, Wm. Randolph Franklin 
%

n=length(vertx);

within=false;
i=1;
j=n;

while i<=n
    
    if ((verty(i)>testy) ~= (verty(j)>testy)) && (testx < (vertx(j)-vertx(i)) * (testy-verty(i)) / (verty(j)-verty(i)) + vertx(i))
        within=~within;        
    end
    
    j=i;
    i=i+1;
end

end