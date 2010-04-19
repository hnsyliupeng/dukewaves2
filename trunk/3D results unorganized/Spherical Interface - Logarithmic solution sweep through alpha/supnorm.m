function [superr,err] = supnorm(disp, uexact, numnod,ifixu,x,y,z)
err = zeros(numnod,1);
for nod=1:numnod
%     rad(nod) = sqrt(x(nod)*x(nod) + y(nod)*y(nod) + z(nod)*z(nod));
    if(ifixu(nod)==0)
%         if(rad(nod)~=0)
            err(nod) = sqrt((uexact(nod)-disp(nod))*(uexact(nod)-disp(nod)));
%         end
    end
end

[superr,nod] = max(err);
nod