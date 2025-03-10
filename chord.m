function z= chord(m,p,x) %mean camber line function
    n=length(x);
    z=zeros(n,1);
    for i=1:n
        if x(i)<=p
            z(i)=m/(p^2)*(2*p*x(i)-(x(i))^2);
        else
            z(i)=m/((1-p)^2)*(1-2*p+2*p*x(i)-(x(i))^2);
        end

    end
end



