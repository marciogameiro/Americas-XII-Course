function [ab] = quadratic_cauchy_product(a,b)

m = length(a);

ab = zeros(m,1);

for n=1:m
    ab(n) = sum(a(1:n).*b(n:-1:1));
end


end

