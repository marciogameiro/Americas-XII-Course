function [s]=quadratic_sum_complex(a,b)

m=(length(a)+1)/2;

B=zeros(2*m-1);

tb=[zeros(m-1,1);flip(b,1);zeros(m-1,1)].';

for k=-m+1:m-1
    B(k+m,:)=tb((-m+1:m-1)-k+(2*m-1));
end

s=B*a;

end