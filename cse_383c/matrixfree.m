classdef matrixfree < handle
  properties 
    epsilon
    O
    s
    b
    m
    n
    N
  end
  
  methods
    function o = matrixfree(epsilon)
       load mf;
       o.n = n;
       o.m = m;
       o.O = O;
       o.s = s;
       o.b = b;
       o.N = N;
       o.epsilon = epsilon;
    end
   
       
   function out = A(o,in)
     [~,k]=size(in);
     out=zeros(o.m,k);
     for j=1:k
       t=dct2(o.v2i(in(:,j)));
       out(:,j) = o.s.*(t(o.O));
     end
   end

   function out = At(o,in)
     [~,k]=size(in);
     out = zeros(o.n,k);
     for j=1:k
       t = o.s .* in(:,j);
       out(o.O,j)=t;
       out(:,j)=o.i2v(idct2(o.v2i(out(:,j))));
     end
   end
   
   function out=vis(o,in)
     imshow(o.v2i(in));
   end
   
   function out=v2i(o,in) 
     out=reshape(in,o.N,o.N);
   end
   
   function out=i2v(o,in)
     out=in(:);
   end
   
  end % methods
end % class
  

