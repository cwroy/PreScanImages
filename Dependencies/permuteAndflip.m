function y = permuteAndflip(x,R)

f=sign(R);
y=permute(x,abs(R));
for d=1:length(f)
   if f(d)==-1
       y=flip(y,d);
   end
end


end