function post_process

fname='anytown_0p99.dat';

d=importdata(fname);

[val,idx]=sort(d.data(:,end),'descend');
idx=idx+1;

fprintf('\n\n Megvaltozasok %%-ban (legnagyobb : %g = 100%%)',val(1));

val(1:8)'/val(1)*100
d.textdata(idx(1:8))'

end