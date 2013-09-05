function f = separate(str, delim)

f={};
idx=[0 find(str==delim) length(str)+1];
for i=1:length(idx)-1
    f{i}=deblank(str(idx(i)+1:idx(i+1)-1));
end;
