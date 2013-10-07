function L=structstrrep(L, orig, target) ;
% L=structstrrep(L, orig, target) ;

fn = fieldnames(L) ;
for k=1:prod(size(L))
  for i=1:length(fn),
    if size(L(k).(fn{i}),1)==1 && ischar(L(k).(fn{i})),
      L(k).(fn{i}) = strrep(L(k).(fn{i}), orig, target) ;
    elseif isstruct(L(k).(fn{i}))
      L(k).(fn{i})=structstrrep(L(k).(fn{i}), orig, target) ;
    end ;
  end ;
end ;


