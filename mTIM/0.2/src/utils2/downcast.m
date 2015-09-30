function baseclass_object=downcast(object, baseclass_name)
% baseclass_object=downcast(object, baseclass_name)

if equal(class(object), baseclass_name),
  baseclass_object=object ;
  return ;
end ;

object=struct(object) ;
names=fieldnames(object) ;

% find iteratively the base_class
for i=length(names),
  name=names{i} ;
  if ~equal(name, 'dummy'),
    baseclass_object=downcast(getfield(object, name), baseclass_name) ;
    if ~isempty(baseclass_object)
      return ;
    end ;
  end ;
end ;
baseclass_object=[] ;