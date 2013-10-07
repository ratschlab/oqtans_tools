function plot_priorities(username)
% plot_priorities(username)
% 
% plots a history of the priorities of all users and contrasts it
% with the priority histogram of the given username (default: whoami)

if nargin<1,
  username=whoami ;
end ;

[ret, all_list]=unix('qstat -pri -u "*" | grep qw') ;
[ret, user_list]=unix(sprintf('qstat -pri -u "%s" | grep qw', username)) ;

[all_pri,username_list]=get_pri(all_list) ;
user_pri=get_pri(user_list) ;

[a,b]=hist(all_pri, 20) ;
[au,bu]=hist(user_pri, b) ;

figure;
bar(b,a,'r'); hold on; bar(bu,au,'b')
legend({'all users', sprintf('user %s', username)}) ;

ul=unique(username_list) ;
for i=1:length(ul),
  idx = find(ismember(username_list, ul{i})) ;
  fprintf('%10s\tmin: %4i\taverage: %4i\tmax: %4i\n', ul{i}, min(all_pri(idx)), mean(all_pri(idx)), max(all_pri(idx))) ;
end ;

return 

function [pri, user_list]=get_pri(list)

lines=separate(list, sprintf('\n')) ;
pri=zeros(1,length(lines)) ;
user_list={} ;
for i=1:length(lines)
  items=separate(lines{i}, ' ', 1) ;
  if isempty(items), pri(i)=nan ; continue ; end ;
  pri(i)=str2num(items{6}) ;
  user_list{i}=items{8} ;
  if length(items)>=13,
    pars=sscanf(items{13}, '%i-%i:%i') ;
    num=pars(2)-pars(1) ;
    pri(end+1:end+num) = pri(i) ; 
  end ;
end ;
pri(isnan(pri))=[] ;



