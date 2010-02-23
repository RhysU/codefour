
a=dir('sol*.txt');
load x.txt
set(gca,'FontSize',18)

for i = 1:length(a)
  u = load(a(i).name);
  plot(x,u,'ko-')
  drawnow
end
