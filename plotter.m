#!/usr/bin/env octave

if nargin < 1
  printf("No HDF5 filename specified")
  exit(1)
endif

set(gca,'FontSize',18)
arg_list = argv()
load("-hdf5", arg_list{1})

for i = 1:size(u,2)
  plot(x,u(:,i),'ko-')
  title(sprintf("t = %g", t(i)))
  xlabel("x")
  ylabel("u(x,t)")
  drawnow
end
