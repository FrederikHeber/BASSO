load rhs_100.m g;
rand("seed", 426);
gp=g-ones(100,1)*0.025+2.*rand(100,1)*0.025;
save rhs_100_delta0.025.m gp;
