f=zeros(100,1);
f(floor(100*0.2):floor(100*0.22))=-5;
f(floor(100*0.3):floor(100*0.32))=5;
f(floor(100*0.5):floor(100*0.52))=7;
f(floor(100*0.8):floor(100*0.84))=3;
save solution_100.m f;
