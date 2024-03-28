fs = fopen('./hmr195_CpG/AF017128_CpG.dat','r');

c = str2num(fgetl(fs));
cgi = str2num(fgetl(fs));
cgi = cgi + ones(1,length(cgi));

fclose(fs);

TR = [0.9997    0.0003;0.0029    0.9971];
EMIS = [0.2690    0.2443    0.2315    0.2553;0.1963    0.3240    0.2002    0.2795];

cgi_e = hmmviterbi(c,TR,EMIS);

plot(cgi_e);hold;plot(cgi);