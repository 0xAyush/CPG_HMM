fs = fopen('./hmr195_CpG/AF017128_CpG.dat','r');

c = str2num(fgetl(fs));
cgi = str2num(fgetl(fs));
cgi = cgi + ones(1,length(cgi));

fclose(fs);

[TR,EMIS] = hmmestimate(c,cgi);

gf = textread('./genes.dat','%s');

gfa = char(gf); % x to char array

[r3, c3] = size(gfa);



for gi = 1:5
fname = gfa(gi,:);
while fname(length(fname)) == ' '
    fname = fname(1:length(fname) - 1);
end
fns = ['./hmr195_CpG/' fname  '_CpG.dat'];
fs = fopen(fns,'r');

c = str2num(fgetl(fs));
cgi = str2num(fgetl(fs));
cgi = cgi + ones(1,length(cgi));

fclose(fs);

[TR,EMIS] = hmmtrain(c,TR,EMIS);
end

disp("TR = ")
disp(TR);
disp("EMIS = ")
disp(EMIS)