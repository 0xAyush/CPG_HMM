gf = textread('./genes.dat','%s');

gfa = char(gf); % x to char array

[r3, c3] = size(gfa);


TP = 0;
TN = 0;
FP = 0;
FN = 0;

for gi = 1:10
fname = gfa(gi,:);
while fname(length(fname)) == ' '
    fname = fname(1:length(fname) - 1);
end


fns = ['./hmr195_CpG/' fname '_CpG.dat'];

fs = fopen(fns,'r');

c = str2num(fgetl(fs));
cgi = str2num(fgetl(fs));
cgi = cgi + ones(1,length(cgi));

[TR,EMIS] = hmmestimate(c,cgi);

fclose(fs);

cgi_e = hmmviterbi(c,TR,EMIS);


for i = 1:length(cgi)
    if cgi_e(i) == cgi(i)
        % True
        if cgi_e(i) == 2
            TP = TP + 1;
        end
        if cgi_e(i) == 1
            TN = TN + 1;
        end
    else
        % False
        if cgi_e(i) == 2
            FP = FP + 1;
        end
        if cgi_e(i) == 1
            FN = FN + 1;
        end
    end
end

end

conf = [TP FP; FN TN];


success_rate = (TP + TN)/ (FP + FN + TP + TN);

disp("TR = ");
disp(TR);
disp("EMIS = ");
disp(EMIS);

disp("Confusion Matrix = ")
disp(conf);
disp("Success rate = ")
disp(success_rate);