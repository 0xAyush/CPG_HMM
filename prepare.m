
gf = textread('./genes.dat','%s');

gfa = char(gf); % x to char array

[r3, c3] = size(gfa);



for gi = 1:r3
fname = gfa(gi,:);
while fname(length(fname)) == ' '
    fname = fname(1:length(fname) - 1);
end
fns = ['./hmr195/' fname  '.dat'];
x = textread(fns,'%s');

a = char(x); % x to char array

[r2, c2] = size(a);

c = 0;

for i = 1:r2
  for j = 1:c2
    c = [c a(i,j)]; % one long string !
  end
end


length_thres = 200;
gc_thres = 0.5;
oer_thres = 0.6;

c = c(2:length(c));

I = [];

for i = 1:length(c)
  code = c(i);

  switch(code)
  case 'A'
    I(:,i) = [1]';
  case 'G'
    I(:,i) = [2]';
  case 'T'
    I(:,i) = [3]';
  case 'C'
    I(:,i) = [4]';
  end
end

wl = length_thres;

cgi = [];
cgc = [];
oer_r = [];

for i = 1:length(I) - wl
    seq = I(i: i + wl);

    CpG = 0;
    C = 0;
    G = 0;

    for j = 1:length(seq) - 1
        if seq(j) == 4 %C
            C = C + 1;
        end

        if seq(j) == 2 %G
            G = G + 1;
        end

        if seq(j) == 4 && seq(j + 1) == 2
            CpG = CpG + 1;
        end
    end

        if seq(length(seq)) == 4
        C = C + 1;
    end

    if seq(length(seq)) == 2
        G = G + 1;
    end

    GC_content = (G + C)/wl;
    OER = (CpG * wl)/(C*G);

    if GC_content >= gc_thres && OER >= oer_thres
        cgi(i) = 1;
    else
        cgi(i) = 0;
    end

    cgc(i) = GC_content;
    oer_r(i) = OER;
end

cgi = [cgi ones(1,wl)];

CpG = 0;
start = 0;
for i = 1:length(cgi)
    if CpG == 0 && cgi(i) == 1
        start = i;
        CpG = 1;
    end
    
    if CpG == 1 && (cgi(i) == 0 || i == length(cgi))
        CpG = 0;
        if i - start < length_thres
            cgi(start:i) = zeros(1,i - start + 1);
        end
    end
end

CpG = 0;
start = 0;
pos = [];
t = 1;
for i = 1:length(cgi)
    if CpG == 0 && cgi(i) == 1
        start = i;
        CpG = 1;
    end
    
    if CpG == 1 && (cgi(i) == 0 || i == length(cgi))
        CpG = 0;
        pos(t,:) = [start,i];
        t = t + 1;
    end
end

disp(["CpG locations for " fname])
disp(pos)
disp("Length = ")
disp(length(c))
fns = ['./hmr195_CpG/' fname  '_CpG.dat'];

fileID = fopen(fns,'w');
formatSpec = '%d ';
fprintf(fileID,formatSpec,I);
fprintf(fileID,'\n');
fprintf(fileID,formatSpec,cgi);
fprintf("CGI : %d, Gene: %d, Gene(2): %d",length(cgi),length(I),length(c));
fclose(fileID);
end