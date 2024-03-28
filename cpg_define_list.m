x = textread('./hmr195/AF017128.dat','%s');

a = char(x); % x to char array

[r2, c2] = size(a);

c = 0;


for i = 1:r2
  for j = 1:c2
    c = [c a(i,j)]; % one long string !
  end
end

length_thres = 200;
gc_thres = 0.55;
oer_thres = 0.65;

c = c(2:length(c));

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

disp("CpG locations = ")
disp(pos)
disp("Length = ")
disp(length(c))
hold
plot(cgi)
plot(cgc)
plot(oer_r)
plot(ones(1,length(cgc))*gc_thres)
plot(ones(1,length(oer_r))*oer_thres)