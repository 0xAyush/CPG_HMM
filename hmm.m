fs = fopen('./hmr195_CpG/AB003306_CpG.dat','r');

c = str2num(fgetl(fs));
cgi = str2num(fgetl(fs));
cgi = cgi + ones(1,length(cgi));

fclose(fs);

[TR,EMIS] = hmmestimate(c,cgi);

cgi_e = hmmviterbi(c,TR,EMIS);




TP = 0;
TN = 0;
FP = 0;
FN = 0;

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

conf = [TP FP; FN TN];

pstates = hmmdecode(c,TR,EMIS);
log_odds = [];

for i = 1:length(pstates(1,:))
    log_odds(i) = log(pstates(2,i)/pstates(1,i));
end

t = min(log_odds);
t_max = max(log_odds);

dt = 0.1;

roc_x = [];
roc_y = [];
j = 1;
while t < t_max
TP = 0;
TN = 0;
FP = 0;
FN = 0;

for i = 1:length(cgi)
    if log_odds(i) > t
        % Positive
        if cgi(i) == 2
            TP = TP + 1;
        else
            FP = FP + 1;
        end
    else
        %Negative
        if cgi(i) == 2
            FN = FN + 1;
        else
            TN = TN + 1;
        end
    end
end

TPR = TP / (TP + FN);
FPR = FP / (FP + TN);
roc_x(j ) = [FPR];
roc_y(j) = [TPR];
j = j + 1;
t = t + dt;
end

roc_area = 0;

for i = 1:length(roc_y)
    roc_area = roc_area + (1/length(roc_y))*roc_y(i);
end

success_rate = (TP + TN)/(TP + TN + FP + FN);

hold on;

subplot(2,2,1);
plot(cgi_e);
hold on;
plot(cgi);

disp("TR = ");
disp(TR);
disp("EMIS = ");
disp(EMIS);

subplot(2,2,2);

plot(log_odds);

disp("Confusion Matrix = ")
disp(conf);
disp("Success rate = ")
disp(success_rate);

subplot(2,2,3);

disp("Advantage over random classifier");
disp(roc_area - 0.5);
plot(roc_x,roc_y);
hold on;
plot((1:length(cgi))/length(cgi),(1:length(cgi))/length(cgi));