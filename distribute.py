dat = ""

with open("./hmrfasta.txt","+r") as fs:
    dat = fs.read()
fs.close()

labels = []
plots = []

lab = False
plot = False
buildLabel = ""
buildPlot = ""

for i in dat:
    if lab and i != '\n':
        buildLabel += i
        continue
    
    elif lab and i == '\n':
        lab = False
        plot = True
        labels.append(buildLabel)
        buildLabel = ""
        continue
    
    if i == '>' and plot:
        lab = True
        plot = False
        if buildPlot[-1] == '\n': buildPlot = buildPlot[:-1]
        plots.append(buildPlot)
        buildPlot = ""
        continue
    elif i == '>':
        lab = True
        continue
    
    if plot and i != '>':
        buildPlot += i
        continue

plots.append(buildPlot)

labelFile = ""

for i in range(len(labels)):
    labelFile += labels[i] + "\n"
    with open("./hmr195/" + labels[i] + ".dat","+w") as fs:
        fs.write(plots[i])
    fs.close()
labelFile = labelFile[:-1]

for i in range(len(labels)):
    with open("./hmr195_CpG/" + labels[i] + "_CpG.dat","+w") as fs:
        fs.write("")
    fs.close()

with open("./genes.dat","+w") as fs:
    fs.write(labelFile)
fs.close()

print("Done")