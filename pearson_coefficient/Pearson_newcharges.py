import pandas as pd
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

experiment = pd.read_excel("tyro_charges.xlsx",sheet_name="Experiment",index_col=0)
print(experiment,"\n")

pearson = {}
spearman = {}
method = ["b3lyp_631g","b3lyp_6311g","b3lyp_def2tzvp","wB97XD_631g","wB97XD_6311g","wB97XD_def2tzvp"]
for i in method:
    qm = pd.read_excel("tyro_charges.xlsx",sheet_name=i,index_col=0)
    qm = qm.drop("AtomQM",axis=1)
    print(qm,"\n")

    for j in qm.columns:
        if j not in pearson:
            pearson[j] = []
        if j not in spearman:
            spearman[j] = []
        
        for k in experiment.columns:
            pearson[j].append(qm[j].corr(experiment[k], method='pearson'))
            spearman[j].append(qm[j].corr(experiment[k], method='spearman'))

# Setting font to Times New Roman and increasing font sizes
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'Times New Roman'
rcParams['font.size'] = 13  # You can adjust this size as needed

# Generating the plot
df = pd.DataFrame(data=pearson)

custom_palette = sns.color_palette(["#0063A6","#0063A6","#0063A6","#0063A6","#0063A6","#0063A6","#0063A6","#0063A6",
                                    "#DD4814","#DD4814","#DD4814","#94C154"]) 
plt.figure(figsize=(10,4))
sns.boxplot(data=df,palette=custom_palette)
plt.ylabel("Pearson coefficient")
plt.vlines(7.5,0.0,1.0)
plt.vlines(10.5,0.0,1.0)
plt.ylim([0.0,1.0])
plt.text(2,0.1, "Electron Density Based")
plt.text(7.7,0.1, "Electrostatic Potential Based")
plt.tight_layout()
plt.savefig("tyrosine_pearson.png")


df = pd.DataFrame(data=spearman)

custom_palette = sns.color_palette(["#0063A6","#0063A6","#0063A6","#0063A6","#0063A6","#0063A6","#0063A6","#0063A6",
                                    "#DD4814","#DD4814","#DD4814","#94C154"]) 
plt.figure(figsize=(10,4))
sns.boxplot(data=df,palette=custom_palette)
plt.ylabel("Pearson coefficient")
plt.vlines(7.5,0.0,1.0)
plt.vlines(10.5,0.0,1.0)
plt.ylim([0.0,1.0])
plt.text(2,0.1, "Electron Density Based")
plt.text(7.7,0.1, "Electrostatic Potential Based")
plt.tight_layout()
plt.savefig("tyrosine_spearman.png")
plt.show()

