import pandas
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import sys
import mating

if len(sys.argv) > 1:
    table = pandas.read_csv(sys.argv[1], sep="\t")

else:
    data = mating.runSim("long")
    headers = data.pop(0)
    table = pandas.DataFrame(data, columns=headers)

#print(table)

plt.figure(figsize=(10,8))
plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["A", "I", "O"], palette=["blue", "green", "red"], data=table)
plot.savefig("output.png")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["M", "F"], palette =["blue", "red"], data=table)
plot.savefig("output2.png")
