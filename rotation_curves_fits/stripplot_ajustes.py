import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set_theme()
df = pd.read_csv('resultados.csv')

# df = df[df['BP'] == True]
df = df[df['BP'] == True]
df = df.drop(columns='BP')

plt.figure(figsize=(8,4), dpi=300)
sns.stripplot(data=df, orient='h', jitter=0.3)
plt.xlim((-10, 200))
plt.xlabel(r'$\chi^2$')
plt.tight_layout()
plt.savefig('stripplot.pdf')