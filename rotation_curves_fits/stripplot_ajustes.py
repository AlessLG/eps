import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# sns.set_theme()
df = pd.read_csv('resultados.csv')

# df = df[df['BP'] == True]
df = df[df['BP'] == True]
x2 = df[['X2red_sim', 'X2red_est', 'X2red_bek']]
rmse = df[['RMSE_sim', 'RMSE_est', 'RMSE_bek']]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,8), dpi=300)

sns.stripplot(data=x2, orient='h', jitter=0.3, ax=ax1, palette=['red', 'blue', 'green'])
sns.stripplot(data=rmse, orient='h', jitter=0.3, ax=ax2, palette=['red', 'blue', 'green'])

ax1.set_xlim((-10, 200))
ax1.set_xlabel(r'$\chi^2_{\text{reducido}}$')
ax2.set_xlabel('RMSE [km/s]')

ax1.set_yticks(['X2red_sim', 'X2red_est', 'X2red_bek'])
ax1.set_yticklabels(['Simple', 'Estándar', 'Bekenstein'])

ax2.set_yticks(['RMSE_sim', 'RMSE_est', 'RMSE_bek'])
ax2.set_yticklabels(['Simple', 'Estándar', 'Bekenstein'])

fig.tight_layout()
fig.savefig('stripplot.pdf')