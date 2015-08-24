
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

df = pd.read_csv('progress.log', parse_dates=['time'])
df.index = df['time']
plt.figure(figsize=(4, 3))
df['words'].plot(grid=False, lw=2)
ax = plt.gca()

today_beginning = df[df['time'].dt.day == df['time'][-1].day]['words'][0]
today_current = df['words'][-1]

plt.ylim(0, max(df['words'])*1.1)
plt.text(0.6, 0.1, 'Written today: {}'.format(today_current - today_beginning),
         horizontalalignment='center',
         verticalalignment='center',
         transform = ax.transAxes,
         fontsize=14)
plt.tight_layout()
plt.savefig('progress.pdf')

