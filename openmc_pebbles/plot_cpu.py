import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("cpu_usage_log.txt")
df.columns = [col.strip() for col in df.columns]  


time = df['Time(s)'] - df['Time(s)'].iloc[0]  
cpu = df['CPU%']
mem = df['Memory(MB)']


fig, ax1 = plt.subplots(figsize=(10, 6))


ax1.plot(time, cpu, color='tab:red', linewidth=2, label='CPU Usage (%)')
#ax1.set_yscale('log')
#ax1.set_xscale('log')
ax1.set_xlabel('Time (s)', fontsize=12)
ax1.set_ylabel('CPU Usage (%)', color='tab:red', fontsize=12)
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.grid(True, which='both', linestyle='--', linewidth=0.5)


ax2 = ax1.twinx()
#ax2.set_yscale('log')
ax2.plot(time, mem, color='tab:blue', linewidth=2, label='Memory Usage (MB)')
ax2.set_ylabel('Memory Usage (MB)', color='tab:blue', fontsize=12)
ax2.tick_params(axis='y', labelcolor='tab:blue')


plt.title('CPU and Memory Usage Over Time', fontsize=14)
fig.tight_layout()
plt.show()

