import pandas as pd
import matplotlib.pyplot as plt
import re
import os


log_dir = "cpu_usages"


volumes = []
max_memories = []


for filename in os.listdir(log_dir):
    if filename.endswith(".txt"):
        
        match = re.search(r"r(\d+)h(\d+)", filename)
        if match:
            r = float(match.group(1))
            h = float(match.group(2))
            volume = 2 * r * h  
            

            df = pd.read_csv(os.path.join(log_dir, filename))
            df.columns = [col.strip() for col in df.columns]
            

            max_mem = df['Memory(MB)'].max()
            

            volumes.append(volume)
            max_memories.append(max_mem)


plt.figure(figsize=(8, 6))
plt.plot(volumes, max_memories, marker='o')
plt.xlabel("2r × h (Bounding Volume, cm³)")
plt.ylabel("Max Memory Usage (MB)")
plt.title("Max Memory Usage vs Bounding Box Volume")
plt.grid(True)
plt.tight_layout()
plt.show()

