import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('results.csv')
plt.figure(figsize=(8, 5))
plt.plot(df['x'], df['u'], 'o-', label='Displacement (m)', color='b', linewidth=2)
plt.title('FEM 1D Bar Displacement')
plt.xlabel('Position X (m)')
plt.ylabel('Displacement U (m)')
plt.grid(True)
plt.legend()
plt.savefig('displacement_plot.png')
print('Plot saved as displacement_plot.png')
