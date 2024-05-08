import pandas as pd
import sqlite_requests2 as sr
import numpy as np
import matplotlib.pyplot as plt


table_ = sr.table_download('Aquasearch_study', 'Spectrums_table_filtered')
table_ = pd.DataFrame(table_, columns = ['mz', 'intensity', 'standard', 'protein', 'sample', 'peptide'])


chicken = table_[table_['protein'] == '1']
chicken_stand = chicken[chicken['standard'] == '1']
chicken_mix = chicken[chicken['standard'] == '0']
plt.scatter(chicken_mix['mz'], chicken_mix['intensity'], s=10)
plt.scatter(list(chicken_stand['mz']),list(chicken_stand['intensity']), s=10)
plt.xlabel("mz")
plt.ylabel("relative int")
plt.title("Chicken albumin")
plt.show()


chicken = table_[table_['protein'] == '5']
chicken_stand = chicken[chicken['standard'] == '1']
chicken_mix = chicken[chicken['standard'] == '0']
plt.scatter(list(chicken_mix['mz']), list(chicken_mix['intensity']), s=10)
plt.scatter(list(chicken_stand['mz']),list(chicken_stand['intensity']), s=10)
plt.xlabel("mz")
plt.ylabel("relative int")
plt.title("Pig albumin")
plt.show()


chicken = table_[table_['protein'] == '3']
chicken_stand = chicken[chicken['standard'] == '1']
chicken_mix = chicken[chicken['standard'] == '0']
plt.scatter(list(chicken_mix['mz']), list(chicken_mix['intensity']), s=10)
plt.scatter(list(chicken_stand['mz']),list(chicken_stand['intensity']), s=10)
plt.xlabel("mz")
plt.ylabel("relative int")
plt.title("Rat ablumin")
plt.show()


chicken = table_[table_['protein'] == '2']
chicken_stand = chicken[chicken['standard'] == '1']
chicken_mix = chicken[chicken['standard'] == '0']
plt.scatter(list(chicken_mix['mz']), list(chicken_mix['intensity']), s=10)
plt.scatter(list(chicken_stand['mz']),list(chicken_stand['intensity']), s=10)
plt.xlabel("mz")
plt.ylabel("relative int")
plt.title("Mounse albumin")
plt.show()


chicken = table_[table_['protein'] == '4']
chicken_stand = chicken[chicken['standard'] == '1']
chicken_mix = chicken[chicken['standard'] == '0']
plt.scatter(list(chicken_mix['mz']), list(chicken_mix['intensity']), s=10)
plt.scatter(list(chicken_stand['mz']),list(chicken_stand['intensity']), s=10)
plt.xlabel("mz")
plt.ylabel("relative int")
plt.title("Human albumin")
plt.show()