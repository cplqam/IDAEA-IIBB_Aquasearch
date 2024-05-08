import sklearn
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import mglearn


col_names = list(dic.values())[0][0]
col_names = list(col_names.keys())
tabla_ = np.array([])

for sample in dic.values():
    score = sample[0]
    l_scores = np.array(list(score.values()))
    
    if tabla_.size == 0:
        tabla_ = np.array([l_scores])
    else:
        tabla_ = np.vstack([tabla_,l_scores])
    
tabla_df = pd.DataFrame(tabla_, columns = col_names)

tabla_df = tabla_df.drop(['P07724', 'P02770', 'P00689', 'P00687', 'Q5XLE4'], axis=1)
for num in range(len(tabla_df.axes[0])):
    row = tabla_df.iloc[num,:]
    max_value = max(row)
    row_new = row/max_value
    
    tabla_df.iloc[num,:] = row_new
tabla_df.to_excel('tabla.xlsx')

tabla_df = tabla_df.iloc[2:12,:]

pca = PCA(n_components=2)
pca.fit(tabla_df)
transformada = pca.transform(tabla_df)

mglearn.discrete_scatter(transformada[:,0],transformada[:,1])