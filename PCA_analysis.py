import pandas as pd
import numpy as np
from sklearn.decomposition import PCA


def PCA_2d(dic_, n_comp):
    
    sample_names = []
    for query, n in zip(dic_.keys(), range(len(dic_.keys()))):
        sample = dic_[query]
        name_sample = query.split('.')[0]
        sample_names.append(name_sample)
        
        if n == 0:
            col_names = list(sample.keys())
            
            lst = []
            for proteins in sample.keys():
                prot = sample[proteins]
                score = prot[0][2]
                lst.append(score)
            arr = np.array([lst])
            
        else:
            lst = []
            for proteins in sample.keys():
                prot = sample[proteins]
                score = prot[0][2]
                lst.append(score)
            arr = np.vstack([arr, lst])
                
        

    df = pd.DataFrame(arr, columns = col_names)
    
    for num in range(len(df.axes[0])):
        row = df.iloc[num,:]
        max_value = max(row)
        row_new = row/max_value
        
        df.iloc[num,:] = row_new
    
    pca = PCA(n_components = n_comp)
    pca.fit(df)
    variance = pca.explained_variance_ratio_
    transformada = pca.transform(df)
    
    return transformada, variance, sample_names
        
def PCA_3d(dic_, n_comp):
    
    
    sample_names = []
    for query, n in zip(dic_.keys(), range(len(dic_.keys()))):
        sample = dic_[query]
        name_sample = query.split('.')[0]
        sample_names.append(name_sample)
        
        if n == 0:
            col_names = list(sample.keys())
            
            lst = []
            for proteins in sample.keys():
                prot = sample[proteins]
                score = prot[0][2]
                lst.append(score)
            arr = np.array([lst])
            
        else:
            lst = []
            for proteins in sample.keys():
                prot = sample[proteins]
                score = prot[0][2]
                lst.append(score)
            arr = np.vstack([arr, lst])
                
        

    df = pd.DataFrame(arr, columns = col_names)
    
    for num in range(len(df.axes[0])):
        row = df.iloc[num,:]
        max_value = max(row)
        row_new = row/max_value
        
        df.iloc[num,:] = row_new
    
    pca = PCA(n_components = n_comp)
    pca.fit(df)
    variance = pca.explained_variance_ratio_
    transformada = pca.transform(df)
    
    return transformada, variance, sample_names

    
    
if __name__ == "__main__":
    a = PCA_2d(dic, 3, 1, 2, [8,6], 70, 'Test')
    # PCA_3d(dic, 3, 1, 2, 3, [8,6], 70, 'Test')
    