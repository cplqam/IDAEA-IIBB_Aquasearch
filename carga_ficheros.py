#He cogido el XML 

import os
import xml.etree.ElementTree as et 
import pandas as pd
import numpy

#Pareamos los datos a partir del fichero XML
carlos = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Prueba/AquaSearch/1'
file = 'peaklist.xml'

path_ = os.path.join(carlos,file)

def parse_xml(path_):
    #Cojo las 2 columnas de mz e intensidad, aunque se podrian coger mas
    df_cols = ["mass","intensity"] 
    rows = []

    xtree = et.parse(path_)
    xroot = xtree.getroot() 

    for node in xroot: 
        s_mass = node.find("mass").text if node is not None else None
        s_int = node.find('absi').text if node is not None else None
    
        rows.append({"mass": s_mass,"intensity": s_int})

    out_df = pd.DataFrame(rows, columns = df_cols)
    out_df = out_df.astype(float)
    
    mz_int = numpy.array([pd.DataFrame.to_numpy(out_df.loc[:,'mass']),pd.DataFrame.to_numpy(out_df.loc[:,'intensity'])])
    mz_int = numpy.transpose(mz_int)
    return out_df,mz_int

df,mz_int = parse_xml(path_)

########TAMBIEN DEJO LA OPCION DE HACERLO A PARTIR DE EXCEL#######

# import os
# import pandas as pd
# import numpy

# carlos = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Prueba/AquaSearch'
# file = 'spectro_test_centroid_deisotoped.xlsx'
# path_ = os.path.join(carlos,file)

# def parse_excel(path_):
#     df = pd.read_excel(path_,header=2)

#     #Cojo la columna de mz, intensidad y area para comparar con la database 
#     df_2 = df[['m/z','Intens.','Area']]

#     mz_int = numpy.array([pd.DataFrame.to_numpy(df_2.loc[:,'m/z']),pd.DataFrame.to_numpy(df_2.loc[:,'Intens.'])])
#     mz_int = numpy.transpose(mz_int)
#     return df_2,mz_int

# df,mz_int = parse_excel(path_)

########TAMBIEN DEJO LA OPCION DE HACERLO A PARTIR DE EXCEL#######

#Aquí habría que hacer la consulta a Live SQL 
