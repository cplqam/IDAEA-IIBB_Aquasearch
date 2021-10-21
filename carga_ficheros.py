#He cogido el XML 

import os
import xml.etree.ElementTree as et 
import pandas as pd

#Pareamos los datos a partir del fichero XML
carlos = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Prueba/AquaSearch/1'
file = 'peaklist.xml'

path_ = os.path.join(carlos,file)

def parse_xml(path_):
    #Cojo las 3 columnas de mz, intensidad y carga, aunque se podrian coger mas
    df_cols = ["mass","intensity","area"]
    rows = []

    xtree = et.parse(path_)
    xroot = xtree.getroot() 

    for node in xroot: 
        s_mass = node.find("mass").text if node is not None else None
        s_int = node.find('absi').text if node is not None else None
        s_area = node.find("area").text if node is not None else None
    
        rows.append({"mass": s_mass,"intensity": s_int, "area": s_area})

    out_df = pd.DataFrame(rows, columns = df_cols)
    return(out_df)

df = parse_xml(path_)

########TAMBIEN DEJO LA OPCION DE HACERLO A PARTIR DE EXCEL#######

# import os
# import pandas as pd

# carlos = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Prueba/AquaSearch'
# file = 'spectro_test_centroid_deisotoped.xlsx'
# path_ = os.path.join(carlos,file)

# def parse_excel(path_):
#     df = pd.read_excel(path_,header=2)

#     #Cojo la columna de mz, intensidad y area para comparar con la database 
#     df_2 = df[['m/z','Intens.','Area']]
#     return(df_2)

# df = parse_excel(path_)

########TAMBIEN DEJO LA OPCION DE HACERLO A PARTIR DE EXCEL#######

#Aquí habría que hacer la consulta a Live SQL 
