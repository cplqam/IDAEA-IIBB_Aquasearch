# Choose XML

import os
import xml.etree.ElementTree as Et
import pandas as pd
import numpy

# Pair the data from an XML file
carlos = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Prueba/AquaSearch/1'
afile = 'peaklist.xml'

path_ = os.path.join(carlos, afile)

def parse_xml(path_):
    # Select mz and intensity columns, could take more if needed
    df_cols = ["mass", "intensity"]
    rows = []

    xtree = Et.parse(path_)
    xroot = xtree.getroot() 

    for node in xroot: 
        s_mass = node.find("mass").text if node is not None else None
        s_int = node.find('absi').text if node is not None else None
    
        rows.append({"mass": s_mass, "intensity": s_int})

    out_df = pd.DataFrame(rows, columns=df_cols)
    out_df = out_df.astype(float)
    
    mz_int = numpy.array([pd.DataFrame.to_numpy(out_df.loc[:, 'mass']),
                          pd.DataFrame.to_numpy(out_df.loc[:, 'intensity'])
                          ])
    mz_int = numpy.transpose(mz_int)
    return out_df, mz_int

if __name__ == '__main__':
    df, mz_int = parse_xml(path_)

# #######ALSO LET THE OPTION TO DO IT FROM EXCEL#######

# import os
# import pandas as pd
# import numpy

# carlos = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Prueba/AquaSearch'
# file = 'spectro_test_centroid_deisotoped.xlsx'
# path_ = os.path.join(carlos,file)

# def parse_excel(path_):
#     df = pd.read_excel(path_,header=2)

#     # Select mz, intensity and area columns to compare with database
#     df_2 = df[['m/z','Intens.','Area']]

#     mz_int = numpy.array([pd.DataFrame.to_numpy(df_2.loc[:,'m/z']),pd.DataFrame.to_numpy(df_2.loc[:,'Intens.'])])
#     mz_int = numpy.transpose(mz_int)
#     return df_2,mz_int

#if __name__ == '__main__':
#   df,mz_int = parse_excel(path_)

# #######ALSO LET THE OPTION TO DO IT FROM EXCEL#######

# Here, we should consult Live SQL
