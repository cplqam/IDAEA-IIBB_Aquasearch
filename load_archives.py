import xml.etree.ElementTree as Et
import pandas as pd
import numpy
import sqlite_requests as sr

# To import the tables from database
def db_table_request(x):
    table = sr.table_download('Aquasearch_study', x)
    table = pd.DataFrame(table, columns=('mz', 'intensity', 'Unique'))
    
    if table.shape[0] > 1:
        mz_s = table['mz'].squeeze()
        int_s = table['intensity'].squeeze()
    else:
        mz_s = [table['mz'].squeeze()]
        int_s = [table['intensity'].squeeze()]
    
    mz_int = numpy.array([mz_s,int_s]).transpose()
    unique_inf = list(table.loc[:,'Unique'])

    return mz_int, unique_inf
    

# To import new requests
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

    mz_s = out_df['mass'].squeeze()
    int_s = out_df['intensity'].squeeze()

    mz_int_ = numpy.array([mz_s, int_s]).transpose()
    return out_df, mz_int_


# To test the function
if __name__ == '__main__':
    # Pair the data from an XML file

    df, mz_int = parse_xml('test_files/mcE61_Figueres.xml')
    table,unique_inf = db_table_request('P19121')

######ALSO LET THE OPTION TO DO IT FROM EXCEL#######

# import pandas as pd
# import numpy

# def parse_excel(path_):
#     df = pd.read_excel(path_,header=2)

#     # Select mz, intensity and area columns to compare with database
#     df_2 = df[['m/z','Intens.','Area']]

#     mz_int = numpy.array([pd.DataFrame.to_numpy(df_2.loc[:,'m/z']),
#                           pd.DataFrame.to_numpy(df_2.loc[:,'Intens.'])])
#     mz_int = numpy.transpose(mz_int)
#     return df_2,mz_int

# if __name__ == '__main__':
#   df,mz_int = parse_excel('test_files/mcE61_Figueres.xml')

# #######ALSO LET THE OPTION TO DO IT FROM EXCEL#######
