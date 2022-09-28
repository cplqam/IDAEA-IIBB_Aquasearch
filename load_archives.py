import xml.etree.ElementTree as Et
import pandas as pd
import numpy
import sqlite_requests as sr

# To import the tables from database
def db_table_request():
    """This function import a table from a SQLite dataset

        >>> res = db_table_request()
        >>> len(res) == list(res.keys())[-1]
        True

    """
    table = sr.table_download('Aquasearch_study', 'Quantitative_information')
    table = pd.DataFrame(table, columns=('protein', 'sequence', 'mz', 'intensity'))
    dic = {}

    n_prot = map(int,table['protein'])
    max_prot = max(n_prot)
    for num in range(int(max_prot)):
        num = num + 1
        dic[num] = table[table['protein'] == str(num)]

    return dic


# To import new requests
def parse_xml(path_):
    """This function import a table from xml file

        >>> df_, mzint = parse_xml('test_files/mcE61_Figueres.xml')
        >>> df_.shape == mzint.shape
        True

    """
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

def parse_txt(path_):
    """This function import a table from txt file

        >>> df_2, mzint_2 = parse_txt('test_files/Standares/pmf_H1.txt')
        >>> df_2.shape == mzint_2.shape
        True

    """

    df_ = pd.read_csv(path_, sep = '\s+')
    df_ = df_.iloc[:,[0,1,2]]

    df_.columns =  ['mass', 'intensity', 'unique']

    mz_s = df_['mass'].squeeze()
    int_s = df_['intensity'].squeeze()
    uniq_s = df_['unique'].squeeze()

    mz_int = numpy.array([mz_s, int_s, uniq_s]).transpose()

    return(df_, mz_int)

def parse_txt_request(path_):
    """This function import a table from txt file

        >>> df_2, mzint_2 = parse_txt('test_files/Standares/pmf_H1.txt')
        >>> df_2.shape == mzint_2.shape
        True

    """

    df_ = pd.read_csv(path_, sep = '\s+')
    df_ = df_.iloc[:,[0,1]]

    df_.columns =  ['mass', 'intensity']

    mz_s = df_['mass'].squeeze()
    int_s = df_['intensity'].squeeze()

    mz_int = numpy.array([mz_s, int_s]).transpose()

    return(df_, mz_int)

# To test the function
if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # Pair the data from an XML file

    prueba, prueba_2 = parse_txt('test_files/Standares/pmf_H1_completed.txt')
    df, mz_int = parse_xml('test_files/mcE61_Figueres.xml')
    try:
        results = db_table_request()
    except ValueError:
        print('no existen las tablas')
