import load_archives
import os
import sqlite3 as sql

def create_db(name):
    """"Creates a new SQLite3 database.

        name: string. Name we want to assign to the DB
    """

    cone = sql.connect(name)
    cone.commit()
    cone.close()
    
# To create a protein dictionary
def create_table_proteins_dic(name, table_n):
    """"Creates a new table in the database.

        name: string. Name of the DB
        table_n: string. Name of the dictionary
    """

    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE IF NOT EXISTS" + table_n + "(Accession text, Protein name text, Organism text)")
    conn.commit()
    conn.close()

def df_2_list_of_tuples_proteins_dict(df):
    """Provides the transformation of the dataframe format (output of parse_zml function from
        load_archives.py) to list of tuples for inserting the information into the table.
    """

    out_list = []
    for i in range(df.shape[0]):
        access = df.iloc[i, 0]
        prot_n = df.iloc[i, 1]
        org = df.iloc[i, 2]
        tup = (access, prot_n, org)
        out_list.append(tup)
    return out_list

def insert_prot_code(name, table_n, df):
    """"Completes a table in a database.

        name: string. The name of the DB.
        table_n: string. The name of the dictionary.
        df: list of tuples. Contains the accession code, protein name and organism.
    """

    conn = sql.connect(name)
    df_def = df_2_list_of_tuples_proteins_dict(df)
    cursor = conn.cursor()
    instruction = f"INSERT INTO " + table_n + " VALUES (?,?,?)" 
    cursor.executemany(instruction, df_def)
    conn.commit()
    conn.close()
    
def table_request_prot_dict(name, table_n, code):
    """"Performs a request to get information from the protein dictionary.

        name: string. The ame of the DB.
        table_n, string. The name of the dictionary.
        code: string. The uniprot code of the protein.
    """

    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("SELECT Protein name, Organism FROM " + table_n + "  WHERE Accession = '" + code + "'")
    datos = cursor.fetchall()
    conn.close()
    return datos


# To introduce a spectrum of reference
def create_table_organism(name, organism, protein):
    """"Creates a new database table.

        name: string. The name of the DB.
        organism: string. The name of the organism the data came from.
        protein: string. The name of the identified protein.
    """

    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE " + organism + "_" + protein + "(mz integer, intensity integer)")
    conn.commit()
    conn.close()
    
def df_2_list_of_tuples(df):
    """Transforms the dataframe format (output of parse_zml function from load_archives.py)
       to lost of tuples for inserting the information into the table.
    """

    out_list = []
    for i in range(df.shape[0]):
        mz = int(df.iloc[i, 0])
        intensity = int(df.iloc[i, 1])
        tup = (mz, intensity)
        out_list.append(tup)
    return out_list


def insert_spectrum(name, organism, protein, df):
    """"Completes a table in a database.

        name: string. The name of the DB.
        organism: string. The name of the organism the data came from.
        protein: string. The name of the identified protein.
        df: list of tuples. Contains the mz and intensity values of the spectrum.
    """

    conn = sql.connect(name)
    df_def = df_2_list_of_tuples(df)
    cursor = conn.cursor()
    instruction = f"INSERT INTO " + organism + "_" + protein + " VALUES (?,?)" 
    cursor.executemany(instruction, df_def)
    conn.commit()
    conn.close()
    
    
def table_request_test(name, table_n, code):
    """"Perform a request to get information from the protein dictionary.

        name: string. The name of the DB.
        table_n: string. The name of the dictionary.
        code: string. Uniprot code of the protein.
    """

    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("SELECT intensity FROM " + table_n + "  WHERE mz = " + code + "")
    datos = cursor.fetchall()
    conn.close()
    return datos

# To test the functions
# if __name__ == '__main__':
#     a = input('Is this a new database? (yes/no): ')
#     if a == 'yes': 
#         create_db('reference_spectra')
        
#     b = input('Is this a new table? (yes/no): ')
#     if b == 'yes': 
#         create_table('reference_spectra', 'pig', 'albumin')
    
#         carlos = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Prueba/AquaSearch/1'
#         afile = 'peaklist.xml'

#         path_ = os.path.join(carlos, afile)
    
#         df, mz_int = carga_ficheros.parse_xml(path_)
    
#     insert_spectrum('reference_spectra', 'pig', 'albumin', df)

    