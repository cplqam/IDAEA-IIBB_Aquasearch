import sqlite3 as sql
import os

def create_db(name):
    """"Creates a new SQLite3 database.

        name: string. Name we want to assign to the DB
    """

    cone = sql.connect(name)
    cone.commit()
    cone.close()
       
def create_db_by_user(name):
    """"Creates a new SQLite3 database.

        name: string. Name we want to assign to the DB
    """
    name2 = 'Aquasearch_study - ' + name
    cone = sql.connect(name2)
    cone.commit()
    cone.close()
    
def db_request():
    """Perform a request to know the datasets
    """
    list_files = os.listdir()
    dbs = []
    for file in list_files:
        if 'Aquasearch_study' in file:
            if '-' in file:
                dbs.append(file.split('- ')[1])
            else:
                dbs.append(file)
    return dbs

# To create a protein dictionary
def create_table_proteins_dic(name, table_n):
    """"Creates a new table in the database.

        name: string. Name of the DB
        table_n: string. Name of the dictionary
    """

    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE IF NOT EXISTS " + table_n + "(Accession text, Protein name text, Organism text)")
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

# To create a table with the protein codes
def protein_codes_table(name):
    
    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("""CREATE TABLE IF NOT EXISTS Protein_codes(
                      id_ integer primary key, 
                      Protein text )""")
    conn.commit()
    conn.close()


def insert_sequence(name, seq, id_):
    """"Completes a table in a database.

        name: string. The name of the DB.
        table_n: string. The name of the table.
        df: list of tuples. Contains the mz and intensity values of the spectrum.
    """

    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("INSERT INTO Protein_codes VALUES (?,?)" , (id_, seq))
    conn.commit()
    conn.close()
    
def insert_protein_code(name, code, id_):
    """"Completes a table in a database.

        name: string. The name of the DB.
        table_n: string. The name of the table.
        df: list of tuples. Contains the mz and intensity values of the spectrum.
    """

    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("INSERT INTO Protein_codes VALUES (?,?)" , (id_, code))
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
    
def spectrums_table(name):
    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("""CREATE TABLE IF NOT EXISTS Spectrums_table (
                      mz real,
                      relative_intensity real, 
                      standard_signal text,
                      protein text,
                      sample text,
                      sequence integer,
                      FOREIGN KEY (protein) REFERENCES Protein_codes(id_),
                      FOREIGN KEY (sequence) REFERENCES Peptide_sequences(id_sequence))""") #########################
    conn.commit()
    conn.close()
    
def df_2_list_of_tuples(df):
    """Transforms the dataframe format (output of parse_zml function from load_archives.py)
        to lost of tuples for inserting the information into the table.
    """

    out_list = []
    for i in range(df.shape[0]):
        mz = float(df.iloc[i, 0])
        stand = df.iloc[i, 2]
        rel_int = float(df.iloc[i, 1])
        seq = float(df.iloc[i, 5])
        samp = df.iloc[i, 4]
        protein = df.iloc[i, 3]
        tup = (mz, rel_int, stand, protein, samp, seq)
        out_list.append(tup)
    return out_list


def insert_new_spectrum(name, signals_i):
    conn = sql.connect(name)
    df_def = df_2_list_of_tuples(signals_i)
    cursor = conn.cursor()
    instruction = "INSERT INTO Spectrums_table VALUES (?,?,?,?,?,?)" 
    cursor.executemany(instruction, df_def)
    conn.commit()
    conn.close()
    
def table_download(name, table_n):
    """"Perform a request to get information from the protein dictionary.

        name: string. The name of the DB.
        table_n: string. The name of the dictionary.
        code: string. Uniprot code of the protein.
    """

    conn = sql.connect(name)
    cur = conn.cursor()
    cur.execute("SELECT * FROM " + table_n)
    datos = cur.fetchall()
    conn.close()
    datos = list(datos)
    return datos 

def eliminate_table(name, table_n):
    """"Eliminates a new database table.

        name: string. The name of the DB.
        table_n: string. The name of the table.
    """
    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE " + table_n)
    conn.commit()
    conn.close()

def table_quant_inf(name):
    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("""CREATE TABLE IF NOT EXISTS Quantitative_information (
                      protein text,
                      sequence text,
                      mz_aver real,
                      relative_int_aver real,
                      FOREIGN KEY (protein) REFERENCES Protein_codes(id_),
                      FOREIGN KEY (sequence) REFERENCES Peptide_sequences(peptide_sequence))""") ##################################
    conn.commit()
    conn.close()

def new_quiant_inf(name, df):
    conn = sql.connect(name)
    df_def = df_2_list_of_tuples_quanti(df)
    cursor = conn.cursor()
    instruction = "INSERT INTO Quantitative_information VALUES (?,?,?,?)" 
    cursor.executemany(instruction, df_def)
    conn.commit()
    conn.close()
    
    
def df_2_list_of_tuples_quanti(df):
    out_list = []
    for i in range(df.shape[0]):
        prot = df.iloc[i, 0]
        pep = df.iloc[i, 1]
        mz = float(df.iloc[i, 2])
        mz = round(mz, 4)
        rel_int = float(df.iloc[i, 3])
        rel_int = round(rel_int, 2)
        tup = (prot, pep, mz, rel_int)
        out_list.append(tup)
    return out_list
    
def table_sequences(name):
    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("""CREATE TABLE IF NOT EXISTS Peptide_sequences (
                      id_sequence integer primary key,
                      peptide_sequence text,
                      uniq text,
                      molecular_mass text)""")

def new_sequence(name, seq, _id, uni, mass):
    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("INSERT INTO Peptide_sequences VALUES (?,?,?,?)" , (_id, seq, uni, mass))
    conn.commit()
    conn.close()
    
def spectrums_table_filtered(name):
    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("""CREATE TABLE IF NOT EXISTS Spectrums_table_filtered (
                      mz real,
                      relative_intensity real, 
                      standard_signal text,
                      protein text,
                      sample text,
                      sequence integer,
                      FOREIGN KEY (protein) REFERENCES Protein_codes(id_),
                      FOREIGN KEY (sequence) REFERENCES Peptide_sequences(id_sequence))""") ###########################
    conn.commit()
    conn.close()

def insert_new_spectrum_filtered(name, signals_i):
    conn = sql.connect(name)
    df_def = df_2_list_of_tuples(signals_i)
    cursor = conn.cursor()
    instruction = "INSERT INTO Spectrums_table_filtered VALUES (?,?,?,?,?,?)" 
    cursor.executemany(instruction, df_def)
    conn.commit()
    conn.close()

def consulta_uniq(name, pep):
    conn = sql.connect(name)
    cur = conn.cursor()
    cur.execute("SELECT uniq FROM Peptide_sequences WHERE id_sequence =" + str(pep) + "")
    datos = cur.fetchall()[0][0]  
    conn.close()
    return datos 

def consulta_pep_seq(name, pep):
    conn = sql.connect(name)
    cur = conn.cursor()
    cur.execute("SELECT peptide_sequence FROM Peptide_sequences WHERE id_sequence =" + str(pep) + "")
    datos = cur.fetchall()[0][0]  
    conn.close()
    return datos

def consulta_prot_id(name, prot):
    conn = sql.connect(name)
    cur = conn.cursor()
    cur.execute("SELECT Protein FROM Protein_codes WHERE id_ =" + str(prot) + "")
    datos = cur.fetchall()[0][0]  
    conn.close()
    return datos

def consulta_prot_id_from_code(name, prot):
    conn = sql.connect(name)
    cur = conn.cursor()
    cur.execute("SELECT id_ FROM Protein_codes WHERE Protein = '" + str(prot) + "'") 
    datos = cur.fetchall()[0][0]
    conn.close()
    return datos

def delete_prot_id(name, prot):
    conn = sql.connect(name)
    cur = conn.cursor()
    cur.execute("DELETE FROM Protein_codes WHERE id_ = '" + str(prot) + "'") 
    conn.commit()
    conn.close()

def delete_spetrums(name, prot):
    conn = sql.connect(name)
    cur = conn.cursor()
    cur.execute("DELETE FROM Spectrums_table WHERE protein = " + str(prot) + "") 
    conn.commit()
    conn.close()
    
def delete_peptides(name, pep):
    conn = sql.connect(name)
    cur = conn.cursor()
    cur.execute("DELETE FROM Peptide_sequences WHERE id_sequence = " + str(pep) + "") 
    conn.commit()
    conn.close()
    
def delete_sample(name, sample_name):
    conn = sql.connect(name)
    cur = conn.cursor()
    cur.execute("DELETE FROM Spectrums_table WHERE sample = '" + str(sample_name) + "'") 
    conn.commit()
    conn.close()
    
def consulta_spect(name, prot):
    conn = sql.connect(name)
    cur = conn.cursor()
    cur.execute("SELECT * FROM Spectrums_table WHERE protein =" + str(prot) + "")
    datos = cur.fetchall() 
    conn.close()
    return datos
    
if __name__ == "__main__":
    # delete_prot_id('Aquasearch_study', 1)
    a = consulta_spect('Aquasearch_study', 3)