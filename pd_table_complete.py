import pandas as pd 
import os
import urllib 
import sqlite_requests as sr

# This function provides a information fom uniprot about a protein
def uniprot_information(accession):
    """Provides the name and organism of origin of a peptide through a Uniprot request.

    INPUT:
        accession: string. the Uniprot peptide ? accession
    """

    data = urllib.request.urlopen("https://www.uniprot.org/uniprot/%s.txt" % accession).read()
    data = str(data)
    data = data.split('\\n')
    organism = ''
    
    for i in data:
        n = i.split('=')
        o = i.split('   ')
        
        if n[0] == 'DE   RecName: Full':
            name = n[1].strip(';')
            name = name.split('{')[0]
        elif o[0] == 'OS':
            o = str(o[1])
            organism = organism + o
            organism = organism.strip('.')
        # TODO: Should be "Eukaryote"
        if o[0] == 'OC' and o[1].split(';')[0] == 'Eukaryota':
            organism = organism.split('(')[0]
         
    return name, organism

# To help to complete the output form Proteome Discoverer
def table_complete(df, d):
    """Helps to complete the output of Proteome Discoverer with the name and the organism of origin of the proteins.

    INPUT:
        df: pandas dataframe. Result of Proteome Discoverer search.
        # TODO: to define. Single char variables are not recommended.
        d: ?
    """
        
    group = df['Protein Group Accessions']
    new_col = []
    for i in group:
        new_col.append(d[i])
    return new_col

# To complete the output form Proteome Discoverer
# TODO: Should be 'protein_dictionary'
def protein_information(df, db='Aquasearch_study', table='protein_diccionary'):
    """Completes the table result of Proteome Discoverer search with information about proteins (protein name
    and organism) with a Uniprot request.

    INPUT:
        d: dictionary. Result of 'uniprot_information' function
        db: string.  Name of the database we are saving the datasets. Default: "Aquasearch_study"
        table: string array. Name of the dable ?? where we are saving the dictionary. Default: "protein_dictionary"
    """
    
    # Try to create a database if it doesn't exist
    if not os.path.exists(db):
        sr.create_db(db)
    
    # Try to create a table if it doesn't exist
    sr.create_table_proteins_dic(db, table)

    group = df['Protein Group Accessions'].unique()
    name = {}
    organ = {}
    # c= 0
    for prot in group:
        # c = c+1
        # print(c)
        protein = prot.split(';')
        l_n = ''
        l_o = ''
            
        for p in protein:
            data = sr.table_request_prot_dict(db, table, p)
            if len(data) == 0:
                n, o = uniprot_information(p)
                df_u = pd.DataFrame({'Accession': [p], 'Protein name': [n], 'Organism': [o]})
                l_n = l_n + '|' + str(n)
                l_o = l_o + '|' + str(o)
                sr.insert_prot_code(db, table, df_u)
            else:
                n = data[0][0]
                o = data[0][1]
                l_n = l_n + '|' + str(n)
                l_o = l_o + '|' + str(o)
            
        name[prot] = l_n[1:]
        organ[prot] = l_o[1:]

    name_p = table_complete(df, name)
    name_o = table_complete(df, organ)
    df['Protein Name'] = name_p
    df['Organism Name'] = name_o

    return df


# To test the functions
if __name__ == '__main__':
    path_ = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Correo de Montse/Identificaciones'
    path_2 = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/IDAEA-IIBB'
    d = 'Protein_dic.xlsx'
    excel = 'mcE61_Figueres_01_peptides_2.xlsx'
        
    afile = os.path.join(path_, excel) 
    df = pd.read_excel(afile)
    
    final = protein_information(df)
