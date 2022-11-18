import pandas as pd
import pd_table_selection
import load_archives
import numpy
pd.options.mode.chained_assignment = None

def split_description(x):
    """Split uniprot protein description in individual items

    description: "Albumin OS=Gallus gallus OX=9031 GN=ALB PE=1 SV=2 - [ALBU_CHICK]"
    returns items: "Albumin", "Gallus gallus"

    """
    description = x['Description']
    items = description.split(' OX')[0]
    protein, species = items.split(" OS=")
    return protein, species


def complete_table_proteins(proteins_file, n=10):
    """Selects and filter the most represented proteins

        INPUT

        proteins_file: string. Path of the file with the identified proteins
        n = integer, Number of proteins selected depending on their representation

    """

    df_full = pd.read_excel(proteins_file)
    df_full = df_full.head(n)
    
    df_full[['Protein Name', 'Organism Name']] = df_full.apply(split_description, axis=1, result_type='expand')   
    df_names = df_full[['Accession', 'Organism Name', 'Protein Name']]
    df_names = df_names.rename(columns={'Accession':'Protein Accession code'})

    return df_full, df_names


def maldi_ident_join(dictionary, maldi):
    """This function joins the result of maldi and the signal identifications

        INPUT
        dictionary: a dictionary with
            keys (int):
                 - mz value
            values (list):
                 - the protein names
                 - species
                 - accessions
                 - Unique/No unique

        maldi: an array variable, it contains the maz and intensity from MALDI spectrum

    """

    mz_maldi = maldi[:, 0]

    size = mz_maldi.shape[0]
    alist = ['-'] * size

    df = pd.DataFrame({'mz': maldi[:, 0], 'intensity': maldi[:, 1], 'Protein': alist,
                       'Organism Name': alist, 'Protein Accession code': alist, 'Unique Pep': alist,
                       'Standard signal':alist, 'Peptide seq': alist, 'Peptide mass': alist
                       })

    for mz_k, value in dictionary.items():
        prot_name, org_name, accession, unique_pep, stand_signal, pep_seq, pep_mm = value[0].split('|')

        pos = numpy.where(mz_maldi == mz_k)[0]
        df.loc[pos, 'Protein'] = prot_name
        df.loc[pos, 'Organism Name'] = org_name
        df.loc[pos, 'Protein Accession code'] = accession
        df.loc[pos, 'Unique Pep'] = unique_pep
        df.loc[pos, 'Standard signal'] = stand_signal
        df.loc[pos, 'Peptide seq'] = pep_seq
        df.loc[pos, 'Peptide mass'] = pep_mm
        
    return df


def txt_complete(txt_, ident_pep, ident_prot, n_=250, ppm=100, unique_=1):
    """Assigns the protein and organism name to the MALDI spectrum signals

       >>> result = txt_complete('test_files/mcE61_Figueres.txt', 'test_files/mcE61_PD14_Figueres_Peptides.xlsx', 'test_files/mcE61_PD14_Figueres_Proteins.xlsx', n_=100, unique_=1)
       >>> result.shape
       (145, 9)

       >>> result.iloc[11,:]
       mz                                    1105.524
       intensity                            3331.8413
       Protein                   Maltase-glucoamylase
       Organism Name                    Homo sapiens 
       Protein Accession code                  O43451
       Unique Pep                                   1
       Standard signal                              0
       Peptide seq                         FAGFPALINR
       Peptide mass                         1106.6198
       Name: 11, dtype: object

        INPUT
        txt_: string. Path of the MALDI txt format file
        ident_pep: string. Path of the file with the identified peptides from PD
        ident_prot: string. Path of the file with the identified proteins from PD
        n_: integer. Number of proteins selected depending on their representation.
            default n_ = 10
        ppm: integer. Maximum error allowed for considering a signal from MALDI and Orbitrap  as the same.
             default ppm = 100
        unique_: integer.
                 unique_ = 1, all the possibilities are considered in non-unique peptides
                 unique_ = 0, 1 peptide is chosen among all non-unique peptides.
                 default unique_ = 1
    """
    target_columns = ['Protein Name', 'Organism Name', 'Protein Group Accessions', 'Unique Pep', 'Standard signal', 'Sequence']
    df, mz_int = load_archives.parse_txt(txt_)
    df_ident, df_names = complete_table_proteins(ident_prot, n_)

    dictionary = {}
    if unique_:
        identifications = pd_table_selection.organism_selection(ident_pep, sel=2)
    else:
        identifications = pd_table_selection.organism_selection(ident_pep, sel=1)

    for ion in mz_int:
        mz_txt = ion[0]
        delta_mz = ppm * mz_txt / 1E6
        list_po = []
        for idx, row in identifications.iterrows():
            # Put in a list all Orbitrap signals, with their identification,
            # which are related with the MALDI signal.
            mz_ident = row['MH+ [Da]']
            if abs(mz_ident - mz_txt) <= delta_mz:
                app = ';'.join([row[value] for value in target_columns])
                app = app.split(';')
                app.append(str(round(mz_ident + 1.00784, 4)))

                options = (len(app)-4) // 3       # Last columns Unique,Satndard,Peptide and mm
                if options == 1:
                    list_po.append(';'.join(app))
                elif options >= 2:
                    for n in range(options):       # Select each possible candidate
                        if n == 0:
                            fin = options * 3 
                        else:
                            fin = options * 3 + n
                                
                        candidate = app[n:fin:options]
                        candidate.append(app[-4])
                        candidate.append(app[-3])
                        candidate.append(app[-2])
                        candidate.append(app[-1])
                        a = ';'.join(candidate)
                        list_po.append(a)
        list_po = list(pd.unique(list_po))
        if len(list_po) >= 2:
            p_protein_name = []
            p_organism_name = []
            p_accession = []
            p_unique = []
            p_stand = []
            p_seq = []
            p_mm = []

            for po in list_po:
                candidate = po.split(';')
                acc = candidate[2]
                for num2 in range(len(df_names)):
                    answer = df_names.iloc[num2, 0]
                    if acc == answer:
                        p_protein_name.append(candidate[0])
                        p_organism_name.append(candidate[1])
                        p_accession.append(candidate[2])
                        p_unique.append('0')
                        p_stand.append(candidate[4])
                        p_seq.append(candidate[5])
                        p_mm.append(candidate[6])
                        

            if len(p_protein_name) >= 1:
                p_protein_name = ';'.join(p_protein_name)
                p_organism_name = ';'.join(p_organism_name)
                p_accession = ';'.join(p_accession)
                p_unique = ';'.join(p_unique)
                p_stand = ';'.join(p_stand)
                p_seq = ';'.join(p_seq)
                p_mm = ';'.join(p_mm)

                list_po = [p_protein_name, p_organism_name, p_accession, p_unique, p_stand, p_seq, p_mm]
                list_po = ['|'.join(list_po)]
                dictionary[mz_txt] = list_po
                    
        elif len(list_po) == 1:

            candidate = list_po[0].split(';')
            acc = candidate[2]
            for num2 in range(len(df_names)):
                answer = df_names.iloc[num2, 0]
                
                if acc == answer:
                    list_po = [candidate[0], candidate[1], candidate[2], candidate[3], candidate[4], candidate[5], candidate[6]]
                    list_po = ['|'.join(list_po)]
                    dictionary[mz_txt] = list_po
    result = maldi_ident_join(dictionary, mz_int)
    return result



if __name__ == '__main__':

    import doctest
    doctest.testmod()


    import common_module as cm
    cm.console_wide_view()
    result_1 = txt_complete('test_files/mcE61_Figueres.txt',
                            'test_files/mcE61_PD14_Figueres_Peptides.xlsx',
                            'test_files/mcE61_PD14_Figueres_Proteins.xlsx', n_=100, ppm = 100, unique_=1)

