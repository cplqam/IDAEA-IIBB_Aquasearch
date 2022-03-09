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

        {1161.503094948756: ['Carcinoembryonic antigen-related cell adhesion;Superoxide dismutase [Cu-Zn]|Homo sapiens ;Gallus gallus |P06731;P80566|No unique;No unique'],
        1213.601431420398: ['Protein AMBP;Immunoglobulin heavy constant alpha;Immunoglobulin alpha-2 heavy chain|Homo sapiens ;Homo sapiens ;Homo sapiens |P02760;P01876;P0DOX2|No unique;No unique;No unique'],
       .................................................................................................................
        2599.235296717235: ['Albumin|Homo sapiens |P02768|Unique'],
        2990.446067568852: ['Alpha-amylase 1A;Alpha-amylase 2B|Homo sapiens (Human);Homo sapiens |P0DUB6;P19961|No unique;No unique']}

        maldi: an array variable, it contains the maz and intensity from MALDI spectrum

    """

    mz_maldi = maldi[:, 0]

    size = mz_maldi.shape[0]
    alist = ['-'] * size

    df = pd.DataFrame({'mz': maldi[:, 0], 'intensity': maldi[:, 1], 'Protein': alist,
                       'Organism Name': alist, 'Protein Accession code': alist, 'Unique Pep': alist
                       })

    for mz_k, value in dictionary.items():
        prot_name, org_name, accession, unique_pep = value[0].split('|')

        pos = numpy.where(mz_maldi == mz_k)[0]
        df.loc[pos, 'Protein'] = prot_name
        df.loc[pos, 'Organism Name'] = org_name
        df.loc[pos, 'Protein Accession code'] = accession
        df.loc[pos, 'Unique Pep'] = unique_pep
    return df


def xml_complete(xml_, ident_pep, ident_prot, n_=250, ppm=100, unique_=1):
    """Assigns the protein and organism name to the MALDI spectrum signals

       >>> result = xml_complete('test_files/mcE61_Figueres.xml', 'test_files/mcE61_PD14_Figueres_Peptides.xlsx', 'test_files/mcE61_PD14_Figueres_Proteins.xlsx', n_=100, unique_=1)
       >>> result.shape
       (69, 6)

       >>> result.iloc[3,:]
       mz                                                              1161.503095
       intensity                                                      16112.227734
       Protein                   Carcinoembryonic antigen-related cell adhesion...
       Organism Name                                  Homo sapiens ;Gallus gallus 
       Protein Accession code                                        P06731;P80566
       Unique Pep                                                              0;0
       Name: 3, dtype: object

        INPUT
        xml_: string. Path of the MALDI xml format file
        ident_pep: string. Path of the file with the identified peptides from PD
        ident_prot: string. Path of the file with the identified proteins from PD
        n_: integer. Number of proteins selected depending on their representation.
            default n_ = 10
        ppm: integer. Maximum error allowed for considering a signal from MALDI and Orbitrap  as the same.
             default ppm = 100
        unique_: integer.
                 unique_ = 1, all the possibilities are considered in non-unique peptides
                 unique_ = 0, the non-unique peptides are included.
                 default unique_ = 1
    """
    target_columns = ['Protein Name', 'Organism Name', 'Protein Group Accessions', 'Unique Pep']
    df, mz_int = load_archives.parse_xml(xml_)
    df_ident, df_names = complete_table_proteins(ident_prot, n_)

    dictionary = {}
    if unique_:
        identifications = pd_table_selection.organism_selection(ident_pep, sel=2)
        for ion in mz_int:
            mz_xml = ion[0]
            delta_mz = ppm * mz_xml / 1E6
            list_po = []
            # TODO: itertuples is faster but needs column names being valid python identifiers
            #       Should replace this with a vectorized method.
            for idx, row in identifications.iterrows():
                # Put in a list all Orbitrap signals, with their identification,
                # which are related with the MALDI signal.
                mz_ident = row['MH+ [Da]']
                if abs(mz_ident - mz_xml) <= delta_mz:
                    app = ';'.join([row[value] for value in target_columns])
                    app = app.split(';')

                    options = (len(app)-1) // 3       # Last column Unique/No Unique
                    if options == 1:
                        list_po.append(';'.join(app))
                    elif options >= 2:
                        for n in range(options):
                            candidate = app[n:-1:options]
                            candidate.append('0')
                            a = ';'.join(candidate)
                            list_po.append(a)

            list_po = list(pd.unique(list_po))
            if len(list_po) >= 2:
                p_protein_name = []
                p_organism_name = []
                p_accession = []
                p_unique = []

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

                if len(p_protein_name) >= 1:
                    p_protein_name = ';'.join(p_protein_name)
                    p_organism_name = ';'.join(p_organism_name)
                    p_accession = ';'.join(p_accession)
                    p_unique = ';'.join(p_unique)

                    list_po = [p_protein_name, p_organism_name, p_accession, p_unique]
                    list_po = ['|'.join(list_po)]
                    dictionary[mz_xml] = list_po

            elif len(list_po) == 1:

                candidate = list_po[0].split(';')
                acc = candidate[2]
                for num2 in range(len(df_names)):
                    answer = df_names.iloc[num2, 0]

                    if acc == answer:
                        list_po = [candidate[0], candidate[1], candidate[2], candidate[3]]
                        list_po = ['|'.join(list_po)]
                        dictionary[mz_xml] = list_po

    else:
        identifications = pd_table_selection.organism_selection(ident_pep, sel=1)
        for i in range(mz_int.shape[0]):
            mz_xml = mz_int[i, 0]
            list_po = []
            for j in range(identifications.shape[0]):
                mz_ident = identifications.loc[j, 'MH+ [Da]']
                if mz_xml > mz_ident:
                    ppm_calculated = (1 - (mz_ident / mz_xml)) * 1E6
                else:
                    ppm_calculated = (1 - (mz_xml / mz_ident)) * 1E6

                if ppm_calculated <= ppm:
                    app = ';'.join([identifications.loc[j, 'Protein Name'],
                                    identifications.loc[j, 'Organism Name'],
                                    identifications.loc[j, 'Protein Group Accessions'],
                                    identifications.loc[j, 'Unique Pep']])
                    list_po.append(app)

            if len(list_po) >= 2:
                positions = []
                for acc in list_po:
                    c = []
                    acc = acc.split(';')[2]

                    for num2 in range(len(df_names)):
                        answer = df_names.iloc[num2, 0]
                        if acc == answer:
                            c.append(list(df_names.iloc[:, 0]).index(acc))

                    if len(c) == 0:
                        c.append(n_+1)

                    positions.append(c[0])

                idx = positions.index(min(positions))
                if positions[idx] < n_ + 1:
                    prot_def = list_po[idx]
                    list_po = prot_def
                    list_po = [list_po]
                    dictionary[mz_xml] = list_po

            elif len(list_po) == 1:
                position = []
                acc = list_po[0]
                acc = acc.split(';')[2]

                for num2, answer in enumerate(df_names):
                    if acc == answer:
                        position.append(list(df_names.iloc[:, 0]).index(acc))

                if len(position) == 0:
                    position.append(n_ + 1)

                if position[0] < n_ + 1:
                    dictionary[mz_xml] = list_po

    result = maldi_ident_join(dictionary, mz_int)
    return result



if __name__ == '__main__':

    import doctest
    doctest.testmod()


    import common_module as cm
    cm.console_wide_view()
    result_1 = xml_complete('test_files/mcE61_Figueres.xml',
                            'test_files/mcE61_PD14_Figueres_Peptides.xlsx',
                            'test_files/mcE61_PD14_Figueres_Proteins.xlsx', n_=100, unique_=1)


# C:\Python38\python.exe C:/Python38/programas/aquasearch/pd_maldi_match.py
#              mz     intensity  ... Protein Accession code           Unique Pep
# 0    832.300294   9016.435667  ...                      -                    -
# 1    927.475515   5074.444521  ...                      -                    -
# 2   1071.546896  10192.549118  ...                      -                    -
# 3   1161.503095  16112.227734  ...          P06731;P80566  No unique;No unique
# 4   1181.560038  21171.262371  ...                      -                    -
# ..          ...           ...  ...                    ...                  ...
# 64  2634.326469  23260.389467  ...                      -                    -
# 65  2812.361993  10658.778184  ...                      -                    -
# 66  2990.446068  17896.442913  ...          P0DUB6;P19961  No unique;No unique
# 67  3003.349389   6511.268556  ...                      -                    -
# 68  3016.452827  30834.559000  ...                      -                    -
#
# [69 rows x 6 columns]
#
# Process finished with exit code 0
