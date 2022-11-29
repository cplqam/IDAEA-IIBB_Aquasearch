from collections import Counter
import pd_table_complete as pdis

def organism_selection(path__, sel=1, db = 'Aquasearch_study'):
    """Completes the peptide output from Proteome Discoverer with the protein and
        organism they belong to, and it selects 1 option among the options of the
        non-unique peptides

        >>> df_sel_2 = organism_selection('test_files/mcE61_PD14_Figueres_Peptides.xlsx', sel = 2)
        >>> df_sel_2.shape
        (1152, 20)

        >>> df_sel_2.loc[11, 'Protein Name']
        'Alpha-amylase 1A;Alpha-amylase 2B;Pancreatic alpha-amylase'

        >>> df_sel_2.loc[11, 'Unique Pep']
        '0'

        INPUT
        path__: string. The path of the peptide output from Proteome Discoverer
        sel = integer
            sel = 1: a selection based on the representation of each specie is performed for
                     non-unique peptides
            sel = 2: all the possibilities are considered for non-unique peptides
    """

    df = pdis.protein_information(path__, db)
    
    if sel == 2:
        uniques = []
        for accessions in df['Protein Group Accessions']:
            if ';' in accessions:
                uniques.append('0')
            else:
                uniques.append('1')
        df.loc[:, 'Unique Pep'] = uniques

    else:
        all_species = []
        for species in df['Organism Name']:
            names = species.split(';')
            all_species.extend(names)

        species_count = Counter(all_species)
        most_common = species_count.most_common()
        freq_organism = dict(most_common)

        protein_selected, accession_selected, organism_selected, uniques = most_abundant_entry_selection(df, freq_organism)

        df.loc[:, 'Protein Group Accessions'] = accession_selected
        df.loc[:, 'Protein Name'] = protein_selected
        df.loc[:, 'Organism Name'] = organism_selected
        df.loc[:, 'Unique Pep'] = uniques

    df['Protein Name'] = df.apply(protein_name_simplification, axis=1)
    df['Standard signal'] = '0'
    return df

def most_abundant_entry_selection(x, freq_organism):
    """Selects a peptide-inferred protein among all non-unique options based on organism prevalence.

        INPUT
        x: DataFrame returned by "pd_table_complete.protein_information.py" with
           the protein options for each identified peptide.
        freq_organism: list of organisms sorted by representation
        """

    final_organism = []
    final_accession = []
    final_protein_name = []
    unique = []

    for n in range(len(x)):
        query = x.loc[n, 'Protein Name']
        accession = x.loc[n, 'Protein Group Accessions']
        nam = x.loc[n, 'Organism Name']

        if ';' not in query:
            final_protein_name.append(query)
            final_accession.append(accession)
            final_organism.append(nam)
            unique.append('1')

        else:
            query_names = nam.split(';')
            species = query.split(';')
            acc_code = accession.split(';')
            reps = []

            for name, acc, name_n in zip(species, acc_code, query_names):
                reps.append((freq_organism[name_n], name, acc, name_n))
            reps.sort(reverse=True)

            higher = 0
            prot_name_chosen = []
            accession_chosen = []
            organism_chosen = []
            for rep, rep_prot, rep_acc, rep_org in reps:
                if rep >= higher:
                    higher = rep       # assign the first entry (the higher in the list) to 'higher'
                    prot_name_chosen.append(rep_prot)  # store all entries with the same frequency
                    accession_chosen.append(rep_acc)
                    organism_chosen.append(rep_org)
                else:
                    break
                
            accession_chosen = list(set(accession_chosen))
            final_accession.append(accession_chosen[0])
            prot_name_chosen = list(set(prot_name_chosen))
            final_protein_name.append(prot_name_chosen[0])
            organism_chosen = list(set(organism_chosen))
            final_organism.append(organism_chosen[0])
            unique.append('0')

    return final_protein_name, final_accession, final_organism, unique

def protein_name_simplification(x):
    """Simplifies protein name (takes only 4 first words max).

        INPUT
        x: a dataframe variable, completed and non-unique peptide selected dataframe output
           from Proteome Discoverer

    """

    proteins = x['Protein Name'].split(';')
    names = []
    for protein in proteins:
        words = protein.split()[:4]
        names.append(' '.join(words))
    return ';'.join(names)


if __name__ == '__main__':

    import doctest
    doctest.testmod()

    df_sel_1 = organism_selection('test_files/mcE61_PD14_Control_Peptides.xlsx', sel=1, db = 'Aquasearch_study')
    df_sel_2 = organism_selection('test_files/mcE61_PD14_Control_Peptides.xlsx', sel=2, db = 'Aquasearch_study')

    print('-----sel 1----------')
    print(df_sel_1)
    print()
    print('-----sel 2----------')
    print(df_sel_2)
