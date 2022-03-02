from collections import Counter
import pd_table_complete as pdis

def organism_selection(path__, sel=1):
    """Completes the peptide output from Proteome Discoverer with the protein and
        organism they belong to, and it selects 1 option among the options of the
        non-unique peptides

        INPUT
        path__: string. The path of the peptide output from Proteome Discoverer
        sel = integer
            sel = 1: a selection based on the representation of each specie is performed for
                     non-unique peptides
            sel = 2: all the possibilities are considered for non-unique peptides
    """

    df = pdis.protein_information(path__)

    if sel == 2:
        uniq = []
        for code_ in df.loc[:, 'Protein Group Accessions']:
            if ';' in code_:
                uniq.append('No unique')
            else:
                uniq.append('Unique')
        df.loc[:, 'Unique Pep'] = uniq

    else:
        all_species = []
        species = list(df['Organism Name'])
        for item in species:
            if ';' in item:
                names = item.split(';')
                all_species.extend(names)
            else:
                all_species.append(item)

        species_count = Counter(all_species)
        most_common = species_count.most_common()
        freq_organism = dict(most_common)

        # TODO: Vectorize all this
        protein_selected, accession_selected, organism_selected, uniq = most_abundant_entry_selection(df, freq_organism)

        df.loc[:, 'Protein Group Accessions'] = accession_selected
        df.loc[:, 'Protein Name'] = protein_selected
        df.loc[:, 'Organism Name'] = organism_selected
        df.loc[:, 'Unique Pep'] = uniq

    df['Protein Name'] = df.apply(protein_name_simplification, axis=1)
    return df

def most_abundant_entry_selection(x, freq_organism):
    """Selects a peptide-inferred protein among all non-unique options based on organism prevalence.

        INPUT
        x: DataFrame returned by "pd_table_complete.protein_information.py" with
           the protein options for each identified peptide.
        freq_organism: list of organisms sorted by representation
        """

    final_n = []
    final_p = []
    final = []
    unique = []

    for n in range(len(x)):
        query = x.loc[n, 'Protein Name']
        prot = x.loc[n, 'Protein Group Accessions']
        nam = x.loc[n, 'Organism Name']

        if ';' not in query:
            final.append(query)
            final_p.append(prot)
            final_n.append(nam)
            unique.append("Unique")

        else:
            n = nam.split(';')
            species = query.split(';')
            species_a = prot.split(';')
            reps = []

            for name, name_a, name_n in zip(species, species_a, n):
                reps.append((freq_organism[name_n], name, name_a, name_n))
            reps.sort(reverse=True)

            higher = 0
            prob = []
            prob_a = []
            prob_n = []
            for rep, item, rep_a, rep_n in reps:
                if rep >= higher:
                    higher = rep       # assign the first entry (the higher in the list) to 'higher'
                    prob.append(item)  # store all entries with the same frequency
                    prob_a.append(rep_a)
                    prob_n.append(rep_n)
                else:
                    break
                prob_a = list(set(prob_a))

            final_p.append(prob_a[0])
            prob = list(set(prob))
            final.append(prob[0])
            prob_n = list(set(prob_n))
            final_n.append(prob_n[0])
            unique.append('No unique')

    return final, final_p, final_n, unique

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
    df_sel_1 = organism_selection('test_files/mcE61_PD14_Figueres_Peptides.xlsx', sel=1)
    df_sel_2 = organism_selection('test_files/mcE61_PD14_Figueres_Peptides.xlsx', sel=2)

    print('-----sel 1----------')
    print(df_sel_1)
    print()
    print('-----sel 2----------')
    print(df_sel_2)