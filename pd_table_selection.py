import os
import pd_table_complete as pdis
from collections import Counter

def organism_selection(path__):
    """This function completes the peptide output from Proteime Discoverer with the protein and organism they belong
        to, and it selects 1 option among the options of the non-unique peptides
        
        INPUT
        path__: string. The path of the peptide output from Proteome Discoverer"""
        
    df = pdis.protein_information(path__)
    
    all_species = []
    species = list(df['Organism Name'])
    for item in species:
        if '|' in item:
            names = item.split('|')
            all_species.extend(names)
        else:
            all_species.append(item)
            
    species_count = Counter(all_species)
    most_common = species_count.most_common()
    frec_organism = dict(most_common)
    
    protein_selected, accession_selected, organism_selected, uniq = most_abundant_entry_selection(df, frec_organism)
    
    df.loc[:, 'Protein Group Accessions'] = accession_selected
    df.loc[:, 'Protein Name'] = protein_selected
    df.loc[:, 'Organism Name'] = organism_selected
    df.loc[:, 'Unique Pep'] = uniq
    
    df['Protein Name'] = df.apply(protein_name_simplification, axis=1)
            
    return df
    
def most_abundant_entry_selection(x, frec_organism):
    """This function selects the peptide chosen among all non-unique options depending on the representation of the organism 
        the peptide belong to
        INPUT
        x: the data frame result from "pd_table_complete.protein_information.py" with the protein options for the peptides
        frec_organism: the sorted list of the organisms depending on their representation
        """
        
    final_n = []
    final_p = []
    final = []
    uni = []

    for n in range(len(x)):
        query = x.loc[n,'Protein Name']
        prot = x.loc[n,'Protein Group Accessions']
        nam = x.loc[n,'Organism Name']
        if '|' not in query:
            final.append(query)
            final_p.append(prot)
            final_n.append(nam)
            uni.append("Unique")
        else:
            n = nam.split('|')
            species = query.split('|')
            species_a = prot.split(';')
            reps = []
            for name,name_a,name_n in zip(species,species_a,n):
                reps.append((frec_organism[name_n], name,name_a,name_n))
            reps.sort(reverse=True)
        
            higher = 0
            prob = []
            prob_a = []
            prob_n = []
            for rep, item, rep_a, rep_n in reps:
                if rep >= higher:
                    higher = rep       # assign the first entry (the higher in the list) to 'higher'.
                    prob.append(item)  # store all entries with the same frecuency
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
            uni.append('No unique')
    
    return final, final_p, final_n, uni

def protein_name_simplification(x):
    """This function perform a protein name simplification
        
        INPUT
        x: a dataframe variable, completed and non-unique peptide selected dataframe output from Proteime Discoverer"""
        
    name = x['Protein Name']
    atoms = name.split()[:4]
    return ' '.join(atoms)
        
    
    
if __name__ == '__main__':
    path_ = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Correo de Montse/Identificaciones'

    excel = 'mcE61_Figueres_01_peptides_2.xlsx'
    afile = os.path.join(path_, excel) 
    
    df_sel = organism_selection(afile)