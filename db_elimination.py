import sqlite_requests as sr
import protein_signals as ps
import os
import pandas as pd

def db_elimination(name):
    if 'Aquasearch_study' in name:
        pass
    else:
        name = 'Aquasearch_study - ' + name
    
    os.remove(name)

def protein_elimination(db, protein):
    id_ = sr.consulta_prot_id_from_code(db, protein)
    sr.delete_prot_id(db, id_)
    
    spectra = sr. consulta_spect(db, id_)
    spectra = pd.DataFrame(spectra, columns=['mz', 'intensity', 'standard', 'protein', 'sample', 'peptide'])
    peptides = spectra.loc[:, 'peptide'].unique()
    
    table_quantitative = sr.table_download(db, 'Quantitative_information')
    table_quantitative = pd.DataFrame(table_quantitative, columns= ['protein', 'sequence', 'mz', 'intensity'])
    peptides_del = []
    
    for pep in peptides:
        qi_p = table_quantitative[table_quantitative['sequence'] == str(pep)]
        if qi_p.shape[0] == 1:
            peptides_del.append(pep)
    
    sr.delete_spetrums(db, id_)
    for peptide in peptides_del:
        sr.delete_peptides(db, peptide)
    
    try:
        ps.quality_filter(db)
        ps.reference_spectrum(db)
        ps.peptide_reevaluation(db)
    except ValueError:
        pass

def sample_elimination(db, sample_):
    sr.delete_sample(db, sample_)
    
    ps.quality_filter(db)
    ps.reference_spectrum(db)
    ps.peptide_reevaluation(db)


# if __name__ == "__main__":
#   protein_elimination('Aquasearch_study - copia', 'P19121')
    # sr.create_db_by_user('prueba')
    # db_elimination('prueba')
    # sample_elimination('Aquasearch_study', 'mcE67_Vic')