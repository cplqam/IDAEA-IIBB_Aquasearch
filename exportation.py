def single_exportation(file, name_sample, name_file):
    dict_keys = list(file.keys())
    
    with open(name_file, 'w') as f:
        f.write('------------------------------------------ \n')
        f.write('Aquasearch results '+name_sample+'\n')
        f.write('------------------------------------------ \n')
        f.write("Protein code\tProtein name\tOrganism\tScore\tNº of peptides\tUnique Peptides \n")
        for key in dict_keys:
            parameters = file[key]
            f.write(key+'\t'+str(parameters[0][0])+'\t'+str(parameters[0][1])+'\t'+
                    str(parameters[0][2])+'\t'+str(parameters[0][3])+'\t'+
                    str(parameters[0][4])+'\n')
        f.write('\n\n\n')
        f.write('INDIVIDUAL PROTEIN IDENTIFICATIONS \n\n')
        
        for key in dict_keys:
            parameters = file[key][1]
            f.write('Protein: '+key+'\n')
            f.write("Peptide sequence\tUnique peptide\tError(ppm) \n")
            try:
                row = parameters.shape[0]
                for n in range(row):
                    f.write(parameters.iloc[n,0]+'\t'+str(parameters.iloc[n,1])+'\t'+
                            str(parameters.iloc[n,2])+'\n')
                    
            except:
                f.write('-'+'\t'+'-'+'\t'+'-'+'\n')
            f.write('\n\n')

    

def multiple_exportation(file, name_file):
    name_samples = list(file.keys())
    
    with open(name_file, 'w') as f:
        for n_sam in name_samples:
            sample = file[n_sam]
            dict_keys = list(sample.keys())
            f.write('------------------------------------------ \n')
            f.write('Aquasearch results of '+n_sam+'\n')
            f.write('------------------------------------------ \n')
            f.write("Protein code\tProtein name\tOrganism\tScore\tNº of peptides\tUnique Peptides \n")
            
            for key in dict_keys:
                parameters = sample[key]
                f.write(key+'\t'+str(parameters[0][0])+'\t'+str(parameters[0][1])+'\t'+
                        str(parameters[0][2])+'\t'+str(parameters[0][3])+'\t'+
                        str(parameters[0][4])+'\n')
            f.write('\n\n\n')
            f.write('INDIVIDUAL PROTEIN IDENTIFICATIONS \n\n')
            
            for key in dict_keys:
                parameters = sample[key][1]
                f.write('Protein: '+key+'\n')
                f.write("Peptide sequence\tUnique peptide\tError(ppm) \n")
                try:
                    row = parameters.shape[0]
                    for n in range(row):
                        f.write(parameters.iloc[n,0]+'\t'+str(parameters.iloc[n,1])+'\t'+
                                str(parameters.iloc[n,2])+'\n')
                        
                except:
                    f.write('-'+'\t'+'-'+'\t'+'-'+'\n')
                f.write('\n\n')
            
    
    
  
    
  
    
  
    
# if __name__ == '__main__':
#     single_exportation(scores2, 'prueba', 'a.txt')
#     multiple_exportation(p, 'a2.txt')