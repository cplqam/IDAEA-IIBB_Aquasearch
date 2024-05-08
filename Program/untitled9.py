import load_archives
import numpy as np



request_xml = 'test_files/mcE61_Figueres.xml'
tolerance = 0.05
shift =  0
mz_power = 0
intensity_power = 1.0
unique_ = 1.5
dat_b='Aquasearch_study'


mz_int_p, u_inf_p = load_archives.db_table_request('P0DUB6')
df, spec1_mz = load_archives.parse_xml('test_files/mcE61_Figueres.xml')

u_inf2 = []
for x in u_inf_p:
    if x == 'Unique':
        u_inf2.append(1)
    elif x == 'No unique':
        u_inf2.append(0)

lowest_idx = 0 
matches = []
for peak1_idx in range(spec1_mz.shape[0]): 
    mz = spec1_mz[peak1_idx, 0] 
    low_bound = mz - tolerance 
    high_bound = mz + tolerance
    for peak2_idx in range(lowest_idx, mz_int_p.shape[0]):  
        mz2 = mz_int_p[peak2_idx, 0] + shift  
        if mz2 > high_bound: 
            break
        if mz2 < low_bound:  
            lowest_idx = peak2_idx
        else:  
            matches.append((peak1_idx, peak2_idx, u_inf2[peak2_idx]))
matches_p = np.array(matches)

#---------------------------------------------#

dat_base = mz_int_p[matches_p[:,1]]
dat_base = np.array([dat_base[:,0], dat_base[:,1], matches_p[:,2]]).transpose()
req = spec1_mz[matches_p[:,0],:]
req[:,1] = (req[:,1]/max(req[:,1]))*100

if len(req) == 0: 
    pass
matching_pairs = [] 

for i, idx in enumerate(req):
    power_prod_spec1 = (req[i, 0] ** mz_power) * (req[i, 1] ** intensity_power) 
    power_prod_spec2 = (dat_base[i, 0] ** mz_power) * (dat_base[i, 1] ** intensity_power)
    if dat_base[i, 2] == 0:
        matching_pairs.append([req[i,0], dat_base[i,0], power_prod_spec1 * power_prod_spec2, 0])
    elif dat_base[i, 2] == 1:
        matching_pairs.append([req[i,0], dat_base[i,0], ((power_prod_spec1 * power_prod_spec2) ** unique_), 1])
        
matching_pairs = np.array(matching_pairs)


#------------------------------------#

score = float(0.0)
used_matches = int(0)
used1 = set()
used2 = set()
unique_used = int(0)
rep = np.empty((0, 4), int)
signals_used = np.empty((0, 3), int)

for i in range(matching_pairs.shape[0]):    
    mz_1 = matching_pairs[i, 0]
    mz_2 = matching_pairs[i, 1]
    s = matching_pairs[i, 2]
    
    # If the peak only appears 1 time, it is added directly to the score
    if np.count_nonzero(matching_pairs[:, 0] == mz_1) == 1 and np.count_nonzero(matching_pairs[:, 1] == mz_2) == 1:
        score += matching_pairs[i, 2] 
        used1.add(matching_pairs[i, 0])  
        used2.add(matching_pairs[i, 1])  
        used_matches += 1
        unique_used += matching_pairs[i, 3]
        
        arr = np.array([matching_pairs[i, 0],matching_pairs[i, 1],matching_pairs[i, 2]]).transpose()
        signals_used = np.concatenate((signals_used, [arr]), axis=0)
        
    # If not, it is added to an array where the repeats will be
    else:
        rep = np.append(rep, np.array([[mz_1, mz_2, s, matching_pairs[i,3]]]), axis=0)
   
rep_def = np.empty((0, 4), int)

# Within the repeated
for i in range(rep.shape[0]):
    mz_1 = rep[i, 0]
    mz_2 = rep[i, 1]
    s = rep[i, 2]
    u = rep[i, 3]
    
    if i == 0:
        # If it is the 1 element, it enters the definitive array
        rep_def = np.append(rep_def, np.array([[mz_1, mz_2, s, u]]), axis=0)
    else: 
        # If it is not the first, it is checked that there is no equal mz in the list, of any of the 2 spectra
        # If there is not the same mz
        if np.count_nonzero(rep_def[:, 0] == mz_1) == 0 and np.count_nonzero(rep_def[:, 1] == mz_2) == 0:
            # They enter the final array
            rep_def = np.append(rep_def, np.array([[mz_1, mz_2, s, u]]), axis=0)
            
        # If there is already an equal mz in the list,
        # I compare the scores and leave the best one in that position of the array
        elif np.count_nonzero(rep_def[:, 0] == mz_1) != 0:  # Compare with first data, if the first is:
            rep1 = np.where(rep_def[:, 0] == mz_1)           # Get the index
            # If its score is better, the values are replaced to get the best score for each signal
            if s > rep_def[rep1, 2]:
                rep_def[rep1, 0] = mz_1
                rep_def[rep1, 1] = mz_2
                rep_def[rep1, 2] = s
                rep_def[rep1, 3] = u
        
        elif np.count_nonzero(rep_def[:, 1] == mz_2) != 0:  # The same is done for the second mz
            # If the two mz match, then it would take the best of the 3 (the two matches and the new one)
            rep2 = np.where(rep_def[:, 1] == mz_2)
            if s > rep_def[rep2, 2]:
                rep_def[rep2, 0] = mz_1
                rep_def[rep2, 1] = mz_2
                rep_def[rep2, 2] = s
                rep_def[rep1, 3] = u

for i in range(rep_def.shape[0]): 
    score += rep_def[i, 2] 
    used1.add(rep_def[i, 0])
    used2.add(rep_def[i, 1])  
    used_matches += 1      
    unique_used += rep_def[i, 3]
    
    arr = np.array([rep_def[i, 0],rep_def[i, 1],rep_def[i, 2]]).transpose()
    signals_used = np.concatenate((signals_used, [arr]), axis=0)

    
# Normalize score:
# spec1_power = (spec1_mz[:, 0]/max(spec1_mz[:, 0])) ** mz_power * (spec1_mz[:, 1]/max(spec1_mz[:, 1])) ** intensity_power
# spec2_power = mz_int_p[:, 0] ** mz_power * mz_int_p[:, 1] ** intensity_power

# score = score/(np.sum(spec1_power ** 2) ** 0.5 * np.sum(spec2_power ** 2) ** 0.5 * used_matches)

i = signals_used[:,0]**mz_power * signals_used[:,2] ** intensity_power * used_matches** (unique_used/used_matches)

score = np.sum(i)/score * 10