###### THIS SCRIP IS A TEST WITH SAMPLES PREPARED FROM A MIX OF OTHER SAMPLES PREVIOUSLY COLLECTED
###### ALL THESE SAMPLES ARE STORED IN PROGRAM_TEST FOLDER

import score_matchms as sm
import os

dic = {}
for f in 'B', 'D', 'F':
    for n in range(12):
        n = n+1

        path_ = 'test_program'
        file_ = 'pmf_' + f + str(n) + '.txt'
        
        dir_ = os.path.join(path_,file_)

        scores= sm.request_scores(dir_, tolerance= 0.05, shift= 0, mz_power= 0.0,
                                            intensity_power = 1.0, dat_b='Aquasearch_study')

        dic[file_] = scores
        print(file_)
