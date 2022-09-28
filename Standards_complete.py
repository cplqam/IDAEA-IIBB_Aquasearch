# simple spcrip to complete the standards information using standard_incorporation.py 

import standard_incorporation as si
import os

path__ = 'test_files/Standares'

unique_file = 'Aquasearch_Proteins_Unique.xlsx'
maldi_H1 = 'pmf_H1.txt'
maldi_H2 = 'pmf_H2.txt'
maldi_H3 = 'pmf_H3.txt'
maldi_H4 = 'pmf_H4.txt'
maldi_H5 = 'pmf_H5.txt'
maldi_H6 = 'pmf_H6.txt'
maldi_H7 = 'pmf_H7.txt'
maldi_H8 = 'pmf_H8.txt'
maldi_H9 = 'pmf_H9.txt'
maldi_H10 = 'pmf_H10.txt'
maldi_H11 = 'pmf_H11.txt'
maldi_H12 = 'pmf_H12.txt'
maldi_H13 = 'pmf_H13.txt'
maldi_H14 = 'pmf_H14.txt'
maldi_H15 = 'pmf_H15.txt'

unique_path = os.path.join(path__, unique_file)
H1_path = os.path.join(path__, maldi_H1)
H2_path = os.path.join(path__, maldi_H2)
H3_path = os.path.join(path__, maldi_H3)
H4_path = os.path.join(path__, maldi_H4)
H5_path = os.path.join(path__, maldi_H5)
H6_path = os.path.join(path__, maldi_H6)
H7_path = os.path.join(path__, maldi_H7)
H8_path = os.path.join(path__, maldi_H8)
H9_path = os.path.join(path__, maldi_H9)
H10_path = os.path.join(path__, maldi_H10)
H11_path = os.path.join(path__, maldi_H11)
H12_path = os.path.join(path__, maldi_H12)
H13_path = os.path.join(path__, maldi_H13)
H14_path = os.path.join(path__, maldi_H14)
H15_path = os.path.join(path__, maldi_H15)


si.standard_complete(H1_path, unique_path,'P07724', 'standard_1', 2.5)
print(1)
si.standard_complete(H2_path, unique_path,'P07724', 'standard_2', 2.5)
print(1)
si.standard_complete(H3_path, unique_path,'P07724', 'standard_3', 2.5)
print(1)
si.standard_complete(H4_path, unique_path,'P02770', 'standard_4', 2.5)
print(1)
si.standard_complete(H5_path, unique_path,'P02770', 'standard_5', 2.5)
print(1)
si.standard_complete(H6_path, unique_path,'P02770', 'standard_6', 2.5)
print(1)
si.standard_complete(H7_path, unique_path,'P02768', 'standard_7', 2.5)
print(1)
si.standard_complete(H8_path, unique_path,'P02768', 'standard_8', 2.5)
print(1)
si.standard_complete(H9_path, unique_path,'P02768', 'standard_9', 2.5)
print(1)
si.standard_complete(H10_path, unique_path,'P08835', 'standard_10', 2.5)
print(1)
si.standard_complete(H11_path, unique_path,'P08835', 'standard_11', 2.5)
print(1)
si.standard_complete(H12_path, unique_path,'P08835', 'standard_12', 2.5)
print(1)
si.standard_complete(H13_path, unique_path,'P19121', 'standard_13', 2.5)
print(1)
si.standard_complete(H14_path, unique_path,'P19121', 'standard_14', 2.5)
print(1)
si.standard_complete(H15_path, unique_path,'P19121', 'standard_15', 2.5)

