import pandas as pd 
import os
import urllib.request

metadata = pd.read_csv('../data/gtex.gtex.BRAIN.MD', delimiter='\t')
urldata = pd.read_csv('../data/recount3_BigWig_list_BRAIN_URLs.csv', index_col = False)

print(metadata.head())
print(metadata.shape)
print(metadata.columns)
print(metadata['DTHHRDY'].describe())
print(metadata[:4].transpose())


print(metadata['SMTSD'].unique())
print(metadata['SMTSD'].value_counts())
print(metadata['AGE'].unique())
print(metadata['AGE'].value_counts())

print(urldata.columns)
print(urldata.describe())
cereb = []
cortex = []
for index, row in metadata.iterrows():
	if row['SMTSD']== 'Brain - Cortex' and row['SEX'] == 2 and row['AGE'] == '40-49':
		cortex.append(row)
	elif (row['SMTSD']== 'Brain - Cerebellum') and row['SEX'] == 2 and row['AGE'] == '40-49':
		cereb.append(row)
all_files = cereb + cortex
print(len(all_files))
print()
for file in all_files:
	#continue
	print(file['rail_id'])
	target = urldata[urldata['rail_id']==file['rail_id']]
	# print(target)
	# print(target['BigWigURL'].shape)
	# print(target['BigWigURL'].item())
	os.system(f'wget {target['BigWigURL'].item()} -P ../data/preprocessing')
	#break
