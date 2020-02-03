import os 
import glob
import pandas as pd
os.chdir("/home/mairin/Documents/GradSchool/PhD/Research/CircuitTheory_Borneo/PRSB_Revision2/Nmix_Models/Oct22_ModelFittingResults/Oct22_2019_ModelFittingResults_AICtables_csv/")

extension = 'csv'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]

# combined all files
combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])
combined_csv.to_csv("Oct22_2019_AllAICTables.csv", index = False, encoding = 'utf-8-sig')
