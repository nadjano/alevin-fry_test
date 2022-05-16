# packages

import matplotlib.pyplot as plt
import pandas as pd
import glob

def parse_table(results_folder):
  

   txtfiles = []
   for file in glob.glob(str(results_folder)+"/*.txt"):
      txtfiles.append(file)


   for file in txtfiles: 
      data=pd.read_csv(file,sep='\t')
   
   print(data)


def bar_plot():
   fig = plt.figure()
   ax = fig.add_axes([0,0,1,1])
   langs = ['C', 'C++', 'Java', 'Python', 'PHP']
   students = [23,17,35,29,12]
   ax.bar(langs,students)
   plt.show()
   plt.savefig('snRNA-mapping-rate/bar_plot.png')


#main 
if __name__ == "__main__":
   parse_table("results")
    
