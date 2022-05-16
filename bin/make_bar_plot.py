# packages

import matplotlib.pyplot as plt
import pandas as pd
import glob
import seaborn as sns

def parse_table(results_folder):
  

   txtfiles = []
   for file in glob.glob(str(results_folder)+"/test*.txt"):
      txtfiles.append(file)

   print(txtfiles)
   list_melt = []
   for file in txtfiles: 
      data=pd.read_csv(file,sep='\t')
      list_melt.append(data.melt( id_vars = ['tool'],var_name="type", 
       value_name="percentage"))
   

   data_concat = pd.concat(list_melt)
   print(data_concat)
   return data_concat

def bar_plot(data_list):
   # fig = plt.figure()
   # ax = fig.add_axes([0,0,1,1])
   # langs = list(data['tool'])
   # students = list(data['MPR1'])
   # ax.bar(langs,students)
   # plt.show()
   # plt.savefig('bar_plot.png')
   
     
   plot = sns.barplot(x = 'tool', y = 'percentage', hue = 'type', data=data_list, ci ="sd")
   fig = plot.figure 
   fig.savefig('bar_plot.png')


#main 
if __name__ == "__main__":
   data=parse_table("results")

   bar_plot(data)
    
