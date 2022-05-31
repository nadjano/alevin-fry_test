# packages

import matplotlib.pyplot as plt
import pandas as pd
import glob
import seaborn as sns

def parse_table(results_folder):
  

   txtfiles = []
   for file in glob.glob(str(results_folder)+"/*.txt"):
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

def bar_plot(data_list1):
   # sns.set(rc={'figure.figsize':(11.7,6)})
   # fig, ax =plt.subplots(1,2,  sharex='col', sharey='row')
     
   plot = sns.barplot(x = 'tool', y = 'percentage', hue = 'type', data=data_list1, capsize=.1)#  ax=ax[0])
   plot.set(title = "scRNA-seq")
   # plot = sns.barplot(x = 'tool', y = 'percentage', hue = 'type', data=data_list2, capsize=.1, ax=ax[1])
   # plot.set(title = "snRNA-seq")
   fig = plot.figure 
   fig.savefig('snRNA-mapping-rate/figures/bar_plot_sc_1.png')


#main 
if __name__ == "__main__":
   data=parse_table("results_sc")

   data_sn = parse_table("results")

   bar_plot(data)
    
