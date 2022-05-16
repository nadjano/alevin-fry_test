# packages

import matplotlib.pyplot as plt

def parse_table(results_folder):
   


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
   bar_plot()
    
