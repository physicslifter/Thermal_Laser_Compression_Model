#%%
import os
import mdutils
from mdutils.mdutils import MdUtils
from mdutils import Html
import markdown
import pickle

folderpath = '../../winhome/Desktop/optimization_data/20220425/LogTest4'
#folderpath = 'C:/Users/Pat/Desktop/optimization_data/20220425/LogTest4'

class SummaryWriter:
    def __init__(self, folderpath):
        self.path = folderpath
        self.run_name = self.path.split('/')[-1]
        self.has_metadata = False
        self.has_optimization_data = False
        self.has_initial_plot = False
        self.has_log = False
        self.has_best_fit = False
        self.has_final_plot = False
        
    def check_files(self):
        for item in os.listdir(self.path):
            #print(item)
            if item == 'metadata.txt':
                self.has_metadata = True
            elif item == 'Best_Fit.png':
                self.has_best_fit = True
            elif item == 'Initial_Fit.png':
                self.has_initial_plot = True
            elif item == 'Optimization_Results.png':
                self.has_final_plot = True
            elif item == 'optimization_output.csv':
                self.has_optimization_data = True
            elif item == 'log.txt':
                self.has_log = True
                
    def write_summary(self):
        #Writing the markdown file
        md_fp = self.path + '/Summary.md'
        html_fp = self.path + '/Summary.html'
        md_title = self.run_name + ': Summary'
        mdFile = MdUtils(file_name = md_fp, title = md_title)
        mdFile.new_header(level = 1, title = 'Overview')
        mdFile.new_header(level = 2, title = 'Log Comments')
        if self.has_log:
            log_path = self.path + '/log.txt'
            with open(log_path, 'r') as f:
                paragraph = f.read()
                mdFile.new_paragraph(paragraph)
        mdFile.new_header(level = 2, title = 'Metadata')
        if self.has_metadata:
            metadata_path = self.path + '/metadata.txt'
            with open(metadata_path, 'r') as f:
                text = f.read()
                text = text.replace('===','').replace('==','')
                text = text.split('\n')
                c=0
                to_print_list = [3, 4, 5, 6, 7, 8] #List containing the lines from metadata that we want printed
                for t in text:
                    printLine = False
                    for l in to_print_list:
                        if l == c:
                            printLine = True
                    if printLine:
                        mdFile.new_paragraph(t)
                    c=c+1
                
        mdFile.new_header(level = 1, title = 'Plots')
        if self.has_final_plot:
            mdFile.new_header(level = 2, title = 'Optimization Results')
            img_path = 'Optimization_Results.png'
            img_line = mdFile.new_inline_image(text = 'Best Fit', path = img_path)
            mdFile.write(text = img_line)
        if self.has_best_fit:
            mdFile.new_header(level = 2, title = 'Best Fit')
            img_path = 'Best_Fit.png'
            img_line = mdFile.new_inline_image(text = 'Best Fit', path = img_path)
            mdFile.write(text = img_line)
            best_fit_fp = self.path + '/Best_Fit.pickle'
            best_fit_file = open(best_fit_fp, 'rb')
            best_fit = pickle.load(best_fit_file)
            if len(best_fit) == 7: #If there is a c param
                mdFile.new_header(level = 3, title = "Best Fit Parameters")
                #mdFile.new_paragraph(text = 'Best fit parameters:  ')
                mdFile.new_paragraph(text = 'iteration: '+str(best_fit[0]))
                mdFile.new_paragraph(text = 'sls: '+str(best_fit[1]))
                mdFile.new_paragraph(text = 'a: '+str(best_fit[2]))
                mdFile.new_paragraph(text = 'b: '+str(best_fit[3]))
                mdFile.new_paragraph(text = 'c: '+str(best_fit[4]))
                mdFile.new_paragraph(text = 'start time: '+str(best_fit[5]))
                mdFile.new_paragraph(text = 'peak temp: '+str(best_fit[6]))
            else:
                mdFile.new_header(level = 3, title = "Best Fit Parameters")
                #mdFile.new_paragraph(text = 'Best fit parameters:  ')
                mdFile.new_paragraph(text = 'iteration: '+str(best_fit[0]))
                mdFile.new_paragraph(text = 'sls: '+str(best_fit[1]))
                mdFile.new_paragraph(text = 'a: '+str(best_fit[2]))
                mdFile.new_paragraph(text = 'b: '+str(best_fit[3]))
                mdFile.new_paragraph(text = 'start time: '+str(best_fit[4]))
                mdFile.new_paragraph(text = 'peak temp: '+str(best_fit[5]))
        if self.has_initial_plot:
            mdFile.new_header(level = 2, title = 'Initial Fit')
            img_path = 'Initial_Fit.png'
            img_code_string = '![Initial Fit]('+img_path+')'
            img_line = mdFile.new_inline_image(text = 'Initial_Fit', path = img_path)
            mdFile.write(text = img_code_string)
        
            #self.b = best_fit
            
        
        mdFile.create_md_file()
        
        
        #Writing the HTML file
        with open(md_fp, 'r') as f:
            text = f.read()
            html = markdown.markdown(text)
        with open(html_fp, 'w') as f:
            f.write(html)
        
a = SummaryWriter(folderpath)
# %%
from mdutils.mdutils import MdUtils
from mdutils import Html

fn = 'Tests/markdowntest.md'
mdFile = MdUtils(file_name = fn, title = 'Testing New Markdown')
mdFile.new_header(level = 1, title = 'Overview')
mdFile.new_paragraph('Howdy doody')

mdFile.new_table_of_contents(table_title = 'Contents', depth = 2)
mdFile.create_md_file()

import markdown
with open(fn, 'r') as f:
    text = f.read()
    html = markdown.markdown(text)
html_name = fn.split('.')[0]+'.html'
with open(html_name, 'w') as f:
    f.write(html)
# %%
