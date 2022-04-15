#%%
import main 


path='../../winhome/Desktop/optimization_data/20220413/Eq3_AllFaces'
q=main.PostOptOperations(path)
r=q.get_best_fit()
#sl=''
#for line in r:
#    sl=sl+line
    
q.group.plot_optimization_group(fname='Test_Save_1.png')
    
#%%
q.get_best_fit()

# %%
import pandas as pd
q1=pd.read_csv(q.fp+'/Clean_Opt_Data.csv')

#%%
q.get_best_fit()
# %%
in_dict=False #Boolean val for whether we are inside a dictionary
is_key=False #Boolean val for whether the current char is part of a dictionary key
looking_for_val=False #Boolean val for whether we are currently looking for the val in a key/val pair
is_val=False #Boolean for telling us whether the current char is part of a string associated w/ a value
is_num=False #Boolean for telling whether 
has_key=False #Boolean indicator for whether or not we currently have a dictionary key
current_vals=[] #Empty list for storing initial vals to first key
my_dict={}
for char in sl: #Iterate through each char in the single line
    
    #=====================================================================================
    #Detecting whether we are in a dictionary
    
    if char=="{": # { is the character denoting the beginning of a dictionary
        in_dict=True
    elif char=="}": # } is the char denoting the closing/ending of a dictionary
        in_dict=False
    #=====================================================================================
    
    #=====================================================================================
    #Detecting whether current char is part of a dictionary key
    
    if in_dict==True: #We only care to find keys if we are currently in a dictionary
        if char=="'": #Quotation marks signal the end/beginning of a key
            if is_key==True: #If we're in a key, this signals the end of the key
                is_key=False
                looking_for_val=True #Now that we've found the key, we're looking for associated values
                has_key=True
            else: #If we are not in a key, this signals the beginning of a key
                is_key=True
                key_string='' #Blank string for saving the key once we have each char associated w it
        else:
            if is_key==True: #If we are currently getting a key, add this character to the key string
                key_string=key_string+char
                
    #=====================================================================================
    
    #=====================================================================================
    #Getting Values for the current key
    
    if looking_for_val==True: #If we're looking for a value
        if is_val==True: #If we've detected a val
            if is_num==False: #If we have detected a val, but not yet a number
                if char.isnumeric()==True: #If our character is a number
                    is_num=True #We are currently at a number, which is part of the val, so set is_num to True
                    string_num=char #Set a string character of our current number so we can access it when we loop over the next character
            else: #otherwise, if we are currently looking at a number
                if char.isnumeric()==True: #There is another char in the number, so append this value
                    string_num=string_num+char
                else: #If the character is not numeric
                    if char=='.': #If the non-numeric character is a decimal, then we are still in a #
                        string_num=string_num+char 
                    else: #If the string is not part of a number...
                        is_num=False
                        current_vals.append(int(string_num)) #write integer form of string_num to current_vals
        if char=='[': #If the character is an opening bracket, this is the beginning of our values
            is_val=True #Change is_val to True so that we know we are in a value string
        if char==']': #This signals the end of a set of values
            is_val=False
            looking_for_val=False
            my_dict[key_string]=current_vals #Write the new entry into the dictionary
            current_vals=[] #Reset current_vals to be an empty list
            has_key=False
                
# %%
