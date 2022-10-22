import re
import operator

"""
extract groupspecific proteins from Subset_X.fasta to dictionary
"""

list_group = []
# subset = open("I:/Methoden/ProteinOrtho/group_specific/results_subsets/subset_3.fasta", "r") #read subset.fasta
subset = open("I:/Methoden/KEGG pathways/KEGG_Pathway_18_10_22/results_subset/subset25.fasta", "r") #read subset.fasta

for line in subset: #search each line in file
    if re.search(">", line): #search after ">" 
        list_group.append(line.split()[1])

# print(list_group)

"""
you need to run isabels script before you can accomplish this step
open G1A-X_ko_pathways.txt and extract those lines (Locustag + KO + which pathway) which are included in the specific groups
"""

pathway = open("I:/Methoden/KEGG pathways/KEGG_Pathway_18_10_22/refs_combined_KO_pathways.txt", "r") #read pathway file

No_of_K0_hits = 0
dictionary = {}

for line in pathway: # search each line in file
    for element in list_group: # search each element of subset
        if re.search(element, line):
            # print(line)
            
            if "\t\t" in line: # if protein go K0-Number
                
                pathways_combined = line.split("\t\t")[-1] # get everything after the K0 Nr  
                pathways_combined = pathways_combined.replace("\n", "") # get rid of \n
                
                
                
                if " // " in pathways_combined: # if pathway contains more than 1 (sep. with //)
                    pathways_combined = pathways_combined.split(" // ")
                    
                    for each_split in pathways_combined:
                        
                        if each_split in dictionary: #if pathway already in  dictionary
                            dictionary[each_split] = dictionary.get(each_split) + 1 #change value to += 1
                        
                        else:
                            dictionary.update({each_split: 1})
                            
                else: #if only 1 pathway 
                    if pathways_combined in dictionary: #if pathway already in  dictionary
                        dictionary[pathways_combined] = dictionary.get(pathways_combined) + 1 #change value to += 1
                        
                    else:
                        dictionary.update({pathways_combined: 1})
                        
                print(pathways_combined)  
                No_of_K0_hits += 1
                

print("\nThere are a total of", No_of_K0_hits, "hits within the groupspecific proteins, that characterise at least one pathway\n")

sorted_dictionary = dict(sorted(dictionary.items(), key=operator.itemgetter(1),reverse=True))
print("Dictionary in descending order by value : ",sorted_dictionary, "\n")

#sollten insagesamt 63 Pathways sein



"""
safe as csv-file
"""
# with open("C:/Users/phh20/Desktop/KEGG pathways/subset_1_pathways.csv", "w") as f:
with open("I:/Methoden/KEGG pathways/KEGG_Pathway_18_10_22/results_subset/subset25_pathways_reduced.csv", "w") as f:
    
    
    f.write("Specific KEGG-Pathway \t Number of assigned K0´s\n")
    
    for key in sorted_dictionary.keys():
        f.write("%s\t %s\n" % (key, sorted_dictionary[key])) #\t tab as seperator

print("Finished, separator of csv-file is \ t")