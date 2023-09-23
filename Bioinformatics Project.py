#!/usr/bin/env python
# coding: utf-8

# In[1]:


import Bio.SeqIO as BSIO
from Bio import AlignIO
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from Bio import SeqIO
#BE.email = 'hodane@utu.fi'



# In[2]:


# Load the alignment from file
alignment = AlignIO.read("eukaryotic delta sample.fasta", "fasta")

columns_of_interest = [875, 877, 980]
MSA_sequences_collections = []
for i in range (0 , len(alignment)):
    if "-" not in (alignment[i].seq[columns_of_interest[0]-1] and alignment[i].seq[columns_of_interest[1]-1] and alignment[i].seq[columns_of_interest[2]-1]):
        MSA_sequences_collections.append(alignment[i])
MSA_sequences_collections
SeqIO.write(MSA_sequences_collections, "hosein_sequences_collection.fasta", "fasta") # a new file



# In[3]:


#alignment = AlignIO.read("hosein_sequences_collection.fasta", "fasta")# renamed by above new file "hosein_sequences_collection.fasta"
align_array = np.array(alignment) 
distributions=[]
for i in range (0 , len(alignment)):
    count = np.count_nonzero(align_array[i:i+1, columns_of_interest[0]+1:columns_of_interest[2]] == "-")
    print('Total occurences of gaps in the sequence between the 2nd & 3rd (H) & ungapped distances:',count,"&",columns_of_interest[2]-columns_of_interest[0]-2-count)
    distributions.append(columns_of_interest[2]-columns_of_interest[0]-2-count)
#print(distributions)
#plt.hist(distributions) # histogram without curve fitting
#plt.show()

mu, std = norm.fit(distributions) 
# Plot the histogram.
plt.hist(distributions, bins=25, density=True, alpha=0.6, color='b')
# Plot the Probability density function curve:
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = "Fit Values: {:.2f} and {:.2f}".format(mu, std)
plt.title(title)
# plt.show()


# ### (3) 

# In[4]:


np.random.seed(10)
data_1 = distributions
smlp_2 = np.random.normal(90, 5, 200) #sample
smlp_3 = np.random.normal(80, 10, 200) #sample
smlp_4 = np.random.normal(70, 15, 100) #sample
data = [data_1,smlp_2,smlp_3,smlp_4]
fig = plt.figure(figsize =(8, 5))
# Creating axes instance
ax = fig.add_axes([0, 0, 1, 1])
# Creating plot
bp = ax.boxplot(data, labels=["delta CAs","smlp_2","smlp_3","smlp_4"])
# show plot
plt.show()




# In[5]:


id_positions=[]
for i in range (0 , len(alignment)):
    cnt = np.count_nonzero(align_array[i:i+1, 0:columns_of_interest[0]+1] == "-")#gaps before the first clmn intrst
#    cnt2= np.count_nonzero(align_array[i:i+1, columns_of_interest[0]+1:columns_of_interest[1]+1] == "-") betweenfirst&second?
    count = np.count_nonzero(align_array[i:i+1, columns_of_interest[0]+1:columns_of_interest[2]] == "-")#gapsbetweensecond&third
    lst1= [alignment[i].id,columns_of_interest[0]-cnt,columns_of_interest[1]-cnt, columns_of_interest[2]-cnt-count]
    id_positions.append(lst1)

letters=[]
for record in alignment:
    amino_acids = ([record.seq[column-1] for column in columns_of_interest])
    lst2=amino_acids
    letters.append(lst2)

for j in range (0 ,len(id_positions)):
    print(f"{id_positions[j]} : {letters[j]}")
#    print(id_positions[j]+ letters[j])      OR
#    print(*indx[j]+ indx2[j], sep = ', ')   OR



# In[6]:


import pandas as pd
counts,values = pd.Series(letters).value_counts().values, pd.Series(letters).value_counts().index
df_results = pd.DataFrame(list(zip(values,counts)),columns=["value","count"])
df_results


# In[7]:


seprated_ids=[]
for j in range (0 ,len(id_positions)):
    seprated_ids.append(id_positions[j]+ letters[j])

for i in range (0 , len(seprated_ids)):
    if seprated_ids[i][4:7]==df_results.value[0]:
        print (df_results.value[0], "=", seprated_ids[i] [:1])
print('─' * 35)
for i in range (0 , len(seprated_ids)):
    if seprated_ids[i][4:7]==df_results.value[1]:
        print (df_results.value[1], "=", seprated_ids[i] [:1])
print('─' * 35)


# In[ ]:




