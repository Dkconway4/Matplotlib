#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# Capomulin and Ramicane reduces the size of tumors the best.
# The 0.84 shows that there is a strong positive correlation, as the weight increases the average tumor volume increases.
# There was almost an equal percentage of male and female mice so you can infer that the data is fairly accurate.

# In[3]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import numpy as np

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
mouse_metadata = pd.read_csv(mouse_metadata_path)
mouse_metadata.head()

# Display the data table for preview


# In[10]:


# Merge data and count mice
combined_df = pd.merge(mouse_metadata, study_results, how='outer', on="Mouse ID")
combined_df.head()

num_mice = combined_df["Mouse ID"].nunique()
num_mice
#Why can't count work?


# In[11]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
dup_mice_ID = combined_df.loc[combined_df.duplicated(subset=['Mouse ID', 'Timepoint']),'Mouse ID'].unique()
dup_mice_ID


# In[12]:


# Optional: Get all the data for the duplicate mouse ID. 
dup_mice_df = combined_df.loc[combined_df["Mouse ID"] == "g989", :]
dup_mice_df


# In[13]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_df = combined_df[combined_df['Mouse ID'].isin(dup_mice_ID)==False]
clean_df.head()


# In[17]:


# Checking the number of mice in the clean DataFrame.
#Why can't .count work
clean_mice = clean_df["Mouse ID"].nunique()

clean_mice


# ## Summary Statistics

# In[21]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
mean = clean_df['Tumor Volume (mm3)'].groupby(clean_df['Drug Regimen']).mean()
median = clean_df['Tumor Volume (mm3)'].groupby(clean_df['Drug Regimen']).median()
var = clean_df['Tumor Volume (mm3)'].groupby(clean_df['Drug Regimen']).var()
std = clean_df['Tumor Volume (mm3)'].groupby(clean_df['Drug Regimen']).std()
sem = clean_df['Tumor Volume (mm3)'].groupby(clean_df['Drug Regimen']).sem()

summary_stat = pd.DataFrame({"Mean Tumor Volume":mean, 
"Median Tumor Volume":median, 
"Tumor Volume Variance":var, 
"Tumor Volume Std. Dev.":std, 
"Tumor Volume Std. Err.":sem})
# Assemble the resulting series into a single summary DataFrame.
print (summary_stat)


# In[22]:


# Generate a summary statistics table of mean, median, variance, standard deviation, 
# and SEM of the tumor volume for each regimen

# Using the aggregation method, produce the same summary statistics in a single line.
summary_aggregation =  clean_df.groupby(['Drug Regimen'])[['Tumor Volume (mm3)']].agg(['mean', 'median', 'var', 'std', 'sem'])
summary_aggregation


# ## Bar and Pie Charts

# In[24]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.
mice_count_df = clean_df["Drug Regimen"].value_counts()
mice_count_df

plot_pandas = mice_count_df.plot.bar(color='r')

plt.xlabel("Drugs")
plt.ylabel("Mice")
plt.title("Mice per Treatment")


# In[26]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using pyplot.

x_axis = mice_count_df.index.values
y_axis = mice_count_df.values

plt.bar(x_axis, y_axis, color='r', alpha=0.8, align='center')


plt.title("Number of Mice Tested per Treatment")
plt.xlabel("Drugs")
plt.ylabel("Mice")
plt.xticks(rotation="vertical")


plt.show()


# In[28]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas
#Can't use "Gender"
gender_data = clean_df["Sex"].value_counts()
plt.title("Female vs. Male Mice")
gender_data.plot.pie(autopct= "%1.1f%%")
plt.show()


# In[29]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot

labels = ['Female', 'Male']
sizes = [49.7999197, 50.200803]
plot = gender_data.plot.pie(y='Total Count', autopct="%1.1f%%")
plt.title('Male vs Female Mouse Population')
plt.ylabel('Sex')
plt.show()


# ## Quartiles, Outliers and Boxplots

# In[35]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

Capomulin_df = clean_df.loc[clean_df["Drug Regimen"] == "Capomulin",:]
Ramicane_df = clean_df.loc[clean_df["Drug Regimen"] == "Ramicane", :]
Infubinol_df = clean_df.loc[clean_df["Drug Regimen"] == "Infubinol", :]
Ceftamin_df = clean_df.loc[clean_df["Drug Regimen"] == "Ceftamin", :]


# In[36]:


# Start by getting the last (greatest) timepoint for each mouse
# Greatest is Cap

Capomulin_last = Capomulin_df.groupby('Mouse ID').max()['Timepoint']
Capomulin_vol = pd.DataFrame(Capomulin_last)
Capomulin_merge = pd.merge(Capomulin_vol, clean_df, on=("Mouse ID","Timepoint"),how="left")
Capomulin_merge.head()

#Cap Info
Capomulin_tumors = Capomulin_merge["Tumor Volume (mm3)"]

quartiles =Capomulin_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


print(f"The lower quartile of Capomulin tumors: {lowerq}")
print(f"The upper quartile of Capomulin tumors: {upperq}")
print(f"The interquartile range of Capomulin tumors: {iqr}")
print(f"The median of Capomulin tumors: {quartiles[0.5]} ")

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)

print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")


# In[37]:


# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
# Ram up next

Ramicane_last = Ramicane_df.groupby('Mouse ID').max()['Timepoint']
Ramicane_vol = pd.DataFrame(Ramicane_last)
Ramicane_merge = pd.merge(Ramicane_vol, clean_df, on=("Mouse ID","Timepoint"),how="left")
Ramicane_merge.head()
Ramicane_merge.to_csv("output.csv")
Ramicane_tumors = Ramicane_merge["Tumor Volume (mm3)"]

quartiles =Ramicane_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


print(f"The lower quartile of Ramicane tumors is: {lowerq}")
print(f"The upper quartile of Ramicane tumors is: {upperq}")
print(f"The interquartile range of Ramicane tumors is: {iqr}")
print(f"The median of Ramicane tumors is: {quartiles[0.5]} ")

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)

print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")


# In[38]:


# Infu Info

Infubinol_last = Infubinol_df.groupby('Mouse ID').max()['Timepoint']
Infubinol_vol = pd.DataFrame(Infubinol_last)
Infubinol_merge = pd.merge(Infubinol_vol, clean_df, on=("Mouse ID","Timepoint"),how="left")
Infubinol_merge.head()


Infubinol_tumors = Infubinol_merge["Tumor Volume (mm3)"]

quartiles =Infubinol_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


print(f"The lower quartile of Infubinol tumors is: {lowerq}")
print(f"The upper quartile of Infubinol tumors is: {upperq}")
print(f"The interquartile range of Infubinol tumors is: {iqr}")
print(f"The median of Infubinol tumors is: {quartiles[0.5]} ")

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)


print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")
Infubinol_merge.to_csv("output.csv")


# In[50]:


# Ceft Info
Ceftamin_last = Ceftamin_df.groupby('Mouse ID').max()['Timepoint']
Ceftamin_vol = pd.DataFrame(Ceftamin_last)
Ceftamin_merge = pd.merge(Ceftamin_vol, clean_df, on=("Mouse ID","Timepoint"),how="left")
Ceftamin_merge.head()


Ceftamin_tumors = Ceftamin_merge["Tumor Volume (mm3)"]

quartiles = Ceftamin_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq

print(f"The lower quartile of treatment is: {lowerq}")
print(f"The upper quartile of temperatures is: {upperq}")
print(f"The interquartile range of temperatures is: {iqr}")
print(f"The the median of temperatures is: {quartiles[0.5]} ")

# Determine outliers using upper and lower bounds

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)

print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")


# In[51]:


# Put treatments into a list for for loop (and later for plot labels)


# Create empty list to fill with tumor vol data (for plotting)


# Calculate the IQR and quantitatively determine if there are any potential outliers. 

    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    
    
    # add subset 
    
    
    # Determine outliers using upper and lower bounds
    
drug_list = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]
drugs = combined_df[combined_df["Drug Regimen"].isin(drug_list)]
drugs.head()

last_timepoint = drugs.groupby(["Drug Regimen", "Mouse ID"]).agg(tumor_size=("Tumor Volume (mm3)", lambda x: x.iloc[-1]))
    
treatment = 0
for drug in drug_list:
    quartiles = last_timepoint[drug].quantile([.25,.5,.75]).round(2)
    lowerq = quartiles[0.25].round(2)
    upperq = quartiles[0.75].round(2)
    iqr = round(upperq-lowerq,2)
    lower_bound = round(lowerq - (1.5*iqr),2)
    upper_bound = round(upperq + (1.5*iqr),2)


    if treatment == 0:
        print(f"------------------------------------------------------------")
    print(f"The lower quartile of {drug} treatments is: {lowerq}")
    print(f"The upper quartile of {drug} treatments is: {upperq}")
    print(f"The interquartile range of {drug} treatments is: {iqr}")
    print(f"Values below {lower_bound} could be {drug} outliers.")
    print(f"Values above {upper_bound} could be {drug} outliers.")
    print(f"------------------------------------------------------------")
    treatment+=1
    


# In[ ]:





# In[69]:


# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.


Final_Tumor_Size = dup_mice_df.drop_duplicates(subset = ("Mouse ID"), keep = "last")
Final_Tumor_by_drug = Final_Tumor_Size[["Mouse ID", "Drug Regimen", "Tumor Volume (mm3)"]]
Final_Tumor_Size_Capomulin = Final_Tumor_by_drug.loc[Final_Tumor_by_drug["Drug Regimen"] == "Capomulin"]
Final_Tumor_Size_Ramicane = Final_Tumor_by_drug.loc[Final_Tumor_by_drug["Drug Regimen"] == "Ramicane"]
Final_Tumor_Size_Infubinol = Final_Tumor_by_drug.loc[Final_Tumor_by_drug["Drug Regimen"] == "Infubinol"]
Final_Tumor_Size_Ceftamin = Final_Tumor_by_drug.loc[Final_Tumor_by_drug["Drug Regimen"] == "Ceftamin"]



box_plot_df = pd.DataFrame({"Capomulin" : Final_Tumor_Size_Capomulin["Tumor Volume (mm3)"],
                            "Ramicane" : Final_Tumor_Size_Ramicane["Tumor Volume (mm3)"],
                           "Infubinol" : Final_Tumor_Size_Infubinol["Tumor Volume (mm3)"],
                        "Ceftamin" : Final_Tumor_Size_Ceftamin["Tumor Volume (mm3)"]})




myboxplot = box_plot_df.boxplot(column = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"], grid = False,
                                flierprops = red_circle,)
                              
plt.title("Tumor Volume at Selected Mouse")
plt.xlabel("Drug")
plt.ylabel("Final Tumor Volume (mm3)")
plt.show()

    


# ## Line and Scatter Plots

# In[73]:


# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin

Line_plot = dup_mice_df.loc[dup_mice_df["Mouse ID"] == "|509"]
x_timepoint = Line_plot["Timepoint"]
y_tumor_volume = Line_plot["Tumor Volume (mm3)"]
plt.plot(x_timepoint, y_tumor_volume, color = "green", marker = "o", markerfacecolor = "white", linewidth = 2)
plt.title(f"Capomulin Treatment for Mouse |509")
plt.xlabel("Time Point in (Days)")
plt.ylabel("Tumor Volume (mm3)")
plt.show()


# In[75]:


# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen

fig1, ax1 = plt.subplots()
avg_cap_vol_df =Capomulin_df.groupby(['Mouse ID']).mean()

marker_size=15
plt.scatter(avg_cap_vol_df['Weight (g)'],avg_cap_vol_df['Tumor Volume (mm3)'], color="blue")
plt.title('Mouse Weight Versus Average Tumor Volume')
plt.xlabel('Weight (g)',fontsize =14)
plt.ylabel('Averag Tumor Volume (mm3)')


# ## Correlation and Regression

# In[18]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen


# In[ ]:




