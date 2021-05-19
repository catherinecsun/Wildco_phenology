# Wildco_phenology
examples of how to extract phenology from CT images

There are 4 scripts for extracting phenology/productivity data from CT timelapse images: one draws the region of inference; a second extracts the vegetation indices; a third fits curves to extract the phenological dates; a fourth organizes all the results from all the CT sites together. The output of each of the first 3 scripts are put into separate folders.

Timelapse images are organized by CT site: each CT site is a folder with subfolders for each deployment period. 

Images should follow this nomenclature: 'Site__Year-Month-Day__Hour-Min-Sec.JPG'; e.g., "CATH08__2019-07-04__12-00-00.JPG"
Deployments Folders should follow this nomeclature: 'Year-Month-Day to Year-Month-Day'; e.g., "2019-07-02 to 2019-08-21"
