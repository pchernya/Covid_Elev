# covid_elev
Code and data for the manuscript "Impact of altitude on COVID-19 infection and death in the United States: a modeling and observational study" submitted to PLOS One. 
Authors: Kenton E Stephens, Pavel Chernyavskiy, and Danielle R Bruns

Unfortunately, due to the size of the spatial object, we could not upload cnty_cases with its polygon geomtery field included as a column. However, we upload cnty_cases as a csv file with a FIPS code, such that it can be matched to any existing continental US shapefile.

The variable RUCC_2013 is a 9-level categorical covariate that reflects the official 2013 USDA Rural-Urban Continuum Codes. 2013 is the most recent year for which data are available. RUCC_2013 = 1 represents the most urban envrionment; RUCC_2013 = 9 represents the most rural environment. For more information and the source file please see (https://www.ers.usda.gov/data-products/rural-urban-continuum-codes/).
