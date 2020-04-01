# Balearic-epi

This repository has been created due to the Covid-19 global health crisis to do research on the impact of epidemics in the Balearic Islands.

The next sections describe the folders in this repository.

## Model

The model that we are using for simulating the epidemics is a compartmental model described by ... **MODEL DESCRIPTION MISSING**

## Procedure for fitting the model to IB data

For this procedure we simulate the model for many combinations of the parameters and find the one that best matches the observed data. 

### Epidemiological data

Time series available are

* Deaths (official data)
* Deaths (press data)
* UCI patients
* Recovered patients
* Infected patients
* Hospitalized patients

### Results of the model

Time series from the model

* Infected
* Recovered

### Fitting model to data

In order to fit the model to the data we proceed in the following way for each combination of parameters:

1. Fix a proportion of recovered patients that actually die (mu_D). 
2. Shift rescaled curve of recovered patients by tau from the model to best fit the curve of deaths from the data. The best fit is assessed by a Chi squared function. *Match with last point of deaths maybe.*
3. With that tau move also the curve of infected from the model and rescale it to best match the curve of UCI patients. Best value of this variable also given by Chi squared (mu_UCI).
4. Assign to this set of parameters the Chi squared that is the sum of the ones in 2. and 3. and save it together with the parameter values, the proportion of deaths and the obtained parameters tau and proportion of infected in UCI.

From all the Chi squared values choose the one that is minimum.

## Folder: Data

The data folder aims to containing data related to

 - Demographic data from [Wikipedia](https://es.wikipedia.org/wiki/Anexo:Municipios_y_comarcas_de_las_Islas_Baleares).
 - Mobility data from the INE regarding the [census 2011](https://www.ine.es/censos2011/tablas/Inicio.do). To obtain these data: crear tabla, Municipio, desglosar para cada unidad administrativa=SI, personas, ocupados de 16 o mas años, filas=municipio de residencia, columnas=datos del lugar de trabajo/lugar de trabajo/municipio de trabajo.
 - Geographic data (shapes of municipalities) from [ESRI](http://opendata.esri.es/datasets/53229f5912e04f1ba6dddb70a5abeb72_0). Data filtered with NUT code ES53 which corresponds to the Balearic Islands.
 - Epidemiological data from the [Spanish Ministry of Health](https://www.mscbs.gob.es/en/profesionales/saludPublica/ccayes/alertasActual/nCov-China/documentos/Actualizacion_61_COVID-19.pdf)
 - Other data related to epidemiological data? **MISSING** 

## Folder: Codes

Codes for 

 - Analysing data (see [Data](#Data))
 - Simulating different types of epidemics and containment strategies
 - Ploting the results.

 ## Folder: Results

Results of the analysis and simulations.

## Folder: Figures

Figures generated with the codes (see [Codes](#Codes)).

## Reports

Reports of the results from the simulations.


```{r, echo=FALSE}
htmltools::includeHTML("./figures/ib_poblacion.html")
```