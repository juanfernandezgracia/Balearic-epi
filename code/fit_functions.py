import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import os
import glob
import sys
from matplotlib.ticker import ScalarFormatter
from matplotlib.dates import date2num, datestr2num

def shift(df, Ndays):
    '''
    shift(df, Ndays) function: 
        shift the dataframe df Ndays to the past
    Inputs:
        >> df: dataframe to be shifted (the one of simulations). Has to have a column 'Día'
        >> Ndays: integer number of days to shift dataframe
    Outputs:
        << dfnew: shifted dataframe.
    '''
    dates = list(df['Día'])
    timedelta = dt.timedelta(days = Ndays)
    newdates = [dates[i] - timedelta for i in range(len(dates))]
    dfnew = df.copy()
    dfnew['Día'] = newdates
    return dfnew

def chi2(df,
         cname1,
         cname2,
         par,
         chi2_type = 'linear',
         minvalue = 1):
    '''
    chi2(df, cname1, cname2, par) function: 
        compute sum of squares of the differences between the column cname1 and the column cnam2*par.
        The sum is restricted to those points where the column cname1 has data (notna) which are 
        over minvalue.
    Inputs:
        >> df: dataframe merged with epi data (cname1) and simu data (cname2)
        >> cname1: str name of the column of epi data.
        >> cname2: str name of the column of simu data.
        >> par: float parameter rescaling the curve in cname2
        >> chi_type: str that tells how to compute chi^2
                Optional, default is 'linear'
                'linear': sum of squares of the differences
                'log': sum of squares of the differences of logs
                'relative': sum of squares of the difference divided by the value of the data
        >> minvalue: int minimum value of the datapoints to fit (minvalue included)
                Optional, default is 1
    Outputs:
        << chi2_: float sum of squares of the differences between the column cname1 and the column cnam2*par.
    '''
    chichi = 0.0
    df2_ = df[df[cname1].notna()] # get rid of nans
    df2_ = df2_[(df2_[[cname1]] >= minvalue).all(axis=1)] # get rid of values below minvalue
    a = np.array(df2_[cname1])
    b = np.array(df2_[cname2])*par
    if chi2_type == 'linear':
        chichi = np.sum(np.square(a - b))
    elif chi2_type == 'log':
        a = np.log(a)
        b = np.log(b)
        chichi = np.sum(np.square(a - b))
    elif chi2_type == 'relative':
        chichi = np.sum(np.square(a - b)/a)
    del df2_
    return chichi

def best_fit(df_epi,
             fname_simu,
             mu_D = 0.01,
             Ndays = [35, 50],
             weight = 0.5,
             chi2_type = 'linear',
             minvalue = 1):
    '''
    best_fit(df_epi, fname_simu [, mu_D = 0.01, Ndays = 50, weight = 0.5, chi2_type = 'linear', minvalue = 1]) function: 
        Given a death rate mu_D and the results of a simulation this function first fits the recovered
        from the simulation multiplied by the death rate mu_D to the deaths from epi data by shifting 
        the recovered*mu_D in the simu by dayshift minimising the chi2_M. Then with that value for dayshift 
        the infected from the simu are shifted and a UCI rate is computed to minimise a chi2 with respect 
        to the epi UCIs (chi2_UCI). In the end one obtains the total chi2 which is the weighted sum of both
        chi2,the dayshift and the UCI rate mu_UCI
    Inputs:
        >> df_epi: dataframe with epi data
        >> fname_simu: str name of the file with the results of the simulation
        >> mu_D: float death rate which rescales the recovered from simu to match deaths from epi. 
                (optional, default is 0.01)
        >> Ndays: list of ints minimum and maximum days tried for the dayshift
                (optional, default is [35, 50])
        >> weight: float how to weight the chi2. 
                The chi2 of the deaths is weighted by weight and the one for UCI with (1-weight)
                (optional, default is 0.5)
        >> chi2_type: str that tells how to compute chi^2
                Optional, default is 'linear'
                'linear': sum of squares of the differences
                'log': sum of squares of the differences of logs
                'relative': sum of squares of the difference divided by the value of the data
        >> minvalue: int minimum value of the datapoints to fit (minvalue included)
                Optional, default is 1
    Outputs:
        << chichi: float weighted sum of chi2 for the best fit. 
                chichi = weight*chi_2_M + (1.0-weight)*chi_2_UCI
        << dayshift: int number of days that the simu has to be shifted for the best fit.
        << mu_UCI: float UCI rate to rescale infected that best fits the UCI from epi data.
                This parameter is scanned in the range (start = 0.00001, end = 0.005, step = 0.00001)
    '''
    df_simu = read_csv_VME(fname_simu)
    chi_2_M = float("inf")
    dayshift = 0
    for iday in range(Ndays[0], Ndays[1]):
        df_shifted = shift(df_simu, iday)
        df_merged = pd.merge(df_shifted, df_epi, on = 'Día')
        chichi = chi2(df_merged, 'M', 'average_recover', mu_D, chi2_type = 'linear', minvalue = 1)
        if chichi < chi_2_M:
            chi_2_M = chichi
            dayshift = iday
    v1=np.arange(0.00001,0.005,0.00001)
    df_shifted = shift(df_simu, dayshift)
    df_merged = pd.merge(df_shifted, df_epi, on = 'Día')
    #print(df_merged.head())
    #sys.exit()
    chi_2_UCI = float("inf")
    mu_UCI = 0.0
    for muci in v1:
        chichi = chi2(df_merged, 'UCI', 'average_infected', muci, chi2_type = 'linear', minvalue = 1)
        if chichi < chi_2_UCI:
            chi_2_UCI = chichi
            mu_UCI = muci
    chichi = weight*chi_2_M + (1.0-weight)*chi_2_UCI
    return chichi, dayshift, mu_UCI

def read_csv_VME(fname):
    '''
    read_csv_VME(fname) function: 
        reads results of simulations provided by the codes of VME. This function also adds a column 'Día'
        with datetime objects having the date.
    Inputs:
        >> fname: str path + filename to access the simulation data
    Outputs:
        << df: dataframe with the simulation data and an added column 'Día' with the dates as datetime objects
    '''
    fname2 = fname.split('.dat')[0]+'_format0.dat'
    os.system("awk '$1=$1' "+fname+" > "+fname2)
    df = pd.read_csv(fname2, sep = ' ')
    os.system('rm '+fname2)
    date0 = dt.datetime(day = 29, month = 2, year = 2020)
    datenums = list(df['#time'])
    dates = [date0 + dt.timedelta(days = datenums[i]) for i in range(len(datenums))]
    df['Día'] = dates
    return df

def limit_tseries(x, y, start_date, end_date):    
    '''
    limit_tseries(x, y, start_date, end_date) function: 
        limits a time series to the dates provided (both dates included)
    Inputs:
        >> x: list of datetime objects that provide dates for the data in y
        >> y: list with the data
        >> start_date: datetime object date of start of the resulting timeseries
        >> end_date: datetime object date of end of the resulting timeseries
    Outputs:
        << xlimit: list of datetime objects that provide dates for the data in y restricted to specified dates
        << ylimit: list with the data restricted to specified dates
    '''
    ddays = (end_date - start_date).days
    xlimit = []
    ylimit = []
    for i in range(len(x)):
        dd = (x[i] - start_date).days
        if dd >=0 and dd <= ddays:
            xlimit.append(x[i])
            ylimit.append(y[i])
    return xlimit, ylimit

def complete_fit(chi2_type = 'log',
                 minvalue = 10,
                 weight = 0.5,
                 percents = [0.01, 0.006, 0.034],
                 Ndays = [35, 50]):
    '''
    complete_fit(chi2_type = 'log',
                 minvalue = 10,
                 weight = 0.5,
                 percents = [0.01, 0.006, 0.034]) function:
            gives the parameters of simulation that best fit the real epi data
    Inputs:
        >> chi2_type: str that tells how to compute chi^2
                Optional, default is 'linear'
                'linear': sum of squares of the differences
                'log': sum of squares of the differences of logs
                'relative': sum of squares of the difference divided by the value of the data
        >> minvalue: int minimum value of the datapoints to fit (minvalue included)
                Optional, default is 10
        >> weight: float how to weight the chi2. 
                The chi2 of the deaths is weighted by weight and the one for UCI with (1-weight)
                Optional, default is 0.5
        >> Ndays: list of ints minimum and maximum days tried for the dayshift
                (optional, default is [35, 50])
        >> percents: list of doubles
            values of percentage of recovered that are deaths. 
            First one is the prediction and then min and max
            Optional, default is [0.01, 0.006, 0.034]
    Outputs:
        << f: dict with keys the percents and values the name of the
            file with results of simu that best fit
        << day_sh: dict with keys the percents and values the number 
            of days to shift the best simulation results
        << mu_U: dict with keys the percents and values the percentage 
            of infected that best fit the UCI cases
    '''
    # read updated epi data HAVE TO IMPLEMENT A DATA UPDATE FUNCTION
    df_epi = pd.read_csv('../data/epi_data_IB.csv')
    df_epi['Día'] = pd.to_datetime(df_epi['Día'], format = '%d/%m/%Y', dayfirst = True)
    # do the fitting of the simulations to the data
    files = glob.glob("../results/time_av_*")
    f = {}
    day_sh = {}
    mu_U = {}
    N = len(files)
    for perc in percents:
        chi_2 = float("inf")
        i = 0
        for fname in files:
            print(perc, N-i)
            i += 1
            chichi, dayshift_, mu_UCI_ = best_fit(df_epi,
                                                  fname_simu = fname,
                                                  mu_D = perc,
                                                  weight = weight,
                                                  chi2_type = chi2_type,
                                                  minvalue = minvalue,
                                                  Ndays = Ndays)
            if chichi < chi_2:
                chi_2 = chichi
                fgood = fname
                dayshift = dayshift_
                mu_UCI = mu_UCI_
        f[perc] = fgood
        day_sh[perc] = dayshift
        mu_U[perc] = mu_UCI
        print(perc, chi_2, fgood, dayshift, mu_UCI)
    return f, day_sh, mu_U

def get_params(fname):
    '''
    params(fname) function:
        gives back the parameters of the simulation that are encoded in the name of the simulation results
    Inputs:
        >> fname: str complete path + filename of simulation results
    Outputs:
        << T_e: int time in state E
        << T_i: int time in state I
        << beta: double infectiousness
    '''
    params = fname.split('.dat')[0].split('_')[2:]
    T_e = int(params[0])
    T_i = int(params[1])
    beta = np.round(float(params[2]),3)
    return T_e, T_i, beta

def plot_deaths(percents,
                f,
                day_sh,
                mu_U,
                chi2_type = 'log',
                minvalue = 10,
                predict = True):
    '''
    plot_deaths(percents,
                f,
                day_sh,
                mu_U,
                chi2_type = 'log',
                minvalue = 10,
                predict = True) function:
            Function to plot the deaths from the data and the adjusted model. 
            The plot is saved on a file '../figures/muertos_day_month_year_chi2type_minvalue.png'
    Inputs:
        >> percents: list of doubles
            values of percentage of recovered that are deaths. 
            First one is the prediction and then min and max
            Optional, default is [0.01, 0.006, 0.034]
        >> f: dict with keys the percents and values the name of the
            file with results of simu that best fit
        >> day_sh: dict with keys the percents and values the number 
            of days to shift the best simulation results
        >> mu_U: dict with keys the percents and values the percentage 
            of infected that best fit the UCI cases
        >> chi2_type: str that tells how to compute chi^2
                Optional, default is 'linear'
                'linear': sum of squares of the differences
                'log': sum of squares of the differences of logs
                'relative': sum of squares of the difference divided by the value of the data
        >> minvalue: int minimum value of the datapoints to fit (minvalue included)
                Optional, default is 10
        >> predict: Boolean tells if we want to plot a prediction of seven days
                True: plot best simu 7 days after last datapoint and last datapoint date
                False: plot best simu from march 15 to 15 days after last datapoint
                Optional, default is True
    Outputs:
        The plot is saved automatically in '../figures/muertos_day_month_year_chi2type_minvalue.png'
        << x_M, y_M: time series of epi data for deaths
        << x2, y2, y2m, y2M: time series with percents[0], percents[1] and percents[2]
    '''
    #datos reales
    df_epi = pd.read_csv('../data/epi_data_IB.csv')
    df_epi['Día'] = pd.to_datetime(df_epi['Día'], format = '%d/%m/%Y', dayfirst = True)
    df2 = df_epi[df_epi['M'].notna()]
    x_M = np.array(df2['Día'])
    y_M = np.array(df2['M'])
    del df2
    #fechas en las que pintar las simulaciones
    if predict:
        start_date = dt.datetime.strptime(str(x_M[-1]), '%Y-%m-%dT%H:%M:%S.%f000')
        end_date = start_date + dt.timedelta(days = 7)
    else:
        start_date = dt.datetime(day = 15, month = 3, year = 2020)
        end_date = dt.datetime.strptime(str(x_M[-1]), '%Y-%m-%dT%H:%M:%S.%f000') + dt.timedelta(days = 15)
    
    #THIS HAS TO BE COMPACTED!! NEED A FUNCTION FOR GETTING THE SERIES TO PLOT
    
    #resultados modelo
        #prediccion con mu_D 1%
    df_simu = read_csv_VME(f[percents[0]])
    df_simu = shift(df_simu, day_sh[percents[0]])
    x2 = list(df_simu['Día'])
    y2 = np.array(df_simu['average_recover'])*percents[0]
    x2, y2 = limit_tseries(x2, y2, start_date, end_date)

        #bandas de predicción

        #min con mu_D 0.6%
    df_simu = read_csv_VME(f[percents[1]])
    df_simu = shift(df_simu, day_sh[percents[1]])
    x2 = list(df_simu['Día'])
    y2m = np.array(df_simu['average_recover'])*percents[1]
    x2, y2m = limit_tseries(x2, y2m, start_date, end_date)

        #max con mu_D 3.4%
    df_simu = read_csv_VME(f[percents[2]])
    df_simu = shift(df_simu, day_sh[percents[2]])
    x2 = list(df_simu['Día'])
    y2M = np.array(df_simu['average_recover'])*percents[2]
    x2, y2M = limit_tseries(x2, y2M, start_date, end_date)
    
    #pintar gráficas de muertes
    plt.figure(figsize = (10, 8))

    plt.plot(x_M, y_M, 'ks', label = 'Muertes oficial', ms = 10)
    plt.plot(x2, y2, 
             c = 'k', 
             label= 'Muertes modelo ('+str(np.round(100*percents[0],2))+'%)',
             lw = 5)
    plt.fill_between(x2, y2m, y2M,
                     color = 'gray', 
                     label= 'Muertes modelo ('+str(np.round(100*percents[1],2))+\
                            '%-'+str(np.round(100*percents[2],2))+'%)',
                     lw = 5,
                     alpha = 0.5)
    ax = plt.gca()
    ax.xaxis.set_tick_params(labelsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    plt.grid()
    plt.yscale('log')
    plt.xlabel('Día', fontsize = 30)
    plt.ylabel('Casos', fontsize = 30)
    ax.yaxis.set_major_formatter(ScalarFormatter())
    plt.xticks(rotation = 45)
    plt.legend(fontsize = 18, loc = 2)
    plt.setp( ax.xaxis.get_majorticklabels(), rotation=45, ha="right" )
    today = dt.datetime.now()
    plt.savefig('../figures/muertes_%02i_%02i_%4i_%s_%i.png' % 
                (today.day, today.month, today.year, chi2_type, minvalue), 
                bbox_inches='tight')
    return x_M, y_M, x2, y2, y2m, y2M

def plot_UCI(percents,
                f,
                day_sh,
                mu_U,
                chi2_type = 'log',
                minvalue = 10,
                predict = True):
    '''
    plot_deaths(percents,
                f,
                day_sh,
                mu_U,
                chi2_type = 'log',
                minvalue = 10,
                predict = True) function:
            Function to plot the UCI cases from the data and the adjusted model. 
            The plot is saved on a file '../figures/UCI_day_month_year_chi2type_minvalue.png'
    Inputs:
        >> percents: list of doubles
            values of percentage of recovered that are deaths. 
            First one is the prediction and then min and max
            Optional, default is [0.01, 0.006, 0.034]
        >> f: dict with keys the percents and values the name of the
            file with results of simu that best fit
        >> day_sh: dict with keys the percents and values the number 
            of days to shift the best simulation results
        >> mu_U: dict with keys the percents and values the percentage 
            of infected that best fit the UCI cases
        >> chi2_type: str that tells how to compute chi^2
                Optional, default is 'linear'
                'linear': sum of squares of the differences
                'log': sum of squares of the differences of logs
                'relative': sum of squares of the difference divided by the value of the data
        >> minvalue: int minimum value of the datapoints to fit (minvalue included)
                Optional, default is 10
        >> predict: Boolean tells if we want to plot a prediction of seven days
                True: plot best simu 7 days after last datapoint and last datapoint date
                False: plot best simu from march 15 to 15 days after last datapoint
                Optional, default is True
    Outputs:
        The plot is saved automatically in '../figures/muertos_day_month_year_chi2type_minvalue.png'
        << x_M, y_M: time series of epi data for deaths
        << x2, y2, y2m, y2M: time series with percents[0], percents[1] and percents[2]
    '''
    #datos reales
    df_epi = pd.read_csv('../data/epi_data_IB.csv')
    df_epi['Día'] = pd.to_datetime(df_epi['Día'], format = '%d/%m/%Y', dayfirst = True)
    df2 = df_epi[df_epi['UCI'].notna()]
    x_M = np.array(df2['Día'])
    y_M = np.array(df2['UCI'])
    del df2
    #fechas en las que pintar las simulaciones
    if predict:
        start_date = dt.datetime.strptime(str(x_M[-1]), '%Y-%m-%dT%H:%M:%S.%f000')
        end_date = start_date + dt.timedelta(days = 7)
    else:
        start_date = dt.datetime(day = 15, month = 3, year = 2020)
        end_date = dt.datetime.strptime(str(x_M[-1]), '%Y-%m-%dT%H:%M:%S.%f000') +\
                   dt.timedelta(days = 15)
    
    #THIS HAS TO BE COMPACTED!! NEED A FUNCTION FOR GETTING THE SERIES TO PLOT
    
    #resultados modelo
        #prediccion con mu_D 1%
    df_simu = read_csv_VME(f[percents[0]])
    df_simu = shift(df_simu, day_sh[percents[0]])
    x2 = list(df_simu['Día'])
    y2 = np.array(df_simu['average_infected'])*mu_U[percents[0]]
    x2, y2 = limit_tseries(x2, y2, start_date, end_date)

        #bandas de predicción

        #min con mu_D 0.6%
    df_simu = read_csv_VME(f[percents[1]])
    df_simu = shift(df_simu, day_sh[percents[1]])
    x2 = list(df_simu['Día'])
    y2m = np.array(df_simu['average_infected'])*mu_U[percents[1]]
    x2, y2m = limit_tseries(x2, y2m, start_date, end_date)

        #max con mu_D 3.4%
    df_simu = read_csv_VME(f[percents[2]])
    df_simu = shift(df_simu, day_sh[percents[2]])
    x2 = list(df_simu['Día'])
    y2M = np.array(df_simu['average_infected'])*mu_U[percents[2]]
    x2, y2M = limit_tseries(x2, y2M, start_date, end_date)
    
    #pintar gráficas de muertes
    plt.figure(figsize = (10, 8))

    plt.plot(x_M, y_M, 'rX', label = 'UCI oficial', ms = 10)
    plt.plot(x2, y2, 
             c = 'r', 
             label= 'UCI modelo',
             lw = 5)
    plt.fill_between(x2, y2m, y2M,
                     color = 'gray', 
                     label= 'Rango UCI',
                     lw = 5,
                     alpha = 0.5)
    ax = plt.gca()
    ax.xaxis.set_tick_params(labelsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    plt.grid()
    plt.yscale('log')
    plt.xlabel('Día', fontsize = 30)
    plt.ylabel('Casos', fontsize = 30)
    ax.yaxis.set_major_formatter(ScalarFormatter())
    plt.xticks(rotation = 45)
    plt.legend(fontsize = 18, loc = 2)
    plt.setp( ax.xaxis.get_majorticklabels(), rotation=45, ha="right" )
    today = dt.datetime.now()
    plt.savefig('../figures/UCI_%02i_%02i_%4i_%s_%i.png' % 
                (today.day, today.month, today.year, chi2_type, minvalue), 
                bbox_inches='tight')
    return x_M, y_M, x2, y2, y2m, y2M

def params_to_df(percents, f, day_sh, mu_U, chi2_type = 'log', minvalue = 10):
    '''
    params_to_df(percents, f, day_sh, mu_U) function:
        Puts all the fitted parameters in a dataframe for easier export and printing.
        It also writes a csv file with the parameters in '../results/fitted_params_day_month_year_chi2type_minvalue.csv'
    Inputs:
        >> percents: list of doubles
            values of percentage of recovered that are deaths. 
            First one is the prediction and then min and max
            Optional, default is [0.01, 0.006, 0.034]
        >> f: dict with keys the percents and values the name of the
            file with results of simu that best fit
        >> day_sh: dict with keys the percents and values the number 
            of days to shift the best simulation results
        >> mu_U: dict with keys the percents and values the percentage 
            of infected that best fit the UCI cases
        >> chi2_type: str that tells how to compute chi^2
                Optional, default is 'linear'
                'linear': sum of squares of the differences
                'log': sum of squares of the differences of logs
                'relative': sum of squares of the difference divided by the value of the data
        >> minvalue: int minimum value of the datapoints to fit (minvalue included)
                Optional, default is 10
    Outputs:
        << df: pandas DataFrame with columns for the parameters adjusted
    '''
    params = dict()
    params['mu_D'] = percents
    params['T_e'] = []
    params['T_i'] = []
    params['beta'] = []
    params['mu_U'] = [mu_U[perc] for perc in percents]
    params['dayshift'] = [day_sh[perc] for perc in percents]
    for perc in percents:
        T_e, T_i, beta = get_params(f[perc])
        params['T_e'].append(T_e)
        params['T_i'].append(T_i)
        params['beta'].append(beta)
    df = pd.DataFrame(data = params)
    today = dt.datetime.now()
    df.to_csv('../results/fitted_params_%02i_%02i_%4i_%s_%i.csv' % (today.day, today.month, today.year, chi2_type, minvalue))
    return df
