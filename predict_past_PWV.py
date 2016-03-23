"""
Goal: To predict past PWV measurements on KPNO based on available data 
from suominet at the sites KITT (KPNO), the nearby desert floor SA48, and 
other nearby sites SA46, P014, and AZAM.
"""
import os.path
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import astropy
from astropy.time import TimeISOT

mypath = os.path.join('c:', os.sep, 'sourcedir')
#gather all files
file_kitt = "C:/Users/Jessica/Desktop/Atmosphere/KITTnrt_20150512.plot.txt"
file_sa48 = "C:/Users/Jessica/Desktop/Atmosphere/SA48nrt_20150701.plot.txt"
file_sa46 = "C:/Users/Jessica/Desktop/Atmosphere/SA46nrt_20150701.plot.txt"
file_p014 = "C:/Users/Jessica/Desktop/Atmosphere/P014nrt_20150701.plot.txt"
file_azam = "C:/Users/Jessica/Desktop/Atmosphere/AZAMnrt_20150701.plot.txt"
file_sa48_2014 = "C:/Users/Jessica/Desktop/Atmosphere/SA48nrt_2014.plot"
file_sa46_2014 = "C:/Users/Jessica/Desktop/Atmosphere/SA46nrt_2014.plot"
file_p014_2014 = "C:/Users/Jessica/Desktop/Atmosphere/P014nrt_2014.plot"
#file_azam_2014 = "C:/Users/Jessica/Desktop/Atmosphere/AZAMnrt_2014.plot"
#no data prior to 2015 available from AZAM
file_sa48_2013 = "C:/Users/Jessica/Desktop/Atmosphere/SA48nrt_2013.plot"
file_sa46_2013 = "C:/Users/Jessica/Desktop/Atmosphere/SA46nrt_2013.plot"
file_p014_2013 = "C:/Users/Jessica/Desktop/Atmosphere/P014nrt_2013.plot"
file_sa48_2012 = "C:/Users/Jessica/Desktop/Atmosphere/SA48nrt_2012.plot"
file_sa46_2012 = "C:/Users/Jessica/Desktop/Atmosphere/SA46nrt_2012.plot"
file_p014_2012 = "C:/Users/Jessica/Desktop/Atmosphere/P014nrt_2012.plot"
file_sa48_2011 = "C:/Users/Jessica/Desktop/Atmosphere/SA48nrt_2011.plot"
file_sa46_2011 = "C:/Users/Jessica/Desktop/Atmosphere/SA46nrt_2011.plot"
file_p014_2011 = "C:/Users/Jessica/Desktop/Atmosphere/P014nrt_2011.plot"
file_sa48_2010 = "C:/Users/Jessica/Desktop/Atmosphere/SA48nrt_2010.plot"
file_sa46_2010 = "C:/Users/Jessica/Desktop/Atmosphere/SA46nrt_2010.plot"
file_p014_2010 = "C:/Users/Jessica/Desktop/Atmosphere/P014nrt_2010.plot"

#data_2015 = (file_kitt, file_sa48, file_sa46, file_p014, file_azam)
data_2014 = (file_sa48_2014, file_sa46_2014, file_p014_2014)
data_2013 = (file_sa48_2013, file_sa46_2013, file_p014_2013)
#data_2012 = (file_sa48_2012, file_sa46_2012, file_p014_2012)
#data_2011 = (file_sa48_2011, file_sa46_2011, file_p014_2011)
#data_2010 = (file_sa48_2010, file_sa46_2010, file_p014_2010)

def get_data(array_files):
    """
    Use np.genfromtxt to retrieve data from an array of files
    
    expects files generated from suominet 
    
    input = array of files, generally one from each location
    output = list of arrays with noted columns 
    """
    orig_data = np.array(array_files)
    data = []
    for x in orig_data:
        data.append(np.genfromtxt(x, usecols=(1,2,7,8,9), names=('date', 'pwv', 'pres', 'temp', 'hum'), dtype=((np.str_, 16), float, float, float, float)))
    return data

def get_date_array(data):
    """
    Construct a sorted list of unique dates from an input list of 
    a data table for each site.
    
    Return date & times expressed in MJD
    """
    
    from astropy.time import Time
    
    datetimes = np.array([site_data['date'] for site_data in data])

    unique_datetimes = []
    for site_datetimes in datetimes:
        site_unique = np.unique(site_datetimes)
        unique_datetimes = np.append(unique_datetimes, site_unique)
    
    unique_datetimes_final = np.unique(unique_datetimes)
    mjd = [Time(t, format = 'isot').mjd for t in unique_datetimes_final]
    
    return sorted(mjd)

def test_get_date_array():
    sample_p014 = get_data(data_2014[2])[0:10]
    return sample_p014


 
def pad_arrays(data):
    """
    Set up arrays, pad values so that all are matching in length 
    and time stamps
    
    data -- array of data for sites
    
    returns array of 0 = not masked, 1 = masked 
    for each possible time at each data site
    """ 
    from astropy.time import Time 
    dates = get_date_array(data)
    
    full_mask = []
    #change for-loops to list comprehensions!
    for site in range(len(data)):
        times = []
        mask = [] 
     
        for time in range(len(data[site])):
            times.append(Time(data[site][time][0], format = 'isot').mjd)
           
        for x in range(len(dates)):
            if dates[x] in times: 
                mask.append(0)
            else:
                mask.append(1)
        full_mask.append([mask])
        #np.concatenate((full_mask[dim], [mask]), axis=0)
        
   
    return full_mask
        
    
    

data = get_data(data_2014)
#print data
#dates = get_date_array(data)

#print data[0]['pwv']
#full_mask = pad_arrays(data)
#print full_mask[1]
#print sum(full_mask)

def mask_bad_vals(data):
    """
    Build on input data mask, by 
    masking times where pwv < 0
    
    Returns updated mask
    """
    times = []
    mask = pad_arrays(data)
    new_mask = []
  
    for site in range(len(data)):
        
        to_mask = mask[site] 
    #want to iterate through all sites bc marking each as good/bad
        #rather than np.intersect1d~
        for time in range(len(mask[site])):
            if (data[site][time]['pwv'] < 0) or (data[site][time]['pres'] == 767.0):
                to_mask[time] = 1
        #np.append(new_mask[site], to_mask)
        #mask[site] = to_mask
        new_mask.append([mask])
    #print new_mask
    return new_mask
    
def final_data(data):
    """
    Builds arrays, uses masking with Astropy table
    to return final data set
    """
    from astropy.table import Table, Column, Row
    
    dates = Column(data=[get_date_array(data)], name='dates')
    mask = mask_bad_vals(data)
    
    t = Table(dates)
    for n in range(len(data)):
        t.add_column(data[n], mask[n])
    #t.mask(mask[0])
    print t

def fit_functions(data):
    """
    Creates a fitting function for each site provided
    
    data -- array of data for sites

    """
    fit_functions = []
    for site in range(1, len(data)):
        fit = (np.polyfit(data[site]['pwv'], data[0]['pwv'], 1))
        fit_functions.append(np.poly1d(fit))
    return fit_functions

#fit = fit_functions(data)
#print fit

    

#pos = pos_dates(data)

final = final_data(data)

def find_dates(data, n):
    """
    Find indices where the recorded PWV is > 0
    for multiple sets of data
    
    data -- padded array of data for sites 
    n    -- number of sites in array
    """
    positive_indices = []
    """
    mask = (data[0]['pwv']>0) & (data[1]['pwv'] > 0)
    data[0]['date'][mask]
    """
    
    good_dates = data[0]['date']
        
    
    for x in data:
        pos_indices = np.where(x['pwv'] > 0)
        positive_indices.append(pos_indices)
    good_dates = np.intersect1d(data[0]['date'][positive_indices[0]], data[1]['date'][positive_indices[1]])
    
    if n <= 2:
        return good_dates

    else:
        final_dates = []
        run = 2
        for x in range(2, n):
            pos_indices = np.where(data[run]['pwv'] > 0)
            positive_indices = []
            positive_indices.append(pos_indices)
            final_dates.append(np.intersect1d(good_dates, data[run]['date'][positive_indices]))
            run +=1
        return final_dates

    
def find_matching_data(data, good_dates, num_sites):
    """
    Find data at matching indices
    Includes masking of certain values
    
    data       -- array of data for sites
    good_dates -- matching dates (array)
    num_sites  -- number of sites in array 
    """
    """
    date_match = data[0]['date']
    for site in range(0, num_sites):
        indices = []
        indices.append(np.searchsorted(data[site]['date'], date_match))
        date_match.append(data[site][indices])
        date_match = data['site'][indices]
        
    go through again, make data match out of final values
    instead of good dates, use date match from last loop 
    
    """
    data_match = []
    data_final = []
    for site in range(0, num_sites):
        indices = []
        indices.append(np.searchsorted(data[site]['date'], good_dates))
        data_match.append(data[site][indices])
    if num_sites >=4:
        masked_pres = ma.masked_equal(data[0]['pres'], 767.0)
        """
        for a given site, step through every value. 
        mask = np.bool_(np.ones(data.shape)) #creates arrays of all trues
        #new??? mask = np.bool_(np.ones(data[x]['date'].shape))
        #big_mask.append(mask) (maybe in else loop???) match w corresponding mask at end 
        ##include mask before you find dates in common!! (we want to get rid of dates)
        for x in range(0, n):
            repeat_num = 0
            for ind, p in enumerate(data[x]['pres']):
                if p == data[x]['pres'][ind - 1]:
                    repeat_num += 1
                else:
                    if repeat_num > 3:
                        mask[x][:][ind-repeat_num:ind] = False ##CHECK THIS!! might not work 
                    
        
        """
        masked_indices = ma.nonzero(masked_pres)
        for site in range(0, num_sites):
            data_final.append(data[site][masked_indices])
            return data_final
    else:
        return data_match

def create_fit_functions(data):
    """
    Creates a fitting function for each site provided
    
    data -- array of data for sites

    """
    fit_functions = []
    for x in range(1, len(data)):
        fit = (np.polyfit(data[x]['pwv'], data[0]['pwv'], 1))
        fit_functions.append(np.poly1d(fit))
    return fit_functions

#still working here
def find_average(data):
    f_avg = []
    for x in range(len(data[0])):
        #data_used = [f_sa48[x],f_sa46[x] , f_p014[x] , f_azam[x] ]
        v = np.median(data)
        s = np.std(data)
        data_new = []
        data_new = [x for x in data_used if (v-(2 * s)) < x < (v+(2 * s))]
        f_avg.append(sum(data_new)/len(data_new))
    return f_avg
   
"""
data__2014 = get_data(data_2014)
good_dates_2014 = find_dates(data__2014, 3)
final_dates_2014 = find_matching_data(data__2014, good_dates_2014, 3)
fit_fns_2014 = create_fit_functions(final_dates_2014, 3)
print fit_fns_2014

"""
