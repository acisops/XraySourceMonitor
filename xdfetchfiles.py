import numpy as np
from astropy.time import Time
import os


def update_BAT_list(fluxmin=20):
    import requests
    import lxml.html as lh
    import pandas as pd
    import numpy as np
    from datetime import date
    from astropy.time import Time
    
    url='https://swift.gsfc.nasa.gov/results/transients/BAT_current.html'
    
    #Create a handle, page, to take in the contents of the website
    page = requests.get(url)
    #Store the contents of the website under doc
    doc = lh.fromstring(page.content)
    #Parse data that are stored between <tr>..</tr> of HTML
    tr_elements = doc.xpath('//tr')
    
    #Create empty list
    col=[]
    i=0
    #For each row, store each first element (header) and an empty list
    for t in tr_elements[0]:
        i+=1
        name=t.text_content()
        print (i,name)
        col.append((name,[]))
        
        #i is the index of our column
        i=0
        
        #Since out first row is the header, data is stored on the second row onwards
    for j in range(1,len(tr_elements)):
        #T is our j'th row
        T=tr_elements[j]
        #If row is not of size 13, the //tr data is not from our table 
        if len(T)!=13:
            break
    
        #i is the index of our column
        i=0
    
        #Iterate through each element of the row
        for t in T.iterchildren():
            data=t.text_content() 
            #Check if row is empty
            if i>0:
            #Convert any numerical value to integers
                try:
                    data=int(data)
                except:
                    pass
            #Append the data to the empty list of the i'th column
            col[i][1].append(data)
            #Increment i for the next column
            i+=1

            
    # Define a dataframe to hold the entries
    Dict={title:column for (title,column) in col}
    df=pd.DataFrame(Dict)


    #### Construct an outputtable array for today's BAT transient list

    #Get the data 
    f1 = np.array(df['Today#'])
    f2 = np.array(f1)  #initialize
    f3 = np.array(f1)  #initialize

    # Pull out entries with fluxes 
    for i in range(len(f1)):
        f=f1[i]
        n=f.find('(')-1
        if n > 0:
            f3[i]='True'
            f2[i]=float(f[0:n])     #Flux lives here
        else:
            f3[i]='False'


    # Screen to retain only those at F > fluxmin (definable)
    for i in range(len(f3)):
        if f3[i] == 'True':
            if float(f2[i]) < fluxmin:    
                f3[i] = 'False'
    
    FF = f2[f3 == 'True']   #   the selected flux list
    src = np.array(df['Source Name'])[f3 == 'True']  #selected source list etc...
    RA  = (np.array(df['RA J2000 Degs'])[f3 == 'True'])
    DEC = (np.array(df['Dec J2000 Degs'])[f3 == 'True'])
    AltName = np.array(df['Alternate Name'])[f3 == 'True']
    stype = np.array(df['Source Type'])[f3 == 'True']

    
    for i in range(len(RA)):
        RA[i] = float(RA[i])
        DEC[i] = float(DEC[i])


    # Cleaning up the name formatting
    for i in range(len(src)):
        src[i]=src[i].replace('        \xa0 ','')
        src[i]=src[i].replace('  ','')
        if src[i][-1] == ' ':
            src[i]=src[i][0:len(src[i])-1]

        src[i]=src[i].replace(' ','_')
    
    for i in range(len(src)):
        AltName[i]=AltName[i].replace('  ','')
    
        if AltName[i] != '':
            if AltName[i][-1] == ' ':
                AltName[i]=AltName[i][0:len(AltName[i])-1]

        AltName[i]=AltName[i].replace(' ','_')


    for i in range(len(src)):
        stype[i]=stype[i].replace('        \xa0 ','')
        stype[i]=stype[i].replace('  ','')
        if stype[i][-1] == ' ':
            stype[i]=stype[i][0:len(stype[i])-1]

            stype[i]=stype[i].replace(' ','_')


    #    Getting today's date in mjd to label the table name 
    dstamp=str(int(Time.now().mjd))
    fname='BAT_record_'+dstamp+'.dat'

    #   Write out today's table
    fil=open(fname,"w+")
    fil.write('%15a    %15a    %10a    %30a    %30a   %20a \n' % ('RA',  'DEC',  'Flux(Crab)',  'Source', 'AltName', 'Type'))

    for i in range(len(src)):
        fil.write('%15.3f    %15.3f    %10.4f    %30a    %30a   %20a \n' % (RA[i],  DEC[i], FF[i]/1000., src[i], AltName[i], stype[i]))

    fil.close()
    


def update_maxi_list():   # INGEST NEW MAXI TOPFLUX TABLES
    import numpy as np
    import subprocess
    from astropy.time import Time
    import glob

    # What is today's MJD?    
    tnow=Time.now()
    mjdnow=int(tnow.mjd)

    # get and sort the present list of fluxtop files
    flist=glob.glob('fluxtop*dat')    
    fsrt = flist.sort()    

    # ID the latest file    
    maxdate=int(flist[-1].replace('fluxtop','').replace('.dat','') )

    #######  Potential troubleshooting here:
    ##  print(['wget fluxtop'+str(foo)+'.dat' for foo in np.arange(maxdate,mjdnow)])
    ##  print('wget fluxtop'+np.arange(maxdate,mjdnow)+'.dat')
    #############################################

    ## start a shell process to grab any newer fluxtop entries
    for el in np.arange(maxdate,mjdnow+1):
        print('fetching fluxtop'+str(el)+'.dat')
        subprocess.Popen('wget http://maxi.riken.jp/fluxtop/fluxtop'+str(el)+'.dat',shell=True)

