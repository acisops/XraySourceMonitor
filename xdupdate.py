import numpy as np

# Routine to send email alerts for common but notable occurences, like a new source popping...
# Aim to make some more common cases like new maxima above (~Crab?) or similar
def _send_email_notice(type='NEW',name='None',sendto='JFS',RA='None',DEC='None',home=False,message='None'):
    import smtplib                                                                                                                                                               
    from email.message import EmailMessage                                                                                                                                       
    msg = EmailMessage() 
    mtxt=''  #initialing null
    if type=='NEW':
        _message=''
        if message != 'None':
            _message = '\n'+message
        msg.set_content('New Source Found! '+str(name)+' at RA='+str(RA)+' & DEC='+str(DEC)+_message)
        mtxt = 'Source'
    else:
        msg.set_content(message)
        mtxt = ''
    fromaddr = 'xraysource.automaton@cfa.harvard.edu'
    if sendto=='JFS':
        toaddr   = 'jsteiner@cfa.harvard.edu'
    else:
        toaddr   = 'acisdude@cfa.harvard.edu'
    #msg['Subject'] = f('X-ray SourceList Alert: '+type+' '+mtxt)
    msg['Subject'] = ('[X-ray Source Alert] '+type+' '+mtxt)
    msg['From'] = fromaddr
    msg['To'] = toaddr
    if home==False:
        s = smtplib.SMTP('localhost') 
        s.send_message(msg) 
        s.quit()
    else:
        print(msg['Subject'])

        
# Converts binary string to text string
def bintostr(inlist):
    return([str(b).replace("b",'').replace("'",'') for b in inlist])
    

# The major routine that updates the master sourcelist, to be run after updating MAXI and BAT daily lists
def update_master_record(masterfile='master_sourcelist.tab',remapfile='override.list', excludelist='blacklist_sources.list',initializing=False,limitdown=30, limitup=80, recentfile='mjd.latest',RecentTimeLimit=100,angthreshold=30,FlickerTimeLimit=50,home=False,mirror=True,maxilist='maxilist.pysav',batlist='batlist.pysav',asmlist='asmlist.pysav'):

    # Keeping the below list handy for debugging run-throughs:
    #masterfile='master_sourcelist.tab';remapfile='override.list'; excludelist='blacklist_sources.list';initializing=False;limitdown=30; limitup=80; recentfile='mjd.latest';RecentTimeLimit=100;angthreshold=30;FlickerTimeLimit=50;mirror=True;home=False;maxilist='maxilist.pysav';batlist='batlist.pysav';asmlist='asmlist.pysav'

    import glob
    from astropy.time import Time
    import os
    import pickle
    import xdreader as mr
    import copy
    import astropy.io.ascii as atab
    import fnmatch
    import string
    from itertools import compress 
    from astropy.table import Table
    
    today=str(int(Time.now().mjd))   # get today in MJD
    today_ = np.float(today)
    
    pfile = open(maxilist, 'rb')     # restore the MAXI data pickle as mres     
    mres= pickle.load(pfile)
    pfile.close()

    pfile = open(batlist, 'rb')      # restore the BAT data pickle as sres     
    sres= pickle.load(pfile)
    pfile.close()


    mtab = atab.read(masterfile)     # read in the master table
    bdata = atab.read(excludelist)   # read in the expressions for names to be excluded
    rdata = atab.read(remapfile)     # read in list of any fixes to be hand supplied
    

    mres_maxiIDstr = np.array(bintostr(mres['maxiID']))  # list of all stored MAXI detections, by MaxiID
    batsourcelist=np.array(bintostr(sres['batID']))      # list of all stored BAT detections, by BATID 


    # SEARCH FOR ALL SOURCES IN BATLIST ABOVE upperlimit and ALL sources above lowerlimit
    recentdata = atab.read(recentfile)    
    MAXIrecentlim = recentdata['MAXILatest'][0]   # set the last time update was run, so will pull from dates past this point...
    BATrecentlim = recentdata['BATLatest'][0]    
   
    lnew_maxi = (np.where(mres['mjd'] >= MAXIrecentlim))[0] #####  indexes of just new dates
    lnew_bat  = (np.where(sres['mjd'] >= BATrecentlim))[0] 


    ###########################   SECTION BELOW SETS FLICKER CONDITION, SOURCES
    ### Activity = 1 --> 0  sources should be checked against these lists to make sure they aren't simply flickering high/low and are steadily below the flux limitdown threshold
    lflick_maxi = (np.where(mres['mjd'] >= (today_-FlickerTimeLimit)))[0] # for sources with noise
    lflick_bat  = (np.where(sres['mjd'] >= (today_-FlickerTimeLimit)))[0]
    
    mrec = mres[lnew_maxi] # recent data only
    srec = sres[lnew_bat]

    mflk1 = mres[lflick_maxi]  # isolating data inside time-window in which to check for flicker activity
    sflk1 = sres[lflick_bat]  
    
    mflk = mflk1[mflk1['flux'] > limitdown]    # any source which has "flickering" at this level stays active under def.
    sflk = sflk1[sflk1['flux'] > limitdown/1000.]  

    uniq_mflicker = list(set(mflk['maxiID']))  # holds list of unique sources that match "flickering" criterion (i.e., recent activity above limitdown)
    uniq_sflicker = list(set(sflk['batID']))
    ##############################3

    
    
    maxihi  = mrec[mrec['flux'] >= limitup]      # Any source in this camp is significant *and* active
    maxilow1 = mrec[mrec['flux'] >= limitdown]   
    maxilow = maxilow1[maxilow1['flux'] < limitup]  # Sources in this camp are only active if they were previously active 
    maxifaint = mrec[mrec['flux'] < limitdown]      # Sources here are below threshold, so inactive unless matching flicker condition

    bathi  = srec[srec['flux'] >= limitup/1000.]
    batlow1 = srec[srec['flux'] >= limitdown/1000.]
    batlow = batlow1[batlow1['flux'] < limitup/1000.]
    batfaint = srec[srec['flux'] < limitdown/1000.]

    uniq_maxihi = list(set(maxihi['maxiID']))     # List of unique sources in maxihi
    uniq_maxilow1 = list(set(maxilow['maxiID']))  
    ulowfilt = [q not in uniq_maxihi for q in uniq_maxilow1]
    uniq_maxilow = list(compress(uniq_maxilow1,ulowfilt))   # List of unique sources in maxilow and *not* in maxihi
        
    uniq_bathi  = list(set(bathi['batID']))   
    uniq_batlow1 = list(set(batlow['batID']))
    ulowfilt = [q not in uniq_bathi for q in uniq_batlow1]
    uniq_batlow = list(compress(uniq_batlow1, ulowfilt))


    uniq_maxifaint1 = list(set(maxifaint['maxiID']))
    ufaintfilt = [q not in uniq_maxihi and q not in uniq_maxilow for q in uniq_maxifaint1]  
    uniq_maxifaint = list(compress(uniq_maxifaint1,ufaintfilt))    # List of unique sources in maxifaint

    _mf_objname = np.array(np.repeat('',len(uniq_maxifaint)),dtype=object)  # initializing to hold maxifaint source name
    _mf_ra      = np.zeros(len(uniq_maxifaint))  # initializing, holds maxifaint RA
    _mf_dec     = np.zeros(len(uniq_maxifaint))  # initializing, holds maxifaint DEC
    
    
    uniq_batfaint1 = list(set(batfaint['batID']))
    ufaintfilt = [q not in uniq_bathi and q not in uniq_batlow for q in uniq_batfaint1]
    uniq_batfaint = list(compress(uniq_batfaint1,ufaintfilt ))

    

    newIDs_maxi = np.array(bintostr(mres[lnew_maxi]['maxiID']))    # *recent* observations only
    newIDs_bat  = np.array(bintostr(sres[lnew_bat]['batID']))
                           
    updated = np.zeros(len(mtab))

    # Define a template for a new table entry 
    blankentry = copy.deepcopy(mtab[0])
    blankentry['OrigObjName'] =  'NONE'
    blankentry['Object']      =  'NONE'
    blankentry['MaxiID']      =  'NONE'
    blankentry['RA']          =  '-99'
    blankentry['DEC']         =  '-99'
    blankentry['Fpeak(Crab)'] =  -1
    blankentry['MJDpeak']     =  -99
    blankentry['RecordInst']  =  'NONE'
    blankentry['Frecent[MAXI]']   = -99
    blankentry['MJDrecent[MAXI]'] = -99
    blankentry['Frecent[BAT]']    = -99
    blankentry['MJDrecent[BAT]']  = -99
    blankentry['Category']        = -99
    blankentry['Activity']        = -99
    blankentry['SigCount']        = -99


    #for each HI source,
    #   (1) Screen against excludelist
    #   (2) Check survivers against master table's sources
    #################################################

    
    lnew = lnew_maxi ## running through maxi first
    for yy in uniq_maxihi:
        xx = str(yy).replace("b'","").replace("'","")  # gets rid of byte formatting
        obj=mr.get_maxi_objname(xx)  # try to get a readable name
        obj0=obj
        tossout=0  # flag whether or not to exclude


        if xx[0]!='J':    # screens for bad maxiID
            tossout=1    
        
        plnewsource = (np.where(mres_maxiIDstr == xx))[0]  # all locations in case source isn't on the table and was faint previously  ** has to have index 0 for some reason, to access the list
        pl = (np.where(newIDs_maxi == xx))[0]   # locations of only *recent* detections  ** has to have index 0 for some python-reason to access the list
            
        mxf   = max(mres['flux'][lnew][pl])/1000.   ### ONLY searching for max in new date range to prevent junk detections from early interceding...
        mxloc = np.argmax(mres['flux'][lnew][pl])
        mxmjd = mres['mjd'][lnew][pl][mxloc]
                                  
        for el in bdata['MATCHTXT']:
            if fnmatch.fnmatch(obj.lower(),(el.lower())):
                tossout = 1
                print('THROWING OUT '+obj+' ::  MATCH to '+el)
        if tossout==0:  # then can continue - not a string expression match for an extended (or bad) source-name 
            obj0=obj
            if obj == '404_Not_Found':
                obj='Unknown'
                obj0=obj
            
            pageworks = mr.test_simbad_page(obj,mirror=mirror)  # test whether simbad page exists
        
        
            if pageworks != 1:
                # try revising the name in some usual ways
                obj=obj.replace('_and_SNR','').replace('_SNR','').replace('Galactic_Center_Region','SgrA*').replace('_Slow-Burster_with_Rapid-Burster','')
                obj=obj.replace('NGC_6814_with_V1432_Aql','V1432_Aql').replace('_Rapid-Burster_with_Slow-Burster','').replace('LS_I_+61_303_with_Swift_J0243.6+6124','Swift_J0243.6+6124')
                obj=obj.replace('_with_GS_0836-429','').replace('_with_Terzan_1','')
                pageworks = mr.test_simbad_page(obj,mirror=mirror)
                if pageworks != 1: 
                    if obj[0:2]=='A_':
                        obj='1'+obj
                        if obj[0:2]=='X_':
                            obj=('foo'+obj).replace('fooX_','3A_')                       
                            pageworks = mr.test_simbad_page(obj,mirror=mirror)

            if obj == '404_Not_Found':  
                obj='Unknown'
                obj0=obj
                coords=['INDEF','INDEF']
            else:
                if pageworks != 1:  # still no joy on the name, so use maxiID to hold coordinates, intended as temporary measure until source is refined
                    if xx[0]=='J':
                        coords=(xx+' '+xx).split()
                    else:
                        coords=['INDEF','INDEF']
                    print('NOOOOOOOOOOO!!!!')
                    _send_email_notice(type='BadEntry',name=obj,sendto='JFS',home=home,message='No location for MAXI source '+xx+' with flux '+np.str(mxf)+'...'+coords[0]+','+coords[1])
                else:
                    coords=mr.get_coordinates_icrs(mr.coord_from_simbad(obj,mirror=mirror)).split()


            ra=coords[0]
            dec=coords[1]

                                  
            frecent=mres['flux'][lnew][pl][-1]/1000.   # most recent flux
            drecent=mres['mjd'][lnew][pl][-1]          # date of most recent detection


            # DEFINE ACTIVITY AND CATEGORY TO BE UPDATED ONLY IN RIGHT CIRCUMSTANCES
            activity=1  #at some point this source was over limitup by selection
            if frecent < limitdown/1000.:
                if (yy not in uniq_mflicker)==True:
                    activity=0   # if no flickering for this *primary* maxiID, can indicate turning off... more checks done later in script to confirm

            if mxf > 10:
                category=5
            else:
                if mxf > 3:
                    category = 4
                else:
                    if mxf > 1:
                        category=3
                    else:
                        if mxf > .3:
                            category=2
                        else:
                            if mxf > .1:
                                category=1
                            else:
                                category=0

            # matches_ holds the index of a matching table entry
            matches_ = []  #initialize
            if xx != 'INDEF':
                matches_ = (np.where(mtab['MaxiID'] == xx))[0]  # first priority - match the maxiID

            if len(matches_) == 0:
                matches_ = (np.where(mtab['Object'] == obj))[0]  # next priority - match the objname

            addaltname=False  # initializing
            if len(matches_) == 0: # final check: search by coordinate
                mfound=0   # match found t/f
                scomplete=0 # search complete t/f
                incr=0
                while (mfound==0 and scomplete==0):   
                    noms = np.str(mtab[incr]['AltNames']).split(',')
                    names = [_.strip() for _ in noms]
                    #matches_ = np.where(names == obj)   ### DISABLING THIS TO AVOID POTENTIAL NAMING MISHAPS... STICKING WITH ANGLE WHICH IS ROBUST UNLESS WANT TO ENABLE THIS AT SOME POINT
                    _RA  = mtab[incr]['RA']
                    _DEC = mtab[incr]['DEC']
                    if len(matches_) > 0:
                        mfound=1
                        matches_ = [incr]
                    if incr == (len(mtab)-1):
                        scomplete=1
                    if mfound==0:
                        # try matching on ra,dec
                        if str(_RA)[0] != 'J' and str(ra)[0] != 'J':    # good coordinates exist for both
                            adif=mr.get_angsep_arcsec(ra,dec,np.float(_RA),np.float(_DEC))
                            if adif < angthreshold:  # then found a match
                                mfound=1 
                                matches_ = [incr]
                                if names[0]=='none' or names[0]=='--' or names[0]=='0':    # screening for nulls
                                    anames = xx
                                    addaltname = True    # update the altnames to include this maxiID
                                else:
                                    if (xx in names) == False:
                                        anames = np.str(names).replace("['","").replace("']","").replace(" ",",").replace("'","")+','+xx
                                        addaltname = True   # update altnames to include this maxiID
                    incr+=1
            if len(matches_) > 0:  # WHEN A MATCH IS FOUND:
                # check brightness, and update source properties
                i=matches_[0]
                mtab[i]['SigCount']+=len(pl)
                mtab[i]['Frecent[MAXI]']=frecent
                mtab[i]['MJDrecent[MAXI]']=drecent
                mtab[i]['Activity']==1  # This update is because *has* been found above hi-threshold at some point in this section
                if (mtab[i]['MJDrecent[BAT]'] < drecent):   #Only if most recent data does this offer a chance to turn it off
                    if activity==0:
                        if obj.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:  # check if name matches a flickering source in the BAT
                            activity=1
                        else:
                            # doublecheck no flickering under a different name
                            noms = np.str(mtab[i]['AltNames']).split(',')
                            names = [_.strip() for _ in noms]
                            for nn in names:
                                if nn[0]=='J':  # likely a MAXI-ID
                                    if nn in bintostr(uniq_mflicker):   # check if matches a flickering ID
                                        activity=1
                                if nn.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:   # check if matches a flickering ID
                                    activity=1                                    
                    mtab[i]['Activity']=activity
                    
                if mxf > mtab[i]['Fpeak(Crab)']:   # if new peak identified, update peak parameters  
                    mtab[i]['Fpeak(Crab)']=mxf
                    mtab[i]['MJDpeak']=mxmjd
                    mtab[i]['RecordInst']='MAXI'
                    mtab[i]['Category']=category
                if addaltname == True:
                    mtab[i]['AltNames']=anames
                updated[i] = 1  # all set
            else:
                # MAKE A NEW ENTRY!
                print('MAKING A NEW ENTRY: name='+obj+' RA = '+np.str(ra)+' & DEC = '+np.str(dec))
                _send_email_notice(type='NEW',name=obj,sendto='JFS',RA=ra,DEC=dec,home=home,message='MAXI detection, peak flux='+np.str(mxf)+' at '+np.str(mxmjd))

                newentry = copy.deepcopy(blankentry)
                                  
                newentry['OrigObjName'] =  obj0
                newentry['Object']      =  obj
                newentry['MaxiID']      =  xx
                newentry['RA']          =  ra
                newentry['DEC']         =  dec
                newentry['Fpeak(Crab)'] =  mxf
                newentry['MJDpeak']     =  mxmjd
                newentry['RecordInst']  =  'MAXI'
                newentry['Frecent[MAXI]']   = frecent
                newentry['MJDrecent[MAXI]'] = drecent
                newentry['Frecent[BAT]']    = -9.999
                newentry['MJDrecent[BAT]']  = -99
                newentry['Category']        = category
                if activity==0:  # odd case - source was active semi-recently but now is off... make a check
                    if obj.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:  # check if name matches a flickering source in the BAT
                            activity=1
                newentry['Activity']        = activity
                newentry['SigCount']        = len(plnewsource)

                mtab=np.append(mtab,newentry)
                updated = np.append(updated,-1)  # flags that update happened but need to search BAT still

                    
    lnew = lnew_bat ## running through BAT now
    for yy in uniq_bathi:
        xx = str(yy).replace("b'","").replace("'","")
        obj=xx
        obj0=obj
        tossout=0

        plnewsource = (np.where(batsourcelist == xx))[0]  # all detections in case source was previously found but not bright enough to register
        pl = (np.where(newIDs_bat == xx))[0]   # only recent detections

                                              
        for el in bdata['MATCHTXT']:
            if fnmatch.fnmatch(obj.lower(),(el.lower())):
                tossout = 1
                print('THROWING OUT '+obj+' ::  MATCH to '+el)
        if tossout==0:  # can continue - not a string match for a bad object            
            ra=sres[lnew][pl][0]['ra']
            dec=sres[lnew][pl][0]['dec']

            mxf   = max(sres['flux'][plnewsource])   # max flux and friends in this block
            mxloc = np.argmax(sres['flux'][plnewsource])
            mxmjd = sres['mjd'][plnewsource][mxloc]
                                  
            frecent=sres['flux'][lnew][pl][-1]
            drecent=sres['mjd'][lnew][pl][-1]


            # DEFINE ACTIVITY AND CATEGORY TO BE UPDATED ONLY IN RIGHT CIRCUMSTANCES
            activity=1  #at some point this source was over limitup by selection
            if frecent < limitdown/1000.:
                if (yy not in uniq_sflicker)==True:
                    activity=0   # if no flickering for this batID, can indicate turning off... more checks done later in script to confirm

            if mxf > 10:
                    category=5
            else:
                if mxf > 3:
                    category = 4
                else:
                    if mxf > 1:
                        category=3
                    else:
                        if mxf > .3:
                            category=2
                        else:
                            if mxf > .1:
                                category=1
                            else:
                                category=0

                                  
            matches_ = (np.where(mtab['Object'] == obj))[0]
            
            addaltname=False  # initializing 
            if len(matches_) == 0:
                mfound=0  # match found
                scomplete=0 # search complete
                incr=0
                while (mfound==0 and scomplete==0):
                    noms = np.str(mtab[incr]['AltNames']).split(',')
                    names = [_.strip() for _ in noms]
                    #matches_ = np.where(names == obj)   ### DISABLING THIS TO AVOID POTENTIAL NAMING MISHAPS... STICKING WITH ANGLE WHICH IS ROBUST UNLESS WANT TO ENABLE THIS AT SOME POINT
                    _RA  = mtab[incr]['RA']
                    _DEC = mtab[incr]['DEC']
                    if len(matches_) > 0:
                        mfound=1
                        matches_ = [incr]
                    if incr == (len(mtab)-1):
                        scomplete=1
                    if mfound==0:
                        #try matching on ra,dec
                        if str(_RA)[0] != 'J':
                            adif=mr.get_angsep_arcsec(ra,dec,np.float(_RA),np.float(_DEC))
                            if adif < angthreshold:  # found a match in the master table
                                mfound=1
                                matches_ = [incr]
                                if names[0]=='none' or names[0]=='--' or names[0]=='0':  # checking if altnames is a null entry
                                    anames = xx
                                    addaltname = True    
                                else:
                                    if (xx in names) == False:
                                        anames = np.str(names).replace("['","").replace("']","").replace(" ",",").replace("'","")+','+xx
                                        addaltname = True

                    incr+=1
            if len(matches_) > 0:
                # check brightness, and update elements
                i=matches_[0]
                mtab[i]['SigCount']+=len(pl)
                mtab[i]['Frecent[BAT]']=np.round(frecent*1000)/1000 # otherwise annoying floating-precision noise keeps showing up
                mtab[i]['MJDrecent[BAT]']=drecent

                mtab[i]['Activity']==1  # This update is because *has* been found above hi-threshold at some point in this section
                if (mtab[i]['MJDrecent[MAXI]'] < drecent):   #Only if most recent data does this offer a chance to turn it off
                    if activity==0:
                        if mtab[i]['MaxiID'] != 'INDEF':                        
                            if mtab[i]['MaxiID'] in bintostr(uniq_mflicker):  # check if name matches a flickering source in MAXI
                                activity=1
                        else:
                            # doublecheck no flickering under a different name
                            noms = np.str(mtab[i]['AltNames']).split(',')
                            names = [_.strip() for _ in noms]
                            for nn in names:
                                if nn[0]=='J':  # likely a MAXI-ID
                                    if nn in bintostr(uniq_mflicker):   # check if matches a flickering ID
                                        activity=1
                                if nn.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:   # check if matches a flickering ID
                                    activity=1                                    
                    mtab[i]['Activity']=activity
                    

                if mxf > mtab[i]['Fpeak(Crab)']:   # if new peak identified, update peak parameters  
                    mtab[i]['Fpeak(Crab)']=mxf
                    mtab[i]['MJDpeak']=mxmjd
                    mtab[i]['RecordInst']='BAT'
                    mtab[i]['Category']=category

                if addaltname == True:
                    mtab[i]['AltNames']=anames

                updated[i] = 1  # all set
            else:
                # MAKE A NEW ENTRY!

                _send_email_notice(type='NEW',name=obj,sendto='JFS',RA=ra,DEC=dec,home=home,message='BAT detection, peak flux='+np.str(mxf)+' at '+np.str(mxmjd))

                newentry = copy.deepcopy(blankentry)
                                  
                newentry['OrigObjName'] =  obj0
                newentry['Object']      =  obj
                newentry['MaxiID']      =  'INDEF'
                newentry['RA']          =  ra
                newentry['DEC']         =  dec
                newentry['Fpeak(Crab)'] =  mxf
                newentry['MJDpeak']     =  mxmjd
                newentry['RecordInst']  =  'BAT'
                newentry['Frecent[BAT]']   = np.round(frecent*1000)/1000
                newentry['MJDrecent[BAT]'] = drecent
                newentry['Frecent[MAXI]']    = -9.999
                newentry['MJDrecent[MAXI]']  = -99
                newentry['Category']        = category
                               
                # IF ACTIVITY = 0 here for a new entry, then this is an odd case - source was active semi-recently but now is off...
                #     until MAXI finds this source and adds it to mtab with a correspodning maxiID, can't crosscheck for mid-brightness
                #     flicker presence (i.e., between "limitup and limitdown"... fortunately, the next section in the code
                #     *should* detect any such case.
                if activity == 0:
                    _send_email_notice(type='FastTransient?',name=obj,sendto='JFS',home=home,message='Transient was bright and then off: '+xx+' with flux '+np.str(mxf)+'...'+coords[0]+','+coords[1]+'\n Verify whether real or screen out if not.')   # not bothering to do this for sources that aren't sufficiently bright#
                newentry['Activity']        = activity

                newentry['SigCount']        = len(plnewsource)
            
                mtab=np.append(mtab,newentry)
                updated = np.append(updated,-1)  # need to search BAT still

                

    #for each LOW detection source, *update* check against active / BRIGHT list (i.e., only remains active if activity is already 1)
    lnew = lnew_maxi ## running through low maxi's first 
    for yy in uniq_maxilow:
        xx = str(yy).replace("b'","").replace("'","")
        obj=mr.get_maxi_objname(  str(xx).replace("b'","").replace("'","") ) 
        obj0=obj
        tossout=0

        if xx[0]!='J':  # screening out bad maxiIDs
            tossout=1

        pl = (np.where(newIDs_maxi == xx))[0]
            
        for el in bdata['MATCHTXT']:
            if fnmatch.fnmatch(obj.lower(),(el.lower())):
                tossout = 1
                print('THROWING OUT '+obj+' ::  MATCH to '+el)
        if tossout==0:  # can continue - not a string match for a bad object
            obj0=obj
            if obj == '404_Not_Found':
                obj='Unknown'
                obj0=obj
            
            pageworks = mr.test_simbad_page(obj,mirror=mirror)        
        
            if pageworks != 1:
                # try out the usual name updates
                obj=obj.replace('_and_SNR','').replace('_SNR','').replace('Galactic_Center_Region','SgrA*').replace('_Slow-Burster_with_Rapid-Burster','')
                obj=obj.replace('NGC_6814_with_V1432_Aql','V1432_Aql').replace('_Rapid-Burster_with_Slow-Burster','').replace('LS_I_+61_303_with_Swift_J0243.6+6124','Swift_J0243.6+6124')
                obj=obj.replace('_with_GS_0836-429','').replace('_with_Terzan_1','')
                pageworks = mr.test_simbad_page(obj,mirror=mirror)
                if pageworks != 1:
                    if obj[0:2]=='A_':
                        obj='1'+obj
                        if obj[0:2]=='X_':
                            obj=('foo'+obj).replace('fooX_','3A_')                       
                            pageworks = mr.test_simbad_page(obj,mirror=mirror)

            if obj == '404_Not_Found':
                obj='Unknown'
                obj0=obj
                coords=['INDEF','INDEF']
            else:
                if pageworks != 1:
                    if xx[0]=='J':
                        coords=(xx+' '+xx).split()
                    else:
                        coords=['INDEF','INDEF']
                    print('NOOOOOOOOOOO!!!!')
#                    _send_email_notice(type='BadEntry',name=obj,sendto='JFS',home=home,message='No location for MAXI source '+xx+' with flux '+np.str(mxf)+'...'+coords[0]+','+coords[1])   # not bothering to do this for sources that aren't sufficiently bright
                else:
                    coords=mr.get_coordinates_icrs(mr.coord_from_simbad(obj,mirror=mirror)).split()

            ra=coords[0]
            dec=coords[1]

            frecent=mres['flux'][lnew][pl][-1]/1000.
            drecent=mres['mjd'][lnew][pl][-1]

            
            # FIND ACTIVITY TO BE UPDATED IN EVENT OF A MATCH 
            activity=0  #initialize under assumption that has not been active
            if frecent < limitdown/1000.:
                if (yy in uniq_mflicker)==True:
                    activity=1   # if flickering for this *primary* maxiID, then activity can be on depending on starting state ... more checks done later in script to confirm

            matches_=[] # initializing
            if xx != 'INDEF':
                matches_ = (np.where(mtab['MaxiID'] == xx))[0]

            if len(matches_) == 0:
                matches_ = (np.where(mtab['Object'] == obj))[0]

            addaltname=False                            
            if len(matches_) == 0:
                mfound=0
                scomplete=0
                incr=0
                while (mfound==0 and scomplete==0):
                    noms = np.str(mtab[incr]['AltNames']).split(',')
                    names = [_.strip() for _ in noms]
                    #matches_ = np.where(names == obj)   ### DISABLING THIS TO AVOID POTENTIAL NAMING MISHAPS... STICKING WITH ANGLE WHICH IS ROBUST UNLESS WANT TO ENABLE THIS AT SOME POINT
                    _RA  = mtab[incr]['RA']
                    _DEC = mtab[incr]['DEC']
                    if len(matches_) > 0:
                        mfound=1
                        matches_ = [incr]
                    if incr == (len(mtab)-1):
                        scomplete=1
                    if mfound==0:
                        #try matching on ra,dec
                        if str(_RA)[0] != 'J' and str(ra)[0] != 'J':
                            #print('fx',incr,obj,ra,dec,_RA,_DEC,mtab[incr]['Object'])
                            adif=mr.get_angsep_arcsec(ra,dec,np.float(_RA),np.float(_DEC),)
                            if adif < angthreshold:
                                mfound=1
                                matches_ = [incr]
                                if names[0]=='none' or names[0]=='--' or names[0]=='0':
                                    anames = xx
                                    addaltname = True    
                                else:
                                    if (xx in names) == False:
                                        anames = np.str(names).replace("['","").replace("']","").replace(" ",",").replace("'","")+','+xx
                                        addaltname = True

                    incr+=1
            if len(matches_) > 0:
                # check brightness, and update elements
                i=matches_[0]
                mtab[i]['SigCount']+=len(pl)
                mtab[i]['Frecent[MAXI]']=frecent
                mtab[i]['MJDrecent[MAXI]']=drecent

                # if mtab Activity has been 1, because it was detected above lowlimit, it remains *on* iff activity=1 (on flickerlist, which is likely for this selection)
                # if mtab Activity has been 0, then should remain 0, since not found on > limitup list.                
                if (mtab[i]['MJDrecent[BAT]'] < drecent):   #Iff most recent data
                    if mtab[i]['Activity']==1:  # case if presently *has* been active
                        if activity==0:  # double-check for flickering
                            if obj.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:  # check if name matches a flickering source in the BAT
                                activity=1
                            else:
                                # doublecheck no flickering under a different name
                                noms = np.str(mtab[i]['AltNames']).split(',')
                                names = [_.strip() for _ in noms]
                                for nn in names:
                                    if nn[0]=='J':  # likely a MAXI-ID
                                        if nn in bintostr(uniq_mflicker):   # check if matches a flickering ID
                                            activity=1
                                    if nn.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:   # check if matches a flickering ID
                                        activity=1                                    
                        mtab[i]['Activity']=activity  ## *Only updates value if current setting (mtab[i]['Activity'])==1
                    
                if addaltname == True:
                    mtab[i]['AltNames']=anames 
                updated[i] = 1  # all set    


    lnew = lnew_bat ## running through low BATs next
    for yy in uniq_batlow:
        xx = str(yy).replace("b'","").replace("'","")
        obj=xx
        obj0=obj
        tossout=0

        pl = (np.where(newIDs_bat == xx))[0]
            
        for el in bdata['MATCHTXT']:
            if fnmatch.fnmatch(obj.lower(),(el.lower())):
                tossout = 1
                print('THROWING OUT '+obj+' ::  MATCH to '+el)
        if tossout==0:  # can continue - not a string match for a bad object            
            ra=sres[lnew][pl][0]['ra']
            dec=sres[lnew][pl][0]['dec']
                                 
            frecent=sres['flux'][lnew][pl][-1]
            drecent=sres['mjd'][lnew][pl][-1]

            
            # FIND ACTIVITY TO BE UPDATED IN EVENT OF A MATCH 
            activity=0  #initialize under assumption that has not been active
            if frecent < limitdown/1000.:
                if (yy in uniq_sflicker)==True:
                    activity=1   # if flickering for this batID, then activity can be on depending on starting state ... more checks done later in script to confirm

                                  
            matches_ = (np.where(mtab['Object'] == obj))[0]
        
            addaltname=False
            if len(matches_) == 0:
                mfound=0
                scomplete=0
                incr=0
                while (mfound==0 and scomplete==0):
                    noms = np.str(mtab[incr]['AltNames']).split(',')
                    names = [_.strip() for _ in noms]
                    #matches_ = np.where(names == obj)   ### DISABLING THIS TO AVOID POTENTIAL NAMING MISHAPS... STICKING WITH ANGLE WHICH IS ROBUST UNLESS WANT TO ENABLE THIS AT SOME POINT
                    _RA  = mtab[incr]['RA']
                    _DEC = mtab[incr]['DEC']
                    if len(matches_) > 0:
                        mfound=1
                        matches_ = [incr]
                    if incr == (len(mtab)-1):
                        scomplete=1
                    if mfound==0:
                        #try matching on ra,dec
                        if str(_RA)[0] != 'J':
                            adif=mr.get_angsep_arcsec(ra,dec,np.float(_RA),np.float(_DEC))
                            if adif < angthreshold:
                                mfound=1
                                matches_ = [incr]
                                if names[0]=='none' or names[0]=='--' or names[0]=='0':
                                    anames = xx
                                    addaltname = True    
                                else:
                                    if (xx in names) == False:
                                        anames = np.str(names).replace("['","").replace("']","").replace(" ",",").replace("'","")+','+xx
                                        addaltname = True

                    incr+=1
            if len(matches_) > 0:
                # check brightness, and update elements
                i=matches_[0]
                mtab[i]['SigCount']+=len(pl)
                mtab[i]['Frecent[BAT]']=np.round(frecent*1000)/1000
                mtab[i]['MJDrecent[BAT]']=drecent
    
                # if mtab Activity has been 1, because it was detected above lowlimit, it remains *on* iff activity=1 (on flickerlist, which is likely for this selection)
                # if mtab Activity has been 0, then should remain 0, since not found on > limitup list.                
                if (mtab[i]['MJDrecent[MAXI]'] < drecent):   #Iff most recent data
                    if mtab[i]['Activity']==1:  # case if presently *has* been active
                        if activity==0:  # double-check for flickering
                            if mtab[i]['MaxiID'] != 'INDEF':                        
                                if mtab[i]['MaxiID'] in bintostr(uniq_mflicker):  # check if name matches a flickering source in MAXI
                                    activity=1
                            else:
                                # doublecheck no flickering under a different name
                                noms = np.str(mtab[i]['AltNames']).split(',')
                                names = [_.strip() for _ in noms]
                                for nn in names:
                                    if nn[0]=='J':  # likely a MAXI-ID
                                        if nn in bintostr(uniq_mflicker):   # check if matches a flickering ID
                                            activity=1
                                    if nn.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:   # check if matches a flickering ID
                                        activity=1                                    
                        mtab[i]['Activity']=activity  ## *Only updates value if current setting (mtab[i]['Activity'])==1
                     
                if addaltname == True:
                    mtab[i]['AltNames']=anames
 
                updated[i] = 1  # all set

                
    # UPDATE ACTIVITY FLAGS FOR ALL TABLE ELEMENTS
    for i in np.arange(len(mtab)):
        el=mtab[i]        
        if el['Activity']==1:
            #see if need to turn activity off due to long absence of detection
            if ((today_-el['MJDrecent[BAT]']) > RecentTimeLimit and (today_-el['MJDrecent[MAXI]']) > RecentTimeLimit):  # no detections by either in [100] days
                mtab[i]['Activity'] = 0
                                                                    

    #########################################
    #   "Hand" editing to adjust unidentified entries; supplied in reference-file
    #########################################
    for u in rdata:
        i=(np.where(mtab['MaxiID']==u['MaxiID']))[0]
        mtab[i][0]['RA']         = u['RA']
        mtab[i][0]['DEC']        = u['DEC']
        mtab[i][0]['Object']     = u['Object']
        mtab[i][0]['OrigObjName']= u['OrigObjName']
    #########################################

    #for non-updated rows, compare against faint points... if low-activity found, can turn activity off &/or update recent flux entries
        
    # SEARCH FOR MISSING SOURCES
    for i in np.arange(len(mtab)):
        if updated[i]==0:    # now look for any new info from the faintlists
            el=mtab[i]
            _RA  = el['RA']
            _DEC = el['DEC']
            _mID = el['MaxiID']
            _src = el['Object']
            #first try to update based on MaxiID, if not INDEF
            if _mID != 'INDEF':
                loc1 = (np.where(newIDs_maxi == _mID))[0]   ## have to index it...                
                if len(loc1) != 0:  # match found... update
                    frecent=mres['flux'][lnew_maxi][loc1][-1]/1000.
                    drecent=mres['mjd'][lnew_maxi][loc1][-1]
                    mtab[i]['SigCount']+=len(loc1)
                    mtab[i]['Frecent[MAXI]']=frecent
                    mtab[i]['MJDrecent[MAXI]']=drecent
            
                    #if Activity already set to 0, do not change, but see about turning off...
                    activity=0 # the default for this case
                    if (mtab[i]['MJDrecent[BAT]'] < drecent):   #Only if most recent data does this offer a chance to turn it off
                        if mtab[i]['Activity']==1:
                            if _src.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:  # check if name matches a flickering source in the BAT
                                activity=1
                            else:
                            # doublecheck no flickering under a different name
                                noms = np.str(mtab[i]['AltNames']).split(',')
                                names = [_.strip() for _ in noms]
                                for nn in names:
                                    if nn[0]=='J':  # likely a MAXI-ID
                                        if nn in bintostr(uniq_mflicker):   # check if matches a flickering ID
                                            activity=1
                                    if nn.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:   # check if matches a flickering ID
                                        activity=1                                    
                            mtab[i]['Activity']=activity   # Only updates for case of initial Activity==1
                        
            else:  # no straightforward match... start a methodical hunt
                for ii in np.arange(len(uniq_maxifaint)):
                    k = str(uniq_maxifaint[ii]).replace("b'","").replace("'","")
                    loc1=(np.where(newIDs_maxi == k))[0]   ## have to index it

                    ra =-99
                    dec=-99
                    if _mf_objname[ii] == '':                        
                        obj=mr.get_maxi_objname(  str(k) )
                        pageworks = mr.test_simbad_page(obj,mirror=mirror)
                        _mf_objname[ii] = obj
                                            
                        if pageworks != 1:
                            # try updating the objectname if a page wasn't found
                            obj=obj.replace('_and_SNR','').replace('_SNR','').replace('Galactic_Center_Region','SgrA*').replace('_Slow-Burster_with_Rapid-Burster','')
                            obj=obj.replace('NGC_6814_with_V1432_Aql','V1432_Aql').replace('_Rapid-Burster_with_Slow-Burster','').replace('LS_I_+61_303_with_Swift_J0243.6+6124','Swift_J0243.6+6124')
                            obj=obj.replace('_with_GS_0836-429','').replace('_with_Terzan_1','')
                            pageworks = mr.test_simbad_page(obj,mirror=mirror)
                            if pageworks != 1:
                                if obj[0:2]=='A_':
                                    obj='1'+obj
                                if obj[0:2]=='X_':
                                    obj=('foo'+obj).replace('fooX_','3A_')                       
                                pageworks = mr.test_simbad_page(obj,mirror=mirror)

                        if pageworks == 1:
                            coords=mr.get_coordinates_icrs(mr.coord_from_simbad(obj,mirror=mirror)).split()
                        
                            ra=coords[0]
                            dec=coords[1]

                        _mf_ra[ii]=ra
                        _mf_dec[ii]=dec

                    ra=_mf_ra[ii]
                    dec=_mf_dec[ii]

                    addaltname=False
                    if dec != -99: # If a real value stored, then check for source offsets
                        adif=mr.get_angsep_arcsec(ra,dec,np.float(_RA),np.float(_DEC))
                        if adif < angthreshold:
                            # Found a match.  Update with the most recent flux data and update the counts.
                            frecent=mres['flux'][lnew_maxi][loc1][-1]/1000.
                            drecent=mres['mjd'][lnew_maxi][loc1][-1]
                            mtab[i]['MaxiID']=_mID
                            mtab[i]['SigCount']+=len(loc1)
                            mtab[i]['Frecent[MAXI]']=frecent
                            mtab[i]['MJDrecent[MAXI]']=drecent

                            noms = np.str(mtab[i]['AltNames']).split(',')
                            names = [_.strip() for _ in noms]
                            if (k not in names):
                                addaltname=True
                                anames = mtab[i]['AltNames']+','+k
                                if names[0] == 'none' or names[0]=='--' or names[0]=='0':
                                    anames=k
                            if addaltname == True:
                                mtab[i]['AltNames']=anames

                            
                            activity=0 # the default for this case
                            if mtab[i]['Activity']==1:
                                if (mtab[i]['MJDrecent[MAXI]'] < drecent):   #Only if most recent data does this offer a chance to turn it off
                                    if k in bintostr(uniq_mflicker):  # check if name matches a flickering source in the BAT
                                        activity=1
                                else:
                                    # doublecheck no flickering under a different name
                                    for nn in names:
                                        if nn[0]=='J':  # likely a MAXI-ID
                                            if nn in bintostr(uniq_mflicker):   # check if matches a flickering ID
                                                activity=1
                                        if nn.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:   # check if matches a flickering ID
                                            activity=1                                    
                                mtab[i]['Activity']=activity   #Only updates for case of initial Activity==1


            if  _src != 'Unknown':
                loc1 = (np.where(newIDs_bat == _src))[0]   ## have to index it...                
                if len(loc1) != 0:  # match found... update
                    frecent=sres['flux'][lnew_bat][loc1][-1]
                    drecent=sres['mjd'][lnew_bat][loc1][-1]
                    mtab[i]['SigCount']+=len(loc1)
                    mtab[i]['Frecent[BAT]']=np.round(frecent*1000)/1000
                    mtab[i]['MJDrecent[BAT]']=drecent
            
                    #if Activity already set to 0, do not change, but see about turning off...
                    activity=0 # the default for this case
                    if (mtab[i]['MJDrecent[MAXI]'] < drecent):   #Only if most recent data does this offer a chance to turn it off
                        if mtab[i]['Activity']==1:
                            if mtab[i]['MaxiID'] != 'INDEF':                        
                                if mtab[i]['MaxiID'] in bintostr(uniq_mflicker):  # check if name matches a flickering source in MAXI
                                    activity=1
                            
                            # doublecheck no flickering under a different name
                            noms = np.str(mtab[i]['AltNames']).split(',')
                            names = [_.strip() for _ in noms]
                            for nn in names:
                                if nn[0]=='J':  # likely a MAXI-ID
                                    if nn in bintostr(uniq_mflicker):   # check if matches a flickering ID
                                        activity=1
                                if nn.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:   # check if matches a flickering ID
                                    activity=1                                    
                        mtab[i]['Activity']=activity   #Only updates for case of initial Activity==1
                        
            else:  # no straightforward match... start a methodical hunt                
                for kk in uniq_batfaint:
                    k = str(kk).replace("b'","").replace("'","")
                    loc1=(np.where(newIDs_bat == k))[0]   

                    addaltname=False
                    if len(loc1) > 0:
                        ra = sres[lnew_bat][loc1][-1]['ra']  # most recent data should have best position
                        dec= sres[lnew_bat][loc1][-1]['dec']
                        if str(_RA)[0] != 'J' and str(ra)[0] != 'J':    # good coordinates exist for both
                            adif=mr.get_angsep_arcsec(ra,dec,np.float(_RA),np.float(_DEC))
                            if adif < angthreshold:  # found a match
                                frecent=sres[lnew_bat][loc1][-1]['flux']
                                drecent=sres[lnew_bat][loc1][-1]['mjd']
                                mtab[i]['SigCount']+=len(loc1)
                                mtab[i]['Frecent[BAT]']=np.round(frecent*1000)/1000
                                mtab[i]['MJDrecent[BAT]']=drecent

                                noms = np.str(mtab[i]['AltNames']).split(',')
                                names = [_.strip() for _ in noms]
                                if (k not in names):
                                    addaltname=True
                                    anames = mtab[i]['AltNames']+','+k
                                    if names[0] == 'none' or names[0]=='--' or names[0]=='0':
                                        anames=k
                                if addaltname == True:
                                    mtab[i]['AltNames']=anames
                                            
                                activity=0 # the default for this case
                                if (mtab[i]['MJDrecent[MAXI]'] < drecent):   #Only if most recent data does this offer a chance to turn it off
                                    if mtab[i]['Activity']==1:
                                        if mtab[i]['MaxiID'] != 'INDEF':                        
                                            if mtab[i]['MaxiID'] in bintostr(uniq_mflicker):  # check if name matches a flickering source in MAXI
                                                activity=1                                        
                                        # doublecheck no flickering under a different name
                                        for nn in names:
                                            if nn[0]=='J':  # likely a MAXI-ID
                                                if nn in bintostr(uniq_mflicker):   # check if matches a flickering ID
                                                    activity=1
                                                if nn.lower() in [q.lower() for q in bintostr(uniq_sflicker)]:   # check if matches a flickering ID
                                                    activity=1                                    
                                        mtab[i]['Activity']=activity   #Only updates for case of initial Activity==1

                        

        #########  Prevent Blank Table AltName Entries from being zerod
        for el in mtab:
            if el['AltNames']=='0' or el['AltNames']=='--':
                el['AltNames']=='none'

                                             
        ##########  SORT THE TABLE & WRITE IT OUT  ###########

        ttmp = Table(mtab)
        scol = copy.deepcopy(ttmp['Fpeak(Crab)'])  
        scol = 1./scol  # to sort by increasing flux
        scol.name = 'sortme'
        ttmp.add_column(scol,index=0)
        
        tout=Table(np.sort(ttmp,order='sortme',axis=0))  
        tout.remove_column('sortme')

        os.system('mv '+masterfile+' '+masterfile.replace('.tab','_backup.tab'))          # Move last file out of the way, and then write update
        atab.write(tout, masterfile, overwrite=True, format='fixed_width',delimiter=' ')

        os.system('mv '+recentfile+' '+recentfile+'.old')   # Move last file out of the way, and then write update
        fil=open(recentfile,"w+")
        fil.write('%15s  %15s\n' % ("MAXILatest", "BATLatest"))
        fil.write('%15.5f %15.5f\n' % (np.float(mres[-1]['mjd']), np.float(sres[-1]['mjd'])))    # keeping some buffer in case of a delayed update... alternatively, could match it to last MAXI entry and last BAT entry....
        fil.close()

