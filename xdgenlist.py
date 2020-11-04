import numpy as np
def gen_init_list(outfile='master_sourcelist.tab',reflist='historical_transients.list',limitup=80,limitdown=30,angthreshold=30,excludelist='blacklist_sources.list',RecentTimeLimit=100,mirror=True, maxilist='maxilist.pysav',batlist='batlist.pysav',asmlist='asmlist.pysav'):
    import glob
    from astropy.time import Time
    import os
    import pickle
    import xdreader as mr
    import copy
    import astropy.io.ascii as atab
    import fnmatch
    import string
    
    #mfils = glob.glob('maxilist_?????.pysav')
    #msrt = mfils.sort()

    
    # Read in table of historical data / historically active sources
    hdata = atab.read(reflist)

    # List of wildcards regexp for excluding extended sources (clusters, etc.) that may trip
    bdata = atab.read(excludelist)
    

    today=str(int(Time.now().mjd))

    # Restore the pickle of MAXI sources over time
    pfile = open(maxilist, 'rb')           
    mres= pickle.load(pfile)
    pfile.close()


    # Restore the pickle of BAT sources over time
    pfile = open(batlist, 'rb')           
    sres= pickle.load(pfile)
    pfile.close()

    # Restore the pickle of (old) ASM files
    pfile = open(asmlist, 'rb')
    ares= pickle.load(pfile)
    pfile.close()

    
    
    # Find all unique instances of MAXI flux (mres) > threshold up (for the rise)
    # Something is empirically screwy for the early data.  So toss out > 10 Crab
    # data for dates < 56250 ... (the 55700-56200 range is problematic)

    #   Find all corresponding unique identifiers in mres
    #   Generate initial list and "dedupe" the list


    # Extract list of significant source detections as "mb"
    #mb = mres[mres['flux'] > limitup and ( mres['mjd'] > 56250 or mres['flux'] < 1e4)]
    ma1 = mres[mres['mjd'] <=  56250]
    ma2 = mres[mres['mjd'] > 56260]
    ma1 = ma1[ma1['flux'] > limitup]
    ma1 = ma1[ma1['flux'] < 500]   ## controlling for spurrious results early on; faint enough to trigger to be on the list are still allowed... just not allowed to be permanent & high
    ma2 = ma2[ma2['flux'] > limitup]
    mb = np.append(ma1,ma2)

    
    # Get the list of IDs for these events
    mbids = mb['maxiID']
    # Find all unique identifiers.  First letter has to be a 'J' to prevent problematic additions... ~40 removed
    s  = list(set(mbids))
    s1 = [q for q in s if q[0:1]==b'J']  # gives the reduced set with ~300 entries

    # Require at least two significant days to count as robust detection for MAXI
    ss = copy.deepcopy(s1)
    for q in s1:
        cnt = sum(1 for i in mbids if i==q)
        if cnt < 2:
            ss = [j for j in ss if j != q]   # screens out anything IDs with < 2 threshold-crossing detections
            

    # Write out a header for the table to be generated
    if outfile=='None':
        outfile='master_sourcelist.tab'
    fil=open(outfile,"w+")
    fil.write('%50a   %25a    %16a    %16a    %16a    %16a   %20a   %16a   %16a   %20a   %16a   %20a   %16a   %16a   %16a   %60a\n' % ('OrigObjName','Object',  'MaxiID',   'RA',   'DEC',  'Fpeak(Crab)',  'MJDpeak', 'RecordInst', 'Frecent[MAXI]', 'MJDrecent[MAXI]', 'Frecent[BAT]','MJDrecent[BAT]','Category', 'Activity','SigCount','AltNames'))


    # Initialize arrays to hold per-source quantities, indexed by the source
    #  (may want to streamline this to have a source class in later iteration...advantage now is flexible types)
    cntr=0
    ll=len(ss)
    _mcnt=np.array(np.zeros(ll),dtype=object)
    _cat =np.array(np.zeros(ll),dtype=object)
    _maxiid   =np.array(np.zeros(ll),dtype=object)
    _activity =np.array(np.zeros(ll),dtype=object)
    
    _origname =np.array(np.zeros(ll),dtype=object)
    _objname  =np.array(np.zeros(ll),dtype=object)
    _ra       =np.array(np.zeros(ll),dtype=object)
    _dec      =np.array(np.zeros(ll),dtype=object)
    _fpeak    =np.array(np.zeros(ll),dtype=object)
    _mjdpeak  =np.array(np.zeros(ll),dtype=object)
    _ipeak    =np.array(np.zeros(ll),dtype=object)
    _frecmaxi =np.array(np.zeros(ll),dtype=object)
    _mjdrecmaxi=np.array(np.zeros(ll),dtype=object)
    _frecbat  =np.array(np.zeros(ll),dtype=object)
    _mjdrecbat=np.array(np.zeros(ll),dtype=object)


    # Populate arrays with per-source characeristics (e.g., maximum flux and date of max, etc.)
    for nom in ss:
        loc=np.where(mb['maxiID'] == nom)

        name=str(nom).replace("b'","").replace("'",'')  # takes care of "byte" formatting
        obj=mr.get_maxi_objname(name)
        obj0=obj
        if obj == '404_Not_Found':
            obj='Unknown'
            obj0=obj

        # Coordinates not supplied by MAXI.  Will come from cross-referencing in simbad
        pageworks = mr.test_simbad_page(obj,mirror=mirror)  # Check to see if simbad page found with identifier name or not (for coordinate extraction)
        
        
        print(cntr,name,obj)

        # Try common or forced fixes to try to update the name to arrive at a simbad match
        if pageworks != 1:
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
    
        # In case no results found, put in placeholding indefinite coords
        if obj == '404_Not_Found':
            obj='Unknown'
            obj0=obj
            coords=['INDEF','INDEF']
        else:
            if pageworks != 1:
                coords = (name+' '+name).split()
                print,'NOOOOOOOOOOO!!!!'  # using the maxiID as placeholder for the coordinates
            else:
                # The default case, will use ICRS coords from simbad
                coords=mr.get_coordinates_icrs(mr.coord_from_simbad(obj,mirror=mirror)).split()

        # Quantities related to the peak detection for the source
        maxflux=max(mb['flux'][loc])
        maxloc = np.argmax(mb['flux'][loc])
        maxdate = mb['mjd'][loc][maxloc]


        # Category used as shorthand for the source risk
        #     Cat5 > 10Crab      <--> HIGH RISK!  Always avoid these
        #     3Crab < Cat 4 <= 10 Crab
        #     1Crab < Cat 3 <= 3 Crab
        # 300 mCrab < Cat 2 <= 1 Crab
        # 100 mCrab < Cat 1 <= 300 mCrab
        #     Cat0  < 100 mCrab  <-->  very low risk
        if maxflux > 10000:
            category=5
        else:
            if maxflux > 3000:
                category = 4
            else:
                if maxflux > 1000:
                    category=3
                else:
                    if maxflux > 300:
                        category=2
                    else:
                        if maxflux > 100:
                            category=1
                        else:
                            category=0
            
        # Current implementation is that 
        #   category 0-2 are subject to active/inactive toggles affecting inclusion on static list
        #   category 3-5 are part of a permanent avoidance list


        #  Now find the most MAXI data on the source
        loc2 = np.where(mres['maxiID'] == nom)
        recentflux = mres['flux'][loc2][-1]

        activity=1  # Start by assuming it is active
        if recentflux < limitdown:  
            activity=0   # if flux is presently low set to inactive 
        
        recentdate = mres['mjd'][loc2][-1]

        # If hasn't been detected for a sufficiently long time, assume inactive and void out recent flux/date data
        if (np.float(today)-np.float(recentdate)) > RecentTimeLimit:    #no detections in last ~year...
            activity=0
            recentdate=-99
            recentflux=-9999
           
        # Update all quantity-storing arrays
        _mcnt[cntr]      =sum(1 for j in mbids if j==nom) # count of significant detections
        _cat[cntr]       =category
        _maxiid[cntr]    =name
        _activity[cntr]  =activity
        _origname[cntr]  =obj0
        _objname[cntr]   =obj
        _ra[cntr]        =coords[0]
        _dec[cntr]       =coords[1]
        _fpeak[cntr]     =maxflux/1e3
        _mjdpeak[cntr]   =maxdate
        _ipeak[cntr]     ='MAXI'
        _frecmaxi[cntr]  =recentflux/1e3
        _mjdrecmaxi[cntr]=recentdate
        _frecbat[cntr]   =-1
        _mjdrecbat[cntr] =-1
        
        cntr=cntr+1
            


    ##  Preventing pathology which can happen if tie in primary sorting column occurs!
    sortie = np.arange(len(_fpeak))/1e6 + _fpeak

    # Sort all arrays by peak flux
    _mcnt=[x for _,x in sorted(zip(sortie,_mcnt))]
    _cat=[x for _,x in sorted(zip(sortie,_cat))]
    _maxiid=[x for _,x in sorted(zip(sortie,_maxiid))]
    _activity=[x for _,x in sorted(zip(sortie,_activity))]
    _origname=[x for _,x in sorted(zip(sortie,_origname))]
    _objname=[x for _,x in sorted(zip(sortie,_objname))]
    _ra=[x for _,x in sorted(zip(sortie,_ra))]
    _dec=[x for _,x in sorted(zip(sortie,_dec))]
    _mjdpeak=[x for _,x in sorted(zip(sortie,_mjdpeak))]
    _ipeak=[x for _,x in sorted(zip(sortie,_ipeak))]
    _frecmaxi=[x for _,x in sorted(zip(sortie,_frecmaxi))]
    _mjdrecmaxi=[x for _,x in sorted(zip(sortie,_mjdrecmaxi))]
    _frecbat=[x for _,x in sorted(zip(sortie,_frecbat))]
    _mjdrecbat=[x for _,x in sorted(zip(sortie,_mjdrecbat))]

    _fpeak=[x for _,x in sorted(zip(sortie,_fpeak))]

    _useful = copy.deepcopy(_cat)    # Screening array: will be set to 0 for duplicates or entries to be discarded for any other condition
    _useful = [1 for q in _useful]
    _alterIDs = copy.deepcopy(_cat)
    _alterIDs = ['' for q in _alterIDs]   #initializing


    # Run through list to check for duplicate entries, and condense / update if found
    for j in np.arange(len(_cat)): 
        cntr=len(_cat)-1-j
        mcntr =  _mcnt[cntr]
        category = _cat[cntr]
        name= _maxiid[cntr]   
        activity= _activity[cntr]
        
        obj0=_origname[cntr]
        obj= _objname[cntr] 
        ra= _ra[cntr]       
        dec= _dec[cntr]       
        fpeak= _fpeak[cntr]     
        maxdate= _mjdpeak[cntr]  
        ipeak  = _ipeak[cntr]    
        frecent= _frecmaxi[cntr]  
        recentdate= _mjdrecmaxi[cntr]
        frecbat= _frecbat[cntr]  
        mjdrecbat= _mjdrecbat[cntr]
        useful   = _useful[cntr]
        alterID= _alterIDs[cntr]
        
        # check for any dupes... and remove / update as needed
        for ii in np.arange(cntr): 
            ra1=_ra[ii]; dec1=_dec[ii]; maxi1=_maxiid[ii] ; obj1=_origname[ii]
            if ra1[0] != 'J' and ra[0] != 'J':    # check first that coordinates are real (from simbad, etc.)
                asep=mr.get_angsep_arcsec(ra1,dec1,ra,dec)  # the angular distance between sources
                if asep < angthreshold:   
                    #a dupe is in the house!
                    #recall, running from brightest to faintest, so the point of comparison is now on the current (brightest) vs the newer (fainter)
                    _useful[ii]=0
                    print('DUPE:',maxi1,name,asep)

                    #check on the frecent, recentdate, and activity flags
                    mcntr +=  _mcnt[ii]
                    _mcnt[cntr] = mcntr
                    
                    if _mjdrecmaxi[ii] > _mjdrecmaxi[cntr]:
                        frecent   = _frecmaxi[ii]
                        recentdate= _mjdrecmaxi[ii]
                        activity  = _activity[ii]
                        
                        _frecmaxi[cntr]   =frecent
                        _mjdrecmaxi[cntr] =recentdate
                        _activity[cntr]   =activity
                        
                    if alterID == '':
                        alterID = maxi1
                    else:
                        alterID = alterID+','+maxi1
                    _alterIDs[cntr] = alterID
                            
    #### CURRENTLY have a deduped (flagged) list of the soures

    #### NEXT bring in the BAT data
    ## List of significant source detections stored as "sb"
    sb = sres[sres['flux'] > limitup/1000.]   # Note the difference here that the BAT data have been stored in Crab units (not mCrab)
    sbids = sb['batID']

    ## Get the unique entries
    ss_  = list(set(sbids))


    
    ### IN THE FUTURE WITH A DEEPER LIBRARY, MAY WANT TO IMPLEMENT THIS BY REQUIRING MULTIPLE DETECTIONS ... FOR NOW, LEAVING ALONE
    #------------------    
    #for q in ss_:
    #    cnt = sum(1 for i in sbids if i==q)
    #    if cnt < 2:   
    #        ss_ = [j for j in ss_ if j != q]   # now should only have elements with multiple detections                                                                          

        

    # For each unique BAT source, get attributes of entry
    # Check and dedupe against the MASTER list, then update MASTER list
    for nom in ss_:
        loc=np.where(sb['batID'] == nom)

        name=str(nom).replace("b'","").replace("'",'')   # takes care of "byte" formatting
        obj=name
        obj0=name

        # Record attributes of maximum 
        maxflux=max(sb['flux'][loc])
        maxloc =np.argmax(sb['flux'][loc])
        maxdate = sb['mjd'][loc][maxloc]


        # Match the category
        if maxflux > 10000/1e3:
            category=5
        else:
            if maxflux > 3000/1e3:
                category = 4
            else:
                if maxflux > 1000/1e3:
                    category=3
                else:
                    if maxflux > 300/1e3:
                        category=2
                    else:
                        if maxflux > 100/1e3:
                            category=1
                        else:
                            category=0


        ra = np.str(np.round(np.float(sb['ra'][loc][-1]),4))    # rounding to prevent an occasional minor annoyance I think is due to a FP precision issue causing the string-length to explode
        dec = np.str(np.round(np.float(sb['dec'][loc][-1]),4))  

        # Determine recent activity
        loc2 = np.where(sres['batID'] == nom)  # Finds any source measurements (including below threshold)
        recentflux = sres['flux'][loc2][-1]  

        activity=1  # Start with assumption that is is active
        if recentflux < limitdown/1000.:
            activity=0     # if flux is presently low set to inactive

        recentdate = sres['mjd'][loc2][-1]
                        
        # If hasn't been detected for a sufficiently long time, assume inactive and void out recent flux/date data
        if (np.float(today)-np.float(recentdate)) > RecentTimeLimit:    #no detections in the last ~year...
            activity=0
            recentdate=-99
            recentflux=-9999

            
            
        ###### Compare against all independent (*_USEFUL=1*) entries to do the deduping
        newadd=1       # Starting assumption that it is a brand new source before comparing
        for kk in np.arange(len(_ra)):
            if _useful[kk]==1 and not _ra[kk][0]=='J':    # Screens out duplicates and entries without firm coordinates
                dang=mr.get_angsep_arcsec(_ra[kk],_dec[kk],ra,dec)
                if dang < angthreshold:
                    # Found a position-based match in the existing catalog... update corresponding entries and do not add a fresh source
                    newadd=0
                    print('...MATCHING BAT '+name+' to '+_objname[kk],dang)
                    if recentdate > _mjdrecbat[kk]:   # use most recent data
                        _frecbat[kk]   = recentflux
                        _mjdrecbat[kk] = recentdate
                        if recentdate > _mjdrecmaxi[kk]:
                            _activity[kk] = activity  # update this if the BAT has more recent data 
                    if _fpeak[kk] < maxflux:
                        # update the instrument, category, peak value, and peak date in the table
                        _fpeak[kk] = maxflux
                        _ipeak[kk] = 'BAT'
                        _mjdpeak[kk] = maxdate
                        _cat[kk]=category
                    # Update AltName entries if needed
                    if _objname[kk] != name:
                        if _alterIDs[kk] == '':
                            _alterIDs[kk] = name
                        else:
                            anames =  [q.strip() for q in _alterIDs[kk].split(',')] 
                            if any(q==name for q in anames):
                                _alterIDs[kk] = _alterIDs[kk]+','+name   #updating this list of "source" names at the reference location +/-XX"
                    _mcnt[kk]+=len(loc2)

        if newadd == 1:   # Add a new entry to the master table
            # First, check whether simbad has better coords:
            pageworks = mr.test_simbad_page(name,mirror=mirror)
            if pageworks == 1:
                coords=mr.get_coordinates_icrs(mr.coord_from_simbad(name,mirror=mirror)).split()
                ra = coords[0]
                dec= coords[1]
                
            # Append attributes to all relevant columns (variables)
            _mcnt    =  np.append(_mcnt,len(loc2))
            _cat     = np.append(_cat,category)
            _maxiid  = np.append(_maxiid,'INDEF')
            _activity= np.append(_activity,activity)
            _origname= np.append(_origname,name)
            _objname = np.append(_objname,name)
            _ra = np.append(_ra,ra)
            _dec = np.append(_dec,dec)
            _fpeak = np.append(_fpeak,maxflux)
            _mjdpeak = np.append(_mjdpeak,maxdate)
            _ipeak   = np.append(_ipeak,'BAT')
            _frecmaxi = np.append(_frecmaxi,-9.999)
            _mjdrecmaxi = np.append(_mjdrecmaxi,-99)
            _frecbat = np.append(_frecbat,recentflux)
            _mjdrecbat = np.append(_mjdrecbat,recentdate)
            _useful    = np.append(_useful,1)
            _alterIDs = np.append(_alterIDs,'')


            # Sort all arrays by peak flux
            sortie = np.arange(len(_fpeak))/1e6 + _fpeak
            
            _mcnt=[x for _,x in sorted(zip(sortie,_mcnt))]
            _cat=[x for _,x in sorted(zip(sortie,_cat))]
            _maxiid=[x for _,x in sorted(zip(sortie,_maxiid))]
            _activity=[x for _,x in sorted(zip(sortie,_activity))]
            _origname=[x for _,x in sorted(zip(sortie,_origname))]
            _objname=[x for _,x in sorted(zip(sortie,_objname))]
            _ra=[x for _,x in sorted(zip(sortie,_ra))]
            _dec=[x for _,x in sorted(zip(sortie,_dec))]
            _mjdpeak=[x for _,x in sorted(zip(sortie,_mjdpeak))]
            _ipeak=[x for _,x in sorted(zip(sortie,_ipeak))]
            _frecmaxi=[x for _,x in sorted(zip(sortie,_frecmaxi))]
            _mjdrecmaxi=[x for _,x in sorted(zip(sortie,_mjdrecmaxi))]
            _frecbat=[x for _,x in sorted(zip(sortie,_frecbat))]
            _mjdrecbat=[x for _,x in sorted(zip(sortie,_mjdrecbat))]
            _useful=[x for _,x in sorted(zip(sortie,_useful))]
            _alterIDs=[x for _,x in sorted(zip(sortie,_alterIDs))]
            
            _fpeak   =[x for _,x in sorted(zip(sortie,_fpeak))]
                    

    ### Next up: Include the ASM data:
    ab = ares[ares['flux'] > limitup/1000.]   # Uses Crab units
    abids = ab['asmID']
    ss_  = list(set(abids))   # The list of unique ASM source-names
 
    for nom in ss_:
        loc=np.where(ab['asmID'] == nom)
                                                                                                                                                                                     
        name=str(nom).replace("b'","").replace("'",'')
        obj=name
        obj0=name

        maxflux=max(ab['flux'][loc])
        maxloc =np.argmax(ab['flux'][loc])
        maxdate = ab['mjd'][loc][maxloc]

        # Determine category-level
        if maxflux > 10000/1e3:
            category=5
        else:
            if maxflux > 3000/1e3:
                category = 4
            else:
                if maxflux > 1000/1e3:
                    category=3
                else:
                    if maxflux > 300/1e3:
                        category=2
                    else:
                        if maxflux > 100/1e3:
                            category=1
                        else:
                            category=0


        ra = np.str(np.round(np.float(ab['ra'][loc][-1]),4))
        dec = np.str(np.round(np.float(ab['dec'][loc][-1]),4))

        loc2 = np.where(ares['asmID'] == nom)
        recentflux = ares['flux'][loc2][-1]

        activity=1
        if recentflux < limitdown/1000.:
            activity=0

        # ASM data should all be old, but no reason not to check in case of overlap with BAT or historical entries pending updates
        recentdate = ares['mjd'][loc2][-1]
        if (np.float(today)-np.float(recentdate)) > RecentTimeLimit:    #no detections in last year...                                                                                           
            activity=0
            recentdate=-99
            recentflux=-9999

            
        # Check against additional entries to see if needs to be added as a new source
        newadd=1  # initial assumption that source is new
        for kk in np.arange(len(_ra)):
            if _useful[kk]==1 and not _ra[kk][0]=='J':    # removes any duplicates and sources wiithout sharp coords
                dang=mr.get_angsep_arcsec(_ra[kk],_dec[kk],ra,dec)   
                if dang < angthreshold:
                    # found a match in the existing master catalog...  update the match and do not add a new entry                                                                                                                                          
                    newadd=0
                    print('...MATCHING ASM '+name+' to '+_objname[kk],dang)
                    if _fpeak[kk] < maxflux:  # check for new maximum
                        # if new max, update the instrument, category, peak value, and peak date in the table                                                                                                               
                        _fpeak[kk] = maxflux
                        _ipeak[kk] = 'ASM'
                        _mjdpeak[kk] = maxdate
                        _cat[kk]=category
                    if _objname[kk] != name:   # Update the AltName list to give the source name
                        if _alterIDs[kk] == '':
                            _alterIDs[kk] = name
                        else:
                            anames =  [q.strip() for q in _alterIDs[kk].split(',')]
                            if any(q==name for q in anames):
                                _alterIDs[kk] = _alterIDs[kk]+','+name   #updating this list of "source" names at the reference location +/-XX
                    _mcnt[kk]+=len(loc2)  # update the detection count
                                                               
        if newadd == 1:  # Adding a new entry
            pageworks = mr.test_simbad_page(name,mirror=mirror)   # test for a simbad page match
            if pageworks == 1:
                coords=mr.get_coordinates_icrs(mr.coord_from_simbad(name,mirror=mirror)).split()  # Use simbad coords if possible
                ra = coords[0]
                dec= coords[1]

            # Add new entry to the table arrays
            _mcnt    = np.append(_mcnt,len(loc2))
            _cat     = np.append(_cat,category)
            _maxiid  = np.append(_maxiid,'INDEF')
            _activity= np.append(_activity,activity)
            _origname= np.append(_origname,name)
            _objname = np.append(_objname,name)
            _ra = np.append(_ra,ra)
            _dec = np.append(_dec,dec)
            _fpeak = np.append(_fpeak,maxflux)
            _mjdpeak = np.append(_mjdpeak,maxdate)
            _ipeak   = np.append(_ipeak,'ASM')
            _frecmaxi = np.append(_frecmaxi,-9.999)
            _mjdrecmaxi = np.append(_mjdrecmaxi,-99)
            _frecbat = np.append(_frecbat,-9.999)
            _mjdrecbat = np.append(_mjdrecbat,-99)
            _useful    = np.append(_useful,1)
            _alterIDs = np.append(_alterIDs,'')

            #REORDER ARRAY
            sortie = np.arange(len(_fpeak))/1e6 + _fpeak
            _mcnt=[x for _,x in sorted(zip(sortie,_mcnt))]
            _cat=[x for _,x in sorted(zip(sortie,_cat))]
            _maxiid=[x for _,x in sorted(zip(sortie,_maxiid))]
            _activity=[x for _,x in sorted(zip(sortie,_activity))]
            _origname=[x for _,x in sorted(zip(sortie,_origname))]
            _objname=[x for _,x in sorted(zip(sortie,_objname))]
            _ra=[x for _,x in sorted(zip(sortie,_ra))]
            _dec=[x for _,x in sorted(zip(sortie,_dec))]
            _mjdpeak=[x for _,x in sorted(zip(sortie,_mjdpeak))]
            _ipeak=[x for _,x in sorted(zip(sortie,_ipeak))]
            _frecmaxi=[x for _,x in sorted(zip(sortie,_frecmaxi))]
            _mjdrecmaxi=[x for _,x in sorted(zip(sortie,_mjdrecmaxi))]
            _frecbat=[x for _,x in sorted(zip(sortie,_frecbat))]
            _mjdrecbat=[x for _,x in sorted(zip(sortie,_mjdrecbat))]
            _useful=[x for _,x in sorted(zip(sortie,_useful))]
            _alterIDs=[x for _,x in sorted(zip(sortie,_alterIDs))]

            _fpeak   =[x for _,x in sorted(zip(sortie,_fpeak))]



    
    ### Final set: Historical Data

    # Run through entries, cross-match against current table, and update
    for el in hdata:
        nom=el['SOURCE']
        name=str(nom).replace("b'","").replace("'",'')
        obj=name
        obj0=name

        maxflux = el['FPEAK']
        maxdate = el['MJDPEAK']

        # Match peak to risk category
        if maxflux > 10000/1e3:
            category=5
        else:
            if maxflux > 3000/1e3:
                category = 4
            else:
                if maxflux > 1000/1e3:
                    category=3
                else:
                    if maxflux > 300/1e3:
                        category=2
                    else:
                        if maxflux > 100/1e3:
                            category=1
                        else:
                            category=0


        ra = el['RA']
        dec = el['DEC']

        
        activity=0   # these sources are nearly all *OLD*

        ###### Now compare against the present table entries (_useful=1)         
        newadd=1  # starting assumption
        for kk in np.arange(len(_ra)):
            if _useful[kk]==1 and not _ra[kk][0]=='J':    # removes duplicate entries and any with bad coordinates
                dang=mr.get_angsep_arcsec(_ra[kk],_dec[kk],ra,dec)
                if dang < angthreshold:
                    # found a match in the existing catalog... check and update the source fields, do not add a new entry                                                                                                         
                    newadd=0
                    print('...MATCHING HISTORICAL-SOURCE '+name+' to '+_objname[kk],dang)
                    if _fpeak[kk] < maxflux:  # See if identified a new maximum
                        # If so, update the instrument, category, peak value, and peak date in the table
                        _fpeak[kk] = maxflux
                        _ipeak[kk] = 'HISTORICAL'
                        _mjdpeak[kk] = maxdate
                        _cat[kk]=category
                    if _objname[kk] != name:   # Update the AltName field if needed
                        if _alterIDs[kk] == '':
                            _alterIDs[kk] = name
                        else:
                            anames =  [q.strip() for q in _alterIDs[kk].split(',')]
                            if any(q==name for q in anames):
                                _alterIDs[kk] = _alterIDs[kk]+','+name   #updating this list of "source" names at the reference location +/-XX                                             
                    _mcnt[kk]+=1  # add 1 to the detection count (though obviously this indicates significantly more data existing as a rule)

        if newadd == 1:   # Add new entry to the master table
            _mcnt    =  np.append(_mcnt,1)
            _cat     = np.append(_cat,category)
            _maxiid  = np.append(_maxiid,'INDEF')
            _activity= np.append(_activity,activity)
            _origname= np.append(_origname,name)
            _objname = np.append(_objname,name)
            _ra = np.append(_ra,ra)
            _dec = np.append(_dec,dec)
            _fpeak = np.append(_fpeak,maxflux)
            _mjdpeak = np.append(_mjdpeak,maxdate)
            _ipeak   = np.append(_ipeak,'HISTORICAL')
            _frecmaxi = np.append(_frecmaxi,-9.999)
            _mjdrecmaxi = np.append(_mjdrecmaxi,-99)
            _frecbat = np.append(_frecbat,-9.999)
            _mjdrecbat = np.append(_mjdrecbat,-99)
            _useful    = np.append(_useful,1)
            _alterIDs = np.append(_alterIDs,'')


            #REORDER ARRAY
            sortie = np.arange(len(_fpeak))/1e6 + _fpeak
            _mcnt=[x for _,x in sorted(zip(sortie,_mcnt))]
            _cat=[x for _,x in sorted(zip(sortie,_cat))]
            _maxiid=[x for _,x in sorted(zip(sortie,_maxiid))]
            _activity=[x for _,x in sorted(zip(sortie,_activity))]
            _origname=[x for _,x in sorted(zip(sortie,_origname))]
            _objname=[x for _,x in sorted(zip(sortie,_objname))]
            _ra=[x for _,x in sorted(zip(sortie,_ra))]
            _dec=[x for _,x in sorted(zip(sortie,_dec))]
            _mjdpeak=[x for _,x in sorted(zip(sortie,_mjdpeak))]
            _ipeak=[x for _,x in sorted(zip(sortie,_ipeak))]
            _frecmaxi=[x for _,x in sorted(zip(sortie,_frecmaxi))]
            _mjdrecmaxi=[x for _,x in sorted(zip(sortie,_mjdrecmaxi))]
            _frecbat=[x for _,x in sorted(zip(sortie,_frecbat))]
            _mjdrecbat=[x for _,x in sorted(zip(sortie,_mjdrecbat))]
            _useful=[x for _,x in sorted(zip(sortie,_useful))]
            _alterIDs=[x for _,x in sorted(zip(sortie,_alterIDs))]

            _fpeak   =[x for _,x in sorted(zip(sortie,_fpeak))]

                    
            
    # As a final step, go through the whole catalog to perform a last check-over against duplicate entries
    for j in np.arange(len(_cat)):

        cntr=len(_cat)-1-j

        mcntr =  _mcnt[cntr]
        category = _cat[cntr]
        name= _maxiid[cntr]
        activity= _activity[cntr]
            
        obj0=_origname[cntr]
        obj= _objname[cntr]
        ra= _ra[cntr]
        dec= _dec[cntr]
        fpeak= _fpeak[cntr]
        maxdate= _mjdpeak[cntr]
        ipeak  = _ipeak[cntr]
        frecent= _frecmaxi[cntr]
        recentdate= _mjdrecmaxi[cntr]
        frecbat= _frecbat[cntr]
        mjdrecbat= _mjdrecbat[cntr]
        useful   = _useful[cntr]
        alterID= _alterIDs[cntr]

        
        ## FINAL DEDUPING CHECK HERE
        if useful:
            for ii in np.arange(cntr): 
                ra1=_ra[ii]; dec1=_dec[ii]; maxi1=_maxiid[ii] ; oobj1=_origname[ii] ; obj1= _objname[ii]
                if ra1[0] != 'J' and ra[0] != 'J': 
                    asep=mr.get_angsep_arcsec(ra1,dec1,ra,dec) 
                    if asep < angthreshold or obj1.lower().strip() == obj.lower().strip() or oobj1.lower().strip() == obj.lower():
                        # a dupe is in the house!
                        # recall, running from brightest to faintest, so the point of comparison is now on the current (brightest) vs the newer (fainter)
                        _useful[ii]=0
                        print('DUPE:',maxi1,name,asep)

                        # check and merge/update the frecent, recentdate, and activity flags
                        mcntr +=  _mcnt[ii]
                        _mcnt[cntr] = mcntr
                            
                        if _mjdrecmaxi[ii] > _mjdrecmaxi[cntr]:                            
                            frecent   = _frecmaxi[ii]
                            recentdate= _mjdrecmaxi[ii]
                            if _mjdrecmaxi[ii]  > _mjdrecbat[cntr]:
                                activity = _activity[ii]
                                    

                        if _mjdrecbat[ii] > _mjdrecbat[cntr]:
                            frecbat   = _frecbat[ii]
                            mjdrecbat =  _mjdrecbat[ii]
                            if _mjdrecbat[ii] > _mjdrecmaxi[cntr]:
                                activity  = _activity[ii]
                        
                        _frecbat[cntr]   =frecbat
                        _mjdrecbat[cntr] =mjdrecbat
                        _frecmaxi[cntr]  =frecent
                        _mjdrecmaxi[cntr]=recentdate
                        _activity[cntr]  =activity
                        
                        if alterID == '':   # Update the AltName field in the table as needed
                            if maxi1 != 'INDEF':
                                alterID = maxi1
                        else:
                            if maxi1 != 'INDEF':                        
                                alterID = alterID+','+maxi1
                        _alterIDs[cntr] = alterID

        


        ## Check against list of excluded name keywords to see if the source is extended and should be discarded
        tossout = 0
        for el in bdata['MATCHTXT']:
            if fnmatch.fnmatch(obj.lower(),(el.lower())):
                tossout = 1
                print('THROWING OUT '+obj+' ::  MATCH to '+el)

        # Write out table with any surviving sources (not excluded, not duplicates)
        if useful==1 and tossout==0:
            fil.write('%50a   %25a    %16a    %16a   %16a    %16.4f   %20.3f   %16a   %16.4f   %20.3f   %16.4f   %20.3f   %16i   %16i   %16i   %60a\n' % (obj0, obj, name,  ra,  dec, fpeak,  maxdate, ipeak,  frecent, recentdate, frecbat, mjdrecbat, category, activity, mcntr, alterID))
            print('%50a   %25a    %16a    %16a   %16a    %16.4f   %20.3f   %16a   %16.4f   %20.3f   %16.4f   %20.3f   %16i   %16i   %16i   %60a\n' % (obj0, obj, name,  ra,  dec, fpeak,  maxdate, ipeak,  frecent, recentdate, frecbat, mjdrecbat, category, activity, mcntr, alterID))
        else:
            print('.....SKIPPING '+obj)
        
    fil.close




# Simple routine to read in the master table and write out a "daily" file of sources to avoid based on the current sky
def output_vtable(reftab='master_sourcelist.tab',output='None'):
    import astropy.io.ascii as atab
    tab=atab.read(reftab)
    c1=tab['Object'][tab['Activity']>0]   # All active sources should be included in the list
    c2=tab['Object'][tab['Category']>3]   # All category 4-5 sources are to be included (peaks > 3 Crab) ## Can easily tweak this if desired


    cc=np.append(np.array(c1),np.array(c2))   # merge both lists
    cs = list(set(cc))   # get the unique entries (since no reason a source couldn't be on both)

    olist = np.ones(len(cs))  # will use this for flux-sorting
    
    for i in np.arange(len(cs)):
        loc=np.where(tab['Object'] == cs[i])
        v = tab[loc][0]
        olist[i]=1./v['Fpeak(Crab)']

    cs=[x for _,x in sorted(zip(olist,cs))]  # the sort step

    # printing/writing out the table as whitespace-separated list: (Name, RA, DEC)
    if output != 'None':
        fil=open(output,"w+")
    for el in cs:
        loc=np.where(tab['Object'] == el)
        v = tab[loc][0]
        
        #print('%-50s  %-16s  %-16s' % (v['Object'], v['RA'], v['DEC']))
        if output != 'None':
            fil.write('%-50s  %-16s  %-16s\n' % (v['Object'], v['RA'], v['DEC']))
        else:
            print('%-50s  %-16s  %-16s' % (v['Object'], v['RA'], v['DEC']))

    if output != 'None':
        fil.close

        
