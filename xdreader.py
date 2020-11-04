import numpy as np


def TableDupCheck(table='master_sourcelist.tab',threshold=30):  # Designed for quick "by-hand" checking of coincident sources when investigating table entries
    xx=ascii.read(table)
    for ii in np.arange(len(xx)): 
        yy=xx[ii]; ra1=yy['RA']; dec1=yy['DEC']; src1=yy['MaxiID'] 
        for jj in np.arange(len(xx)): 
            zz=xx[jj]; ra2=zz['RA']; dec2=zz['DEC']; src2=zz['MaxiID'] 
            if src2 != src1 and ra1[0] != 'J' and ra2[0] != 'J':     
                asep=get_angsep_arcsec(ra1,dec1,ra2,dec2) 
                if asep < threshold: 
                    print(src1,src2,asep) 

    

def get_angsep_arcsec(ra1,dec1,ra2,dec2,c1=-99,c2=-99):
    # Calculates the angular separation of two sources in arcsec
    #
    # Can either supply source coordinates directly (c1,c2 astropy coords)
    #   OR by ra,dec in ICRS
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    if c1==-99:
        c1=SkyCoord(np.float(ra1)*u.degree,np.float(dec1)*u.degree,frame='icrs')
    if c2==-99:
        c2=SkyCoord(np.float(ra2)*u.degree,np.float(dec2)*u.degree,frame='icrs')
    sep = c1.separation(c2).arcsec
    return(sep)


def get_coordinates_icrs(galcoord,precision=5):
    # Computes coordinates in ICRS given galactic coordinates as input (expects l,b input extracted from SIMBAD website)
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    
    c = SkyCoord(l=np.float(galcoord['l'])*u.degree,b=np.float(galcoord['b'])*u.degree,frame='galactic')

    # These lines here as useful examples for common other format manipulations
    #c.transform_to('icrs').ra.hms
    #c.transform_to('icrs').dec.dms
    #return(c.transform_to('icrs').to_string('hmsdms'))
    return(c.transform_to('icrs').to_string('decimal',precision=precision))                                                                                                                                                      


def get_maxi_objname(id):  # Query MAXI site to see if MaxiID linked to a MAXI source-page as back-door to getting Source Name.
    import requests
    a=requests.get('http://maxi.riken.jp/star_data/'+id+'/'+id+'.html')
    name=a.text.split('title')[1].replace('>','').replace('</','').replace(' ','_')
    return(name)

def test_simbad_page(id,mirror=True):  # Tests whether a page for a MaxiID exists (returns 1/0)
    # mirror tells whether to use the CfA mirror or not
    import requests
    if mirror==False:
        page='http://simbad.u-strasbg.fr/simbad/sim-id?Ident='+id.replace('+','%2B').replace('_','+')
    else:
        page='https://simbad.harvard.edu/simbad/sim-id?ident='+id.replace('+','%2B').replace('_','+')          
    print('... testing ... '+page)
    p=requests.get(page)
    txt=str(p.text)
    if txt.find('Gal') != -1 and txt.find('incorrect ') == -1:
        return(1)
    else:
        return(0)
    

def coord_from_simbad(id,page='None',mirror=True):
    # Extracts coordinates from SIMBAD webpage
    # Easiest line to reliably identify looks to be the galactic coords.  So pulling that out as a string
    # Takes standard ID format as input (spaces replaced by underscores), and will reformat slightly in actual query to match their schema
    import requests
    if page == 'None':
        if mirror==False:
            page='http://simbad.u-strasbg.fr/simbad/sim-id?Ident='+id.replace('+','%2B').replace('_','+')
        else:
            page='https://simbad.harvard.edu/simbad/sim-id?ident='+id.replace('+','%2B').replace('_','+')
    p=requests.get(page)
    txt=str(p.text).split('\n')
    endr = False ; i=220 #Starting line pos, actual varies source to source but typically above this by several tens
    while (endr == False):
        i=i+1
        endr=(txt[i][0:3]=='Gal')
    galcrd=txt[i+10]

    gtest = galcrd.split('<')
    if len(gtest) > 6:  # Determined to work from having thrown lots of real and fake queries at it
        galcrd=gtest[0]+gtest[2].split('>')[1]+gtest[4].split('>')[1]+gtest[6].split('>')[1]
    
    GALra = galcrd.split()[0]
    GALdec= galcrd.split()[1]
    coord = {'l':GALra, 'b':GALdec}
    return(coord)



def asmgo(daymod=100,bkdir='bs_pyfile_backups/'):   # Major process designed to read in and cobble together a fresh ASM sourcelist and pickle it
    # Using default settings for subsprocesses
    import gzip
    from astropy.time import Time

    today=str(int(Time.now().mjd))

    asmproc(savefile='asmlist.pysav')

    # Potential option down the road if want to do auto-saving here:
    #if np.mod(np.float(today),daymod)==0:
    #    with open('asmlist.pysav','rb') as f_in, gzip.open('asmlist_'+today+'.pysav.gz','wb') as f_out:
    #        f_out.writelines(f_in)




def batgo(daymod=1,bkdir='bs_pyfile_backups/'): # Major process designed to update the BAT sourcelist and pickle it
    # Using default settings for subprocesses
    import gzip
    import shutil
    from astropy.time import Time
    import xdfetchfiles
    xdfetchfiles.update_BAT_list()

    shutil.move('batlist.pysav','batlist-prev.pysav')
    today=str(int(Time.now().mjd))

    batuproc(restfile='batlist-prev.pysav',savefile='batlist.pysav')

    # Potential option down the road if want to do auto-saving here:
    #if np.mod(np.float(today),daymod)==0:
    #    with open('batlist.pysav','rb') as f_in, gzip.open(bkdir+'batlist_'+today+'.pysav.gz','wb') as f_out:
    #        f_out.writelines(f_in)
                                                            


def maxigo(daymod=1,bkdir='bs_pyfile_backups/'): # Major process designed to update the MAXI sourcelist and pickle it
    # Using default settings for subprocesses
    import gzip
    import shutil
    from astropy.time import Time
    import xdfetchfiles
    xdfetchfiles.update_maxi_list()

    shutil.move('maxilist.pysav','maxilist-prev.pysav')        
    today=str(int(Time.now().mjd))
    
    maxiuproc(restfile='maxilist-prev.pysav',savefile='maxilist.pysav')

    # Potential option down the road if want to do auto-saving here:
    #if np.mod(np.float(today),daymod)==0:
    #    with open('maxilist.pysav','rb') as f_in, gzip.open(bkdir+'maxilist_'+today+'.pysav.gz','wb') as f_out:
    #        f_out.writelines(f_in)

                                                                

def maxiuproc(maxnum=100,restfile='None',savefile='None'):   # Subprocess that updates the MAXI data structures with incremental, new datasets
    # maxnum sets the flux-ordered number of sources that will be considered per day
    # Restores pickle file from restfile, and saves file in savefile
    import numpy as np
    import glob
    import os
    import pickle
    import time
    mfils=(glob.glob('fluxtop?????.dat'))  # stores list of MAXI daily top-flux files
    srtem=mfils.sort()   # glob apparently doesn't order things correctly... this fixes the issues

    d1=os.popen('date').read()    # Tracking the time when run
    init=0
    cntr=0

    # Option to restore from the pickle-save file into "res"
    if restfile != 'None':
        pfile = open(restfile, 'rb')           
        res= pickle.load(pfile)
        pfile.close()
        init=1  # initialized
        cntr=1  # initialized

        vals = [int(g.replace('fluxtop','').replace('.dat','')) for g in mfils]  # dates as values
        maxv = int(res['mjd'][-1])  # finds the most recent date in the pickle-file
        vv=np.array(vals)  
        vv=vv[vv > maxv]  # Contains all dates more recent than the last save
        mfils=['fluxtop'+str(_vv_)+'.dat' for _vv_ in vv]  # Files with most recent dates
        
    #    print(mfils,init,cntr,d1)


    #  Check if each file has data, and if so, ingest
    for f in mfils:
        wc=int(os.popen('wc -l '+f).read().split(  )[0])
        cntr+=1
        print(cntr,f)
        if wc >= 10:  #  contains useful data
            _=maxiread(f,maxnum=maxnum)
            if init == 0:
                res=_
                init=1
            else:
                res=np.append(res,_)
    d2=os.popen('date').read()   # get the date,time of completion
    #time.sleep(0.1)
    print(str(cntr)+' FILES ....')
    print(d1,d2)
    if savefile != 'None':
        pickle.dump(res,open(savefile,"wb"))  # store "res" in a pickle
    return(res)




def maxiproc(maxnum=100,savefile='None'):  # Subprocess to read in a fresh MAXI sourcelist 
    import numpy as np
    import glob
    import os
    import pickle
    mfils=(glob.glob('fluxtop?????.dat'))  # stores list of MAXI daily top-flux files
    srtem=mfils.sort()   # glob apparently doesn't order things correctly... this fixes the issues
    d1=os.popen('date').read()     # Tracking the time when started 
    init=0
    cntr=0
    print(mfils,init,cntr,d1)

    #  Check if each file has data, and if so, ingest
    for f in mfils:
        wc=int(os.popen('wc -l '+f).read().split(  )[0])
        cntr+=1
        print(cntr,f)
        if wc >= 10:  # then contains data
            _=maxiread(f,maxnum=maxnum)
            if init == 0:
                res=_
                init=1
            else:
                res=np.append(res,_)
    d2=os.popen('date').read()      # the time when completed
    print(str(cntr)+' FILES ....')
    print(d1,d2)
    if savefile != 'None':
        pickle.dump(res,open(savefile,"wb"))   # store "res" in a pickle
    return(res)
            

def batuproc(restfile='None',savefile='None'):   # Subprocess that updates the BAT data structures with incremental, new datasets
    # Restores pickle file from restfile, and saves file in savefile
    import numpy as np
    import glob
    import os
    import pickle
    import time
    fils=(glob.glob('BAT_record_?????.dat'))  # stores list of BAT top-flux files
    srtem=fils.sort()  # glob apparently doesn't order things correctly... this fixes the issues

    d1=os.popen('date').read()     # Tracking the time when run
    init=0   
    cntr=0   

    # Option to restore from the pickle-save file into "res"
    if restfile != 'None':
        pfile = open(restfile, 'rb')
        res= pickle.load(pfile)
        pfile.close()
        init=1  # initialized
        cntr=1  # initialized

        vals = [int(g.replace('BAT_record_','').replace('.dat','')) for g in fils]
        maxv = int(res['mjd'][-1])
        vv=np.array(vals)
        vv=vv[vv > maxv]
        fils=['BAT_record_'+str(_vv_)+'.dat' for _vv_ in vv]

    #    print(mfils,init,cntr,d1)                                                                                                                                                  

    
    #  Check if each file has data, and if so, ingest
    for f in fils:
        wc=int(os.popen('wc -l '+f).read().split(  )[0])
        cntr+=1
        print(cntr,f)
        if wc >= 10:
            _=batread(f)
            if init == 0:
                res=_
                init=1
            else:
                res=np.append(res,_)
                d2=os.popen('date').read()
    print(str(cntr)+' FILES ....')
    d2=os.popen('date').read()    
    print(d1,d2)
    if savefile != 'None':
        pickle.dump(res,open(savefile,"wb"))  # store "res" in a pickle
    return(res)


def batproc(savefile='None'):    # Subprocess to read in a fresh BAT sourcelist
    import numpy as np
    import glob
    import os
    import pickle
    fils=(glob.glob('BAT_record_?????.dat'))  # stores list of BAT daily top-flux files
    srtem=fils.sort()    # glob apparently doesn't order things correctly... this fixes the issues

    d1=os.popen('date').read()  # The time when started 
    init=0
    cntr=0

    #  Check if each file has data, and if so, ingest
    for f in fils:
        wc=int(os.popen('wc -l '+f).read().split(  )[0])
        cntr+=1
        print(cntr,f)
        if wc >= 10:   # then contains data
            _=batread(f)
            if init == 0:
                res=_
                init=1
            else:
                res=np.append(res,_)
                d2=os.popen('date').read()  # the time when completed
    print(str(cntr)+' FILES ....')
    d2=os.popen('date').read()
    print(d1,d2)
    if savefile != 'None':
        pickle.dump(res,open(savefile,"wb"))  # store "res" in a pickle
    return(res)



def asmproc(savefile='None'):  # Subprocess to read in the ASM sourcelist 
    import numpy as np
    import glob
    import os
    import pickle
    fils=(glob.glob('ASM_bright_source*py.dat'))  # stores list of ASM top-flux files
    srtem=fils.sort()  # glob apparently doesn't order things correctly... this fixes the issues

    d1=os.popen('date').read()  # Tracking the time when started 
    init=0
    cntr=0

    for f in fils:
        cntr+=1
        print(cntr,f)
        _=asmread(f)
        if init == 0:
            res=_
            init=1
        else:
            res=np.append(res,_)
            d2=os.popen('date').read()     # the time when completed
    print(str(cntr)+' FILES ....')
    d2=os.popen('date').read()
    print(d1,d2)
    if savefile != 'None':
        pickle.dump(res,open(savefile,"wb"))   # store "res" in a pickle
    return(res)



def batread(infile):     # Reads in a BAT daily file
    import astropy.io.ascii as atab
    import copy

    # Per-souce structure to hold in and pass along the table's data 
    dtmplate = np.array([('-1',-1,0,'maxiname',0.,0.,'othernames go here',0.,0.)],dtype=[('dhuman','S12'),('mjd','f4'),('order','i4'),('batID','S20'),('flux','f4'),('err','f4'),('o\
therID','S200'),('ra','f4'),('dec','f4')])

    tab=atab.read(infile)  
    dtmplate['mjd']=str(infile).replace('BAT_record_','').replace('.dat','')  # update the default for the structure template

    init=0
    ii=0
    for l in tab:        
        tres=copy.deepcopy(dtmplate)    # initalize structure for a fresh entry 
        tres['order']=ii
        ii=ii+1
        tres['batID']=l['Source']
        tres['flux']=l['Flux(Crab)']
        tres['err']=0.01
        tres['ra']=l['RA']
        tres['dec']=l['DEC']
        tres['otherID']=l['AltName']

        if init == 0:
            res = tres
            init=1
        else:
            res=np.append(res,tres)  

    return(res)   # Read in to "res"... pass it back


def asmread(infile):   # Reads in a ASM file
    import astropy.io.ascii as atab
    import copy

    # Per-souce structure to hold in and pass along the table's data 
    dtmplate = np.array([('-1',-1,0,'asmname',0.,0.,'othernames go here',0.,0.)],dtype=[('dhuman','S12'),('mjd','f4'),('order','i4'),('asmID','S20'),('flux','f4'),('err','f4'),('otherID','S200'),('ra','f4'),('dec','f4')])

    tab=atab.read(infile,delimiter=' ')
    dtmplate['dhuman']=str(infile).replace('ASM_bright_source_','').replace('_py.dat','')   

    init=0
    ii=0
    for l in tab:
        tres=copy.deepcopy(dtmplate)   # initalize structure for a fresh entry 
        tres['order']=ii
        tres['mjd'] =l['col2']
        ii=ii+1
        tres['asmID']=l['col1']
        tres['flux']=l['col3']
        tres['err']=l['col4']
        tres['ra']=l['col6']
        tres['dec']=l['col7']   

        if init == 0:
            res = tres
            init=1
        else:
            res=np.append(res,tres)

    return(res)   # Read in to "res"... pass it back



def maxiread(infile,maxnum=100):   # Reads in a MAXI daily file
    import numpy as np
    import copy
    from  xdupdate import  _send_email_notice
    
    f = open(infile,'r')    # Read in line by line 
    lines=f.readlines()

    cnd=False  # Initializing
    i=2  # Skip first two header lines

    # Find the first data row
    while cnd==False:
        i=i+1
        cnd=(lines[i][0:4]=='001 ') 
        if i > 30:   # Should not take this long to get there!!
            cnd=True # artificially impose a halt, something wrong 
            _send_email_notice(type='BadMaxiTable',name=infile,sendto='JFS',message='Could no read in '+infile+'... check format for a glitch and see about fixing it.') 
            
    srclist=lines[i:min(len(lines)-3,i+maxnum)]   # read in up to maxnum lines

    # Per-souce structure to hold in and pass along the table's data     
    dtmplate = np.array([('-1',-1,0,'maxiname',0.,0.,'othernames go here')],dtype=[('dhuman','S12'),('mjd','f4'),('order','i4'),('maxiID','S20'),('flux','f4'),('err','f4'),('otherID','S200')]) 

    dtmplate['dhuman']=lines[0].split()[0]  # "human-readable" date, if want to use
    dtmplate['mjd'] = lines[0].split()[1]   #  MJD is in the header (though also in the filename)
    
    init=0
    
    # Build res array of all the entries
    for l in srclist:
        _=(l.strip()).split()
        tres=copy.deepcopy(dtmplate)   # initalize structure for a fresh entry
        print(_)
        tres['order']=_[0]
        tres['maxiID']=_[1]
        tres['flux']=_[2]
        tres['err']=_[3]

        ttxt=''
        for k in _[5:]:   # Additional names stored at the end
            ttxt+=' '+k
        tres['otherID']=ttxt

        if init == 0:
            res = tres
            init=1
        else:
            res=np.append(res,tres)

    return(res)    # Read in to "res"... pass it back

            
        
    
