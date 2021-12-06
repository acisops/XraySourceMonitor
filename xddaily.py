#!/usr/bin/env python3


# Script to execute every day.
# Goals:
#   1.  update BAT & MAXI daily tables
#   2.  update BS Catalogue and extract new static-source list
#   3.  check for errors by making sure no dramatically longer or shorter cat/list. 
#   4.  store daily versions of the cat&list and at longer-intervals, keep a
#            pickle of the stored MAXI & BAT data arrays.

import os 
import numpy as np
import gzip
import xdupdate as uu
import xdreader as rr
import xdgenlist as gg
import shutil
from astropy.time import Time
import glob

from subprocess import check_output

#  Backup files live here
bkdir='bs_pyfile_backups/'

#  Every daymod days, add pickles to the pantry
daymod=5  ## update once more stable

#  Update daily tables  (ignore the daymod keyword here, that is currently not implemented, but may be later...)
rr.maxigo(daymod=1000)
rr.batgo(daymod=1000)

#  Pickles --> Pantry?
today=str(int(Time.now().mjd))
if np.mod(np.float(today),daymod)==0:   
    with open('batlist.pysav','rb') as f_in, gzip.open(bkdir+'batlist_'+today+'.pysav.gz','wb') as f_out:
        f_out.writelines(f_in)
    with open('maxilist.pysav','rb') as f_in, gzip.open(bkdir+'maxilist_'+today+'.pysav.gz','wb') as f_out:
        f_out.writelines(f_in)
    

        
#  The truth is out there ... daily table for potential lr use from "xfile"
xfile='bright_xray.txt'
csvfile='bright_sources.csv'

head=['#############################################################################################',
      '# Creation Date: '+str(Time.now().value)+' UT',
      '#',
      '# bright_sources.csv is generated daily via the xddaily.py and xdgenlist.py routines. This',
      '# is a comma-separated variable standard table, with each entry consisting of the following',
      '# attributes:',
      '#    - Object:      Simbad-recognized standard name for the object, where one is available.',
      '#                     Missing entries are reported as -1.0.',
      '#',
      '#    - RA:          Right ascension in degrees (J2000).',
      '#',
      '#    - DEC:         Declination in degrees (J2000).',
      '#',
      '#    - Fpeak(Crab): The peak reported brightness for the source, in units of the Crab. This',
      '#                     is a reasonable proxy for the potential hazard posed by the source to',
      '#                     the science instruments.',
      '#',
      '#    - Frecent[MAXI]:  The most recent brightness (in Crab) from MAXI monitoring.',
      '#',
      '#    - Drecent[MAXI]:  The most recent date of detection with MAXI.',
      '#',
      '#    - Frecent[BAT]:   The most recent brightness (in Crab) from Swift BAT monitoring.',
      '#',
      '#    - Drecent[BAT]:   The most recent date of detection with Swift BAT.',
      '#',
      '#    - Activity:    A boolean 1/0 indicating whether or not a source has been recently active',
      '#                      (having exceeded 80 mCrab in the last 100 days) or not.',
      '#############################################################################################',
      '#']


# Check present numbers of entries
xf1 = int(check_output(["wc", "-l", xfile]).split()[0])
xm1 = int(check_output(["wc", "-l", 'master_sourcelist.tab']).split()[0])                                    


# Update cat & table
uu.update_master_record(home=False,mirror=True)
gg.output_vtable(output=xfile)
gg.output_csvtab(output=csvfile,header=head,HumanDates=True)


# Create the daily log / backups of cat & table
shutil.copyfile('master_sourcelist.tab',bkdir+'master_sourcelist_'+today+'.tab')
shutil.copyfile(xfile,bkdir+xfile.replace('.txt','_'+today+'.txt'))
shutil.copyfile(csvfile,bkdir+csvfile.replace('.csv','_'+today+'.csv'))
shutil.copyfile(xfile,xfile.replace('bright','static'))   ## hack to make it back-compatible


# Check new numbers of entries
xf2 = int(check_output(["wc", "-l", xfile]).split()[0])
xm2 = int(check_output(["wc", "-l", 'master_sourcelist.tab']).split()[0])


# Email alert if any funny business is afoot
if np.abs(xf1-xf2) > 10:
    uu._send_email_notice(type='BRIGHT-list Issue',message='Bright list len: '+np.str(xf1)+'-->'+np.str(xf2))
if np.abs(xm1-xm2) > 20:
    uu._send_email_notice(type='MASTER-list Issue',message='MASTER list len: '+np.str(xf1)+'-->'+np.str(xf2))


#  TESTING PURPOSES:  can use below to check email functionality working properly        
# uu._send_email_notice(type='MASTER-list Issue',message='Test Only')
    

print(Time.now())
