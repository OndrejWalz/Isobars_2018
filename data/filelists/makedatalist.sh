#!/bin/bash

rm -f isobars_2018.list
get_file_list.pl -keys path,filename -cond trgsetupname=production_isobar_2018,production=P20ic,filetype=daq_reco_picoDst,storage!=HPSS -limit 0 -delim "/" | \
sed 's|^/home/starlib|root://xrdstar.rcf.bnl.gov:1095//home/starlib|g' > isobars_2018.list

