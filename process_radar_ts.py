#!/usr/bin/env python
# coding: utf-8

import sys, getopt

import chilbolton_ts_utils as chilts
import os, getpass, glob
import numpy as np


user = getpass.getuser()


def main(argv):

   inpath = '.'
   outpath = '.'

   try:
      opts, args = getopt.getopt(argv,"hr:d:i:o:p:t:",["radar=","date=","inpath=","outpath=","projectfile=","tag="])
   except getopt.GetoptError:
      print ('process_radar_ts.py -r <radar> -d <date> -i <input_path> -o <output_path> -p <project_file> -t <tag>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('process_radar_ts.py -r <radar> -d <date> -i <input_path> -o <output_path> -p <project_file> -t <tag>')
         sys.exit()
      elif opt in ("-r", "--radar"):
         radar = arg
      elif opt in ("-d", "--date"):
         datestr = arg
      elif opt in ("-i", "--inpath"):
         inpath = arg
      elif opt in ("-o", "--outpath"):
         outpath = arg
      elif opt in ("-p", "--projectfile"):
         project_file = arg
      elif opt in ("-t", "--tag"):
         amof_tag = arg
   print ('Radar is ', radar);
   print ('Input path is ', inpath);
   print ('Output path is ', outpath);
   print ('AMOF tag is ', amof_tag);
   if radar in ("camra", "CAMRa"):
      print("Processing CAMRa data");
      chilts.process_camra_ts(datestr,inpath,outpath,project_file,amof_tag);
   elif radar in ("copernicus", "Copernicus"):
      print("Processing Copernicus data");
      chilts.process_copernicus_ts(datestr,inpath,outpath,project_file,amof_tag);
   elif radar in ("galileo", "Galileo"):
      print("Processing Galileo data");
      chilts.process_galileo_ts(datestr,inpath,outpath,project_file,amof_tag);

if __name__ == "__main__":
   main(sys.argv[1:])



