#!/usr/bin/python

def Run_files():
  file = "Samples_example.txt"
  fh=open(file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,source,dir = l[0],l[1],l[2]
      command1= "python Generate_networks_global_2.0.py "+source+" "+id+" "+dir
      print command1
  fh.close()
  return()

Run_files()





