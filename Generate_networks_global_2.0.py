# Immune-Network-Generation Program developed by Rachael Bashford-Rogers (2018)
# at the University of Cambridge
# E-mail: rb520@cam.ac.uk

# If you find the methods in Immune-Network-Generation useful, please cite the following reference:
# Bashford-Rogers, R. et al. Network properties derived from deep sequencing of the 
# human B-cell receptor repertoires delineates B-cell populations. 
# Genome research, doi:10.1101/gr.154815.113 (2013).

#Copyright (C) 2018  Dr Rachael Bashford-Rogers

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


#!/usr/bin/python
import math
import sys
from collections import defaultdict
import os
import commands
import sys
from operator import itemgetter
sys.path.append('networkx-1.10/')
import networkx as nx

class Tree(defaultdict):
  def __init__(self, value=None):
    super(Tree, self).__init__(Tree)
    self.value = value

def fasta_iterator(fh):
  while True:
    line = fh.readline()
    if line.startswith('>'): break	
  while True:
    header = line[1:-1].rstrip()
    sequence = fh.readline().rstrip()
    while True:
      line = fh.readline()
      if not line: break
      if line.startswith('>'): break
      sequence += line.rstrip()
    yield(header, sequence)
    if not line: return

def Cluster_i (Reduced_file,tmp_file,diff):
  cd_hit_directory = "./cd-hit-v4.6.6-2016-0711/"
  command= cd_hit_directory+"cd-hit -i "+Reduced_file+" -o "+tmp_file+" -c "+str(diff)+" -d 80 -T 10  -M 0 -AL 40 "
  os.system(command)
  return()

def Get_seqs_single (file):
  fh=open(file,"r")
  seqs={}
  for header,sequence in fasta_iterator(fh):
    seqs[header]=sequence
  fh.close()
  return(seqs)

def Get_similarity_single(clust_seqs,file_out,c):
  done=Tree()
  total = len(clust_seqs)
  out='' 
  ind1 = 0 
  mismatch = 1
  for i1 in range(0,total-1):
    read1 = clust_seqs[i1]
    id1 = read1[0]
    seq1= read1[1]
    l1 = read1[2]
    done[id1]=1
    for i2 in range(i1+1,total):
      if (i1<i2):
        read2 = clust_seqs[i2]
        id2 = read2[0]
        seq2= read2[1]
        l2 = read2[2]
        if (id2 not in done):
          (s1, s2, p1)=Trim_sequences(seq1,seq2,l1,l2)
          if (p1==1):
            if(s1==s2):
              out=out+"0\t"+id1+"\t"+id2+"\t"+c+"\t"+str(l1)+"\t"+str(l2)+"\n"
              ind1 =ind1+1
            else:
              (p,mm)=Get_diff(s1,s2,mismatch)
              if (p==1 and mm<=mismatch):
                out=out+str(mm)+"\t"+id1+"\t"+id2+"\t"+c+"\t"+str(l1)+"\t"+str(l2)+"\n"
                ind1 = ind1+1
        if(ind1>=100):
          Write_output(out, file_out)
          ind1=0
          out=''
  Write_output(out, file_out)
  del out
  return()

def Get_diff(s1,s2,mismatch):
  p=1
  (mm)=Do_counting(s1, s2,mismatch)
  if (mm > mismatch):
    p=0
  return (p,mm)

def Do_counting(s1, s2,mismatch):
  mm = sum(1 for a, b in zip(s1,s2) if a != b)
  return (mm)

def Trim_sequences(s1,s2,l1,l2):
  p=0
  i=15
  sample1=s1[i:i+25]
  index=s2.find(sample1)
  if (index!=-1):
    if(index>i):
      if(index-i <=20):s2=s2[index-i:l2]
    else:
      if(i-index <=20):s1=s1[i-index:l1]
    min_len=min([len(s1),len(s2)])
    if ((max([len(s1),len(s2)]) - min_len) <25):
      s1=s1[0:min_len]
      s2=s2[0:min_len]
      p=1
  else:
    i=l1-50
    sample1=s1[i:i+25]
    index=s2.find(sample1)
    if (index!=-1):
      if(index>i):
        if(index-i <=20):s2=s2[index-i:l2]
      else:
        if(i-index <=20):s1=s1[i-index:l1]
      min_len=min([len(s1),len(s2)])
      if ((max([len(s1),len(s2)]) - min_len) <25):
        s1=s1[0:min_len]
        s2=s2[0:min_len]
        p=1
      else:
        p=0
    else:
      p=0
  return (s1, s2, p)

def Write_output(out, file):
  fh=open(file, "a")
  fh.write(out)
  fh.close()
  return ()

def Get_cluster_similarities_single(seqs,coclust, cluster,file_out,inv):
  fh=open(file_out,"w")
  fh.close()
  ind=0
  total=len(cluster)
  t=0 
  for c in cluster:
    if (c not in inv):
      ind=ind+1
      clust_seqs=[]
      t=0 
      for id in cluster[c]:
        clust_seqs.append((id, seqs[id], len(seqs[id])))
        t=t+1
      for c1 in coclust[c]:
        for id in cluster[c1]:
          clust_seqs.append((id, seqs[id], len(seqs[id])))
          t=t+1
      clust_seqs=sorted(clust_seqs, key=itemgetter(2), reverse=True)
      if(len(clust_seqs)>1):
        if(len(clust_seqs)>500):print len(clust_seqs), total, ind
        Get_similarity_single(clust_seqs,file_out,c)
  return()

def Get_clusters (file):
  fh=open(file,"r")
  cluster=Tree()
  clust =''
  for l in fh:
    l=l.strip().split()
    if(l[0][0]==">"):
      clust = l[1]
    else:
      id = l[2].replace("...","").replace(">","")
      cluster[clust][id].value=1
  fh.close()
  return(cluster)

def Get_cluster_sizes_single (file):
  (cluster)= Get_clusters(file+".clstr")
  print len(cluster)
  sizes=[] 
  for c in cluster:
    tot=len(cluster[c])
    t=(c, tot)
    if (tot>50):
      sizes.append(t)
  s=sorted(sizes, key=itemgetter(1), reverse=True)
  return(s, cluster)

def Get_vaguely_similar_seqs(s1,s2,mis):
  l1=len(s1)
  l2=len(s2)
  trim=15
  seg_length = (l1-(2*trim))/mis
  p=0
  for i in range(0,mis):
    pos=trim+(i*seg_length)
    seg = s1[pos:pos+seg_length]
    index=s2.find(seg)
    if (index!=-1):
      if(index>pos):
        s2=s2[index-pos:l2]
      else:
        s1=s1[pos-index:l1]
      min_len=min([len(s1),len(s2)])
      s1=s1[0:min_len]
      s2=s2[0:min_len]
      p=1
      break
  return (s1, s2, p)

def Count_diffs (s1, s2, mis):
  i1=0
  i2=0
  mm=0
  p=1
  for i in range(0,len(s2)-1):
    if (s1[i1]==s2[i2]):
      i1=i1+1
      i2=i2+1
    else:
      if (s1[i1+1]==s2[i2]):
        i1=i1+1
        mm=mm+1
      else:
        if (s1[i1]==s2[i2+1]):
          i2=i2+1
          mm=mm+1
        else:
          mm=mm+1
    if (mm>mis):
      p=0
      break
  return (mm,p)

def Get_similar_clusters(s_sizes, cluster, seqs,tmp_file):
  coclust=Tree()
  inv={}
  mis = 5
  out = ''
  indw = 0
  fh1=open (tmp_file, "w")
  fh1.close()
  comp = 5
  for i in range(0,len(s_sizes)):
    print len(s_sizes), i
    clust=s_sizes[i]
    c1=clust[0]
    if (c1 not in inv):
      s1=''
      coclust[c1][c1].value=1
      inv[c1]=c1
      seqs1= []
      ind = 0
      for id in cluster[c1]:
        if (id in seqs):
          s1=seqs[id]
          seqs1.append(s1)
          ind = ind+1
          if(ind>comp):
            break
      for i2 in range(i+1,len(s_sizes)):
        clust2=s_sizes[i2]
        c2=clust2[0]
        ind = 0
        found = 0
        for id in cluster[c2]:
          if (id in seqs):
            ### Get difference with c1 sequences
            ### If sequence is different > X% then move onto next cluster
            ### If sequence is less different than Y% then co-cluster the clusters
            p=0
            ind = ind+1
            for s1 in seqs1:
              (s1, s2,p)=Get_vaguely_similar_seqs(s1, seqs[id], mis)
              if (p==1):
                (mm, q)=Count_diffs(s1, s2,mis)
                if (q==1):
                  out=out+c1+"\t"+c2+"\n"
                  indw = indw+1
                  if (indw >100):
                    Write_output(out, tmp_file)
                    out = ''
                    indw=0
                  found = 1
                  break
            if (ind > comp or found ==1):
              break
  Write_output(out, tmp_file)
  print out
  return()

def Get_coclustered (file):
  inv={}
  coclust = Tree()
  fh=open (file,"r")
  for l in fh:
    l=l.strip()
    l=l.split()
    inv[l[1]]=l[0]
    coclust[l[0]][l[1]].value = 1
  fh.close()
  return(inv, coclust)

def Decon_edges(att_file,file_seqs,file_edges):
  (inverse)=Get_inverse_ids(file_seqs)
  edges,edges23 = Tree(), Tree()
  fh1 = open (att_file , "r")
  for l in fh1:
    l=l.strip()
    l1=l.split()
    if (int(l1[0])==1 or int(l1[0])==2):
      edges[l1[1]][l1[2]].value=1
  fh1.close()
  Print_single_edges(file_edges, inverse, edges,tmp_file1)
  del inverse, edges, edges23
  return()

def Get_inverse_ids(file_seqs):
  fh = open (file_seqs, "r")
  inverse = {}
  ind =  0
  for l in fh:
    if (l[0]==">"):
      l=l.strip()
      l=l.replace(">","")
      l=l.split("|")
      inverse[l[0]]=l[1]
      ind = ind+1
  fh.close()
  return(inverse)
    
def Decon_identical (seq_file, att_file, file_vertex, file_seqs, read_number_division):
  fh = open (seq_file , "r")
  all,seqs={},{}
  for header,sequence in fasta_iterator(fh):
    seqs[header]=sequence
    all[header]=sequence
  fh.close()
  fh1 = open (att_file , "r")
  same1 = Tree()
  for l in fh1:
    l=l.strip()
    l1=l.split()
    cluster = l1[3]
    if (l1[0]=="0"):
      same1[cluster][l1[2]][l1[1]].value = 1
      same1[cluster][l1[1]][l1[2]].value =1 
  fh1.close()
  same,inverse, inv, out, ind, length = {},{},{},'',0,{}
  fh=open (file_seqs, "w")
  fh1.close()
  for c in same1:
    sub_same = Tree()
    for id1 in same1[c]:
      for id2 in same1[c][id1]:
        sub_same[id1][id2].value=1
        sub_same[id2][id1].value=1
    (sub_same, sub_inv)=Deconvolute_same_array (sub_same)
    for i in sub_same:
      s=seqs[i]
      total = 0
      mins=s
      for j in sub_same[i]:
        freq = j.split(read_number_division)
        total = total+ int(freq[1])
        inverse[j]=i
        out = out+">"+j+"|"+i+"\n"+seqs[j]+"\n"
        s=seqs[j]
        if(len(s)<len(mins)):
          mins=s
        ind = ind+1
        if(ind>100):
          Write_output(out, file_seqs)
          out = ""
          ind = 0
      same[i]=total
      length[i]=mins
  Write_output(out, file_seqs)
  del seqs
  Print_vertices(all, inverse, same, file_vertex, length,read_number_division)
  return()

def Print_single_edges(file_edges, inverse, edges,tmp_file):
  fh = open (tmp_file, "w")
  fh.close()
  edge,ind = '',0
  for id1 in edges:
    ida=id1 
    if (id1 in inverse):
      ida = inverse[id1]
    for id2 in edges[id1]:
      idb = id2
      if (id2 in inverse):
        idb = inverse[id2]
      if (ida!=idb):
        edge=edge+ida+"\t"+idb+"\t"+str(1)+"\t"+id1+"\t"+id2+"\n"
        ind = ind+1
        if (ind>300):
          Write_output(edge, tmp_file)
          edge = ''
          ind = 0
  Write_output(edge, tmp_file)
  del edges
  return()

def Print_vertices(all, inverse, same, file_vertex,length,read_number_division):
  out =''
  total = 0
  fh = open (file_vertex, "w")
  fh.close()
  ind =0
  for id in all:
    ind = ind+1
    if (id in inverse):
      if (id in same):
        if(id not in length):
          out = out +id+"\t"+str(same[id])+"\t"+all[id]+"\n"
        else:
          out = out +id+"\t"+str(same[id])+"\t"+length[id]+"\n"
    else:
      freq =id.split(read_number_division)
      out = out +id+"\t"+str(freq[1])+"\t"+all[id]+"\n"
    if (ind>300):
      Write_output(out, file_vertex)
      ind =0
      out = ''
  Write_output(out, file_vertex)
  return()

def Deconvolute_same_array (tree):
  decon = Tree()
  inv = {}
  index = 0
  for i in tree:
    if (i not in inv):
      index = index+1
      decon[i][i].value = 1
      inv[i]=i
      array = []
      for j in tree[i]:
        decon[i][j].value = 1
        inv[j]=i
        array.append(j)
      array_new=[]
      found = 0
      if (len(array)>0):
        found = 1
      while (found == 1):
        array_new =[]
        for k in array:
          for l in tree[k]:
            if (l not in inv):
              array_new.append(l)
              inv[l]=i
              decon[i][l].value = 1
        array = array_new
        if (len(array_new)==0):
          found = 0
  return (decon, inv)

def Reduce_edges(file_in, file_out):
  done,out, ind = Tree(),'',0
  fh = open (file_in, "r")
  fh1 = open (file_out, "w")
  fh1.close()
  for l in fh:
    l=l.strip()
    l=l.split()
    if (l[0] not in done[l[1]] and l[1] not in done[l[0]]):
      if (l[0]!=l[1]):
        out = out+l[0]+"\t"+l[1]+"\t"+l[2]+"\n"
        done[l[0]][l[1]].value =1
        ind = ind+1
        if (ind>300):
          ind =0
          Write_output(out, file_out)
          out=''
  Write_output(out, file_out)
  del done
  return()

def Deconvolute_edges(seq_file,  att_file, file_vertex, file_seqs, file_edges,read_number_division):
  Decon_identical (seq_file, att_file, file_vertex, file_seqs,read_number_division)
  Decon_edges(att_file,file_seqs,file_edges)
  Reduce_edges (tmp_file1, file_edges)
  return()

def Get_network_clusters(file_vertex, file_edges, cluster_file):
  (G,scale)=Read_graphical_inputs(file_vertex, file_edges)
  Output_cluster_file (G,cluster_file)
  return()

def Read_graphical_inputs(file_vertex, file_edges):
  fh = open (file_vertex, "r")
  freq = {}
  size=[]
  ind = 0
  G=nx.Graph()
  G.rtt={}
  for l in fh:
    l=l.strip()
    l=l.split()
    freq[l[0]]=int(l[1])
    size.append(int(l[1]))
    G.add_node(l[0])
    G.rtt[l[0]] = int(l[1])
    ind = ind+1
  scale1 = max(size)
  fh.close()
  fh = open (file_edges, "r")
  for l in fh:
    l=l.strip()
    l=l.split()
    if (l[0] in freq and l[1] in freq):
      G.add_edge(l[0],l[1])
  fh.close()
  return(G,scale1)

def Output_cluster_file (G, cluster_file):
  con= nx.connected_components(G)
  ind = 0
  ind1 = 0
  ind2 = 0
  fh = open(cluster_file, "w")
  fh.close()
  out = '# Connected_components\n'
  max_f,t,nvertmax = 0,0,0
  for i in con:
    ind = ind+1
    tc = 0
    nvert = 0
    for j in i:
      ind1=ind1+1
      ind2 = ind2+1
      out = out+str(ind1)+"\t"+str(ind)+"\t"+j+"\t"+str(G.rtt[j])+"\n"
      tc,t = tc+G.rtt[j],t+G.rtt[j]
      nvert = nvert+1
      if (ind2>100):
        Write_output(out, cluster_file)
        out = ""
        ind2=0
    if(tc>max_f):max_f,nvertmax = tc,nvert
  Write_output(out, cluster_file)
  print file_vertex,"Maximum cluster:",max_f*100.0/t,"%", nvertmax , "vertices"
  return()

def Get_network_input(file_cluster, outfile,edge_file):
  fh=open(outfile,"w")
  fh.close()
  fh=open(file_cluster, "r")
  cluster=Tree()
  ids = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cluster[l[1]][l[1]+"\t"+l[2]+"\t"+l[3]].value=1
      ids[l[2]]=1
  fh.close()
  (out, ind)=('',0)
  for c in cluster:
    if(len(cluster[c])>1):
      for l in cluster[c]:
        out=out+l+"\n"
        ind = ind+1
        if(int(l.split("\t")[2])>1000):print l
    else:
      for l in cluster[c]:
        if(int(l.split("\t")[2])>1):
          out=out+l+"\n"
          ind = ind+1
    if(ind>300):
      Write_output(out, outfile)
      (out, ind)=('',0)
  Write_output(out, outfile)
  return()

def Reduce_identical_sequences(Reduced_file,file_vertex,read_number_division):
  fh=open(Reduced_file,"w")
  fh.close()
  fh=open(file_vertex,"r")
  (ind, out) = (0,'')
  for l in fh:
    l=l.strip().split()
    out=out+">"+l[0].split(read_number_division)[0]+read_number_division+l[1]+"\n"+l[2]+"\n"
    ind = ind+1
    if(ind>500):
      Write_output(out, Reduced_file)
      (ind, out) = (0,'')
  fh.close()
  Write_output(out, Reduced_file)
  return()

def Check_same(ids, seqs_sub, same,inv):
  l=len(ids)
  ind = 0
  for i in range(0,l):
    for j in range(i,l):
      if(seqs_sub[i]==seqs_sub[j]):
        same [ids[i]][ids[j]].value = 1
        inv[ids[j]]=ids[i]
        ind = ind+1
      else:
        s1 = seqs_sub[i]
        s2 = seqs_sub[j]
        (sa,sb,p)=Trim_sequences(s1,s2,len(s1),len(s2))
        if(sa==sb):
          same [ids[i]][ids[j]].value = 1
          inv[ids[j]]=ids[i]
          ind = ind+1
    if(ind>l/2):break
  return(same,inv)

def Get_sequence_same_single(seqs, tmp_file1,reduced_sequences, tmp_pre, read_number_division):
  fh=open(reduced_sequences, "w")
  fh.close()
  (clust)= Get_clusters(tmp_file1+".clstr")
  (out, ind)=('',0)
  same = Tree()
  inv = {}
  for c in clust:
    if(len(clust[c])>1):
      ids=[]
      seqs_sub=[]
      for id in clust[c]:
        ids.append(id)
        seqs_sub.append(seqs[id])
      (same,inv)=Check_same(ids, seqs_sub,same,inv)
  for id in seqs:
    if(id in inv):
      if(id in same):
        f = 0
        for id1 in same[id]:
          id1 = id1.split(read_number_division)
          f = f + int(id1[len(id1)-1])
        out=out+">"+id.split(read_number_division)[0]+read_number_division+str(f)+"\n"+seqs[id]+"\n"
        ind = ind+1
    else:
      s = seqs[id]
      if(id.count(read_number_division)>1):
        id = id.split(read_number_division)
        id = id[0]+read_number_division+id[len(id)-1]
      out=out+">"+id+"\n"+s+"\n"
      ind = ind+1
    if(ind>100):
      Write_output(out, reduced_sequences)
      (out, ind)=('',0)
  Write_output(out, reduced_sequences)
  return()

def Generate_networks_pre(Sequence_file,tmp_pre,reduced_sequences,edge_lengths,read_number_division):
  Cluster_i(Sequence_file, tmp_pre,edge_lengths)
  (seqs)=Get_seqs_single (Sequence_file)
  Get_sequence_same_single(seqs, tmp_pre,reduced_sequences, tmp_pre,read_number_division)
  return()

def Generate_networks(Sequence_file, tmp_file1,edge_lengths, tmp_file, file_out):
  Cluster_i(Sequence_file, tmp_file1,edge_lengths)
  (s_sizes, cluster)=Get_cluster_sizes_single (tmp_file1)
  (seqs)=Get_seqs_single (Sequence_file)
  Get_similar_clusters(s_sizes, cluster, seqs,tmp_file+"_coclustered")
  (inv,coclust)=Get_coclustered (tmp_file+"_coclustered")
  Get_cluster_similarities_single(seqs,coclust, cluster,file_out,inv)
  return()

def Cleanup(tmp_reduced_sequences,tmp_pre,tmp_file1,tmp_file,att_file,file_seqs):
  for f in [tmp_reduced_sequences, tmp_pre, tmp_pre+".clstr",tmp_file1,tmp_file1+".clstr",tmp_file+"_coclustered",att_file,file_seqs]:
    os.system("rm "+f)
  return()

###########################
source = sys.argv[1]
id = sys.argv[2]
dir = sys.argv[3]
########################## Files for clustering
att_file = dir+"Vertex_relations_"+id+".txt"
file_vertex = dir+"Att_"+id+".txt"
file_edges = dir+"Edges_"+id+".txt"
cluster_file = dir+"Cluster_identities_"+id+".txt"
Reduced_file = dir+"Fully_reduced_"+id+".fasta"
######## checked_edges = dir+"Checked_edges_"+id+".txt"
plot_ids_file = dir+"Plot_ids_"+id+".txt"
file_seqs = dir+"Tmp_sequences_"+id+".txt"
tmp_reduced_sequences = dir+"Tmp_primary_reduced_sequences_"+id+".fasta"
tmp_pre = dir+"Tmp_Pre_tmp_"+id
tmp_file1 = dir+"Tmp_cluster_"+id+".1"
edge_lengths=0.85
tmp_file = dir+"Tmp_cluster_"+id+"."
read_number_division = "__" ## symbol indicating how the frequency of the reads is separated from the id

######################### Commands
          #### Clustering: 1. Reduce large datasets into unique sequences 
Generate_networks_pre(source,tmp_pre,tmp_reduced_sequences,1,read_number_division)
         #### Clustering: 2. Cluster related sequences together
Generate_networks(tmp_reduced_sequences, tmp_file1,edge_lengths, tmp_file, att_file)
Deconvolute_edges(tmp_reduced_sequences,  att_file, file_vertex, file_seqs, file_edges,read_number_division)
Get_network_clusters(file_vertex, file_edges, cluster_file)
Reduce_identical_sequences(Reduced_file,file_vertex,read_number_division)
Get_network_input(cluster_file,plot_ids_file,file_edges)
Cleanup(tmp_reduced_sequences,tmp_pre,tmp_file1,tmp_file,att_file,file_seqs)


