import os
from Bio import SeqIO
import subprocess
import pandas as pd
from collections import Counter
def extract_prot_name(path):
 protls=[]
 for record in SeqIO.parse(path+"/raw_data/clbrener.gb", "genbank"):
        #protein=(record.description).split("[Trypanosoma cruzi]")
        if "RecName: Full=" in record.description:
            prot=(record.description).split(";")
            print record.organisum
            protname=prot[0].lstrip("RecName: Full=")
            protls.append(protname)
        elif "[Trypanosoma cruzi strain CL Brener]" in record.description:
            prot=(record.description).split("[Trypanosoma cruzi strain CL Brener]") 
            protls.append(prot[0])
        else:
            protls.append(record.description)

 with open(path+"/clbrener_gene.txt","w") as outfl:
     for item in protls:
       outfl.write(item+"\n")
def remove_hypothetical(path):
    with open("clbrener_rmv_hypo.txt","w") as rmv_hypo:
      with open(path+"/clbrener_gene.txt","r") as read_fil:
          for line in read_fil:
             line=line.strip()
             if "hypothetical protein" not in line:
                rmv_hypo.write(line+"\n")
def remove_putative_partial(path):
    create_key=open("clbrener_putative.txt","w")
    create_key1=open("clbrener_partial.txt","w") 
    with open("clbrener_unique.txt","w") as rmv_hypo:
      with open(path+"/clbrener_rmv_hypo.txt","r") as read_fil:
          protls=[]
          protein=[]
         
          for line in read_fil:
             
             line=line.strip()
             protein.append(line)
             put=0
             put2=0
             
             par=0
            
             if ", putative" in line.lower(): 
                put+=1 
                create_key.write(line+"\n") 
             if "putative" in line.lower():
                put2+=1
                create_key.write(line+"\n") 
             if ", partial" in line.lower():  
                par+=1 
                create_key1.write(line+"\n")  
             if put==1:
               line=(line.lower()).replace(", putative","")
               protls.append(line.strip()) 
             if put2==1:
               line=(line.lower()).replace("putative ","")
               protls.append(line.strip()) 
             
               line=(line.lower()).replace(", partial","")
               protls.append(line.strip()) 
            
             if put==0 and put2==0 and par==0 :
               protls.append(line.strip())               
      
             
          df=pd.DataFrame()
          df1=pd.DataFrame()
          count={}
          prot=[]
          Freq=[]
          count=Counter(protls)
          for key,val in count.items():
             prot.append(key)
             Freq.append(val)
          df["Protein"]=prot
          df["Frequency"]=Freq
          df.to_csv("clbrener.txt",sep="\t",index=False) 
          prot_seq={} 
          prot1=[]
          Freq1=[]
          prot_seq=Counter(protein)
          for key,val in prot_seq.items():
             prot1.append(key)
             Freq1.append(val)
          df1["Protein"]=prot1
          df1["Frequency"]=Freq1
          df1.to_csv("clbrener_partial_putative.txt",sep="\t",index=False)              
          for item in set(protls):
                rmv_hypo.write(item+"\n")
          

                
if __name__=='__main__':
    path=os.getcwd()
    extract_prot_name(path)
    remove_hypothetical(path)
    remove_putative_partial(path) 
   
    
   
