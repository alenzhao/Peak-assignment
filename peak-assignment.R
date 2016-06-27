ss=library(calibrate)
ss=library(base)
path<-"~/Dropbox/condmatbiophysics/R/data/"
file<-paste(path,"giant.csv",sep="")
s=read.csv(file, header=T)
#giant.csv is the file used to store all of the database files

#Read in the file as a data frame with first line as column names
s$Monoisotopic.mass[s$Monoisotopic.mass[]==-999999]<-NA
#Replace -999999 with NA

s$Molecular_weight[s$Molecular_weight[]==""]<-NA
#Replace blank with NA in Molecular_weight

ion_mass<-matrix(nrow=dim(s)[1],ncol=3)
colnames(ion_mass)<-c("[M+H]+","[M+K]+","[M+Na]+")
ion_mass_neg<-matrix(nrow=dim(s)[1],ncol=3)
colnames(ion_mass_neg)<-c("[M-H]-","[M+K-2H]-","[M+Na-2H]-")
#ion_mass is used to store ions in positive mode
ion_mass[,1]=s$Monoisotopic.mass[]+1.00782-0.000548
ion_mass[,2]=s$Monoisotopic.mass[]+38.96371-0.000548
ion_mass[,3]=s$Monoisotopic.mass[]+22.98980-0.000548
#ion_mass_neg is used to store ions in negative mode
ion_mass_neg[,1]=s$Monoisotopic.mass[]-1.00782+0.000548
ion_mass_neg[,2]=s$Monoisotopic.mass[]+38.96371-2*(1.00782)+0.000548
ion_mass_neg[,3]=s$Monoisotopic.mass[]+22.98980-2*(1.00782)+0.000548

file_pos<-paste(path,"peak_pos.csv",sep="")
file_neg<-paste(path,"peak_neg.csv",sep="")
pos=read.csv(file_pos)
neg=read.csv(file_neg)
#peak_pos.csv and peak_neg.csv are used to store peak values from positive mode and negative mode

output_pos <- c()
output_neg <- c()
peak_assign_pos<-c()
peak_assign_neg<-c()
pos_da<-c()
neg_da<-c()
sub_pos<-c()
sub_neg<-c()
peak_num_pos=0
x_pos<-c()
peak_num_neg=0
x_neg<-c()

for (i in 1:dim(pos)[1]){
t=abs(ion_mass[,1:3]-pos[i,1]*ones(dim(s)[1],3))
#Use ion mass minus the peak value, t is the abosolute value of the difference.
y<-which(t<0.055,arr.ind=T)
#return row number and column number of the matching ones
if (dim(y)[1]>0){
        result<-matrix(nrow=dim(y)[1],ncol=8)
        colnames(result)<-c("Compound name","Ions","Chemical_formula","Measured_mass","Caculated_mass","Mass_difference","Reaction_equation","Pathway")
        peak_num_pos=peak_num_pos+1
        x1_pos<-c(peak_num_pos,pos[i,1])
        x_pos<-rbind(x_pos,x1_pos,y)   
        #peak_num_pos is to strore the number of peaks actually found matching peaks
        #x_pos is to store the row and column value from giant.csv
                for (j in 1:dim(y)[1]){
                    result[j,1]<-toString(s$Compound_common_name[y[j,1]])
                    result[j,2]<-toString(colnames(ion_mass)[y[j,2]])
                    result[j,3]<-toString(s$Chemical_formula[y[j,1]]) 
                    result[j,4]=pos[i,1]
                    result[j,5]=ion_mass[y[j,1],y[j,2]]
                    result[j,6]=t[y[j,1],y[j,2]]
                    result[j,7]=toString(s$Reaction_equation[y[j,1]])
                    result[j,8]=toString(s$Pathway[y[j,1]])
                }
                output_pos = rbind(output_pos, result)
                rm(result)
            
}  
}
write.csv(output_pos,"output_pos.csv")
#output_pos.csv is all the possible matching metabolites
pos_da<-data.frame(output_pos)
peak_assign_pos<-unique(output_pos[!duplicated(pos_da$Compound.name),])
write.csv(peak_assign_pos,"peak_assign_pos.csv")
#peak_assign_pos.csv is all of the matching metabolites without redundancy
sub_pos<-unique(pos_da[c("Compound.name","Ions")])
write.csv(pos_da[row.names(sub_pos),],"peak_assign_pos1.csv")
#peak_assign_pos1.csv stores all of the matching metabolites with different possible ions
write.csv(x_pos,"assigned_peak_pos.csv")
#assigned_peak_pos.csv stores the assigned peaks and the matching row and column numbers in giant.csv
for (i in 1:dim(neg)[1]){
  t1=abs(ion_mass_neg[,1:3]-neg[i,1]*ones(dim(s)[1],3))
  y1<-which(t1<0.055,arr.ind=T)
  if (dim(y1)[1]>0){
    result<-matrix(nrow=dim(y1)[1],ncol=8)
    colnames(result)<-c("Compound name","Ions","Chemical_formula","Measured_mass","Caculated_mass","Mass_difference","Reaction_equation","Pathway")
    peak_num_neg=peak_num_neg+1
    x1_neg<-c(peak_num_neg,neg[i,1])
    x_neg<-rbind(x_neg,x1_neg,y1)
    for (j in 1:dim(y1)[1]){
      result[j,1]<-toString(s$Compound_common_name[y1[j,1]])
      result[j,2]<-toString(colnames(ion_mass_neg)[y1[j,2]])
      result[j,3]<-toString(s$Chemical_formula[y1[j,1]]) 
      result[j,4]=neg[i,1]
      result[j,5]=ion_mass_neg[y1[j,1],y1[j,2]]
      result[j,6]=t1[y1[j,1],y1[j,2]]
      result[j,7]=toString(s$Reaction_equation[y1[j,1]])
      result[j,8]=toString(s$Pathway[y1[j,1]])
    }
    output_neg = rbind(output_neg, result)
    rm(result)    
  }  
}
write.csv(output_neg,"output_neg.csv")
neg_da<-data.frame(output_neg)
peak_assign_neg<-unique(output_neg[!duplicated(neg_da$Compound.name),])
write.csv(peak_assign_neg,"peak_assign_neg.csv")
sub_neg<-unique(neg_da[c("Compound.name","Ions")])
write.csv(neg_da[row.names(sub_neg),],"peak_assign_neg1.csv")
write.csv(x_neg,"assigned_peak_neg.csv")