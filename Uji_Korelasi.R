IHS=read.table("E:/1/Data & Syntax/Data_IHS_JII.txt",header=TRUE)
x1=IHS[,3]
x2=IHS[,4]
y1=IHS[,1]
y2=IHS[,2]
n=length(x1)
##Hasil Uji Korelasi Pearson##
r=cor(y1,y2, method = c("pearson"))
thitung<-(r*sqrt(n-2))/sqrt(1-(r^2))
v<-n-2
ttabel<-qt(1-(0.05/2),v)
print(r)
print(thitung)
print(ttabel)

