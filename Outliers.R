##############
# libraries

library(robustbase)
library(astro)
#library(mclust)
library(cosmoFns)
#library(lawstat)
#library(weights)
library(nortest)
library(outliers)
library(matrixStats)

######################

#setwd("/")
#setwd("/home/silas/Mestrado/Projeto/Remoção de Outliers/Gap fixo e variavel")
options(warn=-1)

######################
cspeed=299792.458
H0=100
pi=3.141593
G=4.3*10^(-9.0)
####################

#cls=read.table("ask78.txt",head=T)

##############################

###### initial cutoff ###########

Vmax<-function(ra, dec, redshift){
	
	cls<-data.frame(ra, dec, redshift)
	zgal=cls$redshift

	vgal=cspeed*(((zgal+1)^2 -1)/((zgal+1)^2 +1))

	mgal=median(vgal)

	vgmax=mgal + 4000
	vgmin=mgal - 4000

	cls=cbind(cls,vgal)

	#cls=subset(cls, cls$vgal <= vgmax & cls$vgal >= vgmin)

	######### finding the field baricenter  #########

	RA=cls$ra
	DEC=cls$dec
	Z_SPEC=cls$redshift

	CenRA=median(RA)
	CenDEC=median(DEC)


	#######################################

	MedianZ=median(Z_SPEC)
	sigma_r=1E-9
	sigma_v=50

	q = 25

	Vel=cspeed*(((Z_SPEC+1)^2 -1)/((Z_SPEC+1)^2 +1))

	vall=median(Vel)

	sigall=s_Qn(Vel)

	vlos=(Vel - vall)/(1+MedianZ)


	# Work out projected distance from group centre for group and nearby galaxies. Use CosDist for this.
	# hc is defined as 0.15 *(AngDist/320 h^-1 Mpc) in Serra et al. 2011.

		groupdist = CosDist(MedianZ,H0=100,OmegaM=0.3,OmegaL=0.7)[1,'AngDist']
		hc = 0.15*(groupdist/320)

			xcen = groupdist*cos(CenRA*pi/180)*cos(CenDEC*pi/180)
			ycen = groupdist*sin(CenRA*pi/180)*cos(CenDEC*pi/180)
			zcen = groupdist*sin(CenDEC*pi/180)

		projdist = {}
		xgal = {}
		ygal = {}
		zgal = {}
			xgal = groupdist*cos(RA*pi/180)*cos(DEC*pi/180)
			ygal = groupdist*sin(RA*pi/180)*cos(DEC*pi/180)
			zgal = groupdist*sin(DEC*pi/180)
			dist = sqrt((xgal-xcen)^2 + (ygal-ycen)^2 + (zgal-zcen)^2)

			dx=xgal-xcen
			dy=ygal-ycen
			dz=zgal-zcen


	############# projdist ################


	theta=acos(sin(DEC*pi/180)*sin(CenDEC*pi/180) + cos(DEC*pi/180)*cos(CenDEC*pi/180))*cos((RA - CenRA)*pi/180)
	projdist=groupdist*theta


	###########################################

	ratio=projdist/dist

	phi=asin(ratio)
	phi[!is.finite(phi)] <- 1.570796

	#########################

	vec=data.frame(cbind(cls,Vel,vlos,dist,projdist,phi,dx,dy,dz))

	vec=vec[order(vec$projdist),]


	######################## first method ---> VMAX ###################


	##########################  vmax ####################

	cspeed=299792.458
	H0=100
	pi=3.141593
	G=4.3*10^(-9.0)


	H=H0*(sqrt(0.3*(1+MedianZ)^3 + 0.7))

	GG=4.3*10^(-9.0)

	NN=nrow(vec)

	rr=vector()
	mr=vector()
	vc=vector()
	vi=vector()
	vm=vector()
	sig=vector()
	rp=vector()
	vl=vector()
	class1=vector()
	nphi=vector()

	clo=0

	ij=0

	lim=NN

	for (ii in 1:lim){

	ij=ij+1

	temp=vec[ii:NN,]

	ngal=nrow(temp)

	rp[ij]=vec$projdist[ii]
	vl[ij]=vec$vlos[ii]

	sig[ij]=s_Qn(temp$vlos)
	nphi[ij]=vec$phi[ii]

	if (sig[ij]==0) sig[ij]=mean(sig[ij -1],sig[ij -2],sig[ij -3])

	rr[ij]=sqrt(3)*sig[ij]/(10*H)

	mr[ij]=log10((3*rr[ij]*sig[ij]^2/GG))

	vc[ij]=sqrt((GG*10^mr[ij])/rr[ij])

	vi[ij]=sqrt(2)*vc[ij]

	vm[ij]=max(vc[ij]*sin(nphi[ij]),vi[ij]*cos(nphi[ij]))

	if (abs(vl[ij]) >= vm[ij]) clo=1
	if (abs(vl[ij]) < vm[ij]) clo=0

	class1[ij]=clo

	}


	### class points out a member or not ###########

	nvec=cbind(vec,vl,vm,class1)
	#############################################

	out1=subset(nvec, nvec$class1==1)
	return (out1)

}


	############## Second method ---> ########################
	#cls=cbind(cls,vlos)
	#cls=cbind(cls,projdist)
	#cls<-cbind(cls,flag=sample(T,length(cls$vlos), replace=TRUE))




############ shiftting gapper -->  fixed gap ############

GapFix<-function(RA, DEC, Z_SPEC, projdist, vlos){
	cls<-data.frame(RA, DEC, Z_SPEC, projdist, vlos)
	cls<-data.frame(cls,flag=sample(T,length(cls$vlos), replace=TRUE))
	cont2<-2
	cspeed=299792.458
	H0=100
	pi=3.141593
	G=4.3*10^(-9.0)
	while(cont2!=0){
		cls<-cls[order(cls$projdist),]
		clsbin<-c()
		tam<-0
		maxprojdist<-0
		cont<-0
		cont2<-0
		cls<-cls[which(cls$flag == T),]
		while(length(cls$vlos)>tam){
			cont=cont+1
			print(maxprojdist)
			clsbin[[cont]]<-cls[which(cls$projdist <= (0.4+maxprojdist) & cls$projdist > maxprojdist),]
			if(length(clsbin[[cont]]$projdist)<15){
				aux<-length(clsbin[[cont]]$projdist)
				if((tam+15)<=length(cls$vlos)){
					clsbin[[cont]]<-rbind(clsbin[[cont]], cls[(tam+aux+1):(tam+15),])
				}else{
					if((tam+aux)<length(cls$vlos)){
						clsbin[[cont]]<-rbind(clsbin[[cont]], cls[(tam+aux+1):length(cls$vlos),])
					}
				}
				maxprojdist<-max(clsbin[[cont]]$projdist)
			}else{
				maxprojdist<-maxprojdist+0.4
			}
			tam<-tam+length(clsbin[[cont]]$projdist)
			clsbin[[cont]]<-clsbin[[cont]][order(clsbin[[cont]]$vlos),]
		}
		cls<-c()
		for(i in 1:cont){
			for(j in 1:(length(clsbin[[i]]$vlos)-1)){
				if(length(clsbin[[i]]$vlos)>1){
				
					if( abs((clsbin[[i]][j+1,]$vlos - clsbin[[i]][j,]$vlos))<350){
						clsbin[[i]][j+1,]$flag<-T
					}else{
						if(clsbin[[i]][j,]$vlos<0){
							clsbin[[i]][j,]$flag<-F
						}else{
							clsbin[[i]][j+1,]$flag<-F
						}
						cont2<-cont2+1
					}
				}
			}
			cls<-rbind(cls,clsbin[[i]])
		}
		print(cont2)
		#print(summary(clsbin[[i]]))
	}

	clsgapfix<-cls
	#plot(clsante$projdist,clsante$vlos)
	#points(clsgapfix$projdist,clsgapfix$vlos,phc=20,col="green")
	return (clsgapfix)
}





############## variable gap ############################
GapVar<-function(RA, DEC, Z_SPEC, projdist, vlos){
	cls<-data.frame(RA, DEC, Z_SPEC, projdist, vlos)
	cls<-cbind(cls,flag=sample(T,length(cls$vlos), replace=TRUE))
	cont2<-2
	cspeed=299792.458
	H0=100
	pi=3.141593
	G=4.3*10^(-9.0)
	while(cont2!=0){
		cls<-cls[order(cls$projdist),]
		clsbin<-c()
		tam<-0
		maxprojdist<-0
		cont<-0
		cont2<-0
		cls<-cls[which(cls$flag == T),]
		while(length(cls$vlos)>tam){
			print(maxprojdist)
			cont=cont+1
			clsbin[[cont]]<-cls[which(cls$projdist <= (0.4+maxprojdist) & cls$projdist > maxprojdist),]
			if(length(clsbin[[cont]]$projdist)<15){
				aux<-length(clsbin[[cont]]$projdist)
				if((tam+15)<=length(cls$vlos)){
					clsbin[[cont]]<-rbind(clsbin[[cont]], cls[(tam+aux+1):(tam+15),])
				}else{
					if((tam+aux)<length(cls$vlos)){
						clsbin[[cont]]<-rbind(clsbin[[cont]], cls[(tam+aux+1):length(cls$vlos),])
					}
				}
				maxprojdist<-max(clsbin[[cont]]$projdist)
			}else{
				maxprojdist<-maxprojdist+0.4
			}
			tam<-tam+length(clsbin[[cont]]$projdist)
			clsbin[[cont]]<-clsbin[[cont]][order(clsbin[[cont]]$vlos),]
		}
		cls<-c()
		for(i in 1:cont){
			quantis=data.frame(quantile(clsbin[[i]]$vlos,probs=c(0.25,0.75)))
			fl=quantis[1,]
			fu=quantis[2,]
			sf = (fu - fl)/1.349

			print(sf)
			for(j in 1:(length(clsbin[[i]]$vlos)-1)){
				if(length(clsbin[[i]]$vlos)>1){
					if( ((clsbin[[i]][j+1,]$vlos - clsbin[[i]][j,]$vlos))<sf){
						clsbin[[i]][j+1,]$flag<-T
					}else{
						if(clsbin[[i]][j,]$vlos<0){
							clsbin[[i]][j,]$flag<-F
						}else{
							clsbin[[i]][j+1,]$flag<-F
						}
						cont2<-cont2+1
					}
				}
			}
			cls<-rbind(cls,clsbin[[i]])
		}
		print("A")
		print(cont2)
	}

	clsgapvar<-cls
	#points(clsgapvar$projdist,clsgapvar$vlos,phc=20,col="green")

	#length(clsbin[[1]]$vlos)+length(clsbin[[2]]$vlos)+length(clsbin[[3]]$vlos)+length(clsbin[[4]]$vlos)+length(clsbin[[5]]$vlos)
	return (clsgapvar)
}
