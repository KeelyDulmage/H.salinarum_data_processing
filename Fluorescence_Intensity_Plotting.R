#Image J multiplot data must be loaded into R without a header.

norm_ave_plot3 = function(plot_data){
	par(mfrow=c(2,2))
	plot_data2=plot_data
	plot_data3=data.frame(matrix(nrow=501,ncol=dim(plot_data2)[2]/2))
	plot_data4=plot_data
	a_seq=seq(from=0, to=1, by=.002)
	
	for(i in 1:(length(plot_data[1,])/2)){
		cell_length= max(plot_data[,(i*2)-1], na.rm=T)
		plot_data2[,(i*2)-1] = plot_data[,(i*2)-1]/cell_length
		plot_data2[,(i*2)]= (plot_data[,(i*2)]/max(plot_data[,(i*2)], na.rm=T))
		}
		
		#approximate so that all cells have same X coordinates
	for(i in 1:(length(plot_data2[1,])/2)){
		new.coord=approx(plot_data2[,(i*2)-1],plot_data2[,(i*2)], xout=a_seq)
		plot_data3[,i] = new.coord$y
				}
	ave_y=rowMeans(plot_data3)
	x_val=a_seq
	plot(x_val,ave_y,ylim=c(0,1), xlim=c(0,1), type='l', lwd=3, xlab= 'Fraction of Cell Length', ylab='Fluorescence Intensity (A.U.)', cex.axis=1.5, cex.lab=1.5, col='blue')
	
	for(i in 2:(length(plot_data[1,])/2)){
	points(plot_data2[,(i*2)-1], plot_data2[,(i*2)], type='l', lwd=0.5, col='gray')
				}
				
	points(x_val,ave_y,ylim=c(0,1), xlim=c(0,1), type='l', lwd=3, cex.axis=1.5, cex.lab=1.5, col='blue')
	norm_data= data.frame(x_val, ave_y)
	
	#plot where max points are
	max_table=data.frame(matrix(nrow=length(plot_data3[,1]), ncol=2))
	colnames(max_table)=c('Position','Intensity')
	
	for(i in 1:length(plot_data3[1,])){
		max_i=max(plot_data3[,i])
		max_table[i,2]=max_i
		m_index=which(plot_data3[,i]== max_i)
		m_index2=m_index[1]
		max_table[i,1]=x_val[m_index2]
	}
	hist(max_table$Position, xlab='Fraction of cell length',cex.lab=1.5, cex.axis=1.5, col='gray', main='')
	
	#let's also plot max intensity versus length of cell.
	ml=data.frame(matrix(nrow=length(plot_data3[1,]), ncol=2))
	colnames(ml)=c('Cell.Length','Max.Intensity')
	for(i in 1:(length(plot_data3[1,]))){
		max_l=max(plot_data4[,(i*2)-1],na.rm=T)
		max_i=max(plot_data4[,(i*2)],na.rm=T)
		ml[i,1]=max_l
		ml[i,2]=max_i
		}
	plot(ml$Cell.Length, ml$Max.Intensity, xlab='Cell Length (microns)', ylab= 'Max Intensity (A.U)', pch=16, cex.lab=1.5, cex.axis=1.5)
	ml.lm=lm(ml$Max.Intensity~ml$Cell.Length)
	abline(ml.lm, col='red', lwd=1.5, lty=2)
	
	datalist=list(norm_data,max_table,ml)	
	datalist
}

# data can be accessed like this: object[[1]], etc.
