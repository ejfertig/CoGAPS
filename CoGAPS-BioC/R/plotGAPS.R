# plotGAPS: function to plot the decomposed A and P matrices
# History: EJF - original CoGAPS 

# Inputs: A - A matrix
#         P - P matrix
#         outputPDF - string for PDF output, if null uses screen

# Output: PDF file if indicated, otherwise none

plotGAPS <- function(A, P, outputPDF="") {
	if (outputPDF != "") {
	  pdf(file=paste(outputPDF,"_Patterns",".pdf",sep=""))
	} else {
	  dev.new()
	}

	arrayIdx <- 1:ncol(P)
	matplot(arrayIdx, t(P), type='l', lwd=10, xlim=c(1,ncol(P)), ylim=c(0,1))
	title(main='Inferred patterns')
  
	if (outputPDF == "") {
	  dev.new()
	} else {
	  dev.off()
	  pdf(file=paste(outputPDF,"_Amplitude",".pdf",sep=""))
	}

	heatmap(A, Rowv=NA, Colv=NA) 
	
	if (outputPDF != "") {
		dev.off()
	}
  
}
