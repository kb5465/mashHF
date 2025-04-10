library(lattice)
library(colorRamps)

# Compute the correlations from covariance matrix
k <- 1
covmat <- U_ed
x <- cov2cor(covmat[[k]])
colnames(x) <- names
rownames(x) <- names

# Trait names
names  <- gsub("BETA_", "", colnames(Bhat))

# Indices of outcome traits (all in this case)
h <- rep(TRUE, length(names))

# Generate heatmap of Uk covariance matrix
clrs <- colorRampPalette(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                               "#E0F3F8","#91BFDB","#4575B4"))(64)
lat <- x[rev(h),rev(h)]
lat[lower.tri(lat)] <- NA
n <- nrow(lat)
print(levelplot(lat[n:1, ],col.regions = clrs, xlab = "", ylab = "", main = paste0("Covariance Matrix Sigma_", k),
      colorkey = TRUE, scales=list(x=list(rot=60))))

# Plot eigenvectors capturing predominant patterns
for(g in 1:3){
  v <- svd(covmat[[k]])$v[h,]
  d <- svd(covmat[[k]])$d
  par(mar=c(8,4.1,4.1,2.1))
  barplot(v[,g] / v[which.max(abs(v[ ,g])), g], las = 2,
          main=paste("Eigenvector",g,"of Sigma_",k),
          cex.names = 0.5,names=names[h])
}
