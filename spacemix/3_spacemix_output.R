library(SpaceMix)

prefix="SNP_sol_B_n415_LDprune"
run="92458"

prefix="SNP_sol_A_n109_LDprune"
run="383525"

setwd(paste0("/local/workdir/yc2644/Spisula/spacemix/grp/run_",run,
             "/",prefix,"_LongRun"))

spmix.data <- paste0(prefix,"_spacemix.data.Robj")
mcn.freq.list <- paste0(prefix,"_MCN.frequencies.list.Robj")
mcmc.output <- paste0(prefix, "_space_MCMC_output1.Robj")

load(spmix.data)
load(mcn.freq.list)
load(mcmc.output)

#check the data
str(MCN.frequencies.list)
str(spacemix.data)
ls()

###################### Evaluate MCMC Performance ##############################

#Plot MCMC trace 
plot(Prob,xlab="MCMC iterations",ylab="value",
     main="Posterior probability trace plot",type='l')
# If the chain is mixing well, the trace plot will resemble a 
# “fuzzy caterpillar.” If the trace plot has not plateaued, 
# it is an indication that the chain has not converged on 
# the stationary posterior distribution, and must be run longer. 

#### Plot Trace of Nugget parameters
matplot(t(nugget),type='l',
        xlab="MCMC iterations",ylab="Parameter value",
        main="Trace plot of nugget parameters")  
# If the trace plot of a parameter exhibits high autocorrelation, 
# the user may wish to “thin” the chain by decreasing the frequency 
# with which the chain is sampled

#Joint Marginal Plots
plot(a0,a1,xlab="a0",ylab="a1",
     main="Joint marginal of a0 and a1",pch=20,
     col=adjustcolor(rainbow(1000,start=4/6,end=6/6),0.3))
legend(x="bottomright",pch=19,cex=0.8,
       col=rainbow(1000,start=4/6,end=6/6)[c(1,500,1000)],
       legend=c("Sampled MCMC iteration 1",
                "Sampled MCMC iteration 500",
                "Sampled MCMC iteration 1000"))
# If the trace plot of a parameter exhibits high autocorrelation, 
# the user may wish to “thin” the chain by decreasing the frequency 
# with which the chain is sampled

#Acceptance Rate
plot(accept_rates$a0_accept_rate,
     xlab="MCMC iterations",ylab="Acceptance rate",
     main="Acceptance rate of a0",type='l',
     ylim=c(0,1))
abline(h=0.44,col="gray",lty=2)
# If the acceptance rate has not plateaued by the end of an analysis, 
# it may be an indication that the chain is still “going somewhere” 
# in parameter space, and the user may wish to perform subsequent analyses.

#Acceptance Rates of the nugget
matplot(t(accept_rates$nugget_accept_rate),
        xlab="MCMC iterations",ylab="Acceptance rate",
        main="Acceptance rates of nuggets",type='l',
        ylim=c(0,0.7))
abline(h=0.44,col="gray",lty=2)

######################Measurement of Model Adequacy ##############################
#calculate the sample covariance from the mean centered and normalized sample allele frequencies.
sample.covariance <- cov(t(MCN.frequencies.list$mean.centered.normalized.sample.frequencies),
                         use="pairwise.complete.obs")

k <- nrow(MCN.frequencies.list$mean.centered.normalized.sample.frequencies)
MC.matrix <- diag(k) - matrix(1/last.params$inv.mean.sample.sizes / 
                                (sum(1/last.params$inv.mean.sample.sizes)),
                              nrow=k,ncol=k,byrow=TRUE)

MC.parametric.covariance <- (MC.matrix) %*% last.params$admixed.covariance %*%  t(MC.matrix) #

index.matrix <- upper.tri(sample.covariance,diag=TRUE)

plot(sample.covariance[index.matrix], 
     MC.parametric.covariance[index.matrix],
     col=adjustcolor("black",0.3),pch=20,
     xlab="Sample covariance",
     ylab="Parametric covariance",
     main="Model adequacy:\n matrix comparison")
abline(0,1,col="red")    
# Finally, compare the sample covariance to the parametric
#   covariance.  Ideally, there will be a very tight correspondence 
#   between the data and the model.  If there is not, it may 
#   be an indication either that the MCMC has not converged on 
#   the stationary distribution or that the process that generated 
#   the data is only poorly approximated by SpaceMix's model.

plot(last.params$D[1:k,1:k][index.matrix], 
     sample.covariance[index.matrix],
     pch=19,col="black",
     xlab="geogenetic distance",
     ylab="covariance",
     main="Model adequacy:\n IBD patterns")
points(last.params$D[1:k,1:k][index.matrix], 
       MC.parametric.covariance[index.matrix],col="red",pch=20)
legend(x="topright",pch=19,col=c(1,2),
       legend=c("observed","model estimate"))

###################### Making geogenetic maps ##############################
sample.names <- unlist(sort(unique(read_tsv(paste0("../../",prefix,".grp_lon_lat"),col_names = T)$Region)))
coord <- read_tsv(paste0("../../coord.grp.",prefix),col_names = F)
coord <- matrix(unlist(coord), ncol =2, nrow=dim(coord)[1])

colors_A <- unlist(c("#366bfd", #SCC1
                   rep("#50ff35",4), #SLI1-4
                   rep("#a5e12f",2)) #SLI5-6
                   )

colors_B <- unlist(c(rep("#4ac0ff",3), #CCB1-3
                     rep("#f9731a",5), #DEL1-5
                     rep("#e20ffb",2), #GBE1-2
                     "#8505ff", #NAT
                     rep("#f4a261",2), #NJ1-2
                     rep("#4ac0ff",3), #PLY,PT1-2
                     "#366bfd", #SCC2
                     rep("#50ff35",3), #SLI2,3,4
                     rep("#a5e12f",2) #SLI5-6
                     )) #

spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = mcmc.output,
                                                    geographic.locations = coord,
                                                    name.vector = sample.names,
                                                    color.vector = colors_B,
                                                    quantile=0.95,
                                                    burnin=0)

make.spacemix.map(spacemix.map.list,
                          ellipses=T,
                  source.option=T,
                  text=TRUE,xlim = c(-75,-69), ylim=c(37.5,44))

plot.admix.arrows(admix.source.coords=spacemix.map.list$MAPP.admix.source.coords,
                  geogen.coords=spacemix.map.list$MAPP.geogen.coords,
                  admix.proportions = spacemix.map.list[["MCMC.output"]][["admix.proportions"]],
                  colors=NULL,
                  length=0.2)

arrows(x0 = spacemix.map.list$MAPP.admix.source.coords[,1],
       y0 = spacemix.map.list$MAPP.admix.source.coords[,2],
       x1 = spacemix.map.list$MAPP.geogen.coords[,1],
       y1 = spacemix.map.list$MAPP.geogen.coords[,2],
       lwd = spacemix.map.list[["MCMC.output"]][["admix.proportions"]]*100, #Increase the visibility
       col = spacemix.map.list$color.vector,
       length=0.2)

make.spacemix.map(spacemix.map.list,
                  ellipses=T,
                  source.option=T,
                  text=F,xlim = c(-74.5,-70.5), ylim=c(38.5,43))

plot.admix.arrows(admix.source.coords=spacemix.map.list$MAPP.admix.source.coords,
                  geogen.coords=spacemix.map.list$MAPP.geogen.coords,
                  admix.proportions = spacemix.map.list[["MCMC.output"]][["admix.proportions"]],
                  colors=NULL,
                  length=0.2)

arrows(x0 = spacemix.map.list$MAPP.admix.source.coords[,1],
       y0 = spacemix.map.list$MAPP.admix.source.coords[,2],
       x1 = spacemix.map.list$MAPP.geogen.coords[,1],
       y1 = spacemix.map.list$MAPP.geogen.coords[,2],
       lwd = spacemix.map.list[["MCMC.output"]][["admix.proportions"]]*100, #Increase the visibility
       col = spacemix.map.list$color.vector,
       length=0.2)


