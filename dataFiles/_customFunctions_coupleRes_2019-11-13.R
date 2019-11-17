
###################################################################
### SELECT PCs according to Broken Stick Model (MacArthur 1957) ###
###################################################################

# Takes output from plotTangentSpace or prcomp as argument
# Returns selected PCs and saves number of selected PCs in variable called "n.PCs.[argument name]"
# Adapted from function "evplot" (Francois Gillet, http://adn.biol.umontreal.ca/~numericalecology/numecolR/)

selectPCs <- function(PCA.output) {
  ev<-PCA.output$sdev^2
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  
  test<-cbind(100*ev/sum(ev), bsm$p[n:1])
  n.PCs<-sum(test[,1] >= test[,2])
  
  # Save number of principal components
  arg_name <- deparse(substitute(PCA.output))
  var_name <- paste("n.PCs", arg_name, sep=".")
  assign(var_name, n.PCs, .GlobalEnv)
  
  # Print PCs and variance explained
  if (!is.null(PCA.output$pc.summary$importance)) {
    return(PCA.output$pc.summary$importance[,1:n.PCs])
  } else {
    temp<-summary(PCA.output)
    return(temp$importance[,1:n.PCs])
  }
}

##########################
### SPLIT VIOLIN PLOTS ###
##########################

# from https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})
geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

###################################
### RAINCLOUD/FLAT VIOLIN PLOTS ###
###################################

# from https://micahallen.org/2018/03/15/introducing-raincloud-plots/

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

raincloud_theme = theme(
  text = element_text(size = 10),
  axis.text = element_text(size = 10),
  axis.title.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "none",
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')
)

#####################################
### Calculate Procrustes distance ###
#####################################

# Input arguments:
# - coords.matrix: matrix of shape coordinates, first column needs to be ID
# - pairs.table: table that has pairs of IDs as rows (column1 = individual 1, column 2 = individual 2)

calcProcD <- function(coords.matrix, pairs.table) {
  
  R <- nrow(pairs.table)
  procDist <- numeric(R)
  
  for (i in 1:R) {
    
    # find index of line in coords.matrix that has same ID as individual on ith line of pairs.table
    ind1.id<-which(rownames(coords.matrix) %in% (pairs.table[i, 1]))
    ind2.id<-which(rownames(coords.matrix) %in% (pairs.table[i, 2]))

    # extract the coordinates for that particular woman/man
    ind1 <- as.matrix(coords.matrix[ind1.id,])
    ind2 <- as.matrix(coords.matrix[ind2.id,])
    
    # calculate Euclidean distance, save as ith element in proc.dist vector
    procDist[i] <- sqrt(sum((ind1 - ind2)**2))
  }
  
  procDist.table<-tibble(
    "ind1" = pairs.table[, 1],
    "ind2" = pairs.table[, 2],
    "procDist" = procDist
  )
  
  return(procDist.table)
  
}

###############################################################
### Create k x r cross-validated principal component scores ###
###############################################################

cv.prcomp<-function(coords.matrix, output.path, k = 10, r = 10){
  # output.path: string without trailing "/"
  # k = number of folds
  # r = number of repeats
  
  # Create matrix listing all possible combinations of k-1 folds
  fold.combs <- combn(seq_len(k), k - 1)
  
  # Empty dataframe to store number of PCs extracted in each resample
  extractedPCs <- data_frame("extractedPCs" = numeric(k*r)) 
  
  # Counter for inner for-loop
  counter <- 0
  
  # Create output directory if it does not exist yet
  ifelse(!dir.exists(output.path), dir.create(output.path), FALSE)
  
  # Load data on sex of faces [to keep number of men and women in folds balanced]
  sex<-read_csv("./dataFiles/sexOfFaces_anon.csv", col_types = list("c","c")) %>%
    as.data.frame() %>%
    arrange(id)
  
  for (j in 1:r) {
    
    # Set seed for each repeat
    set.seed(100 + j)
    
    # Create the k folds for each repeat and split data
    folds <- createFolds(sex$sex, k)
    splits <- lapply(folds, function(ind, dat) dat[ind,], dat = coords.matrix)
    
    # For each repeat, re-run inner loop for all possible combinations of k-1 training folds
    for (i in 1:ncol(fold.combs)) {
      
      # bind training folds into matrix
      temp.data <- do.call(rbind, splits[fold.combs[,i]])
      
      # run PCA on training folds
      PCA <- prcomp(temp.data)
      
      # select PCs
      selectPCs(PCA)
      
      # save number of extracted PCs
      extractedPCs[counter+i,]<-n.PCs.PCA
      
      # predict and save PC scores for *all* observations
      predict(PCA,coords.matrix)[,1:n.PCs.PCA] %>% 
        as.data.frame() %>%
        write_csv(paste0(output.path,"/PCScores_repeat_",j,".fold_", i,".csv"))
    }
    counter<-counter+k
  }
  #write_csv(extractedPCs,paste0(output.path,"/numberOfExtractedPCs.csv"))
}

##########################################
### Re-create 3D mesh from point cloud ###
##########################################

# Using package alphashape3d
# Adapted from nat::as.mesh3d.ashape3d (Gregory S X E Jefferis and James D Manton (2014). NeuroAnatomy Toolbox v1.5.2. ZENODO. 10.5281/zenodo.10171)

convertPointsToMesh <- function(point.matrix, alpha, showme = TRUE) {
  
  shape.alpha <- ashape3d(point.matrix, alpha = alpha)
  selrows <- shape.alpha$triang[,9] >= 2
  tr <- shape.alpha$triang[selrows, c("tr1", "tr2", "tr3")]
  
  mesh<-rgl::tmesh3d(
    vertices = t(shape.alpha$x),
    indices = t(tr),
    homogeneous = FALSE
  )
  
  if (showme == TRUE){
    plot3d(mesh, type="wire", col="gray", box=FALSE, axes=FALSE, xlab="", ylab="", zlab="")
  }
  
  return(mesh)
}

####################################
### Plot 3D principal components ###
####################################

# Adapted from geomorph::plotTangentSpace

make3Dpcs <- function(pca.output, numberOfPCs, visSD, mean_shape, plot=TRUE, return=TRUE) {
  
  # pca.output = prcomp output
  # numberOfPCs = how many PCs to visualize
  # visSD = how many SDs to visualize, e.g. 3
  # mean_shape = reference coordinates which to manipulate
  # plot or not.
  
  shapes <- shape.names <- NULL
  
  for (i in 1:(numberOfPCs)) {
    
    multiplier<-pca.output$sdev[i]*visSD
    
    minus<-as.vector(t(mean_shape))-as.matrix(multiplier * (t(pca.output$rotation[,i])))
    plus<-as.vector(t(mean_shape))+as.matrix(multiplier * (t(pca.output$rotation[,i])))
    
    shapes <- rbind(shapes, minus, plus)
    shape.names <- c(shape.names,
                     paste("PC", i, "_minus",visSD,"SD",sep=""),
                     paste("PC", i, "_plus",visSD,"SD",sep="")
    )
  }
  shapes <- arrayspecs(shapes, dim(mean_shape)[1], dim(mean_shape)[2])
  shapes <- lapply(seq(dim(shapes)[3]), function(x) shapes[, , x])
  names(shapes)<-shape.names
  
  mesh.list<-list()
  
  for (i in 1:length(shapes)){
    mesh.list[[i]]<-convertPointsToMesh(shapes[[names(shapes)[i]]],1.2,FALSE)
  }
  
  names(mesh.list)<-shape.names
  
  if (plot==TRUE){
    open3d()
    mfrow3d(numberOfPCs, 2, sharedMouse = TRUE)
    for (i in 1:length(shapes)){
      plot3d(mesh.list[[i]],type="wire",col="gray",box=FALSE,axes=FALSE,xlab="",ylab="",zlab="",specular="black")
    }
    rglwidget()
  }
  
  if (return==TRUE){
  return(mesh.list)
  }
}

###################################
### CUSTOM SUMMARY FOR LMERTEST ###
###################################
lmer.summary <- function(model, statsname = "stats") {
  
  ci<-confint(model, method = "Wald") %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    filter(!is.na(`2.5 %`))
  
  lmer.summary<-summary(model)
  
  # tidy up summary output table
  coefs <- lmer.summary$coefficients %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    left_join(ci,by="rowname") %>%
    rename(
      effect = rowname,
      estimate = "Estimate",
      p = "Pr(>|t|)",
      SE = "Std. Error",
      t= "t value"
    ) %>%
    select(1:2,7:8,3:6) %>%  
    mutate_at(vars(estimate:t), round, 2) %>%
    mutate(p=format(round(p,3),nsmall=3),
           effect=ifelse(str_detect(effect,":"),gsub(':','\\_',effect),effect)
    )
  
  if (ncol(coefs)>5) {
    coefs <- coefs %>%
      mutate(
        sig = ifelse(p<.001, "***",
                     ifelse(p<.01,   "**",
                            ifelse(p<.05,   "*",
                                   ifelse(p<.10,   "+", "")))))
  }
  
  
  # export individual rows as object for access in markdown
  stats <- list()
  for (i in 1:nrow(coefs)){
    s <- coefs[i,]
    p <- ifelse(s$p < .001, "p<.001", paste0("p=", s$p))
    ef <- str_replace_all(coefs$effect[i], ":", "_")
    stats[ef] = paste0("(estimate=",  s$estimate, 
                       ", 95% CI = [", s$'2.5 %', ", ", s$'97.5 %', 
                       "], ", p, ")" )
    
  }
  assign(statsname, stats, envir = .GlobalEnv)
  
  coefs <- coefs %>% kable() %>% kable_styling()
  
  # print tidied table
  return(list(lmer.summary$ngrps, coefs))
}


##############################
### CUSTOM SUMMARY FOR LM ###
#############################

custom_lm.summary <- function(model, statsname = "stats") {
  
  ci<-confint(model, method = "Wald") %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    filter(!is.na(`2.5 %`))
  
  lm.summary<-summary(model)
  
  # tidy up summary output table
  coefs <- lm.summary$coefficients %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    left_join(ci,by="rowname") %>%
    rename(
      effect = rowname,
      estimate = "Estimate",
      p = "Pr(>|t|)",
      SE = "Std. Error",
      t= "t value"
    ) %>%
    select(1:2,6:7,3:5) %>%  
    mutate_at(vars(estimate:t), round, 2) %>%
    mutate(p=format(round(p,3),nsmall=3),
           effect=ifelse(str_detect(effect,":"),gsub(':','\\_',effect),effect)
    )
  
  if (ncol(coefs)>5) {
    coefs <- coefs %>%
      mutate(
        sig = ifelse(p<.001, "***",
                     ifelse(p<.01,   "**",
                            ifelse(p<.05,   "*",
                                   ifelse(p<.10,   "+", "")))))
  }

  # print tidied table
  return(coefs %>% kable() %>% kable_styling())
}
