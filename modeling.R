####---- required packages ----#####
library(DBI)
library(RMySQL)
library(gplots)
library(ROCR)
library(lattice)
library(ggplot2)
library(randomForest)
library(caret)
library(dplyr)
library(GPArotation)
library(tidyr)
library(reshape2)
library(dplyr)
library(tidyverse)
library(h2o)
library(sqldf)
library(foreach)
library(iterators)
library(parallel)

#### ---- load data ---- ####
### ---- start MySQL ---- #####
con <- dbConnect( MySQL(), host = "localhost", dbname = "cell_line", user = "root" , password = "root" )  ##

DEMETER_score <- dbGetQuery(con, "select * from DEMETER_score;" )
cell_line <- dbGetQuery( con, "select distinct(cell_line) from DEMETER_score;" )
expression_data <- dbGetQuery( con, "select * from rpkm_meandata;" ) 
mutation_data <- dbGetQuery( con, "select gene, cell_line, FATHMM_score from mutaion_data;" )
methylation <- dbGetQuery(con , "select a.gene, a.site, b.score, b.cell_line from meth_genes as a, meth_score as b  where a.probe = b.probe & b.pvalue <= 0.05;" )

dbDisconnect(con)  ## close MySQL

#### ---- data transformation ---- ####
expression_data <- data.frame(expression_data, expr="expr")
expression_data <- unite( expression_data, "gene", gene, expr, sep = "_" )
expression_data2 <- dcast( expression_data, cell_line~gene, value.var = "value", fun.aggregate = mean )

mutation_data <- data.frame( mutation_data, mut = 'mut' )
mutation_data <- unite( mutation_data, "gene", gene, mut, sep = "_")
mutation_data2 <- dcast( mutation_data, cell_line ~ gene, value.var = "FATHMM_score" )

methylation_data <- methylation_data[which(methylation_data$site != "Body"), ]
methylation_data <- methylation_data[which(methylation_data$site != "3UTR"), ]
methylation_data <- unite( methylation_data, "gene", gene, site, sep = "_" )
methylation_data2 <- dcast( methylation_data, cell_line ~ gene, value.var = "score" )

df <- merge(expression_data2, mutation_data2, by="cell_line", all.x = T )
df[ is.na(df) ] <- 0
df <- merge( df, methylation_data2, by="cell_line", all.x = T )
df[ is.na(df) ] <- 0.5

infile <- commandArgs(TRUE)

gene_name = "gene_name" 

#### ---- feature selection ---- ####
DEMETER_gene <- DEMETER_score[ which(DEMETER_score$gene == gene_name), c(2, 4) ]
dat <- merge( df, DEMETER_gene, by = "cell_line" )
dat1 <- dat[ which(dat$eg == "Y"), ]  
dat2 <- dat[ which(dat$eg == "N"), ]   
x <- sample( 1:nrow(dat2), nrow(dat1), replace = F ) 
dat3 <-  dat2[x, ] 
data <- rbind( dat1, dat3 )
 
for ( i in 2:( ncol(data)-1 ) ){
  data[ ,i] <- as.numeric( data[ ,i] )
}

data <- cbind( data[ ,which(colSums(data[2:(ncol(data)-1) ]) != 0)],data[,ncol(data)] )
data <- cbind( data[,which(colSums(data[2:(ncol(data)-1) ]) != nrow(data) * 235)], data[ ,ncol(data) ] )
names( data )[ ncol(data) ] <- "eg"

### --- data filter ---

zerovar=nearZeroVar( data[,2:(ncol(data)-1)] )
data2 <- data[ ,2:(ncol(data)-1) ][ ,-zerovar ]

data$eg <- as.factor(data$eg)

col_num <- ncol(data)
ctrl= rfeControl(functions = rfFuncs, method = "cv",number = 10)
results <- rfe(data2, data[ ,col_num ], sizes = c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000), rfeControl=ctrl) 

selected_features <- data.frame( results$fit$importance )
selected_features <- data[ ,c(colnames(data)[which(selected_features$MeanDecreaseGini >= 0) +1], "eg")]
selected_features <- cbind(data[ ,1], selected_features)
names(selected_features)[1] <- "cell_line"

selected_features <- cbind(selected_features[ ,which(colSums(selected_features[2:(ncol(selected_features)-1)]) >  0 )], selected_features[ ,ncol(selected_features)])
names(selected_features)[ncol(selected_features)] <- "eg"

write.table(selected_feature, "feature_path",col.names = T , row.names = F, quote = F, sep = "\t")

#### ---- machine learning ----
## --- start h2o ---
h2o.shutdown( prompt = F ); 
localH2O <- h2o.init(ip = "localhost", port = 50003, nthreads = 10,  startH2O = TRUE, max_mem_size = "100G")


dat5 <- selected_features[ ,-1]

### --- chose training set and testing set ---
dat5.train <- dat5 %>% sample_n( nrow(dat5) * 0.66 );
table(dat5.train$eg);

dat5.test <- subset( dat5, !rownames(dat5) %in% rownames( dat5.train ) );
table(dat5.test$eg);

dat5.train.h2o <- as.h2o(dat5.train)
dat5.train.h2o$eg <- as.factor(dat5.train.h2o$eg)

dat5.test.h2o  <-  as.h2o(dat5.test);


## multi-thread
library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

model.5 <-
  h2o.randomForest(x = colnames(dat5.train.h2o)[ colnames(dat5.train.h2o) != "eg" ],  # column IDs for the features
                   y = c("eg"),   # column ID for the label
                   training_frame = dat5.train.h2o, # data in H2O format
                   balance_classes = T);

stopCluster(cl) ## close multi-thread

## -- predict using training data --
pre.5 <- predict(object = model.5, newdata = dat5.train.h2o);
df.pre.5 <- as.data.frame(pre.5);

## -- performance of the prediction --
pre.5.eval <- data.frame( rel = dat5.train$eg, pre = df.pre.5$predict, good = dat5.train$eg == df.pre.5$predict );
( model.5.performance <- table( pre.5.eval[, c("rel", "good")] ) );
( perf_train <- round(model.5.performance / rowSums( model.5.performance ) * 100,3 ));

#############################################################
## --- predict using test data --
pre.5.test <- predict(object = model.5, newdata = dat5.test.h2o);
df.pre.5.test <- as.data.frame(pre.5.test);

## -- performance of the prediction --
pre.5.test.eval <- data.frame( rel = dat5.test$eg, pre = df.pre.5.test$predict, good = dat5.test$eg == df.pre.5.test$predict );
( model.pre.5.test.eval.performance <- table( pre.5.test.eval[, c("rel", "good")] ) );
( perf.5.test <- round( model.pre.5.test.eval.performance / rowSums( model.pre.5.test.eval.performance ) * 100 , 3) );


### --- performance evaluation --
perf.train <- h2o.performance(model.5, newdata = dat5.train.h2o);
perf.test <- h2o.performance(model.5, newdata = dat5.test.h2o);

dat.plot <-
  list( perf.train, perf.test ) %>%
  map(function(x) x  %>%
        # from all these 'paths' in the object
        .@metrics %>% .$thresholds_and_metric_scores %>%
        # extracting true positive rate and false positive rate
        .[c('tpr','fpr')] %>%
        # add (0,0) and (1,1) for the start and end point of ROC curve
        add_row(tpr=0,fpr=0,.before=T) %>%
        add_row(tpr=0,fpr=0,.before=F)) %>%
  # add a column of model name for future grouping in ggplot2
  map2(c('TRAIN', 'TEST'),
       function(x,y) x %>% add_column(Model=y)) %>%
  # reduce four data.frame to one
  reduce(rbind);

## -- plot --
pdf_name = paste( gene_name, ".pdf", sep= "" )
pdf(pdf_name)
p=ggplot(dat.plot, aes(fpr,tpr,col=Model))+
  geom_line( size = 1 )+
  geom_segment(aes(x=0,y=0,xend = 1, yend = 1),linetype = 2,col='grey')+
  xlab('False Positive Rate')+
  ylab('True Positive Rate')+
  ggtitle('ROC Curve for two models') +
  theme(legend.position="top") +
  scale_color_manual( values = c( "red", "blue"),
                      name = "Model", breaks = c("TRAIN", "TEST"),
                      labels = c( paste("on training data (AUC=", round( perf.train@metrics$AUC, 3), ")" ),
                                  paste("on test data (AUC=", round( perf.test@metrics$AUC, 3), ")" ))
  );

print(p)
dev.off()
### --- result 
features <- model.5@model$variable_importancescho
features <- features[ order( -features$percentage )[1:3], 1 ] ## choose top3 features

result_data <- list( gene_name,perf.5.test[1, -1], perf.5.test[2,-1], round( perf.train@metrics$AUC, 3),
                    round( perf.test@metrics$AUC, 3), ncol(dat5)-1, nrow(dat5), nrow( dat5[ which( dat5$eg == 'Y'), ] ), features[1],
                    features[2], features[3] )

names(result_data) <- c( "gene_name", "test_N", "test_Y", "train_AUC", "test_AUC", "features_num", "cell_line_num", "eg_num",
                        "feature1", "feature2", "feature3" )

result_data <- data.frame( result_data )

### -- save the predictive model ---
if ( round( perf.test@metrics$AUC, 3) >= 0.8 ){
    model_path <- h2o.saveModel(model.5,path = "model_path",force = FALSE) ### save models 
    model_name <- paste(pwd2,gene_name,".path",sep = "")
    write.csv(model_path,model_name,col.names = F,row.names = F)
  }
  features <- model.5@model$variable_importances
  features <- features[order(-features$percentage)[1:10],1]
  result_data <- list(gene_name,perf.5.test[1,-1],perf.5.test[2,-1],round( perf.train@metrics$AUC, 3),
                      round( perf.test@metrics$AUC, 3),ncol(dat5)-1,nrow(dat5),nrow(dat5[which(dat5$eg == 'Y'),]),features[1],features[2],
                      features[3],features[4],features[5],features[6],features[7],features[8],features[9],features[10])
  names(result_data) <- c("gene_name","test_N","test_Y","train_AUC","test_AUC","features_num","cell_line_num","eg_num",
                          "feature1","feature2","feature3","feature4","feature5","feature6","feature7","feature8","feature9","feature10")
  result_data <- data.frame(result_data)
  file_name = paste(pwd2,gene_name,".list",sep = "")
  write.table(result_data,file_name,col.names = T,row.names = F,quote = F,sep = "\t")
  
}
