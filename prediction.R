#### ---- required package ---- ####
library(h2o)

#### ---- load data ---- ### 

df <- read.table("all_features.txt", header = T,sep = "\t")

#### ---- start h2o ----
h2o.shutdown( prompt = F ); 
localH2O <- h2o.init(ip = "localhost", port = 50013, nthreads = 10,  startH2O = TRUE, max_mem_size = "100G")

gene_name = "gene_name" ### input the name from the saved models 

model_file_name <- paste("/mnt/raid5/haoxw/machine-learning/randomForest/eg_50/", gene_name, ".path", sep = "")
model_name <- read.csv(model_file_name, header = T)
model_path <- as.character(model_name[1,1])
  
feature_file <- paste("feature_path",gene_name, ".features", sep = "")
features <- read.table(feature_file, header = T, sep = "\t")
feature_names <- colnames(features)
  
selected_features <- df[,which(colnames(df) %in% feature_names)]

### --- load the model ---
loaded_model <- h2o.loadModel(model_path)
  
test_result <- predict(loaded_model, newdata = as.h2o(selected_features))
result_performance <- cbind(selected_features[,1],as.data.frame(test_result))
names(result_performance)[1] <- "cell_line"
write.table(result_performance,paste("predicton_path",gene_name,".txt", sep = ""), col.names = T, row.names = F, quote = F,sep = "\t")
