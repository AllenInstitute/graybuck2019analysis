map_sampling2 <- function (train.dat, train.cl, test.dat, markers, markers.perc = 0.8, 
          iter = 100, method = "median") 
{
  map.result = sapply(1:iter, function(i) {
    print(paste0("Running round ",i))
    tmp.markers = sample(markers, round(length(markers) * 
                                          markers.perc))
    map_by_cor(train.dat[tmp.markers, ], train.cl, test.dat[tmp.markers, 
                                                            ], method = method)
  }, simplify = F)
  map.cl = sapply(map.result, function(x) x$pred.df$pred.cl)
  row.names(map.cl) = colnames(test.dat)
  map = as.data.frame(as.table(as.matrix(map.cl)))
  map.freq <- table(map$Var1, map$Freq)
  map.df = data.frame(pred.cl = setNames(colnames(map.freq)[apply(map.freq, 
                                                                  1, which.max)], row.names(map.freq)), prob = rowMaxs(map.freq)/iter)
  return(list(map.df = map.df, map.freq = map.freq))
}
