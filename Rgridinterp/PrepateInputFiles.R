library(plotly)
library(pracma)
#=============== TEST 1D With constant X ===========

{ # 1D interpolations
  # Test function
  x <- seq(40.234,64.38,0.643)
  y <- sin(x/5.2)^3
}

{# write Axis
  orig <- 40.234
  dx <- 0.643
  N <- 38
  writeAxis("test1_xaxis_const.tmp", c(orig, dx, N), "CONST")
}

{# write data
  writeData("test1_cnst_data.tmp", data = y, "LINEAR", axisFiles = c("Rgridinterp/test1_xaxis_const.tmp"))
}


{ # Read the results
  test_l <- read.table("res_cnst_test1_l.tmp", header = F, col.names = c("X","Y"))
  test_n <- read.table("res_cnst_test1_n.tmp", header = F, col.names = c("X","Y"))
}


{ # plot n compare
  p <- plot_ly()
  p <- add_trace(p, x = x, y = y , type = 'scatter', mode = 'lines', name = 'Data')
  p <- add_trace(p, x = test_l$X, y = test_l$Y , type = 'scatter', mode = 'lines', name = "Linear")
  p <- add_trace(p, x = test_n$X, y = test_n$Y , type = 'scatter', mode = 'lines', name = "Nearest")
  p %>%
    layout(title = "Evenly spaced data")
}

#=============== TEST 1D With variable X ===========
{# Test function
  xvar = 38.425 + cumsum(runif(n = 38,min = 0.3, max = 2))
  y <- sin(x/5.2)^3
}

{# write Axis
  writeAxis("test1_xaxis_var.tmp", xvar, "VAR")
}

{# write data
  writeData("test1_var_data.tmp", data = y, "LINEAR", axisFiles = c("Rgridinterp/test1_xaxis_var.tmp"))
}

{ # Read the results
  test_l <- read.table("res_var_test1_l.tmp", header = F, col.names = c("X","Y"))
  test_n <- read.table("res_var_test1_n.tmp", header = F, col.names = c("X","Y"))
}

{
  p <- plot_ly()
  p <- add_trace(p, x = xvar, y = y , type = 'scatter', mode = 'lines+markers', name = "Data")
  p <- add_trace(p, x = test_l$X, y = test_l$Y , type = 'scatter', mode = 'lines', name = "Linear")
  p <- add_trace(p, x = test_n$X, y = test_n$Y , type = 'scatter', mode = 'lines', name = "Nearest")
  p %>%
    layout(title = "Variably spaced data")
}

#=============== TEST 2D With constant X, Y ===========
{# test function
  N = 20
  P <- peaks(v = N)
  x <- 10*P$X[1, ]
  y <- 10*P$Y[, 1]
}

{ # write axis
  orig <- -30
  dx <- diff(x)[1]
  # For the peak example we will use the same ticks in both axis so we print one file
  writeAxis("peakAxis_cnst.tmp", c(orig, dx, N), "CONST")
}



{ # write data
  writeData("peak_data_cnst.tmp", data = P$Z, "LINEAR", 
            axisFiles = c("Rgridinterp/peakAxis_cnst.tmp", "Rgridinterp/peakAxis_cnst.tmp"))
}

{# read the results
  peaks_l <- read.table("res_cnst_peak_l.tmp", header = F, col.names = c("X", "Y", "Z"))
  peaks_n <- read.table("res_cnst_peak_n.tmp", header = F, col.names = c("X", "Y", "Z"))
  
}

{
  p <- plot_ly()
  p <- add_surface(p, x=x, y = y, z = P$Z, name = 'Data',showscale = FALSE)
  p <- add_markers(p, x = peaks_l$X, y = peaks_l$Y, z = peaks_l$Z, marker = list(size = 1), name = 'Linear')
  p <- add_markers(p, x = peaks_n$X, y = peaks_n$Y, z = peaks_n$Z, marker = list(size = 1), name = 'Nearest')
  p %>%
    layout(title = "Evenly spaced data")
    
  
}

#=============== TEST 2D With variable X, Y ===========
{# Test function
  xvar <- -2 + cumsum(runif(n = 38,min = 0.03, max = 0.17))
  yvar <- -1 + cumsum(runif(n = 20,min = 0.03, max = 0.15))
  XYgrid <- meshgrid(xvar, yvar)
  z <- (4 - 2.1*XYgrid$X^2 + (XYgrid$X^4)/3)*XYgrid$X^2 + XYgrid$X*XYgrid$Y + (-4 + 4*XYgrid$Y^2)*XYgrid$Y^2
  
  xvar <- 10*xvar
  yvar <- 10*yvar
}

{ # write axis n data
  writeAxis("test2D_Xvar.tmp", xvar, "VAR")
  writeAxis("test2D_Yvar.tmp", yvar, "VAR")
  writeData("test2D_data_var.tmp", data = z, "LINEAR", 
            axisFiles = c("Rgridinterp/test2D_Xvar.tmp", "Rgridinterp/test2D_Yvar.tmp"))
}

{# read the results
  test2_l <- read.table("res_var_test2D_l.tmp", header = F, col.names = c("X", "Y", "Z"))
  test2_n <- read.table("res_var_test2D_n.tmp", header = F, col.names = c("X", "Y", "Z"))
  
}


{
  p <- plot_ly()
  p <- add_surface(p, x = xvar, y = yvar, z = z, name = 'Data')
  p <- add_markers(p, x = test2_l$X, y = test2_l$Y, z = test2_l$Z, marker = list(size = 1), name = 'Linear')
  p <- add_markers(p, x = test2_n$X, y = test2_n$Y, z = test2_n$Z, marker = list(size = 1), name = 'Nearest')
  p %>%
    layout(title = "Variably spaced data")
  
}
{# logo image
  ax <- list(title = "", zeroline = FALSE, showline = F, showticklabels = F, showgrid = F)
  p <- plot_ly()
  p <- add_markers(p, x = test2_l$X, y = test2_l$Y, z = test2_l$Z, marker = list(size = 1))
  p <- add_markers(p, x = test2_n$X, y = test2_n$Y, z = test2_n$Z, marker = list(size = 1))
  p %>%
    layout(scene = list(xaxis = ax, yaxis = ax, zaxis = ax),showlegend = FALSE)
}


######### WRITE FUNCTIONS =======================

writeAxis <- function(filename, x, spacing){
  if (spacing == "CONST"){
    N <- x[3]
  }
  else{
    N <- length(x)
  }
  write(paste(spacing,N), file = filename, append = F)
  if (spacing == "CONST"){
    write(x[1:2], file = filename, append = T)
  }
  else{
    write.table(x, file = filename, append = T, row.names = F, col.names = F)
  }
}

writeData <- function(filename, data,method, axisFiles){
  if (is.null(dim(data))){
    write(paste(method, length(data)), file = filename, append = F)
    write(axisFiles, file = filename, append = T, ncolumns = length(axisFiles))
    write.table(data, file = filename, append = T, row.names = F, col.names = F)
  }
  else{
    if (length(dim(data)) == 2){
      write(paste(method, dim(data)[2], dim(data)[1]), file = filename, append = F)
      write(axisFiles, file = filename, append = T, ncolumns = length(axisFiles))
      write.table(data, file = filename, append = T, row.names = F, col.names = F)
    }
    else if (length(dim(data)) == 3){
      write(paste(method, dim(data)[2], dim(data)[1], dim(data)[3]), file = filename, append = F)
      write(axisFiles, file = filename, append = T, ncolumns = length(axisFiles))
      for (i in 1:dim(data)[3]) {
        write.table(data[,,i], file = filename, append = T, row.names = F, col.names = F)
      }
    }
  }
}
  