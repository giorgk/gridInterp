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
  writeData("test1_data_p.tmp", data = y, mode = "POINT", "LINEAR", axisFiles = c("Rgridinterp/test1_xaxis_const.tmp"))
  writeData("test1_data_c.tmp", data = y[1:length(y)-1], mode = "CELL", "LINEAR", axisFiles = c("Rgridinterp/test1_xaxis_const.tmp"))
}


{ # Read the results
  tstpl <- read.table("res_cnst_test1D_pl.tmp", header = F, col.names = c("X","Y"))
  tstpn <- read.table("res_cnst_test1D_pn.tmp", header = F, col.names = c("X","Y"))
  tstcl <- read.table("res_cnst_test1D_cl.tmp", header = F, col.names = c("X","Y"))
  tstcn <- read.table("res_cnst_test1D_cn.tmp", header = F, col.names = c("X","Y"))
}


{ # plot n compare
  p <- plot_ly()
  p <- add_trace(p, x = x, y = y , type = 'scatter', mode = 'lines')
  p <- add_trace(p, x = tstpl$X, y = tstpl$Y , type = 'scatter', mode = 'lines')
  p <- add_trace(p, x = tstpn$X, y = tstpn$Y , type = 'scatter', mode = 'lines')
  p <- add_trace(p, x = tstcl$X, y = tstcl$Y , type = 'scatter', mode = 'lines')
  p <- add_trace(p, x = tstcn$X, y = tstcn$Y , type = 'scatter', mode = 'lines')
  p
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
  writeData("test1_data_var_p.tmp", data = y, mode = "POINT", "LINEAR", axisFiles = c("Rgridinterp/test1_xaxis_var.tmp"))
  writeData("test1_data_var_c.tmp", data = y[1:length(y)-1], mode = "CELL", "LINEAR", axisFiles = c("Rgridinterp/test1_xaxis_var.tmp"))
}

{ # Read the results
  tstpl <- read.table("res_var_test1D_pl.tmp", header = F, col.names = c("X","Y"))
  tstpn <- read.table("res_var_test1D_pn.tmp", header = F, col.names = c("X","Y"))
  tstcl <- read.table("res_var_test1D_cl.tmp", header = F, col.names = c("X","Y"))
  tstcn <- read.table("res_var_test1D_cn.tmp", header = F, col.names = c("X","Y"))
}

{
  p <- plot_ly()
  p <- add_trace(p, x = xvar, y = y , type = 'scatter', mode = 'lines+markers')
  p <- add_trace(p, x = tstpl$X, y = tstpl$Y , type = 'scatter', mode = 'lines')
  p <- add_trace(p, x = tstpn$X, y = tstpn$Y , type = 'scatter', mode = 'lines')
  p <- add_trace(p, x = tstcl$X, y = tstcl$Y , type = 'scatter', mode = 'lines')
  p <- add_trace(p, x = tstcn$X, y = tstcn$Y , type = 'scatter', mode = 'lines')
  p
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
  writeData("peak_data_cnst_p.tmp", data = P$Z, mode = "POINT", "LINEAR", 
            axisFiles = c("Rgridinterp/peakAxis_cnst.tmp", "Rgridinterp/peakAxis_cnst.tmp"))
}

{# read the results
  peakspl <- read.table("res_cnst_peak_pl.tmp", header = F, col.names = c("X", "Y", "Z"))
  peakspn <- read.table("res_cnst_peak_pn.tmp", header = F, col.names = c("X", "Y", "Z"))
  
}

{
  p <- plot_ly()
  p <- add_surface(p, x=x, y = y, z = P$Z)
  p <- add_markers(p, x = peakspl$X, y = peakspl$Y, z = peakspl$Z, marker = list(size = 1))
  p <- add_markers(p, x = peakspn$X, y = peakspn$Y, z = peakspn$Z, marker = list(size = 1))
  p
  
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
  writeData("test2D_data_var_p.tmp", data = z, mode = "POINT", "LINEAR", 
            axisFiles = c("Rgridinterp/test2D_Xvar.tmp", "Rgridinterp/test2D_Yvar.tmp"))
  
  writeData("test2D_data_var_c.tmp", data = z[1:19,1:37], mode = "CELL", "LINEAR", 
            axisFiles = c("Rgridinterp/test2D_Xvar.tmp", "Rgridinterp/test2D_Yvar.tmp"))
}

{# read the results
  test2pl <- read.table("res_var_test2D_pl.tmp", header = F, col.names = c("X", "Y", "Z"))
  test2pn <- read.table("res_var_test2D_pn.tmp", header = F, col.names = c("X", "Y", "Z"))
  test2cl <- read.table("res_var_test2D_cl.tmp", header = F, col.names = c("X", "Y", "Z"))
  test2cn <- read.table("res_var_test2D_cn.tmp", header = F, col.names = c("X", "Y", "Z"))
  
}


{
  p <- plot_ly()
  p <- add_surface(p, x = xvar, y = yvar, z = z)
  p <- add_markers(p, x = test2pl$X, y = test2pl$Y, z = test2pl$Z, marker = list(size = 1))
  p <- add_markers(p, x = test2pn$X, y = test2pn$Y, z = test2pn$Z, marker = list(size = 1))
  p <- add_markers(p, x = test2cl$X, y = test2cl$Y, z = test2cl$Z, marker = list(size = 1))
  p <- add_markers(p, x = test2cn$X, y = test2cn$Y, z = test2cn$Z, marker = list(size = 1))
  p
  
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

writeData <- function(filename, data, mode, method, axisFiles){
  if (is.null(dim(data))){
    write(paste(mode, method, length(data)), file = filename, append = F)
    write(axisFiles, file = filename, append = T, ncolumns = length(axisFiles))
    write.table(data, file = filename, append = T, row.names = F, col.names = F)
  }
  else{
    if (length(dim(data)) == 2){
      write(paste(mode, method, dim(data)[2], dim(data)[1]), file = filename, append = F)
      write(axisFiles, file = filename, append = T, ncolumns = length(axisFiles))
      write.table(data, file = filename, append = T, row.names = F, col.names = F)
    }
    else if (length(dim(data)) == 3){
      write(paste(mode, method, dim(data)[2], dim(data)[1], dim(data)[3]), file = filename, append = F)
      write(axisFiles, file = filename, append = T, ncolumns = length(axisFiles))
      for (i in 1:dim(data)[3]) {
        write.table(data[,,i], file = filename, append = T, row.names = F, col.names = F)
      }
    }
  }
}


  writegridInterp1D <- function(filename, x, spacing, v = NA, mode = "POINT", method = "LINEAR"){
    if (is.na(v)){
      has_v = 0  
    }
    else{
      has_v = 1
    }
    write(paste(mode, method, spacing, has_v), file = filename, append = F)
    if (spacing == "CONST"){
      N = x[3]
      write(x, file = filename, append = T)
      if (!is.na(v[1])){
        write.table(v, file = filename, append = T, row.names = F, col.names = F)
      }
    }
    else{
      if (mode == "POINT"){
        N = length(x)
        write(N, file = filename, append = T)
        if (!is.na(v[1])){
          xv <- cbind(x,v)
        }
        else{
          xv <- x
        }
        write.table(xv, file = filename, append = T, row.names = F, col.names = F)
      }
      else if (mode == "CELL"){
        N = length(x) - 1
        write(N, file = filename, append = T)
        if (!is.na(v[1])){
          xv <- cbind(x[1:N],v)
        }
        else{
          xv <- x[1:N]
        }
        write.table(xv, file = filename, append = T, row.names = F, col.names = F)
        write(x[N+1], file = filename, append = T)
      }
    }
  }
  
  
  writegridInterp2D <- function(filename, xaxis, yaxis, v, mode = "POINT", method = "LINEAR"){
    write(paste(mode, method, dim(v)[1], dim(v)[2]), file = filename, append = F)
    write(paste("XAXIS", xaxis), file = filename, append = T)
    write(paste("YAXIS", yaxis), file = filename, append = T)
    write.table(v, file = filename, append = T, row.names = F, col.names = F)
  }
  
  