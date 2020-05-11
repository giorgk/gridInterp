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
  
  xvar = 38.425 + cumsum(runif(n = 38,min = 0.3, max = 2))
  writeAxis("test1_xaxis_var.tmp", xvar, "VAR")
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
  x = 38.425 + cumsum(runif(n = 35,min = 0.3, max = 2))
  y <- sin(x/5.2)^3
}

{
  writegridInterp1D(filename = "test1D_var_pl.tmp", x = x, v = y, 
                    spacing = "VAR", mode = "POINT", method = "LINEAR")
  writegridInterp1D(filename = "test1D_var_pn.tmp", x = x, v = y, 
                    spacing = "VAR", mode = "POINT", method = "NEAREST")
  # for CELL MODE the length of x is length(v) + 1
  xx <- c(x, x[length(x)] + runif(n = 1,min = 0.3, max = 2) )
  writegridInterp1D(filename = "test1D_var_cl.tmp", x = xx, v = y, 
                    spacing = "VAR", mode = "CELL", method = "LINEAR")
  writegridInterp1D(filename = "test1D_var_cn.tmp", x = xx, v = y, 
                    spacing = "VAR", mode = "CELL", method = "NEAREST")
}

{ # Read the results
  tstpl <- read.table("res_var_test1D_pl.tmp", header = F, col.names = c("X","Y"))
  tstpn <- read.table("res_var_test1D_pn.tmp", header = F, col.names = c("X","Y"))
  tstcl <- read.table("res_var_test1D_cl.tmp", header = F, col.names = c("X","Y"))
  tstcn <- read.table("res_var_test1D_cn.tmp", header = F, col.names = c("X","Y"))
}

{
  p <- plot_ly()
  p <- add_trace(p, x = tstpl$X, y = tstpl$Y , type = 'scatter', mode = 'lines')
  p <- add_trace(p, x = tstpn$X, y = tstpn$Y , type = 'scatter', mode = 'lines')
  p <- add_trace(p, x = tstcl$X, y = tstcl$Y , type = 'scatter', mode = 'lines')
  p <- add_trace(p, x = tstcn$X, y = tstcn$Y , type = 'scatter', mode = 'lines')
  p <- add_trace(p, x = x, y = y , type = 'scatter', mode = 'lines+markers')
  p
}

#=============== TEST 2D With constant X, Y ===========
{# test function
  P <- peaks()
  x <- P$X[1, ]
  y <- P$Y[, 1]
}

{ # write files
  # write the axis files
  writegridInterp1D("XaxisTest2D.tmp", c(-3, 0.125, 49), "CONST", mode = "POINT", method = "LINEAR")
  writegridInterp2D("Test2D_pl.tmp", xaxis = "XaxisTest2D.tmp", yaxis = "XaxisTest2D.tmp", v = P$Z,
                    mode = "POINT", method = "LINEAR")
  
}


{
  p <- plot_ly()
  p <- add_surface(p, x=x, y = y, z = P$Z)
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
    write(axisFiles, file = filename, append = T)
    write.table(data, file = filename, append = T, row.names = F, col.names = F)
  }
  else{
    if (length(dim(data)) == 2){
      write(paste(mode, method, dim(data)[1], dim(data)[2]), file = filename, append = F)
      write(axisFiles, file = filename, append = T)
      write.table(data, file = filename, append = T, row.names = F, col.names = F)
    }
    else if (length(dim(data)) == 3){
      write(paste(mode, method, dim(data)[1], dim(data)[2], dim(data)[3]), file = filename, append = F)
      write(axisFiles, file = filename, append = T)
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
  
  