# from https://gist.githubusercontent.com/mvarela/e3b868fe36d232a39bf5d87796c8ac28/raw/b4c4d31f562179573f58b7d8a48ef6d9286aba32/binVis2d.r

library(tidyverse)

# binviz Veles-like binary visualizaiton
binViz2d <- function(filename, alpha = 1/100, maxsize = 5000000,
                    save = TRUE, polar = FALSE, sample = FALSE,
                    sample_size = 2000000, do_density = FALSE){

  # setting dens_plot as NA simplifies the logic below a bit
  dens_plot = NA

  # we read the file as a stream of bytes, and prepare our tibble
  # We'll add a column indexing the trigram position in the file
  # This will come in handy later if we want to facet the plot by position
  # as done in the Veles article. We'll just mutate binViz here, to save memory.
  rawdata          <- readBin(filename, integer(), n=maxsize, size = 1, signed = FALSE)
  size             <- rawdata %>% as.tibble %>% nrow
  binViz           <- cbind(0:(size - 1),rawdata, lead(rawdata), lead(rawdata,n=2L))
  colnames(binViz) <- c('idx', 'x', 'y', 'z')
  

  # We then remove any missing values from the dataset
  toplot <- binViz %>% as.tibble %>% na.omit

  # If sampling is required, we do it now. Sampling is important
  # if doing the density plots, as going beyond 1M points gets SLOW
  if(sample){
   toplot <- toplot %>% sample_n(min(count(toplot), sample_size))
  }

  # The actual plotting
  theplot <- binViz2d_do_plot(toplot, alpha, polar) +
    ggtitle(title_spec(filename, sample, sample_size))
  if(do_density){
    dens_plot <- binViz2d_do_density_plot(toplot, polar)
  }

  # Saving the plots
  if(save){
    namespec <- name_spec(filename, sample, sample_size, polar)
    binViz2d_save(namespec, theplot, dens_plot)
  }
  return(list(binViz_plot = theplot, dens_plot = dens_plot))
}

binViz2d_do_plot <- function(data, alpha, polar){
  theplot <- data  %>% ggplot(mapping = aes(x,y)) +
    geom_point(mapping = aes(color=z), alpha = alpha, size = 0.75) +
    scale_color_gradient(low="blue", high="orange") +
    coord_fixed(ratio = 1)+
    labs(x="i", y="i+1", z="i+2")

  if(polar){
    theplot <- theplot + coord_polar()
  }
  return(theplot)
}

binViz2d_do_density_plot <- function(toplot, polar){
  dens_plot <- toplot  %>% ggplot(mapping = aes(x,y)) +
    stat_density2d(aes(fill = ..density..), geom="raster", contour = FALSE) +
    scale_fill_gradient(low="steelblue4", high="sienna2") +
    coord_fixed(ratio = 1)+
    labs(x="i", y="i+1")
  return(dens_plot)
}

title_spec <- function(name, sampled, nsamples){
  if(sampled){
    title <- paste(name, "-", nsamples, "samples.")
  }else{
    title <- name
  }
  return(title)
}

# We create a name separated by underscores, this simplifies later parsing
# of file names, if needed, to automate e.g., reports creation
name_spec <- function(name, sampled, nsamples, polar){
  polar_str     <- ""
  if(polar){
    polar_str   <- "polar"
  }
  sampled_str   <- ""
  if(sampled){
    sampled_str <- paste("sampled", nsamples, sep="_")
  }
  basename <- chartr('/.', '::',
                     paste("plot", polar_str, sampled_str, name, sep = "_"))
  return(paste(basename, ".png", sep=""))
}

binViz2d_save <- function(namespec, binViz_plot, dens_plot){
    png(namespec, width = 15, height = 15, units = "cm", res = 300)
    print(binViz_plot)
    dev.off()
    if(!is.na(dens_plot)){
      png(paste("density",namespec,sep="_"), width = 15, height = 15,
          units = "cm", res = 300)
      print(dens_plot)
      dev.off()
    }
}
