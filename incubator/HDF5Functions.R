#'Write a COMPASSContainer to an HDF5 file
#'
#'This function implement h5write from the rhdf5 package 
#'to support writing of COMPASSContainer objects to HDF5 files.
#'
#'@param x is a \code{COMPASSContainer} object to be written to hdf5
#'@param filename is a \code{character} name of the file to write to.
#'@param group is a \code{character} parent group name / path to the dataset
#'in the HDF5 file where the object will be written. If it doesn't exist, it will be created.
#'@param overwrite is a \code{boolean} specifying whether to overwrite the file. Defaults to \code{TRUE}
#'
#'@details h5write.COMPASSContainer will write the object \code{x} to file at the group defined by \code{group}.
#'The function will overwrite any existing file with name \code{filename} by default because deleting and replacing HDF5
#'groups and datasets is apparently not trivially supported by rhdf5 at the moment. 
#'
#'@note Need to implement adding to an existing HDF5 file under a new group name. For example if we have 
#'a set of CD4 matrices written in an HDF5 file, we may want to add CD8 matrices to the same file
#'rather than creating a new HDF5 file.
#'@importFrom rhdf5 h5write
#'@S3method h5write COMPASSContainer
h5write.COMPASSContainer<-function(x=NA,filename=NA,group=NA,overwrite=TRUE){
  if(any(c(is.na(x),is.na(filename),is.na(group)))){
    stop("All arguments must be specified")
  }
  #remove if overwrite = TRUE
  if(overwrite&file.exists(filename)){
    file.remove(filename)
  }
  
  #Create the h5 file if it doesn't exist
  if(!file.exists(filename)){
    h5createFile(filename)
  }
  #TODO handle the case where overwrite=FALSE but the file exists
  #How do we check if a group exists?
  
  top <- strsplit(group,"/")[[1]][1]
  
  loc <- H5Fopen(filename)
  group_exists <- !H5Lexists(h5loc=loc,name=group)
  H5Fclose(loc)

  if(group_exists){
    h5createGroup(filename,top)
    h5createGroup(filename,group)
    #list the empty matrices and fill them with a zero row so that we can store them in hdf5
    #also store the names of the zero samples
    zeros<-names(which(do.call(c,lapply(x$data,function(y)nrow(y)==0))))
    h5write(zeros,file=filename,sprintf("%s/empty_matrices",group))
    for(i in (zeros)){
      x$data[[i]]<-matrix(0,nrow=1,ncol=ncol(x$data[[i]]),dimnames=dimnames(x$data[[i]]))
    }
    h5write(x$data,file=filename,sprintf("%s/data",group))
    h5write(as.data.frame(x$meta),file=filename,sprintf("%s/metadata",group))
    h5write(x$counts,file=filename,sprintf("%s/counts",group))
    h5write(x$individual_id,file=filename,sprintf("%s/individual_id",group))
    h5write(x$sample_id,file=filename,sprintf("%s/sample_id",group))
  }else{
    message("Group exists, can't replace, you must set overwrite=TRUE.\n")
  }
}
