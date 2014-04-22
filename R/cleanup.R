## This function:
## 1. Concatenates all .c files into <pkg>_c.c,
## 2. Concatenates all .cpp files into <pkg>_cpp.cpp,
## 3. But leaves RcppExports.cpp and <pkg>_init.c unmodified.
## This greatly speeds compilation time.
cleanup <- function() {

  ## Bail if this is an R CMD INSTALL
  if (Sys.getenv("R_INSTALL_PKG") != "") {
    message("It looks like this isn't being called from R CMD build; bailing out")
    return(invisible(NULL))
  }

  files <- list.files()
  stopifnot("DESCRIPTION" %in% files)

  ## get the package name from the DESCRIPTION file
  DESCRIPTION <- as.list( read.dcf("DESCRIPTION")[1, ] )
  pkg_name <- DESCRIPTION$Package
  pkg_version <- DESCRIPTION$Version

  ## copy the files to a 'build' directory
  buildDir <- getwd()

  ## in the build directory, 'cat' all the src files together
  ## have a separate file for .c, .cpp files
  src_files <- list.files( full.names=TRUE,
    file.path( buildDir, "src" )
  )

  ## but don't concatenate init.c; copy it separately
  src_files <- grep("init.c", src_files, value=TRUE, fixed=TRUE, invert=TRUE)

  ## regex: the regex to match for picking up files
  ## ext: the file extension to use on the outputted file
  concatenate_src <- function(regex, ext) {
    files <- shQuote(grep(regex, src_files, value=TRUE))
    final <- shQuote( paste( sep='', buildDir, "/src/", pkg_name, "_", gsub("\\.", "", ext), ext ) )
    system( paste("touch", final) )
    files <- files[ files != final ]
    for( file in files ) {
      system( paste("cat", file, ">>", final) )
      system( paste("echo '' >>", final ) ) ## adds a newline just in case
      system( paste("rm", file) )
    }
  }

  cat("Concatenating source files...\n")
  concatenate_src("\\.c$", ".c")
  concatenate_src("\\.cpp$", ".cpp")

  ## remove all the .o, .so, .dll files
  for( file in grep("o$|so$|dll$", list.files( full.names=TRUE, file.path( buildDir, "src" ) ), value=TRUE ) ) {
    file.remove(file)
  }

}
