
 foo <- system("cd ../glmm/src; grep -E 'Calloc[(]|Free[(]' *.c",
     intern = TRUE)
 foo
 bar <- sub(": *|:\\t", ":", foo)
 bar
 baz <- strsplit(bar, ":")
 baz

 file <- sapply(baz, function(x) x[1])
 code <- sapply(baz, function(x) x[2])

 ufile <- unique(file)
 ufile

 code

 is.calloc <- grepl("Calloc", code)
 is.free <- grepl("Free", code)
 identical(is.calloc, ! is.free)

 get.calloc.pointer <- function(x)
     sapply(strsplit(x, "\\*| *= *"), function(x) x[2])
 get.calloc.pointer(code)

 get.free.pointer <- function(x)
     sub(");", "", sub("Free(", "", x, fixed = TRUE), fixed = TRUE)
 get.free.pointer(code)

 pointer <- ifelse(is.calloc, get.calloc.pointer(code), get.free.pointer(code))
 pointer

 woof <- list()
 for (u in ufile) {
     i <- file == u
     ipointer <- pointer[i]
     icalloc <- is.calloc[i]
     woof[[u]] <- list(calloc = ipointer[icalloc], free = ipointer[! icalloc])
 }

 woof

 sapply(woof, function(x) sort(x$calloc) == sort(x$free))

 sapply(woof, function(x)
     length(x$calloc) == length(x$free) && all(sort(x$calloc) == sort(x$free)))

 #### Oops !!!!
 woof$hess.c
 setdiff(woof$hess.c$calloc, woof$hess.c$free)
 setdiff(woof$hess.c$free, woof$hess.c$calloc)
