---
title: "class_7"
author: "Gary Le"
date: "April 25, 2018"
output: 
  html_document: 
    code_folding: show
    keep_md: yes
---

##Function warning() and stop() for errorhandling
*warning() will warn users when something is weird
*stop() will halt the function when encountered and display a message

Our rescale function from last class but modified to have error handling

```r
rescale <- function(x, na.rm=TRUE, plot=FALSE, ...) {
  # Our rescale function from the end of lecture 9
  
  #call. = FALSE means that the stop message won't show where the function failed in the error message.
  if(!is.numeric(x)){
    stop("The variable x should be numeric", call.= FALSE)
  }
  
  if(na.rm) {
    rng <-range(x, na.rm=TRUE)
  } else {
    rng <-range(x)
  }

  answer <- (x - rng[1]) / (rng[2] - rng[1])
  if(plot) { 
    plot(answer, ...) 
  }

  return(answer)
}
```


When "call. = FALSE"

```r
rescale("a")
```

```
## Error: The variable x should be numeric
```

When "call. = TRUE"

```r
rescale("a")
```

```
## Error: The variable x should be numeric
```

###Designing More Functions


```r
x <- c(1,2,NA,3,NA)
y<- c(NA,3,NA,3,4)

#checking for NA values in a vector
is.na(x)
```

```
## [1] FALSE FALSE  TRUE FALSE  TRUE
```

```r
is.na(y)
```

```
## [1]  TRUE FALSE  TRUE FALSE FALSE
```

```r
#1st attempt to see if we can find when both indexes are NA
#doesn't work, shows when both indexs are TRUE or when both are FALSE
is.na(x) == is.na(y)
```

```
## [1] FALSE  TRUE  TRUE  TRUE FALSE
```

```r
#This logical operator worked well
is.na(x) & is.na(y)
```

```
## [1] FALSE FALSE  TRUE FALSE FALSE
```

```r
#Which() shows which index is TRUE so...
which(is.na(x) & is.na(y))
```

```
## [1] 3
```

```r
#and it works! the line above returns the index where two vectors are NA

#how do we tell how many indexes are both missing info?
sum(is.na(x & is.na(y)))
```

```
## [1] 1
```

With our code chunks complete, let's write the function

```r
both_na <- function(x,y) {
  return(sum(is.na(x) & is.na(y)))
}
```

Testing function:

```r
both_na(x,y)
```

```
## [1] 1
```


```r
x1 <- c(NA,NA,NA)
y1 <- c(1,NA,NA)
y2 <- c(1, NA, NA, NA)

both_na(x1,y2)
```

```
## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
## shorter object length
```

```
## [1] 3
```
We see that the function doesn't handle mismatching vectors well. And we see that it found 3 mismatches between the vectors, meaning that it used the 1st index of x1 to compare with the 4th of y2. This is 'recycling'


```r
both_na2 <- function(x, y) {
  ## Check for NA elements in both input vectors and don't allow re-cycling 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  }
  sum( is.na(x) & is.na(y) )
}
```

Testing the new version:

```r
both_na2(x1,y2)
```

```
## Error: Input x and y should be vectors of the same length
```

Making the function more complex: Telling user both how many and where the vectors are missing values. Sample of more complex function structure.

```r
both_na3 <- function(x, y) {
  ## Print some info on where NA's are as well as the number of them 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  }
  
  # Calculating the truth vector for where both vectors are missing values
  na.in.both <- ( is.na(x) & is.na(y) )
  
  # Counting up number of missing values
  na.number  <- sum(na.in.both)
  # Showing where the values are missing in both
  na.which   <- which(na.in.both)

  message("Found ", na.number, " NA's at position(s):", 
          paste(na.which, collapse=", ") ) 
  
  # Returns a list of answers. Useful for searching through!
  return( list(number=na.number, which=na.which) )
}
```

Testing function v3

```r
ans <- both_na3(x,y)
```

```
## Found 1 NA's at position(s):3
```

```r
ans$which
```

```
## [1] 3
```

```r
ans$number
```

```
## [1] 1
```

###Yet an even more complex function

Finding the intersection between two data frames and accessing their values

Start with hardcoded example data

```r
## Find common genes in two lists

df1 <- data.frame(IDs=c("gene1", "gene2", "gene3"),
                  exp=c(2,1,1),
                  stringsAsFactors=FALSE)

df2 <- data.frame(IDs=c("gene2", "gene4", "gene3", "gene5"),
                  exp=c(-2, NA, 1, 2),
                  stringsAsFactors=FALSE)

df3 <- data.frame(IDs=c("gene2", "gene2", "gene5", "gene5"),
                  exp=c(-2, NA, 1, 2),
                  stringsAsFactors=FALSE)

#Start with a simpler proble: finding common IDs
a <- df1$IDs
b <- df2$IDs

#As our friend Google says, there is an existing function that can find the intersections
intersect(a,b)
```

```
## [1] "gene2" "gene3"
```

```r
#So how do we return the index of our intersectons? Google again...
#x %in% y shows if each indeces of x is in y
a %in% b
```

```
## [1] FALSE  TRUE  TRUE
```

```r
b %in% a
```

```
## [1]  TRUE FALSE  TRUE FALSE
```

Now that we have indices, we can retrieve the data from each vector

```r
a[a %in% b]
```

```
## [1] "gene2" "gene3"
```

```r
b[b %in% a]
```

```
## [1] "gene2" "gene3"
```

```r
#put them together...
#cbind() can also turn rows into columns
cbind(a[a %in% b], b[b %in% a])
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```


```r
cbind( c("Hello", "Help"), c("Please","Me"))
```

```
##      [,1]    [,2]    
## [1,] "Hello" "Please"
## [2,] "Help"  "Me"
```

Time to wrap the code chunk in a function:

```r
gene_intersect <- function(x, y) { 
   cbind( x[ x %in% y ], y[ y %in% x ] )
}
```

Test:

```r
gene_intersect(a,b)
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```

Testing again with data frames instead of vectors:

```r
gene_intersect2 <- function(df1, df2) { 
   cbind( df1[ df1$IDs %in% df2$IDs, ], 
          df2[ df2$IDs %in% df1$IDs, "exp"] )
}
```


```r
gene_intersect2(df1, df2)
```

```
##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
## 2 gene2   1                               -2
## 3 gene3   1                                1
```

This works, but the code looks ugly and the output is strangely formatted. The column retrieved is also hard-coded to "exp"

Let's add some flexibility.

```r
gene.colname = "IDs"
```


```r
gene_intersect3 <- function(df1, df2, gene.colname="IDs") { 
   cbind( df1[ df1[,gene.colname] %in% df2[,gene.colname], ], 
          exp2=df2[ df2[,gene.colname] %in% df1[,gene.colname], "exp"] )
}
```
The code should work now! But... it's still ugly.


```r
gene_intersect3(df1,df2)
```

```
##     IDs exp exp2
## 2 gene2   1   -2
## 3 gene3   1    1
```


Making it look better by adding local variables with nicer names

```r
gene_intersect4 <- function(df1, df2, gene.colname="IDs") { 

  df1.name <- df1[,gene.colname]
  df2.name <- df2[,gene.colname]

  df1.inds <- df1.name %in% df2.name
  df2.inds <- df2.name %in% df1.name

   cbind( df1[ df1.inds, ], 
          exp2=df2[ df2.inds, "exp"] )
}
```

Great! Now we can test and break the function!


```r
merge(df1,df2, by="IDs")
```

```
##     IDs exp.x exp.y
## 1 gene2     1    -2
## 2 gene3     1     1
```

#Packages

##Installing Bioconductor Packages
Bioconductor packages are tightly orchestrated together, with a lot of interdependencies.
To install them:


```r
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("GenomicFeatures")
```

