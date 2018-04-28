## Version 1.3.3

## Test environments
* local OS X install, R 3.3.1
* win-builder (devel and release)
* ubuntu 12.04 on travis-ci, (devel)

## R CMD check results

Status: OK

R CMD check results
0 errors | 0 warnings | 1 note 
checking dependencies in R code ... NOTE
Missing or unexported objects:
  ‘lmerTest::anova’ ‘lmerTest::summary’

This update introduces compatibility with current and previous version of lmerTest. These two functions were renamed in the new version of that package hence the note. 
  

## Resubmission

This is an minor update to fix error caused by new version of lmerTest.

## Downstream dependencies
There are currently no downstream dependencies for this package



## Version 1.3.2

## Test environments
* local OS X install, R 3.3.1
* win-builder (devel and release)
* ubuntu 12.04 on travis-ci, (devel)

## R CMD check results

Status: OK

R CMD check results
0 errors | 0 warnings | 0 notes

## Resubmission

This is an update 

## Downstream dependencies
There are currently no downstream dependencies for this package




## Version 1.2.1

## Test environments
* local OS X install, R 3.3.1
* win-builder (devel and release)
* ubuntu 12.04 on travis-ci, (oldrel, release and devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
  
  * checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Angel Garcia de la Garza <agarciadlg@gmail.com>’

New submission

This is my first submission

## Resubmission
This is a resubmission. In this version I have:

* Explained the meaning of NIfTI and included a reference in the description.

## Downstream dependencies
There are currently no downstream dependencies for this package
